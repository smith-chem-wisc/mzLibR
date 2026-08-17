# Label-free quantification with FlashLFQ, backed by mzLib.
#
# This is the most valuable tranche in mzLibR and the most dangerous, because almost everything
# that goes wrong here produces a plausible number rather than an error. The bake-off arm that
# did this task without the library reported 284 peptides, 173 proteins and 8 MBR transfers
# against a truth of 257, 140 and 2 - wrong in the believable direction, which is the one nobody
# checks.
#
# Three facts drive every design decision in this file:
#
#   1. **The peptide roll-up drops most MBR transfers.** Read `peaks`. On mzLib's own K562 pair
#      there are 140 true transfers and the peptide table shows 52 - a 63% under-count - and a
#      whole run's transfers can vanish entirely (run_3: 62 from peaks, 0 from the roll-up).
#   2. **A peptide intensity of 0 and a protein intensity of NA mean different things**, and R
#      is the only one of the three bindings whose type system can say so. See below.
#   3. **max_threads other than 1 makes results non-reproducible** (smith-chem-wisc/mzLib#1111).

# The wire's spelling of a number, independent of the R session's locale.
#
# The bridge parses with InvariantCulture, so a comma-decimal locale must never reach it.
# `format()` honours `options(OutDec)`, which a user in a European locale may well have set, so
# the decimal mark is stated explicitly rather than assumed.
wire_number <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(mzlib_usage_error(paste0(
      name, " must be a single finite number; got ", paste(deparse(value), collapse = " "), "."
    )))
  }
  # trimws because passing `decimal.mark` routes formatC through prettyNum, which pads the
  # result out to the digit width - "              10" rather than "10".
  trimws(formatC(value, format = "g", digits = 15L, decimal.mark = ".", big.mark = ""))
}

# ---------------------------------------------------------------- the spectra design

# Render the spectra argument as the bridge's tab-separated stdin lines.
#
# One run per line: `path[\tcondition[\tbiorep[\ttechrep[\tfraction]]]]`. Omitted trailing design
# fields are simply not written, and the bridge then applies the same defaults MetaMorpheus uses
# with no experimental-design file: blank condition, each run its own biological replicate,
# fraction 0, technical replicate 0.
#
# On stdin rather than argv because a real experiment's worth of paths goes past argv's ~32 KB
# ceiling.
flashlfq_spectra_stdin <- function(spectra) {
  if (is.character(spectra)) {
    spectra <- data.frame(path = spectra, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(spectra)) {
    stop(mzlib_usage_error(paste0(
      "spectra must be a character vector of mzML paths, or a data.frame with a 'path' column ",
      "and optional condition, biological_replicate, technical_replicate and fraction columns."
    )))
  }
  if (!"path" %in% names(spectra)) {
    stop(mzlib_usage_error("spectra is a data.frame with no 'path' column."))
  }
  if (nrow(spectra) == 0L) {
    stop(mzlib_usage_error("At least one spectra file is required."))
  }

  paths <- as.character(spectra$path)
  if (any(is.na(paths) | !nzchar(trimws(paths)))) {
    stop(mzlib_usage_error("Every spectra path must be a non-empty string."))
  }
  if (any(grepl("[\t\r\n]", paths))) {
    stop(mzlib_usage_error(
      "A spectra path may not contain a tab or a newline; those are the wire's field and record separators."
    ))
  }

  # mzML only. Rejecting here rather than letting mzLib fail deep inside a reader means the
  # message names the file and the reason.
  wrong_format <- paths[!grepl("\\.mzML$", paths, ignore.case = TRUE)]
  if (length(wrong_format) > 0L) {
    stop(mzlib_usage_error(paste0(
      "FlashLFQ reads mzML only; got '", wrong_format[1L], "'. Convert .raw or .d first ",
      "(MSConvert, or ThermoRawFileParser for .raw)."
    )))
  }

  design_field <- function(name) {
    if (!name %in% names(spectra)) {
      return(rep("", nrow(spectra)))
    }
    values <- spectra[[name]]
    rendered <- character(length(values))
    for (index in seq_along(values)) {
      value <- values[index]
      if (is.na(value)) {
        rendered[index] <- ""
        next
      }
      if (!is.numeric(value) || value != round(value) || value < 0) {
        stop(mzlib_usage_error(paste0(
          "spectra$", name, "[", index, "] must be a non-negative whole number; got ",
          paste(deparse(value), collapse = " "), "."
        )))
      }
      rendered[index] <- formatC(value, format = "d")
    }
    rendered
  }

  condition <- if ("condition" %in% names(spectra)) {
    ifelse(is.na(spectra$condition), "", as.character(spectra$condition))
  } else {
    rep("", nrow(spectra))
  }

  fields <- list(
    paths, condition,
    design_field("biological_replicate"),
    design_field("technical_replicate"),
    design_field("fraction")
  )

  vapply(seq_len(nrow(spectra)), function(row) {
    line <- vapply(fields, function(column) column[row], character(1L))
    # Drop trailing empties so a bare path stays a bare path on the wire.
    while (length(line) > 1L && !nzchar(line[length(line)])) {
      line <- line[-length(line)]
    }
    paste(line, collapse = "\t")
  }, character(1L))
}

# ---------------------------------------------------------------- argument assembly

flashlfq_build_args <- function(psms, normalize, ppm_tolerance, isotope_ppm_tolerance, integrate,
                                match_between_runs, mbr_ppm_tolerance, mbr_q_value_threshold,
                                use_shared_peptides_for_protein_quant, bayesian_protein_quant,
                                use_pep_q_value, max_threads, output_directory) {
  if (!is.character(psms) || length(psms) != 1L || is.na(psms) || !nzchar(trimws(psms))) {
    stop(mzlib_usage_error(
      "psms must be a path to a PSM result file, e.g. a MetaMorpheus AllPSMs.psmtsv."
    ))
  }

  flag <- function(name, value) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(mzlib_usage_error(paste0(name, " must be TRUE or FALSE.")))
    }
    isTRUE(value)
  }

  if (!is.numeric(max_threads) || length(max_threads) != 1L || is.na(max_threads) ||
    max_threads != round(max_threads) || max_threads == 0 || max_threads < -1) {
    stop(mzlib_usage_error(paste0(
      "max_threads must be a positive whole number, or -1 for every core; got ",
      paste(deparse(max_threads), collapse = " "), "."
    )))
  }

  args <- c("quant", "flashlfq", "--psms", trimws(psms))
  if (flag("normalize", normalize)) {
    args <- c(args, "--normalize")
  }
  args <- c(args, "--ppm", wire_number(ppm_tolerance, "ppm_tolerance"))
  args <- c(args, "--isotope-ppm", wire_number(isotope_ppm_tolerance, "isotope_ppm_tolerance"))
  if (flag("integrate", integrate)) {
    args <- c(args, "--integrate")
  }
  if (flag("match_between_runs", match_between_runs)) {
    args <- c(args, "--mbr")
  }
  args <- c(args, "--mbr-ppm", wire_number(mbr_ppm_tolerance, "mbr_ppm_tolerance"))
  args <- c(args, "--mbr-q", wire_number(mbr_q_value_threshold, "mbr_q_value_threshold"))
  if (flag("use_shared_peptides_for_protein_quant", use_shared_peptides_for_protein_quant)) {
    args <- c(args, "--shared-peptides")
  }
  if (flag("bayesian_protein_quant", bayesian_protein_quant)) {
    args <- c(args, "--bayesian")
  }
  if (flag("use_pep_q_value", use_pep_q_value)) {
    args <- c(args, "--use-pep-q")
  }
  args <- c(args, "--threads", formatC(max_threads, format = "d"))

  if (!is.null(output_directory)) {
    if (!is.character(output_directory) || length(output_directory) != 1L ||
      is.na(output_directory) || !nzchar(trimws(output_directory))) {
      stop(mzlib_usage_error(
        "output_directory must be a non-empty path, or NULL to write nothing."
      ))
    }
    args <- c(args, "--out", trimws(output_directory))
  }

  args
}

# ---------------------------------------------------------------- parsing

# Turn a `{run: value}` map into long rows, with an explicit answer for `null`.
#
# `missing` is the crux of this whole module and differs by table - 0 for a peptide, NA for a
# protein. See `flashlfq_parse()`.
flashlfq_long_map <- function(entry, field, missing) {
  values <- entry[[field]]
  if (!is.list(values) || length(values) == 0L) {
    return(list(file_name = character(0), value = numeric(0)))
  }
  list(
    file_name = names(values),
    value = vapply(values, function(value) {
      if (is.null(value) || (length(value) == 1L && is.logical(value) && is.na(value))) {
        missing
      } else {
        as.numeric(value)
      }
    }, numeric(1L), USE.NAMES = FALSE)
  )
}

flashlfq_parse_peptides <- function(entries) {
  empty <- data.frame(
    sequence = character(0), base_sequence = character(0), protein_groups = character(0),
    file_name = character(0), intensity = numeric(0), detection_type = character(0),
    stringsAsFactors = FALSE
  )
  if (length(entries) == 0L) {
    return(empty)
  }

  parts <- lapply(entries, function(entry) {
    # A peptide intensity of 0 means "not measured in this run" - it is never "could not be
    # resolved", which is a protein-level outcome only. Reading a null here as NA would put an
    # NA into a column where every arithmetic operation then propagates it, and mean() of a
    # peptide column would go from a number to NA across the board.
    intensities <- flashlfq_long_map(entry, "intensities", 0)
    if (length(intensities$file_name) == 0L) {
      return(NULL)
    }
    detection <- entry[["detection_types"]]
    data.frame(
      sequence = wire_field(entry, "sequence", "character", NA_character_),
      base_sequence = wire_field(entry, "base_sequence", "character", NA_character_),
      protein_groups = wire_field(entry, "protein_groups", "character", ""),
      file_name = intensities$file_name,
      intensity = intensities$value,
      detection_type = vapply(intensities$file_name, function(run) {
        value <- if (is.list(detection)) detection[[run]] else NULL
        if (is.null(value) || (length(value) == 1L && is.logical(value) && is.na(value))) {
          NA_character_
        } else {
          as.character(value)
        }
      }, character(1L), USE.NAMES = FALSE),
      stringsAsFactors = FALSE
    )
  })

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(empty)
  }
  do.call(rbind, parts)
}

flashlfq_parse_proteins <- function(entries) {
  empty <- data.frame(
    protein_group = character(0), gene_name = character(0), organism = character(0),
    file_name = character(0), intensity = numeric(0), stringsAsFactors = FALSE
  )
  if (length(entries) == 0L) {
    return(empty)
  }

  parts <- lapply(entries, function(entry) {
    # NA, not 0. A protein intensity of `null` on the wire means FlashLFQ's median-polish
    # produced NaN - it could not resolve a number at all - which is categorically different
    # from a measured zero. R is the only one of the three bindings whose types can say this
    # without prose: arithmetic propagates NA, so mean(proteins$intensity) returns NA rather
    # than a confidently wrong number, and na.rm = TRUE becomes a visible choice the analyst
    # makes rather than a default they never saw.
    #
    # mzLibR never passes na.rm = TRUE on the user's behalf, anywhere.
    #
    # The keys of `intensities` are sample labels, not file names: protein quant groups runs by
    # condition and biological replicate first. flashlfq_quantify() gives each run its own sample,
    # so the two coincide and the column is named file_name to match the peptides frame. A verb
    # that grouped runs into replicates would break that coincidence.
    intensities <- flashlfq_long_map(entry, "intensities", NA_real_)
    if (length(intensities$file_name) == 0L) {
      return(NULL)
    }
    data.frame(
      protein_group = wire_field(entry, "protein_group", "character", NA_character_),
      gene_name = wire_field(entry, "gene_name", "character", ""),
      organism = wire_field(entry, "organism", "character", ""),
      file_name = intensities$file_name,
      intensity = intensities$value,
      stringsAsFactors = FALSE
    )
  })

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(empty)
  }
  do.call(rbind, parts)
}

flashlfq_parse_peaks <- function(entries) {
  if (length(entries) == 0L) {
    return(data.frame(
      file_name = character(0), sequence = character(0), base_sequence = character(0),
      intensity = numeric(0), detection_type = character(0), retention_time = numeric(0),
      num_identifications = numeric(0), protein_groups = character(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    file_name = vapply(entries, wire_field, character(1L), "file_name", "character", NA_character_),
    sequence = vapply(entries, wire_field, character(1L), "sequence", "character", NA_character_),
    base_sequence = vapply(entries, wire_field, character(1L), "base_sequence", "character", NA_character_),
    intensity = vapply(entries, wire_field, numeric(1L), "intensity", "numeric", NA_real_),
    detection_type = vapply(entries, wire_field, character(1L), "detection_type", "character", NA_character_),
    retention_time = vapply(entries, wire_field, numeric(1L), "retention_time", "numeric", NA_real_),
    num_identifications = vapply(entries, wire_field, numeric(1L), "num_identifications", "numeric", NA_real_),
    protein_groups = vapply(entries, wire_field, character(1L), "protein_groups", "character", ""),
    stringsAsFactors = FALSE
  )
}

flashlfq_parse_spectra_files <- function(entries) {
  if (length(entries) == 0L) {
    return(data.frame(
      file_name = character(0), full_path = character(0), condition = character(0),
      biological_replicate = numeric(0), technical_replicate = numeric(0), fraction = numeric(0),
      peak_count = numeric(0), mbr_peak_count = numeric(0), stringsAsFactors = FALSE
    ))
  }
  numeric_field <- function(name) {
    vapply(entries, wire_field, numeric(1L), name, "numeric", NA_real_)
  }
  data.frame(
    file_name = vapply(entries, wire_field, character(1L), "file_name", "character", NA_character_),
    full_path = vapply(entries, wire_field, character(1L), "full_path", "character", ""),
    condition = vapply(entries, wire_field, character(1L), "condition", "character", ""),
    biological_replicate = numeric_field("biological_replicate"),
    technical_replicate = numeric_field("technical_replicate"),
    fraction = numeric_field("fraction"),
    peak_count = numeric_field("peak_count"),
    mbr_peak_count = numeric_field("mbr_peak_count"),
    stringsAsFactors = FALSE
  )
}

flashlfq_parse <- function(data) {
  entries <- function(name) {
    value <- data[[name]]
    if (is.list(value)) value else list()
  }

  parameters <- data[["parameters"]]
  if (!is.list(parameters)) {
    parameters <- list()
  }

  output <- data[["output_directory"]]
  if (is.null(output) || (length(output) == 1L && is.logical(output) && is.na(output))) {
    output <- NA_character_
  }

  structure(
    list(
      psm_file = as.character(wire_field(data, "psm_file", "character", NA_character_)),
      identification_count = as.numeric(wire_field(data, "identification_count", "numeric", NA_real_)),
      parameters = parameters,
      output_directory = as.character(output),
      spectra_files = flashlfq_parse_spectra_files(entries("spectra_files")),
      peptides = flashlfq_parse_peptides(entries("peptides")),
      proteins = flashlfq_parse_proteins(entries("proteins")),
      peaks = flashlfq_parse_peaks(entries("peaks"))
    ),
    class = "mzlibr_quant"
  )
}

# ---------------------------------------------------------------- the public surface

#' Quantify a search's peptides across mzML runs with FlashLFQ
#'
#' @param psms Path to a PSM result file. **Use a MetaMorpheus `.psmtsv` or `.osmtsv`.** Every
#'   run named in the file must have a matching mzML in `spectra`; FlashLFQ matches
#'   identifications to runs by base file name.
#' @param spectra The mzML runs. Either a character vector of paths, or a data.frame with a
#'   `path` column and optional `condition`, `biological_replicate`, `technical_replicate` and
#'   `fraction` columns. Omitted design fields get the defaults MetaMorpheus uses with no
#'   experimental-design file.
#'
#'   **mzML only.** `.raw` and `.d` are rejected up front; convert them first.
#' @param normalize Whether to normalise between runs.
#' @param ppm_tolerance Mass tolerance for peak-finding, in ppm.
#' @param isotope_ppm_tolerance Mass tolerance for isotope-envelope matching, in ppm.
#' @param integrate Whether to integrate peak areas rather than take apex intensity.
#' @param match_between_runs Whether to transfer identifications between runs.
#'
#'   MBR needs a **complete, balanced design** to work properly, and `mbr_q_value_threshold` is
#'   its FDR control - without it roughly 80% of transfers are false. Set `condition` and
#'   `biological_replicate` on `spectra` so FlashLFQ knows which runs are comparable.
#'
#'   **Whatever you do, count transfers from `peaks`, never from `peptides`.** See the return
#'   value.
#' @param mbr_ppm_tolerance Mass tolerance for match-between-runs, in ppm.
#' @param mbr_q_value_threshold FDR threshold for accepting a transfer.
#' @param use_shared_peptides_for_protein_quant Whether peptides shared between protein groups
#'   contribute to protein quantification.
#'
#'   `FALSE` is the default and it is the main reason protein intensities come back as **0**: on
#'   mzLib's K562 pair, **847** of 943 protein groups are 0 in both runs, mostly because their
#'   only evidence is shared peptides.
#' @param bayesian_protein_quant Whether to use the Bayesian protein quantification model.
#' @param use_pep_q_value Whether to use PEP q-values rather than q-values for filtering.
#' @param max_threads Worker threads. **Defaults to 1 here, which differs from pyMzLib's -1, and
#'   deliberately.**
#'
#'   With more than one thread the peptide roll-up nondeterministically drops MBR intensities,
#'   so identical inputs give different protein-level answers roughly **1 run in 6** - a
#'   borderline protein was unresolvable in 5 of 6 repeats. Any figure produced multithreaded
#'   may not reproduce (smith-chem-wisc/mzLib#1111). mzLibR warns if you set anything else.
#' @param output_directory Where FlashLFQ should write its TSV output, or `NULL` to write none.
#' @param timeout Seconds to allow, or `NULL` to wait as long as it takes.
#'
#' @return An `mzlibr_quant`: `psm_file`, `identification_count`, `parameters`,
#'   `output_directory`, and four tidy data.frames - `spectra_files`, `peptides`, `proteins`
#'   and `peaks`.
#'
#' @section Read peaks, not the peptide roll-up:
#'
#' `peptides` is FlashLFQ's roll-up and **it drops most match-between-runs transfers.** On
#' mzLib's own K562 pair there are **140** true transfers in `peaks` and the peptide table shows
#' **52** - a 63% under-count. Worse, it is not evenly spread: per run the peaks give run_3
#' **62** and run_4 **78**, while the roll-up gives run_3 **0** and run_4 52. Read only the
#' roll-up and MBR appears not to have worked at all in half the experiment.
#'
#' "Peptides quantified in both runs" is **257** from `peaks` and **169** from `peptides`.
#'
#' Nor can you reproduce the roll-up by pivoting `peaks` yourself: where a run has several peaks
#' for one peptide the roll-up reports one rather than their sum.
#'
#' @section `proteins$file_name` is a sample, not a file:
#'
#' FlashLFQ measures peptides in **files** but resolves proteins across **samples**, grouping runs
#' by condition and biological replicate before the median-polish roll-up. `flashlfq_quantify()`
#' gives every run its own sample, so `proteins$file_name` and `peptides$file_name` carry the same
#' run base names and the two frames join cleanly on it.
#'
#' The distinction only bites if runs are ever grouped into replicates: `proteins$file_name` would
#' then hold a sample label (`"condition_replicate"`) covering several runs, while
#' `peptides$file_name` stayed per run. The column name is kept for symmetry between the two
#' frames; read it as "the thing this intensity was measured over".
#'
#' @section What 0 and NA mean, and which is rare:
#'
#' `peptides$intensity` is **0** when the peptide was not measured in that run. It is never NA.
#'
#' `proteins$intensity` is **NA** when FlashLFQ could not resolve a number at all - its
#' median-polish produced NaN. Arithmetic propagates it, so `mean()` on a protein column returns
#' NA rather than a confidently wrong number, and `na.rm = TRUE` is a choice you make visibly.
#' mzLibR never applies it on your behalf.
#'
#' **NA is the rare outcome and 0 is the common one.** On the K562 pair, **2** proteins are NA
#' and **847** are 0 in both runs. "No usable number" is 849; "could not be resolved" is 2.
#'
#' @seealso [flashlfq_mbr_peaks()], [flashlfq_mbr_rescued_peptide_count()]
#' @export
flashlfq_quantify <- function(psms, spectra, normalize = FALSE, ppm_tolerance = 10,
                              isotope_ppm_tolerance = 5, integrate = FALSE,
                              match_between_runs = FALSE, mbr_ppm_tolerance = 10,
                              mbr_q_value_threshold = 0.05,
                              use_shared_peptides_for_protein_quant = FALSE,
                              bayesian_protein_quant = FALSE, use_pep_q_value = FALSE,
                              max_threads = 1, output_directory = NULL, timeout = NULL) {
  args <- flashlfq_build_args(
    psms, normalize, ppm_tolerance, isotope_ppm_tolerance, integrate, match_between_runs,
    mbr_ppm_tolerance, mbr_q_value_threshold, use_shared_peptides_for_protein_quant,
    bayesian_protein_quant, use_pep_q_value, max_threads, output_directory
  )
  stdin <- flashlfq_spectra_stdin(spectra)

  # The warning goes here, at the call that will produce the unreproducible answer, rather than
  # only in the help a user may never open.
  if (!identical(as.numeric(max_threads), 1)) {
    warning(
      "max_threads = ", max_threads, ": FlashLFQ's peptide roll-up nondeterministically drops ",
      "MBR intensities when multithreaded, so identical inputs give different protein-level ",
      "answers roughly 1 run in 6 (smith-chem-wisc/mzLib#1111). Use max_threads = 1 for any ",
      "result you intend to publish.",
      call. = FALSE
    )
  }

  flashlfq_parse(bridge_invoke(args, stdin = stdin, timeout = timeout))
}

#' Number of quantified peptides
#'
#' Distinct peptide sequences, not rows: `peptides` is long, with one row per peptide per run.
#'
#' @param results A [flashlfq_quantify()] result.
#' @return A single count.
#' @export
flashlfq_peptide_count <- function(results) {
  stopifnot(inherits(results, "mzlibr_quant"))
  length(unique(results$peptides$sequence))
}

#' Number of quantified protein groups
#'
#' @param results A [flashlfq_quantify()] result.
#' @return A single count.
#' @export
flashlfq_protein_count <- function(results) {
  stopifnot(inherits(results, "mzlibr_quant"))
  length(unique(results$proteins$protein_group))
}

#' Total match-between-runs peaks across every run
#'
#' Counts transferred **peaks**, not distinct peptides: one peptide rescued in two runs is two
#' peaks here. For "how many peptides did MBR rescue", use
#' [flashlfq_mbr_rescued_peptide_count()].
#'
#' Either way, do not count MBR from `peptides` - it under-counts, badly and unevenly. Zero
#' unless `match_between_runs` was on.
#'
#' @param results A [flashlfq_quantify()] result.
#' @return A single count.
#' @export
flashlfq_mbr_peak_count <- function(results) {
  stopifnot(inherits(results, "mzlibr_quant"))
  sum(results$spectra_files$mbr_peak_count, na.rm = TRUE)
}

#' Exactly the peaks transferred by match-between-runs
#'
#' @param results A [flashlfq_quantify()] result.
#' @return The rows of `results$peaks` whose `detection_type` is `"MBR"`.
#' @export
flashlfq_mbr_peaks <- function(results) {
  stopifnot(inherits(results, "mzlibr_quant"))
  results$peaks[!is.na(results$peaks$detection_type) & results$peaks$detection_type == "MBR", ]
}

#' Distinct peptides rescued by match-between-runs
#'
#' **What this actually computes: the number of distinct modified sequences among the MBR
#' peaks.** That is stated in code terms on purpose, because the prose definition a reader
#' supplies - "peptides quantified in at least one run *only* by MBR" - is subtly different and
#' gives a different number. On the K562 pair this is **140** while the strict reading is
#' **135**; the five that differ have both an MBR peak and a zero-intensity MSMS peak in the
#' same run.
#'
#' It equals [flashlfq_mbr_peak_count()] only when no peptide was rescued in more than one run.
#'
#' @param results A [flashlfq_quantify()] result.
#' @return A single count.
#' @export
flashlfq_mbr_rescued_peptide_count <- function(results) {
  length(unique(flashlfq_mbr_peaks(results)$sequence))
}

#' Print a quantification result
#'
#' A compact summary that names the two things most likely to be misread: how many
#' match-between-runs transfers the peptide roll-up is hiding, and how many protein intensities
#' are `NA` rather than `0`.
#'
#' @param x A [flashlfq_quantify()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_quant <- function(x, ...) {
  cat("<mzlibr_quant> ", basename(x$psm_file), ", ",
    format(x$identification_count), " identifications\n",
    sep = ""
  )
  cat("  ", nrow(x$spectra_files), " runs, ",
    flashlfq_peptide_count(x), " peptides, ",
    flashlfq_protein_count(x), " protein groups, ",
    nrow(x$peaks), " peaks\n",
    sep = ""
  )

  from_peaks <- nrow(flashlfq_mbr_peaks(x))
  if (from_peaks > 0L) {
    from_rollup <- sum(!is.na(x$peptides$detection_type) & x$peptides$detection_type == "MBR")
    cat("  MBR: ", from_peaks, " transfers in peaks, ",
      flashlfq_mbr_rescued_peptide_count(x), " distinct peptides\n",
      sep = ""
    )
    if (from_rollup < from_peaks) {
      # Where the mistake is made: someone printing a result and reading off an MBR number.
      cat("  ! the peptide roll-up shows only ", from_rollup,
        " of those ", from_peaks, " transfers - read peaks, not peptides.\n",
        "    See ?flashlfq_quantify (smith-chem-wisc/mzLib#1111 is a separate issue).\n",
        sep = ""
      )
    }
  }

  unresolved <- sum(is.na(x$proteins$intensity))
  zeroed <- sum(!is.na(x$proteins$intensity) & x$proteins$intensity == 0)
  if (unresolved > 0L || zeroed > 0L) {
    cat("  proteins: ", unresolved, " NA (could not be resolved), ",
      zeroed, " zero (not measured)\n",
      sep = ""
    )
  }

  threads <- x$parameters[["max_threads"]]
  if (!is.null(threads) && !identical(as.numeric(threads), 1)) {
    cat("  ! max_threads = ", format(threads),
      " - this result may not reproduce (mzLib#1111).\n",
      sep = ""
    )
  }
  invisible(x)
}
