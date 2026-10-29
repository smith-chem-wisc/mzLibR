# Spectral libraries: mzLib's binary `.msl` format, read with mzLib's own MslLibrary.
#
#   library_read()        one library, one row per precursor, or per fragment ion on request
#   library_read_many()   many libraries into one long table, in ONE bridge call
#   library_calibrate_rt() each run's minutes onto a library's iRT scale, from anchors, and back
#
# Ported from pyMzLib's `library.py` (read, read_many, LibraryRecords, LibraryBatch), which decided
# the verb, the wire fields and the caveats. What changes here is only the projection into R's
# idiom, exactly as for the readers it sits beside: the table is a data.frame, `records`; `NA` where
# pyMzLib has `None`; `precursor_index` and `source_index` 1-based, as every position in R is; and
# the per-library facts of a batch a `files` data.frame (R/bulk.R).
#
# A library holds, for each precursor - a peptide or proteoform at one charge - the m/z, retention
# time and ion mobility a search expects to see it at, and the fragment ions it expects in the MS2
# spectrum. mzLib's "not available" values (an ion mobility of 0.0, a NaN q-value, an NCE of 0, an
# empty protein accession) arrive as NA from the bridge; nothing here converts or fills anything.

# The pyMzLib release whose bridge first dispatches `library read`.
LIBRARY_SINCE <- "0.6.0"

# ---------------------------------------------------------------- argument shaping

# `--fragments` and `--index-only-threshold`, validated before the bridge is spawned.
library_options <- function(fragments, index_only_threshold) {
  if (!is.logical(fragments) || length(fragments) != 1L || is.na(fragments)) {
    stop(mzlib_usage_error(paste0(
      "fragments must be TRUE or FALSE; got ", paste(deparse(fragments), collapse = " "), "."
    )))
  }
  args <- if (fragments) "--fragments" else character(0)
  if (!is.null(index_only_threshold)) {
    if (!is.numeric(index_only_threshold) || length(index_only_threshold) != 1L ||
      is.na(index_only_threshold) || !is.finite(index_only_threshold) ||
      index_only_threshold != round(index_only_threshold) || index_only_threshold < 0) {
      stop(mzlib_usage_error(paste0(
        "index_only_threshold must be a whole number of bytes, 0 or more, or NULL for mzLib's ",
        "default; got ", paste(deparse(index_only_threshold), collapse = " "), "."
      )))
    }
    # sprintf, not formatC(format = "d"): a threshold past 2 GiB is past an R integer.
    args <- c(args, "--index-only-threshold", sprintf("%.0f", index_only_threshold))
  }
  args
}

# The one-library request: the shared readers window, under the library module, then the options.
library_build_read_args <- function(path, fragments, index_only_threshold, limit, offset, out) {
  if (is.character(path) && length(path) > 1L) {
    stop(mzlib_usage_error(paste0(
      "library_read() reads one library; got ", length(path), " paths. Read several in one call with ",
      "library_read_many()."
    )))
  }
  extra <- library_options(fragments, index_only_threshold)
  c(readers_build_read_args(path, limit, offset, out, "read", module = "library"), extra)
}

# ---------------------------------------------------------------- parsing

# `precursor_index` 1-based, as every position in R is, in place.
library_one_based <- function(table) {
  if (!is.null(table) && "precursor_index" %in% names(table)) {
    table$precursor_index <- as.numeric(table$precursor_index) + 1
  }
  table
}

library_parse_records <- function(data) {
  out <- c(
    readers_parse_common(data),
    readers_parse_block(data),
    list(
      format_version = as.numeric(wire_field(data, "format_version", "numeric", NA_real_)),
      is_compressed = wire_field(data, "is_compressed", "logical", NA),
      is_index_only = wire_field(data, "is_index_only", "logical", NA),
      target_count = as.numeric(wire_field(data, "target_count", "numeric", NA_real_)),
      decoy_count = as.numeric(wire_field(data, "decoy_count", "numeric", NA_real_)),
      proteoform_count = as.numeric(wire_field(data, "proteoform_count", "numeric", NA_real_)),
      min_precursor_mz = as.numeric(wire_field(data, "min_precursor_mz", "numeric", NA_real_)),
      max_precursor_mz = as.numeric(wire_field(data, "max_precursor_mz", "numeric", NA_real_)),
      fragments_included = isTRUE(data[["fragments_included"]]),
      index_only_threshold = as.numeric(wire_field(data, "index_only_threshold", "numeric", NA_real_)),
      row_count = as.numeric(wire_field(data, "row_count", "numeric", NA_real_))
    )
  )
  out$records <- library_one_based(out$records)
  structure(out, class = "mzlibr_library_records")
}

library_parse_batch <- function(data) {
  batch <- bulk_parse_batch(data, "library read", "mzlibr_library_batch")
  batch$records <- library_one_based(batch$records)
  batch
}

# ---------------------------------------------------------------- the public surface

#' Read an mzLib spectral library (`.msl`): one row per precursor, or per fragment ion
#'
#' Answers **what is in this spectral library?** For each precursor - a peptide or proteoform at
#' one charge - the m/z, retention time and ion mobility a search expects to see it at, its labels
#' and provenance, and on request the fragment ions it expects in the MS2 spectrum. MetaMorpheus and
#' mzLib write libraries in mzLib's binary `.msl` format; this reads one with mzLib's own reader,
#' `MslFileTypeHandler.Open` and `MslLibrary.GetEntry`.
#'
#' Rows come in the library's stored order (`precursor_index`), not by m/z: sort on `precursor_mz`
#' for that. mzLib's "not available" values - an ion mobility of 0.0, a NaN q-value, an NCE of 0,
#' an empty protein accession - arrive as `NA`.
#'
#' Large libraries open **index-only**: precursors in memory, fragment blocks left on disk. mzLib
#' decides by file size (1 GiB by default in the pinned build); pass `index_only_threshold` to move
#' the line. The rows are the same either way.
#'
#' @param path Path to a `.msl` file. Any other extension is an `mzlib_usage_error`: text libraries
#'   (`.msp`) are not read here.
#' @param fragments One row per **fragment ion** instead of one per precursor, with the precursor's
#'   columns repeated on each row. `FALSE` by default: a library has tens of fragments per
#'   precursor, so pair it with `limit` or `out` for a large one.
#' @param index_only_threshold A library this many **bytes** or larger opens index-only; a smaller
#'   one loads whole. `NULL` (the default) sends nothing, so mzLib's own default applies:
#'   `MslFileTypeHandler.DefaultIndexOnlyThresholdBytes`, 1073741824 bytes (1 GiB) in the pinned
#'   build. `0` opens every library index-only. The rows are the same either way.
#' @param limit Maximum **precursors** to return; with `fragments` each gives one row per fragment
#'   ion. `NULL`, the default, returns all of them.
#' @param offset **Precursors** to skip, in `precursor_index` order.
#' @param out Write a tab-separated table here and return only a summary. The path must end in
#'   `.tsv` (any case); any other is an `mzlib_usage_error` before a read. The file keeps the wire's
#'   0-based `precursor_index`.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait.
#'
#' @return An `mzlibr_library_records`. `record_count` counts the precursors in the whole library,
#'   before `offset` and `limit`; `returned_count` the precursors returned (precursors, not rows: the
#'   unit `limit` and `offset` count in), 0 with `out`; `row_count` the rows of `records`, one per
#'   fragment ion with `fragments`; `offset` the precursors skipped; `truncated` whether precursors
#'   were left behind.
#'
#'   The library's facts: `reader` (`"MslLibrary"`); `format_version`, the file's `.msl` format
#'   version; `is_compressed` (a compressed library always loads whole) and `is_index_only`;
#'   `target_count`, precursors that are not decoys, entrapment targets included; `decoy_count`,
#'   decoy precursors, entrapment decoys included; `proteoform_count`, precursors whose
#'   `molecule_type` is `"Proteoform"`; and `min_precursor_mz` and `max_precursor_mz`, in m/z, `NA`
#'   when the library has no precursors. `fragments_included` says whether the rows are fragment
#'   ions, and `index_only_threshold` the threshold applied, in bytes: yours, or mzLib's default.
#'
#'   `file_type` is always `NA` (`.msl` is not one of mzLib's `SupportedFileType` members);
#'   `rows_not_read` is always `NA`: rows are not counted for a binary library, because mzLib checks
#'   the whole file before returning any precursor;
#'   `retention_time_unit` is `NA` at the top level, because the unit varies by precursor - it is a
#'   column. `absent_fields` names the columns this library has no source for, so every value is
#'   `NA` (`ion_mobility` when no precursor has one, `protein_accession` and `protein_name` without
#'   a protein table, `gene_name` without gene names); `failed_fields` the columns whose read threw on
#'   some rows; `excluded_fields` is always empty; `caveats` says what the table cannot be trusted
#'   to mean. `output` says where the table was written, or is `NULL`.
#'
#'   `records` is a data.frame, or `NULL` with `out`. One row per precursor: `precursor_index`
#'   (**1-based in R**, the position in the library and the row order), `full_sequence` (mzLib
#'   bracket notation), `base_sequence`, `precursor_mz` in m/z, `charge_state` (a charge),
#'   `retention_time` in iRT or minutes **per row** and `retention_time_unit` (`"iRT"` or
#'   `"minutes"`), `ion_mobility` in 1/K0 (V s/cm2; `NA` where mzLib stores 0.0), `is_decoy`,
#'   `is_entrapment`, `decoy_contaminant_target` (`"T"`, `"D"`, `"ET"` or `"ED"`),
#'   `protein_accession`, `protein_name` and `gene_name` (`NA` where mzLib stores an empty string),
#'   `is_proteotypic`, `q_value` (a fraction, 0 to 1; `NA` where mzLib stores NaN, a predicted
#'   library), `elution_group_id`, `molecule_type`, `source_type` (`"Predicted"`, `"Empirical"` or
#'   `"EmpiricalRefined"`), `dissociation_type` and `nce` (an NCE, a whole number; `NA` where mzLib
#'   stores 0).
#'
#'   With `fragments = TRUE`, one row per fragment ion: those columns, then `fragment_count` (fragment
#'   ions of the precursor), `fragment_mz` in m/z, `fragment_intensity` (relative intensity within
#'   the precursor, 0 to 1), `product_type`, `secondary_product_type` and
#'   `secondary_fragment_number` (`NA` for a terminal ion), `fragment_number`, `residue_position`
#'   (`NA` for a custom-loss fragment, whose stored position mzLib reads back as 0),
#'   `fragment_charge` (a charge), `neutral_loss` in Da (0 is no loss), `is_internal`,
#'   `is_diagnostic`, `exclude_from_quant` and `annotation`. A precursor with no fragment ions is one
#'   row with those cells `NA`.
#'
#' @section Three things to know before you use the numbers:
#'
#' **Retention time has a unit per row.** A library stores either an indexed retention time (iRT)
#' or run minutes for each precursor (`MslLibraryEntry.RtIsCalibrated`), and one library can hold
#' both. Read `retention_time_unit` on each row; there is no conversion between them without a
#' calibration against a run, which [library_calibrate_rt()] fits.
#'
#' **"Not available" is `NA`, not a number.** A text library converted to `.msl` (from `.msp`) has
#' no NCE, dissociation type, q-value, ion mobility or protein, so those are `NA` (and
#' `dissociation_type` is `"Unknown"`).
#'
#' **The library is unfiltered.** Decoys (`D`, `ED`) and entrapment precursors (`ET`, `ED`) are
#' rows; filter on `is_decoy` or `decoy_contaminant_target` before counting targets.
#'
#' @seealso [library_read_many()], [readers_read_spectra()] for the runs a library is searched
#'   against
#' @spec library.read
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' lib <- library_read("synthetic_library.msl")
#' lib
#' lib$records[, c("precursor_index", "full_sequence", "decoy_contaminant_target",
#'                 "retention_time", "retention_time_unit")]
#' # A predicted precursor: mzLib's "not available" values are NA.
#' lib$records[5, c("ion_mobility", "q_value", "nce", "protein_accession")]
#'
#' # One precursor's fragment ions:
#' ions <- library_read("synthetic_library.msl", fragments = TRUE, limit = 1)
#' c(precursors = ions$returned_count, rows = ions$row_count)
#' ions$records[, c("annotation", "fragment_mz", "fragment_intensity", "neutral_loss")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
library_read <- function(path, fragments = FALSE, index_only_threshold = NULL, limit = NULL, offset = 0,
                         out = NULL, timeout = NULL) {
  # Every argument is checked before any process is started, the verb check included.
  args <- library_build_read_args(path, fragments, index_only_threshold, limit, offset, out)
  bridge_require_verb("library read", LIBRARY_SINCE)
  library_parse_records(bridge_invoke(args, timeout = timeout))
}

#' Read many spectral libraries into one long table, in one bridge call
#'
#' The many-files form of [library_read()]: a predicted library and the empirical one refined from
#' it, say. ONE bridge process reads every library, `threads` at a time, into ONE table, in input
#' order; each library's facts (counts, load mode, format version) are its row of `files`.
#'
#' @param paths A character vector of `.msl` paths. Read in this order; a file's position is its
#'   `source_index`. A repeated path is refused rather than read twice.
#' @param fragments One row per fragment ion, in every library. `FALSE` by default.
#' @param index_only_threshold The size in **bytes** at and above which a library opens index-only,
#'   for every library. `NULL` (the default) leaves it to mzLib (1 GiB in the pinned build).
#' @param out Write the long table here as tab-separated text, one library at a time, in order, and
#'   return only a summary. A batch that stops on an error removes its partial table. The path must
#'   end in `.tsv` (any case).
#' @param threads Libraries read at once, a count of files: `1` (the default), more, or `-1` for one
#'   per core. The result is **byte-identical at any value**, so this trades memory for speed and
#'   never changes an answer. The default is 1 because a library loaded whole is held in memory while
#'   its rows are written.
#' @param on_error `"fail"` (the default) stops at the first library that cannot be read and
#'   raises, naming it; `"skip"` records the failure in that library's row of `files` and reads the
#'   rest.
#' @param timeout Seconds to allow for the whole batch, or `NULL` (the default) to wait.
#'
#' @return An `mzlibr_library_batch` (an `mzlibr_batch`). `records` begins with `source_index`
#'   (1-based, into `paths`) and `source_path`, then the columns of [library_read()], rows grouped
#'   by library in input order, each library's in `precursor_index` order (1-based) - **whatever**
#'   `threads` **was**. Their units are as there: `precursor_mz` and `fragment_mz` in m/z,
#'   `charge_state` and `fragment_charge` a charge, `retention_time` in iRT or minutes per row (read
#'   `retention_time_unit`), `ion_mobility` in 1/K0, `q_value` a fraction, `nce` an NCE,
#'   `fragment_count` in fragment ions, `fragment_intensity` a relative intensity (0 to 1) and
#'   `neutral_loss` in Da. `file_count`, `read_count` and `failed_count` count files; `record_count`
#'   counts precursors summed over the libraries read, `returned_count` the precursors returned and
#'   `row_count` the rows (one per fragment ion with `fragments`), both 0 with `out`.
#'   `fragments_included` and `index_only_threshold` (in bytes) are as for [library_read()].
#'
#'   `files` has one row per input: the facts [library_read()] reports at its top level, under the
#'   same names - `reader`, `format_version`, `is_compressed`, `is_index_only`, `target_count`,
#'   `decoy_count`, `proteoform_count`, `min_precursor_mz`, `max_precursor_mz` (m/z), `record_count`
#'   (precursors), `caveats`, `absent_fields` - and `error_kind`, `error_type` and `error_message`,
#'   `NA` for a library that was read. **For a library that could not be read, the facts it could
#'   not establish are `NA`, not zero**: nothing was counted.
#'
#' @seealso [library_read()]
#' @spec library.read bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' batch <- library_read_many(c("synthetic_library.msl", "myPrositLib.msl"))
#' batch
#' batch$files[, c("record_count", "target_count", "decoy_count", "format_version")]
#' batch$files$absent_fields[[2]]
#' table(batch$records$source_index)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
library_read_many <- function(paths, fragments = FALSE, index_only_threshold = NULL, out = NULL, threads = 1,
                              on_error = "fail", timeout = NULL) {
  extra <- library_options(fragments, index_only_threshold)
  request <- bulk_request(c("library", "read"), paths, threads, on_error, out, "library_read_many")
  bridge_require_verb("library read", LIBRARY_SINCE)
  library_parse_batch(bridge_invoke(c(request$args, extra), stdin = request$stdin, timeout = timeout))
}

# ---------------------------------------------------------------- calibrate-rt (pyMzLib 0.7.0)

# Ported from pyMzLib's `library.py` (calibrate_rt, RtCalibration, RtCalibrationGroup, FitFailure),
# which decided the verb, the wire fields and the caveats. In R the anchors and queries are
# data.frames (or lists of `list(group, rt_minutes, irt)`), `NA` is a time not given, and the
# per-group facts are a `groups` data.frame whose `error` is split into error_kind, error_type and
# error_message, as a batch's `files` is (R/bulk.R).

# The first pyMzLib whose bridge has `library calibrate-rt`: its spec's since.pymzlib.
LIBRARY_CALIBRATE_SINCE <- "0.7.0"

# The values of `kind` in library_calibrate_rt(): mzLib's IrtCalibrationKind.
LIBRARY_CALIBRATION_KINDS <- c("linear", "lowess")

# One time on a stdin line: NA and NaN are blank, a number its exact text.
library_time_cell <- function(value, what) {
  if (length(value) != 1L || !(is.numeric(value) || (is.logical(value) && is.na(value)))) {
    stop(mzlib_usage_error(paste0(what, " must be a number or NA; got ", paste(deparse(value), collapse = " "), ".")))
  }
  if (is.na(value)) "" else stats_number_text(as.numeric(value))
}

# The entries of `anchors` or `queries` as list(group, rt_minutes, irt): a data.frame's rows, by
# column name (a missing time column is NA throughout), or a list's elements, by position or name.
library_calibration_entries <- function(entries, what) {
  if (is.data.frame(entries)) {
    if (!"group" %in% names(entries)) {
      stop(mzlib_usage_error(paste0(
        what, " needs a `group` column; got ", paste(names(entries), collapse = ", "), "."
      )))
    }
    column <- function(name) if (name %in% names(entries)) entries[[name]] else rep(NA_real_, nrow(entries))
    rt <- column("rt_minutes")
    irt <- column("irt")
    return(lapply(seq_len(nrow(entries)), function(i) list(entries$group[[i]], rt[[i]], irt[[i]])))
  }
  if (!is.list(entries)) {
    stop(mzlib_usage_error(paste0(
      what, " must be a data.frame with group, rt_minutes and irt columns, or a list of ",
      "list(group, rt_minutes, irt); got ", class(entries)[1L], "."
    )))
  }
  lapply(seq_along(entries), function(i) {
    entry <- entries[[i]]
    if (is.list(entry) && !is.null(names(entry)) && "group" %in% names(entry)) {
      return(list(entry[["group"]],
        if (is.null(entry[["rt_minutes"]])) NA_real_ else entry[["rt_minutes"]],
        if (is.null(entry[["irt"]])) NA_real_ else entry[["irt"]]))
    }
    if (!is.list(entry) || length(entry) != 3L) {
      stop(mzlib_usage_error(paste0(
        what, "[[", i, "]] must be list(group, rt_minutes, irt); got ", paste(deparse(entry), collapse = " "), "."
      )))
    }
    entry
  })
}

# One stdin line per anchor or query: group, rt_minutes, irt, tab-separated.
library_calibration_lines <- function(entries, what, query) {
  if (is.null(entries)) {
    return(character(0))
  }
  parsed <- library_calibration_entries(entries, what)
  vapply(seq_along(parsed), function(i) {
    values <- parsed[[i]]
    group <- values[[1L]]
    if (is.factor(group)) group <- as.character(group)
    if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(trimws(group))) {
      stop(mzlib_usage_error(paste0(what, " entry ", i, " needs a non-empty group name; got ",
        paste(deparse(group), collapse = " "), ".")))
    }
    if (grepl("[\t\r\n]", group)) {
      stop(mzlib_usage_error(paste0(what, " entry ", i, "'s group name holds a tab or newline: '", group, "'.")))
    }
    if (identical(trimws(group), PROTEINS_SECTION)) {
      stop(mzlib_usage_error(paste0(
        what, " entry ", i, "'s group name is '--', which separates anchors from queries on the wire."
      )))
    }
    cells <- c(
      library_time_cell(values[[2L]], paste0(what, " entry ", i, " rt_minutes")),
      library_time_cell(values[[3L]], paste0(what, " entry ", i, " irt"))
    )
    if (query && (nzchar(cells[[1L]]) == nzchar(cells[[2L]]))) {
      stop(mzlib_usage_error(paste0(
        what, " entry ", i, " must give exactly one of rt_minutes and irt, the other NA: the missing ",
        "one is what the calibration computes."
      )))
    }
    paste(c(group, cells), collapse = "\t")
  }, character(1L))
}

# `--name value` when given; nothing when NULL, so mzLib's default applies.
library_calibration_option <- function(name, value, integer) {
  if (is.null(value)) {
    return(character(0))
  }
  argument <- gsub("-", "_", name, fixed = TRUE)
  ok <- is.numeric(value) && length(value) == 1L && !is.na(value) && (!integer || value == round(value))
  if (!ok) {
    stop(mzlib_usage_error(paste0(
      argument, " must be ", if (integer) "a whole number" else "a number", ", or NULL for mzLib's default; got ",
      paste(deparse(value), collapse = " "), "."
    )))
  }
  c(paste0("--", name), if (integer) sprintf("%.0f", value) else stats_number_text(value))
}

library_build_calibrate_request <- function(anchors, queries, kind, minimum_anchors, bandwidth, knots,
                                            robustness_iterations, on_error) {
  if (!is.null(kind) && (!is.character(kind) || length(kind) != 1L || !kind %in% LIBRARY_CALIBRATION_KINDS)) {
    stop(mzlib_usage_error(paste0(
      "kind must be one of ", paste(LIBRARY_CALIBRATION_KINDS, collapse = ", "),
      " (mzLib's IrtCalibrationKind), or NULL for mzLib's default; got ", paste(deparse(kind), collapse = " "), "."
    )))
  }
  if (!is.character(on_error) || length(on_error) != 1L || !on_error %in% c("fail", "skip")) {
    stop(mzlib_usage_error(paste0("on_error must be 'fail' or 'skip'; got ", paste(deparse(on_error), collapse = " "), ".")))
  }
  args <- c("library", "calibrate-rt")
  if (!is.null(kind)) {
    args <- c(args, "--kind", kind)
  }
  args <- c(
    args,
    library_calibration_option("minimum-anchors", minimum_anchors, integer = TRUE),
    library_calibration_option("bandwidth", bandwidth, integer = FALSE),
    library_calibration_option("knots", knots, integer = TRUE),
    library_calibration_option("robustness-iterations", robustness_iterations, integer = TRUE)
  )
  if (identical(on_error, "skip")) {
    args <- c(args, "--on-error", "skip")
  }
  anchor_lines <- library_calibration_lines(anchors, "anchors", query = FALSE)
  if (length(anchor_lines) == 0L) {
    stop(mzlib_usage_error("library_calibrate_rt() needs at least one anchor."))
  }
  query_lines <- library_calibration_lines(queries, "queries", query = TRUE)
  lines <- c(anchor_lines, if (length(query_lines) > 0L) c(PROTEINS_SECTION, query_lines))
  list(args = args, stdin = paste0(paste(lines, collapse = "\n"), "\n"))
}

# The per-group facts as a data.frame, one row per group, `error` split in three.
library_parse_calibration_groups <- function(value) {
  groups <- if (is.list(value)) value else list()
  text <- function(entry, name) wire_field(entry, name, "character", NA_character_)
  number <- function(entry, name) as.numeric(wire_field(entry, name, "numeric", NA_real_))
  error <- function(entry, name) {
    failure <- entry[["error"]]
    if (is.list(failure)) wire_field(failure, name, "character", NA_character_) else NA_character_
  }
  data.frame(
    group = vapply(groups, text, character(1L), "group"),
    kind = vapply(groups, text, character(1L), "kind"),
    anchor_count = vapply(groups, number, numeric(1L), "anchor_count"),
    query_count = vapply(groups, number, numeric(1L), "query_count"),
    residual_sd = vapply(groups, number, numeric(1L), "residual_sd"),
    first_anchor_minutes = vapply(groups, number, numeric(1L), "first_anchor_minutes"),
    last_anchor_minutes = vapply(groups, number, numeric(1L), "last_anchor_minutes"),
    error_kind = vapply(groups, error, character(1L), "kind"),
    error_type = vapply(groups, error, character(1L), "type"),
    error_message = vapply(groups, error, character(1L), "message"),
    stringsAsFactors = FALSE
  )
}

library_parse_calibration <- function(data) {
  structure(
    list(
      kind = wire_field(data, "kind", "character", NA_character_),
      minimum_anchors = stats_count(data, "minimum_anchors"),
      bandwidth = stats_count(data, "bandwidth"),
      knots = stats_count(data, "knots"),
      robustness_iterations = stats_count(data, "robustness_iterations"),
      on_error = wire_field(data, "on_error", "character", "fail"),
      group_count = stats_count(data, "group_count"),
      fitted_count = stats_count(data, "fitted_count"),
      failed_count = stats_count(data, "failed_count"),
      anchor_count = stats_count(data, "anchor_count"),
      query_count = stats_count(data, "query_count"),
      row_count = stats_count(data, "row_count"),
      groups = library_parse_calibration_groups(data[["groups"]]),
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_rt_calibration"
  )
}

#' Fit each run's map from minutes onto a library's iRT scale, and convert times either way
#'
#' Answers **where should this library precursor elute in my run?** and **what is this run's minute
#' on the library's scale?** A spectral library stores retention time on a run-independent scale,
#' iRT (Escher et al. 2012); a run reports minutes. mzLib's `IrtCalibration.Fit` maps one onto the
#' other from **anchors**: precursors you identified confidently in the run, each with its observed
#' apex time in minutes and its library iRT. The fit is robust: it starts from binned medians and
#' downweights anchors far from the curve (Tukey bisquare), so some wrong identifications do not
#' bend it. It is strictly increasing, so it inverts. By default it is a LOWESS curve (Cleveland
#' 1979) that follows a gradient that bends; `kind = "linear"` fits one straight line.
#'
#' Each **group** (a run, typically) is fitted on its own anchors, and the **queries** you send for
#' it are converted with that fit: give a time in minutes to get its iRT, or an iRT to get the
#' minutes at which that precursor should elute in this run.
#'
#' **You choose the anchors.** mzLib applies no q-value cutoff, and mzLibR adds none. **Send the
#' queries in the same call as their anchors**: mzLib keeps the fitted curve's points private, so a
#' calibration cannot be returned and used later. Refitting the same anchors gives the same map,
#' because mzLib sorts them first.
#'
#' @param anchors One row per anchor: a data.frame with columns `group` (any name: the run),
#'   `rt_minutes` (the observed retention time, in **minutes**) and `irt` (the library's indexed
#'   retention time, in **iRT**), or a list of `list(group, rt_minutes, irt)`. They travel on stdin,
#'   one line each. mzLib refuses an `NA` or `NaN` time.
#' @param queries Times to convert, in the same shape: each row gives its `group` and exactly one of
#'   `rt_minutes` (minutes, converted to iRT) or `irt` (iRT, converted to minutes), the other `NA`.
#'   A data.frame may leave out the column it never gives. Each group's queries use that group's
#'   fit. `NULL` (the default) converts nothing: you get each group's fit and its anchors' residuals.
#' @param kind `"linear"` or `"lowess"`, mzLib's `IrtCalibrationKind`. `NULL` (the default) sends
#'   nothing, so mzLib's default, LOWESS, applies.
#' @param minimum_anchors A group with fewer anchors than this is refused, not fitted (anchors).
#'   `NULL` (the default) leaves mzLib's default, 20 anchors at the pin. mzLib requires at least 3.
#' @param bandwidth LOWESS only: the **fraction of anchors** in each local fit, in (0, 1]. Smaller
#'   follows the gradient, and its noise, more closely. `NULL` (the default) leaves mzLib's default,
#'   0.3.
#' @param knots LOWESS only: how many evenly spaced times between the first and last anchor the
#'   curve is fitted at (points); between them it is a straight line. `NULL` (the default) leaves
#'   mzLib's default, 100 points. At least 3.
#' @param robustness_iterations Rounds of Tukey bisquare reweighting after the first fit
#'   (iterations). `NULL` (the default) leaves mzLib's default, 3 iterations. 0 or more.
#' @param on_error `"fail"` (the default) stops at the first group mzLib refuses, with an
#'   `mzlib_usage_error` naming it; `"skip"` records the refusal in that group's row of `groups`
#'   and fits the rest. An option mzLib refuses fails the call either way.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_rt_calibration`.
#'
#'   `records` is a data.frame with one row per anchor, then one per query, in the order sent:
#'   `group`, as sent; `role`, `"anchor"`, `"to_irt"` (you gave minutes) or `"to_minutes"` (you gave
#'   iRT); `rt_minutes` (min), as given, or the converted time for a `"to_minutes"` query, `NA` when
#'   the group was not fitted and you gave none; `irt` (iRT), as given, or the converted value for a
#'   `"to_irt"` query, `NA` likewise; `fitted_irt` (iRT), anchors only, the curve's iRT at the
#'   anchor's time, so `irt - fitted_irt` is its residual, `NA` on a query or when the group was not
#'   fitted; and `is_extrapolated`, whether `rt_minutes` lies outside the group's anchors, where the
#'   map is a straight-line continuation rather than a fit, `NA` when the group was not fitted.
#'
#'   `groups` is a data.frame, one row per group in order of first appearance: `group`; `kind`, the
#'   shape mzLib fitted (`NA`: not fitted); `anchor_count`, the anchors it was fitted from (anchors;
#'   the anchors sent, if refused); `query_count`, the times converted (queries); `residual_sd`, how
#'   far the anchors that fit sit from the curve, in **iRT units**, mzLib's robust residual SD with
#'   wrong anchors set aside (divide by the curve's slope, iRT per minute, to read it in minutes;
#'   `NA`: not fitted); `first_anchor_minutes` and `last_anchor_minutes`, the span the map is fitted
#'   over, in minutes (`NA`: not fitted); and `error_kind`, `error_type` and `error_message`, why
#'   mzLib refused the group under `on_error = "skip"` (`"usage"`, `"ArgumentException"` and mzLib's
#'   message), `NA` when it was fitted. `groups[groups$group == name, ]` is one group;
#'   `groups[!is.na(groups$error_kind), ]` the groups mzLib refused.
#'
#'   `kind` is the kind applied, yours or mzLib's default; `minimum_anchors` (anchors), `bandwidth`
#'   (fraction of anchors), `knots` (points) and `robustness_iterations` (iterations) the options
#'   applied, yours or mzLib's; `on_error` as requested. `group_count` counts the groups sent (groups),
#'   `fitted_count` those mzLib fitted (groups), `failed_count` those it refused (groups; non-zero only
#'   under `"skip"`), `anchor_count` the anchors sent (anchors), `query_count` the queries sent
#'   (queries), and `row_count` the rows of `records` (rows): `anchor_count + query_count`.
#'   `column_names` lists the columns. `caveats` says what a calibration can and cannot tell you,
#'   plus any that apply to this call (queries outside the anchors, groups not fitted).
#'
#' @section The traps:
#'
#' **iRT is the scale of one library.** A calibration maps one run onto the scale of the library the
#' anchors' iRT values came from, so iRT from two libraries is not comparable. **`residual_sd` is
#' robust**: it understates the spread when many anchors are wrong, and a LOWESS fit to few anchors
#' can pass through nearly all of them, so a `residual_sd` near 0 there means an overfit, not a
#' perfect run. With few anchors, prefer `kind = "linear"`. **Outside the first and last anchor the
#' map is a straight line**, not a fit: read `is_extrapolated` before trusting a converted time
#' there.
#'
#' @seealso [library_read()] for a library's iRT values, [readers_read_records()] for a search's
#'   observed retention times
#' @spec library.calibrate-rt
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # Three 24-minute DIA runs of PBMCs, from a DIA-NN report in mzLib's test data. One precursor,
#' # AAAAAAAAGAFAGR at charge 2, is held out of the anchors, and each run predicts where it should
#' # elute from its library iRT:
#' pbmc <- read.delim(system.file("extdata", "calibration", "diann_pbmc_rt.tsv", package = "mzLibR"))
#' held <- pbmc$precursor_id == "AAAAAAAAGAFAGR2"
#' anchors <- data.frame(group = pbmc$run[!held], rt_minutes = pbmc$retention_time[!held],
#'                       irt = pbmc$indexed_retention_time[!held])
#' queries <- data.frame(group = pbmc$run[held], irt = 11.886096)
#' calibration <- library_calibrate_rt(anchors, queries, kind = "linear", minimum_anchors = 10)
#' calibration
#' predicted <- calibration$records[calibration$records$role == "to_minutes", ]
#' # Each prediction is within a quarter of a minute of where the precursor was seen:
#' data.frame(run = substr(predicted$group, 1, 8), predicted = round(predicted$rt_minutes, 2),
#'            seen = round(pbmc$retention_time[held], 2))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
library_calibrate_rt <- function(anchors, queries = NULL, kind = NULL, minimum_anchors = NULL, bandwidth = NULL,
                                 knots = NULL, robustness_iterations = NULL, on_error = "fail", timeout = 60) {
  request <- library_build_calibrate_request(
    if (missing(anchors)) NULL else anchors, queries, kind, minimum_anchors, bandwidth, knots,
    robustness_iterations, on_error
  )
  bridge_require_verb("library calibrate-rt", LIBRARY_CALIBRATE_SINCE)
  library_parse_calibration(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Print a retention-time calibration
#'
#' @param x A [library_calibrate_rt()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_rt_calibration <- function(x, ...) {
  cat("<mzlibr_rt_calibration> ", x$kind, ": ", format(x$fitted_count), " of ", format(x$group_count),
    " group(s) fitted from ", format(x$anchor_count), " anchors, ", format(x$query_count), " queries converted\n",
    sep = ""
  )
  for (i in seq_len(nrow(x$groups))) {
    g <- x$groups[i, ]
    cat("  ", g$group, ": ",
      if (is.na(g$error_kind)) {
        paste0(g$kind, ", ", format(g$anchor_count), " anchors over ", format(g$first_anchor_minutes, digits = 4),
          " to ", format(g$last_anchor_minutes, digits = 4), " min, residual SD ", format(g$residual_sd, digits = 3), " iRT")
      } else {
        paste0("not fitted - ", g$error_message)
      }, "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Print a spectral library read
#'
#' @param x A [library_read()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_library_records <- function(x, ...) {
  cat("<mzlibr_library_records> ", basename(x$path), ": ", format(x$record_count), " precursors (",
    format(x$target_count), " target, ", format(x$decoy_count), " decoy), .msl format ",
    format(x$format_version), if (isTRUE(x$is_index_only)) ", opened index-only" else "", "\n",
    sep = ""
  )
  if (!is.na(x$min_precursor_mz)) {
    cat("  precursor m/z ", format(x$min_precursor_mz, digits = 6), " to ", format(x$max_precursor_mz, digits = 6),
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$records)) {
    cat("  ", format(x$returned_count), " precursor(s) returned",
      if (isTRUE(x$fragments_included)) paste0(" as ", format(x$row_count), " fragment-ion rows") else "",
      if (isTRUE(x$offset > 0)) paste0(" from offset ", format(x$offset)) else "", "\n",
      sep = ""
    )
    units <- unique(x$records$retention_time_unit)
    units <- units[!is.na(units)]
    if (length(units) > 1L) {
      cat("  ! retention_time is in ", paste(units, collapse = " and "), ": read retention_time_unit on each row\n",
        sep = ""
      )
    }
  }
  if (length(x$absent_fields) > 0L) {
    cat("  no source for: ", paste(x$absent_fields, collapse = ", "), " (every value NA)\n", sep = "")
  }
  if (!is.null(x$output)) {
    cat("  written to ", x$output$path, " (", x$output$format, ", ", format(x$output$row_count), " rows)\n", sep = "")
  }
  invisible(x)
}
