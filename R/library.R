# Spectral libraries: mzLib's binary `.msl` format, read with mzLib's own MslLibrary.
#
#   library_read()        one library, one row per precursor, or per fragment ion on request
#   library_read_many()   many libraries into one long table, in ONE bridge call
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
#' calibration against a run.
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
