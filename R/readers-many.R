# The `_many` forms of the readers: many files in, one long table out, in ONE bridge call.
#
# See R/bulk.R for the rule and the shape. Each function here only names its verb and passes its
# options: the list travels on stdin, the bridge reads `threads` files at once, and the rows come
# back grouped by file in input order whatever `threads` was - the same list always gives the same
# table. None of them loops over files, and none may.

readers_read_many <- function(verb, paths, out, threads, on_error, timeout, fn, extra = character(0)) {
  request <- bulk_request(c("readers", verb), paths, threads, on_error, out, fn)
  bulk_parse_batch(
    bridge_invoke(c(request$args, extra), stdin = request$stdin, timeout = timeout),
    paste("readers", verb), "mzlibr_read_batch"
  )
}

#' Identify many files in one bridge call
#'
#' The many-files form of [readers_identify()]: one bridge process resolves every path, `threads`
#' at a time, and answers in input order.
#'
#' @param paths A character vector of paths. A Bruker `.d` directory is accepted.
#' @param threads Paths identified at once. `1`, the default, or `-1` for one per core. The answer
#'   does not depend on it.
#' @param on_error `"fail"`, the default, stops at the first path that cannot be identified.
#'   `"skip"` records why in that path's row and identifies the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_identify_batch`. `files` is a data.frame with one row per path, in input
#'   order: `path`, `file_type`, `extension`, `reader`, a `views` list column, `is_quantifiable`,
#'   and `error_kind`, `error_type` and `error_message`, which are `NA` for a path that was
#'   identified. `file_count`, `read_count` and `failed_count` count paths.
#'
#' @seealso [readers_identify()], [readers_formats()]
#' @spec readers.identify bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' found <- readers_identify_many(
#'   c("PXD078927_msgf_1_1_0.mzid", "no-such-run.mzML",
#'     "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"),
#'   on_error = "skip"
#' )
#' found$files[, c("file_type", "is_quantifiable", "error_message")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
readers_identify_many <- function(paths, threads = 1, on_error = "fail", timeout = 60) {
  request <- bulk_request(c("readers", "identify"), paths, threads, on_error, NULL, "readers_identify_many")
  batch <- bulk_parse_batch(
    bridge_invoke(request$args, stdin = request$stdin, timeout = timeout),
    "readers identify", "mzlibr_identify_batch"
  )
  views <- batch$files$views
  batch$files$is_quantifiable <- if (is.null(views)) {
    logical(nrow(batch$files))
  } else {
    vapply(views, function(v) READERS_QUANTIFIABLE %in% v, logical(1L))
  }
  batch
}

#' Read many result files into the uniform record view, in one bridge call
#'
#' The many-files form of [readers_read_results()].
#'
#' @param paths A character vector of paths to files offering the `"quantifiable"` view.
#' @param out Write the long table here as tab-separated text, one file at a time, and return only
#'   a summary.
#' @param threads Files read at once. `1`, the default, holds one whole file in memory at a time;
#'   `-1` uses one per core. The table does not depend on it.
#' @param on_error `"fail"`, the default, stops at the first file that cannot be read. `"skip"`
#'   records the failure in that file's row of `files` and reads the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`. `records` begins with `source_index` (1-based, into `paths`) and
#'   `source_path`, then the columns of [readers_read_results()]: `retention_time` in each file's
#'   own `retention_time_unit` - minutes for all four formats at mzLib 1.0.592 - `charge_state` a
#'   charge, `monoisotopic_mass` in Da. Convert times with [readers_retention_time_in_minutes()],
#'   which applies each file's unit to its rows. `files` has one row per input: `record_count` in
#'   records, `rows_not_read` in rows, `retention_time_unit`, `caveats`, `absent_fields` and the
#'   error columns. `file_count`, `read_count` and `failed_count` count files; `record_count` and
#'   `returned_count` count records over the files read; `row_count` counts the rows of `records`.
#'
#' @seealso [readers_read_results()]
#' @spec readers.read-results bulk
#' @examples
#' \dontrun{
#' # Your own files: no recording of a many-file read of this view exists to replay.
#' batch <- readers_read_results_many(c("run_1/AllPSMs.psmtsv", "run_2/AllPSMs.psmtsv"),
#'   threads = 2)
#' table(batch$records$source_path)
#' }
#' @export
readers_read_results_many <- function(paths, out = NULL, threads = 1, on_error = "fail",
                                      timeout = NULL) {
  readers_read_many("read-results", paths, out, threads, on_error, timeout, "readers_read_results_many")
}

#' Read many files, each into its format's own fields, in one bridge call
#'
#' The many-files form of [readers_read_records()]. Every file must have the same mzLib record type,
#' because one table has one column set: a list mixing record types is refused, naming the groups,
#' before any file is parsed. Split such a list by `file_type` from [readers_identify_many()]
#' first.
#'
#' @param paths A character vector of paths to any files mzLib recognises.
#' @param out Write the long table here as tab-separated text and return only a summary.
#' @param threads Files read at once. `1`, the default, or `-1` for one per core. The table does
#'   not depend on it.
#' @param on_error `"fail"`, the default, or `"skip"` to record an unreadable file in `files` and
#'   read the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`. `records` begins with `source_index` (1-based, into `paths`) and
#'   `source_path`, then the format's own columns. `files` has one row per input: `record_type`,
#'   `record_count` in records, `skipped_count` in items for mzIdentML, `caveats`,
#'   `excluded_fields`, `failed_fields` and the error columns. `file_count`, `read_count` and
#'   `failed_count` count files; `record_count` and `returned_count` count records; `row_count`
#'   counts the rows of `records`.
#'
#' @seealso [readers_read_records()]
#' @spec readers.read-records bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' batch <- readers_read_records_many(
#'   c("PXD078927_msgf_1_1_0.mzid", "PXD078927_msgf_1_1_0.mzid.gz"),
#'   threads = 2
#' )
#' batch
#' table(batch$records$source_index)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
readers_read_records_many <- function(paths, out = NULL, threads = 1, on_error = "fail",
                                      timeout = NULL) {
  readers_read_many("read-records", paths, out, threads, on_error, timeout, "readers_read_records_many")
}

#' Read MS1 features from many files, in one bridge call
#'
#' The many-files form of [readers_read_features()].
#'
#' @param paths A character vector of paths to `_ms1.feature` or Dinosaur `.feature.tsv` files.
#' @param out Write the long table here as tab-separated text and return only a summary.
#' @param threads Files read at once. `1`, the default, or `-1` for one per core. The table does
#'   not depend on it.
#' @param on_error `"fail"`, the default, or `"skip"` to record an unreadable file in `files` and
#'   read the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`. `records` begins with `source_index` (1-based, into `paths`) and
#'   `source_path`, then the columns of [readers_read_features()]: `mz` in m/z, `charge`,
#'   `retention_time_start` and `retention_time_end` in each file's own `retention_time_unit`, and
#'   `intensity` in the instrument's intensity units. `files` has one row per input, with its
#'   `retention_time_unit` - `"unknown"` for `_ms1.feature`. `file_count`, `read_count` and
#'   `failed_count` count files; `record_count` and `returned_count` count features; `row_count`
#'   counts the rows of `records`.
#'
#' @seealso [readers_read_features()]
#' @spec readers.read-features bulk
#' @examples
#' \dontrun{
#' # Your own files: no recording of a many-file read of this view exists to replay.
#' batch <- readers_read_features_many(c("a.feature.tsv", "b.feature.tsv"), threads = 2)
#' batch$files[, c("path", "retention_time_unit")]
#' }
#' @export
readers_read_features_many <- function(paths, out = NULL, threads = 1, on_error = "fail",
                                       timeout = NULL) {
  readers_read_many("read-features", paths, out, threads, on_error, timeout, "readers_read_features_many")
}

#' Read identifications from many files, in one bridge call
#'
#' The many-files form of [readers_read_matches()].
#'
#' @param paths A character vector of paths to MsPathFinderT, Casanovo or mzIdentML files.
#' @param scores Make the table long by score - one row per match and engine score - as
#'   [readers_read_matches()] does.
#' @param out Write the long table here as tab-separated text and return only a summary.
#' @param threads Files read at once. `1`, the default, or `-1` for one per core. The table does
#'   not depend on it.
#' @param on_error `"fail"`, the default, or `"skip"` to record an unreadable file in `files` and
#'   read the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`. `records` begins with `source_index` (1-based, into `paths`) and
#'   `source_path`, then the columns of [readers_read_matches()]: `q_value` a fraction from 0 to 1,
#'   and with `scores = TRUE` `score_value` in the engine's own units, named by `score_name`.
#'   `files` has one row per input, with `skipped_count` in items for mzIdentML. `file_count`,
#'   `read_count` and `failed_count` count files; `record_count` and `returned_count` count
#'   matches; `row_count` counts rows, more than matches with `scores = TRUE`.
#'
#' @seealso [readers_read_matches()]
#' @spec readers.read-matches bulk
#' @examples
#' \dontrun{
#' # Your own files: no recording of a many-file read of this view exists to replay.
#' batch <- readers_read_matches_many(c("a.mzid", "b.mzid"), threads = 2)
#' batch$records[batch$records$pass_threshold %in% TRUE, ]
#' }
#' @export
readers_read_matches_many <- function(paths, scores = FALSE, out = NULL, threads = 1,
                                      on_error = "fail", timeout = NULL) {
  readers_read_many(
    "read-matches", paths, out, threads, on_error, timeout, "readers_read_matches_many",
    readers_scores_args(scores)
  )
}

#' Read the scans of many spectra files, in one bridge call
#'
#' The many-files form of [readers_read_spectra()]: every file read by one bridge process,
#' `threads` at a time, into one scan table - so a study's runs can be compared by instrument and
#' acquisition time from `files` without a loop.
#'
#' @param paths A character vector of paths to spectra files. Bruker `.d` directories are accepted.
#' @param ms_order Keep only scans at this MS level in every file, or `NULL` for every level.
#' @param peaks Include the `mz` and `intensity` arrays. `FALSE` by default.
#' @param out Write the long table here as tab-separated text and return only a summary.
#' @param threads Files read at once. `1`, the default, holds one whole file in memory at a time;
#'   `-1` uses one per core. The table does not depend on it.
#' @param on_error `"fail"`, the default, or `"skip"` to record an unreadable file in `files` and
#'   read the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`. `records` begins with `source_index` (1-based, into `paths`) and
#'   `source_path`, then the columns of [readers_read_spectra()]: `retention_time` in minutes;
#'   `injection_time` in ms; `total_ion_current` and `selected_ion_intensity` in the instrument's
#'   intensity units; `peak_count` in peaks; `compensation_voltage` in volts;
#'   `selected_ion_charge_state_guess` a charge; `scan_window_lower_mz`, `scan_window_upper_mz`,
#'   `isolation_mz`, `isolation_width`, `selected_ion_mz` and `selected_ion_monoisotopic_guess_mz`
#'   in m/z; and with `peaks = TRUE` the list columns `mz` (m/z) and `intensity`.
#'   `files` has one row per input: `scan_count` in scans, `record_count` in scans after the
#'   `ms_order` filter, the run's `instrument_model`, `instrument_model_accession`,
#'   `instrument_serial_number`, `acquisition_start_time` and `acquisition_start_time_is_utc`
#'   (flattened from each file's `source`; `NA` where the file does not record them), `caveats`
#'   and the error columns. `file_count`, `read_count` and `failed_count` count files;
#'   `record_count` and `returned_count` count scans over the files read; `row_count` counts the
#'   rows of `records`.
#'
#' @seealso [readers_read_spectra()]
#' @spec readers.read-spectra bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' runs <- readers_read_spectra_many(
#'   c("sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"),
#'   threads = 2, on_error = "skip"
#' )
#' runs
#' runs$files[, c("file_type", "scan_count", "instrument_model", "error_kind")]
#' table(runs$records$source_index)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
readers_read_spectra_many <- function(paths, ms_order = NULL, peaks = FALSE, out = NULL, threads = 1,
                                      on_error = "fail", timeout = NULL) {
  readers_read_many(
    "read-spectra", paths, out, threads, on_error, timeout, "readers_read_spectra_many",
    readers_spectra_args(ms_order, peaks)
  )
}
