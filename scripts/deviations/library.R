# Spectral libraries: library read. See scripts/spec-facts.R for what each list means.

R_TABLE["library read"] <- "records"
R_DEVIATIONS[["library read"]] <- list(
  "field.columns" = RECORDS,
  "field.error" = ONE_PATH_ERROR,
  "field.precursor_index" = as_r("precursor_index", "1-based in R, like every position")
)
# The many-libraries page, as every readers `_many` page (scripts/deviations/readers.R, which is
# sourced after this file, so its list is not reused here).
R_DEVIATIONS[["library read (bulk)"]] <- list(
  "param.paths-stdin" = as_r("paths", "the character vector of paths; it travels on stdin, one per line"),
  "field.columns" = RECORDS,
  "field.files" = as_r("files", paste(
    "a data.frame, one row per input: `error` is split into error_kind, error_type and",
    "error_message, and a list becomes a list column"
  ))
)

PARENT_MAP[c("library_read", "library_read_many")] <- c("library.read", "library.read_many")
PARENT_OMISSIONS[c("LibraryFileReport.ok", "LibraryBatch.failed_files")] <- c(
  "a column test on `files`: `is.na(error_kind)`",
  "a row filter on `files`: `batch$files[!is.na(batch$files$error_kind), ]`"
)

FIELD_CHECKS[["library read"]] <- list("library_read_synthetic.json", function(d) mz$library_parse_records(d))

# The spectral-library article: one precursor's fragment ions, and an index-only read.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "library read", fixture = "library_read_fragments.json")
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "library read", fixture = "library_read_index_only.json")

# Not yet ported: the pyMzLib 0.7.0 bridge's library calibrate-rt (mzLib #1387), a later mzLibR port.
PARENT_OMISSIONS[c("library.calibrate_rt", "RtCalibration.group", "RtCalibration.failed_groups",
                   "RtCalibrationGroup.ok")] <- c(
  "not ported yet: wire `library calibrate-rt` (mzLib #1387), new in the pyMzLib 0.7.0 bridge",
  "not ported yet: part of `library.calibrate_rt`'s result",
  "not ported yet: part of `library.calibrate_rt`'s result",
  "not ported yet: part of `library.calibrate_rt`'s result"
)
