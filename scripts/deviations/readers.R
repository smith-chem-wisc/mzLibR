# Reading files. See scripts/spec-facts.R for what each list means.

R_TABLE["readers formats"] <- "(the returned data.frame itself)"
R_TABLE[c(
  "readers read-results", "readers read-records", "readers read-features", "readers read-matches",
  "readers read-spectra"
)] <- "records"

R_DEVIATIONS[["readers formats"]] <- list(
  "field.format_count" = not_here("the result is the formats data.frame; this is nrow() of it"),
  "field.formats" = not_here("the result is this list itself, as a data.frame, one row per format")
)
R_DEVIATIONS[["readers identify"]] <- list(
  "field.error" = ONE_PATH_ERROR
)
R_DEVIATIONS[["readers read-results"]] <- list(
  "field.columns" = RECORDS,
  "field.reader" = pending(),
  "field.absent_fields" = pending(),
  "field.failed_fields" = pending(),
  "field.excluded_fields" = pending(),
  "field.error" = ONE_PATH_ERROR
)
R_DEVIATIONS[["readers read-records"]] <- list(
  "field.columns" = RECORDS,
  "field.retention_time_unit" = pending(),
  "field.caveats" = pending(),
  "field.absent_fields" = pending(),
  "field.skipped_count" = pending(),
  "field.skipped" = pending(),
  "field.rows_not_read" = pending(),
  "field.error" = ONE_PATH_ERROR
)
R_DEVIATIONS[["readers read-features"]] <- list(
  "field.columns" = RECORDS,
  "field.reader" = pending(),
  "field.rows_not_read" = pending(),
  "field.absent_fields" = pending(),
  "field.failed_fields" = pending(),
  "field.excluded_fields" = pending(),
  "field.error" = ONE_PATH_ERROR
)
R_DEVIATIONS[["readers read-matches"]] <- list(
  "param.scores" = pending(),
  "field.columns" = RECORDS,
  "field.reader" = pending(),
  "field.skipped_count" = pending(),
  "field.skipped" = pending(),
  "field.rows_not_read" = pending(),
  "field.retention_time_unit" = pending(),
  "field.absent_fields" = pending(),
  "field.failed_fields" = pending(),
  "field.excluded_fields" = pending(),
  "field.scores_included" = pending(),
  "field.row_count" = pending(),
  "field.q_value" = pending(),
  "field.rank" = pending(),
  "field.pass_threshold" = pending(),
  "field.match_index" = pending(),
  "field.score_name" = pending(),
  "field.score_value" = pending(),
  "field.error" = ONE_PATH_ERROR
)
R_DEVIATIONS[["readers read-spectra"]] <- list(
  "field.columns" = RECORDS,
  "field.source" = pending(),
  "field.rows_not_read" = pending(),
  "field.absent_fields" = pending(),
  "field.failed_fields" = pending(),
  "field.excluded_fields" = pending(),
  "field.error" = ONE_PATH_ERROR
)

PARENT_MAP[c(
  "readers_formats", "readers_identify", "readers_read_results", "readers_read_records",
  "readers_read_features", "readers_read_matches", "readers_read_spectra"
)] <- c(
  "readers.formats", "readers.identify", "readers.read_results", "readers.read_records",
  "readers.read_features", "readers.read_matches", "readers.read_spectra"
)
PARENT_ADDITIONS["readers_retention_time_in_minutes"] <- paste(
  "one function with a `column` argument where pyMzLib has one property per column",
  "(ResultRecords.retention_time_in_minutes, FeatureRecords.retention_time_start_in_minutes and",
  "_end_in_minutes), because R dispatches on the object's class and Python does not"
)
PARENT_OMISSIONS[c(
  "Format.is_quantifiable", "FileInfo.is_quantifiable", "ResultRecords.retention_time_in_minutes",
  "FeatureRecords.retention_time_start_in_minutes", "FeatureRecords.retention_time_end_in_minutes",
  "ScanRecords.total_ion_current"
)] <- c(
  "a column of readers_formats(), `is_quantifiable`",
  "an element of readers_identify(), `is_quantifiable`",
  "readers_retention_time_in_minutes()",
  "readers_retention_time_in_minutes(column = \"retention_time_start\")",
  "readers_retention_time_in_minutes(column = \"retention_time_end\")",
  "a column of `records`, `total_ion_current`"
)

# readers read-matches: ?readers_read_matches shows Casanovo, whose missing decoy label is the
# point of its example; the spec's own example is mzIdentML.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "readers read-matches", fixture = "readers_matches_casanovo.json")

FIELD_CHECKS[["readers formats"]] <- list("readers_formats.json", function(d) mz$readers_parse_formats(d))
FIELD_CHECKS[["readers identify"]] <- list("readers_identify_mzid.json", function(d) mz$readers_parse_file_info(d))
FIELD_CHECKS[["readers read-results"]] <- list("readers_results_fragger.json", function(d) mz$readers_parse_records(d))
FIELD_CHECKS[["readers read-records"]] <- list("readers_records_toppic.json", function(d) mz$readers_parse_native_records(d))
FIELD_CHECKS[["readers read-features"]] <- list("readers_features_topfd.json", function(d) mz$readers_parse_feature_records(d))
FIELD_CHECKS[["readers read-matches"]] <- list("readers_matches_casanovo.json", function(d) mz$readers_parse_match_records(d))
FIELD_CHECKS[["readers read-spectra"]] <- list("readers_spectra_mzml.json", function(d) mz$readers_parse_scan_records(d))
