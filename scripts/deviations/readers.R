# Reading files. See scripts/spec-facts.R for what each list means.

R_TABLE["readers formats"] <- "(the returned data.frame itself)"
R_TABLE[c(
  "readers read-results", "readers read-records", "readers read-features", "readers read-matches",
  "readers read-spectra", "readers read-protein-groups", "readers read-quantified-peptides",
  "readers read-occupancy"
)] <- "records"

# Every `_many` page: the list of paths is an argument, the long table is `records`, and each
# input's facts are a row of the `files` data.frame.
BULK_READ <- list(
  "param.paths-stdin" = as_r("paths", "the character vector of paths; it travels on stdin, one per line"),
  "field.columns" = RECORDS,
  "field.files" = as_r("files", paste(
    "a data.frame, one row per input: `error` is split into error_kind, error_type and",
    "error_message, a list becomes a list column, and a nested object (a spectra file's `source`)",
    "is flattened into its own columns"
  ))
)

R_DEVIATIONS[["readers formats"]] <- list(
  "field.format_count" = not_here("the result is the formats data.frame; this is nrow() of it"),
  "field.formats" = not_here("the result is this list itself, as a data.frame, one row per format")
)
R_DEVIATIONS[["readers identify"]] <- list(
  "field.error" = ONE_PATH_ERROR
)
for (verb in c(
  "readers read-results", "readers read-records", "readers read-features", "readers read-matches",
  "readers read-spectra", "readers read-protein-groups", "readers read-quantified-peptides",
  "readers read-occupancy"
)) {
  R_DEVIATIONS[[verb]] <- list("field.columns" = RECORDS, "field.error" = ONE_PATH_ERROR)
  R_DEVIATIONS[[paste0(verb, " (bulk)")]] <- BULK_READ
}
rm(verb)

PARENT_MAP[c(
  "readers_formats", "readers_identify", "readers_read_results", "readers_read_records",
  "readers_read_features", "readers_read_matches", "readers_read_spectra",
  "readers_read_protein_groups", "readers_read_quantified_peptides", "readers_read_occupancy",
  "readers_identify_many", "readers_read_results_many", "readers_read_records_many",
  "readers_read_features_many", "readers_read_matches_many", "readers_read_spectra_many",
  "readers_read_protein_groups_many", "readers_read_quantified_peptides_many",
  "readers_read_occupancy_many"
)] <- c(
  "readers.formats", "readers.identify", "readers.read_results", "readers.read_records",
  "readers.read_features", "readers.read_matches", "readers.read_spectra",
  "readers.read_protein_groups", "readers.read_quantified_peptides", "readers.read_occupancy",
  "readers.identify_many", "readers.read_results_many", "readers.read_records_many",
  "readers.read_features_many", "readers.read_matches_many", "readers.read_spectra_many",
  "readers.read_protein_groups_many", "readers.read_quantified_peptides_many",
  "readers.read_occupancy_many"
)
PARENT_ADDITIONS["readers_retention_time_in_minutes"] <- paste(
  "one function with a `column` argument where pyMzLib has one property per column",
  "(ResultRecords.retention_time_in_minutes, FeatureRecords.retention_time_start_in_minutes and",
  "_end_in_minutes, ReadBatch.in_minutes), because R dispatches on the object's class and Python",
  "does not"
)
PARENT_OMISSIONS[c(
  "Format.is_quantifiable", "FileInfo.is_quantifiable", "ResultRecords.retention_time_in_minutes",
  "FeatureRecords.retention_time_start_in_minutes", "FeatureRecords.retention_time_end_in_minutes",
  "ScanRecords.total_ion_current", "ReadBatch.in_minutes", "ReadBatch.failed_files",
  "IdentifyBatch.failed_files", "FileReport.ok", "SpectraSource.acquired_at"
)] <- c(
  "a column of readers_formats(), `is_quantifiable`",
  "an element of readers_identify(), `is_quantifiable`",
  "readers_retention_time_in_minutes()",
  "readers_retention_time_in_minutes(column = \"retention_time_start\")",
  "readers_retention_time_in_minutes(column = \"retention_time_end\")",
  "a column of `records`, `total_ion_current`",
  "readers_retention_time_in_minutes() of the batch, which converts each file's rows by its unit",
  "a row filter on `files`: `batch$files[!is.na(batch$files$error_kind), ]`",
  "a row filter on `files`: `batch$files[!is.na(batch$files$error_kind), ]`",
  "a column test on `files`: `is.na(error_kind)`",
  "`as.POSIXct()` of `source$acquisition_start_time`, with `tz = \"UTC\"` only when `acquisition_start_time_is_utc`"
)

# readers read-matches: ?readers_read_matches also shows Casanovo, whose missing decoy label is
# the point of that part of its example, and a --scores read; the spec's own example is mzIdentML.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "readers read-matches", fixture = "readers_matches_casanovo.json")
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "readers read-matches", fixture = "readers_matches_mzid_scores.json")
# readers identify: the many-files example; the spec lists only the one-path recording.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "readers identify", fixture = "readers_many_identify.json")
# The verb guard (bridge_require_verb()) asks the stand-in which verbs it has before the three
# quantification readers run; bridge_version.json is already staged for `version` by its spec.

FIELD_CHECKS[["readers formats"]] <- list("readers_formats.json", function(d) mz$readers_parse_formats(d))
FIELD_CHECKS[["readers identify"]] <- list("readers_identify_mzid.json", function(d) mz$readers_parse_file_info(d))
FIELD_CHECKS[["readers read-results"]] <- list("readers_results_fragger.json", function(d) mz$readers_parse_records(d))
FIELD_CHECKS[["readers read-records"]] <- list("readers_records_mzid_gz.json", function(d) mz$readers_parse_native_records(d))
FIELD_CHECKS[["readers read-features"]] <- list("readers_features_topfd.json", function(d) mz$readers_parse_feature_records(d))
FIELD_CHECKS[["readers read-matches"]] <- list("readers_matches_mzid.json", function(d) mz$readers_parse_match_records(d))
FIELD_CHECKS[["readers read-spectra"]] <- list("readers_spectra_mzml.json", function(d) mz$readers_parse_scan_records(d))
FIELD_CHECKS[["readers read-protein-groups"]] <- list(
  "readers_protein_groups.json",
  function(d) mz$readers_parse_quant_records(d, "mzlibr_protein_group_records")
)
FIELD_CHECKS[["readers read-quantified-peptides"]] <- list(
  "readers_quantified_peptides.json",
  function(d) mz$readers_parse_quant_records(d, "mzlibr_quantified_peptide_records")
)
FIELD_CHECKS[["readers read-occupancy"]] <- list(
  "readers_occupancy.json",
  function(d) mz$readers_parse_quant_records(d, "mzlibr_occupancy_records")
)
