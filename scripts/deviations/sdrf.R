# SDRF-Proteomics. See scripts/spec-facts.R for what each list means.

# SDRF's column names repeat, so a document keeps its header as `columns` - a character vector
# that may hold duplicates - and its cells as `rows`, rather than as a data.frame.
R_DEVIATIONS[["sdrf read"]] <- list(
  "field.column_names" = as_r("columns", "the header, duplicates kept, is `columns`")
)
R_DEVIATIONS[["sdrf pool"]] <- list(
  "param.stdin" = as_r("documents", "one stdin line per element of `documents`; its name is the label"),
  "field.column_names" = as_r("columns", "the header, duplicates kept, is `columns`")
)

PARENT_MAP[c(
  "sdrf_read", "sdrf_pool", "sdrf_index_of", "sdrf_indexes_of", "sdrf_value", "sdrf_all",
  "sdrf_records", "sdrf_has_repeated_columns", "sdrf_ragged_row_count", "sdrf_source_documents"
)] <- c(
  "sdrf.read", "sdrf.pool", "SdrfDocument.index_of", "SdrfDocument.indexes_of",
  "SdrfDocument.value", "SdrfDocument.all", "SdrfDocument.records",
  "SdrfDocument.has_repeated_columns", "SdrfDocument.ragged_row_count",
  "PooledSdrf.source_documents"
)

FIELD_CHECKS[["sdrf read"]] <- list("sdrf_read_PXD000070.json", function(d) mz$sdrf_parse_document(d))
FIELD_CHECKS[["sdrf pool"]] <- list("sdrf_pool_two.json", function(d) mz$sdrf_parse_pooled(d))

# ---------------------------------------------------------------- validate, lint, assess, samples, ages

R_TABLE[c("sdrf validate", "sdrf lint", "sdrf assess", "sdrf samples", "sdrf parse-age")] <- "records"

# One document and many are two functions (sdrf_validate / sdrf_validate_many), the cross-binding
# decision in the bridge's BULK.md: the bulk options exist only on the `_many` form.
sdrf_one_document <- function(bulk) {
  list(
    "param.paths-stdin" = not_here(paste0(bulk, "() is the many-documents form")),
    "param.threads" = not_here(paste0("only ", bulk, "() takes threads")),
    "param.on-error" = not_here(paste0("only ", bulk, "() takes on_error"))
  )
}
SDRF_BULK_PAGE <- list(
  "param.paths-stdin" = as_r("paths", "a character vector, sent one path per stdin line"),
  "param.threads" = as_r("threads", "the same option, spelled for R"),
  "param.on-error" = as_r("on_error", "the same option, spelled for R"),
  "field.files" = as_r("files", "a data.frame, one row per input; `error` becomes error_kind and error_message")
)
SDRF_SOURCE_INDEX <- as_r("source_index", "1-based in R, so it indexes `files` directly")

R_DEVIATIONS[["sdrf validate"]] <- c(sdrf_one_document("sdrf_validate_many"), list(
  "field.columns" = RECORDS,
  "field.row_index" = as_r("row_index", "1-based in R, like every position in this module"),
  "field.source_index" = SDRF_SOURCE_INDEX
))
R_DEVIATIONS[["sdrf validate (bulk)"]] <- SDRF_BULK_PAGE
R_DEVIATIONS[["sdrf lint"]] <- list(
  "param.stdin" = as_r("documents", "one stdin line per element of `documents`; its name is the label"),
  "field.columns" = RECORDS,
  "field.finding_index" = as_r("finding_index", "1-based in R"),
  "field.variant_rank" = as_r("variant_rank", "1-based in R: 1 is the majority spelling"),
  "field.documents" = as_r("documents", "a list column of character vectors")
)
R_DEVIATIONS[["sdrf assess"]] <- c(sdrf_one_document("sdrf_assess_many"), list(
  "field.columns" = RECORDS,
  "field.source_index" = SDRF_SOURCE_INDEX
))
R_DEVIATIONS[["sdrf assess (bulk)"]] <- c(SDRF_BULK_PAGE, list(
  "field.verdict_counts" = as_r("verdict_counts", "a named list: informative, partial, skeleton")
))
R_DEVIATIONS[["sdrf samples"]] <- c(sdrf_one_document("sdrf_samples_many"), list(
  "field.columns" = RECORDS,
  "field.position" = as_r("position", "1-based in R, like every position in this module"),
  "field.source_index" = SDRF_SOURCE_INDEX
))
R_DEVIATIONS[["sdrf samples (bulk)"]] <- SDRF_BULK_PAGE
R_DEVIATIONS[["sdrf parse-age"]] <- list(
  "param.stdin" = as_r("cells", "one stdin line per element of `cells`; NA is sent as an empty cell"),
  "field.columns" = RECORDS
)

PARENT_MAP[c(
  "sdrf_validate", "sdrf_validate_many", "sdrf_lint", "sdrf_assess", "sdrf_assess_many",
  "sdrf_samples", "sdrf_samples_many", "sdrf_parse_ages"
)] <- c(
  "sdrf.validate", "sdrf.validate_many", "sdrf.lint", "sdrf.assess", "sdrf.assess_many",
  "sdrf.samples", "sdrf.samples_many", "sdrf.parse_ages"
)
PARENT_OMISSIONS[c(
  "SdrfValidation.messages", "SdrfValidation.errors", "SdrfValidation.warnings",
  "SdrfValidationBatch.messages", "SdrfDrift.findings", "SdrfAssessmentBatch.paths_with",
  "SdrfSamples.ages", "SdrfSamples.conflicts"
)] <- c(
  "the findings are `records`, a data.frame; its rows are the messages",
  "`records[records$severity == \"Error\", ]`",
  "`records[records$severity == \"Warning\", ]`",
  "the findings are `records`, a data.frame, with source_index and source_path",
  "`split(records, records$finding_index)`",
  "`files$path[files$verdict %in% c(...)]`",
  "`records[records$column_name == \"characteristics[age]\", ]`",
  "`records[records$status == \"conflicting\", c(\"source_name\", \"column_name\")]`"
)

FIELD_CHECKS[["sdrf validate"]] <- list("sdrf_validate_cohort.json", function(d) mz$sdrf_parse_validation(d))
FIELD_CHECKS[["sdrf lint"]] <- list("sdrf_lint_cohort.json", function(d) mz$sdrf_parse_drift(d))
FIELD_CHECKS[["sdrf assess"]] <- list("sdrf_assess_cohort.json", function(d) mz$sdrf_parse_assessment(d))
FIELD_CHECKS[["sdrf samples"]] <- list("sdrf_samples_cohort.json", function(d) mz$sdrf_parse_samples(d))
FIELD_CHECKS[["sdrf parse-age"]] <- list("sdrf_parse_age.json", function(d) mz$sdrf_parse_ages_result(d))
