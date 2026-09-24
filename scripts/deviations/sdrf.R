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
