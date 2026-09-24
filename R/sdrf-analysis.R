# SDRF-Proteomics analysis: is a document well-formed, do several agree, does it describe its
# samples, and what are its samples' ages.
#
# sdrf.R reads and pools documents. This file asks mzLib's four judgements about them - all
# public in mzLib since 1.0.592 and exposed by the bridge from pyMzLib 0.2.0:
#
#   sdrf_validate()   SdrfValidator: the specification's structural rules, for one document.
#   sdrf_lint()       SdrfDriftLint: concepts a SET of documents wrote inconsistently.
#   sdrf_assess()     SdrfSampleInformativeness: whether the sample half describes a design, or is
#                     a valid skeleton that says nothing.
#   sdrf_samples()    SdrfSampleBlock.BySourceName: each sample's characteristics and factor
#                     values, merged over its rows, with ages parsed.
#   sdrf_parse_ages() SdrfAge.TryParse: age cells into years, refusing anything that needs a guess.
#
# The three per-document verbs also have a `_many` form that reads a whole corpus in ONE bridge
# process (R/bulk.R); nothing here loops over files.
#
# Ported from pyMzLib's sdrf.py. What changes is only the projection into R: each result's table
# is a data.frame, `records`, rather than a map of columns; NA rather than None; and every index
# or position is 1-based, as elsewhere in this module - `row_index`, `finding_index`,
# `variant_rank` and `position` are one more than on the wire, and a batch's `source_index` indexes
# its `files` data.frame directly.

# The three verdicts sdrf_assess() can return, best first, as mzLib spells them.
SDRF_VERDICTS <- c("Informative", "Partial", "Skeleton")

# ---------------------------------------------------------------- parsing

# The wire's 0-based index columns, made 1-based in place. NA stays NA.
sdrf_one_based <- function(table, columns) {
  if (is.null(table)) {
    return(table)
  }
  for (column in intersect(columns, names(table))) {
    table[[column]] <- as.numeric(table[[column]]) + 1
  }
  table
}

# The findings / evidence / samples table of a one-document result. An empty result still has its
# columns, so `nrow()` and `$` work without a NULL check.
sdrf_analysis_table <- function(data, one_based = character(0)) {
  table <- readers_parse_records_table(data)
  if (is.null(table)) {
    names <- readers_parse_strings(data, "column_names")
    table <- as.data.frame(
      structure(rep(list(logical(0)), length(names)), names = names),
      stringsAsFactors = FALSE, optional = TRUE
    )
  }
  sdrf_one_based(table, one_based)
}

# Every scalar and string list of a payload, under its wire name, plus `records`.
sdrf_analysis_fields <- function(data, table) {
  out <- list()
  for (key in names(data)) {
    if (key %in% c("columns")) next
    value <- data[[key]]
    out[[key]] <- if (wire_null(value)) {
      NA
    } else if (is.list(value) && is.null(names(value))) {
      wire_strings(value)
    } else {
      value
    }
  }
  if (!is.null(out$column_names)) out$column_names <- readers_parse_strings(data, "column_names")
  if (!is.null(out$caveats)) out$caveats <- readers_parse_strings(data, "caveats")
  out$records <- table
  out
}

sdrf_parse_validation <- function(data) {
  structure(
    sdrf_analysis_fields(data, sdrf_analysis_table(data, "row_index")),
    class = "mzlibr_sdrf_validation"
  )
}

sdrf_parse_drift <- function(data) {
  table <- sdrf_analysis_table(data, c("finding_index", "variant_rank"))
  # Each variant is used by one or more documents, so `documents` is a list column even when every
  # row happens to name one - a character column would change type with the data.
  if ("documents" %in% names(table)) {
    raw <- data[["columns"]][["documents"]]
    table$documents <- I(lapply(raw, wire_strings))
  }
  structure(sdrf_analysis_fields(data, table), class = "mzlibr_sdrf_drift")
}

sdrf_parse_assessment <- function(data) {
  structure(
    sdrf_analysis_fields(data, sdrf_analysis_table(data)),
    class = "mzlibr_sdrf_assessment"
  )
}

sdrf_parse_samples <- function(data) {
  fields <- sdrf_analysis_fields(data, sdrf_analysis_table(data, "position"))
  fields$problems <- readers_parse_strings(data, "problems")
  structure(fields, class = "mzlibr_sdrf_samples")
}

sdrf_parse_ages_result <- function(data) {
  structure(
    sdrf_analysis_fields(data, sdrf_analysis_table(data)),
    class = "mzlibr_sdrf_ages"
  )
}

# A `_many` result: the shared batch shape, with this module's 0-based columns made 1-based.
sdrf_parse_batch <- function(data, verb, class, one_based = character(0)) {
  batch <- bulk_parse_batch(data, verb, class)
  batch$records <- sdrf_one_based(batch$records, one_based)
  batch
}

# ---------------------------------------------------------------- arguments

sdrf_one_path <- function(verb, path) {
  c("sdrf", verb, "--path", readers_normalise_path(path))
}

# ---------------------------------------------------------------- validate

#' Check one SDRF file against the specification's structural rules
#'
#' Calls mzLib's `SdrfValidator.Validate`. Structure only: required and recommended columns,
#' column order, casing, malformed names, ragged rows, integer replicate and fraction columns,
#' reserved-word casing, and the one hard row rule - `source name` + `assay name` +
#' the `comment` column for the label must be unique. Controlled-vocabulary accessions are **not**
#' resolved.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_validation`: `path`, `is_valid`, and the counts - `error_count` and
#'   `warning_count` findings of each severity, `message_count` findings in all, `row_count` data
#'   rows in the document (the header is not a row) - plus `caveats`, `column_names`, and
#'   `records`, a data.frame with one row per finding: `severity` (`"Error"` or `"Warning"`),
#'   `rule` (a stable name such as `"RequiredColumn"`), `message`, `row_index` (the data row,
#'   1-based), `line_number` (the file line, 1-based, the header being line 1) and
#'   `column_name`. `row_index` and `line_number` are `NA` for a finding about the whole
#'   document; `column_name` is `NA` when the finding is not about one column.
#'
#' @section Warnings never make a document invalid:
#'
#' `is_valid` is `TRUE` whenever there is no `Error`. mzLib calibrated every severity against the
#' 1,236-file curated corpus, and a rule that fired on most curated files was judged wrong rather
#' than the files. And a valid document may still say nothing about its samples: that is what
#' [sdrf_assess()] asks.
#'
#' @seealso [sdrf_validate_many()], [sdrf_assess()], [sdrf_lint()]
#' @spec sdrf.validate
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' result <- sdrf_validate("sdrf_skeleton.sdrf.tsv")
#' result
#' result$records[result$records$severity == "Error", c("rule", "line_number", "column_name")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_validate <- function(path, timeout = 60) {
  # Built first, so a bad path is refused before the bridge is even looked for.
  args <- sdrf_one_path("validate", path)
  sdrf_parse_validation(bridge_invoke(args, timeout = timeout))
}

#' Validate many SDRF files in one bridge call
#'
#' One bridge process for the whole corpus instead of one per file, with the parallelism inside
#' the bridge where it belongs. The result does not depend on `threads`: rows are always in input
#' order.
#'
#' @param paths A character vector of `.sdrf.tsv` paths, in the order to report them. Each may be
#'   given once.
#' @param threads Documents validated at once. `1`, the default, holds one document in memory at a
#'   time; `-1` uses every core.
#' @param on_error `"fail"`, the default, raises on the first unreadable document in input order.
#'   `"skip"` records the failure in that document's row of `files` and carries on.
#' @param timeout Seconds to allow for the whole batch, or `NULL`, the default, to wait
#'   indefinitely.
#'
#' @return An `mzlibr_sdrf_validation_batch`: `file_count` documents given, `read_count`
#'   documents validated, `failed_count` documents that could not be read, `valid_count` documents
#'   read with no `Error`, `message_count` findings over every document read, `record_count` rows
#'   in `records`, and `caveats`.
#'
#'   `records` has the columns of [sdrf_validate()]'s, preceded by `source_index` - the document's
#'   1-based position in `paths`, so it indexes `files` - and `source_path`. `files` has one row per
#'   input: `path`, `is_valid`, `error_count`, `warning_count`, `message_count`, `row_count`, and
#'   `error_kind` / `error_message`, which are `NA` for a document that was read. A document that
#'   was not read has `NA` for every count.
#'
#' @seealso [sdrf_validate()]
#' @spec sdrf.validate bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' batch <- sdrf_validate_many(
#'   c("sdrf_skeleton.sdrf.tsv", "sdrf_cohort.sdrf.tsv", "missing.sdrf.tsv"),
#'   on_error = "skip"
#' )
#' batch
#' batch$files[, c("path", "is_valid", "error_count", "error_kind")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_validate_many <- function(paths, threads = 1, on_error = "fail", timeout = NULL) {
  request <- bulk_request(c("sdrf", "validate"), paths, threads, on_error, fn = "sdrf_validate_many")
  sdrf_parse_batch(
    bridge_invoke(request$args, stdin = request$stdin, timeout = timeout),
    "sdrf validate", "mzlibr_sdrf_validation_batch", "row_index"
  )
}

# ---------------------------------------------------------------- lint

#' Find the concepts a set of SDRF documents annotated inconsistently
#'
#' Validity is a property of one file; comparability is a property of the relationship between
#' files. Every document can pass [sdrf_validate()] and the pooled table still be unusable because
#' one writes `"Homo sapiens"` and another `"homo sapiens"`. This is mzLib's `SdrfDriftLint`.
#'
#' Columns unique per row by construction are never compared - `source name`, `assay name`,
#' `comment[data file]`, `comment[searched data file]`, `comment[file uri]` and
#' `comment[source document]` - so differing file names are not drift. Reserved words and empty
#' cells are skipped too, so documents that say nothing lint clean: use [sdrf_assess()] for that.
#'
#' @param documents A character vector of paths. **Name it** to choose each document's label, as
#'   for [sdrf_pool()]: `c(cohort = "a.sdrf.tsv", partner = "b.sdrf.tsv")`. Name every element or
#'   none; mzLib's default label depends on where the files sit.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_drift`: `document_count` documents linted, `paths`, `labels`,
#'   `finding_count` distinct findings (not rows), `caveats`, `column_names`, and `records`, a
#'   long data.frame with **one row per finding x variant**: `finding_index` (1-based, most
#'   impactful first), `kind`, `concept`, `column_name` (`NA` for a `ColumnNameVariant`, which is
#'   about a name rather than a column's values), `variant_rank` (1 for the majority spelling),
#'   `value` (this spelling, as written), `occurrences` (how many documents used it), and
#'   `documents`, a list column of the labels of those documents.
#'
#' @section Kinds of drift:
#'
#' `AccessionNameConflict` - one accession written with several names. `NameAccessionConflict` -
#' one name with several accessions, the serious one. `MixedTermAndFreeText` - a column mixing
#' controlled-vocabulary terms and free text. `ColumnNameVariant` - column names differing only by
#' case or spacing. `ValueCaseVariant` - values differing only by case. The majority spelling is
#' not advice: it is only the most common.
#'
#' @seealso [sdrf_pool()], [sdrf_validate()]
#' @spec sdrf.lint
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' drift <- sdrf_lint(c(cohort = "sdrf_cohort.sdrf.tsv", partner = "sdrf_cohort_partner.sdrf.tsv"))
#' drift
#' drift$records[, c("finding_index", "kind", "value", "occurrences")]
#' drift$records$documents[[1]]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_lint <- function(documents, timeout = 60) {
  request <- sdrf_build_pool_request(documents, NULL, NULL, 0)
  sdrf_parse_drift(bridge_invoke(c("sdrf", "lint"), stdin = request$stdin, timeout = timeout))
}

# ---------------------------------------------------------------- assess

#' Decide whether an SDRF file describes its samples, or is a valid skeleton
#'
#' A file generated from a list of data files - reserved words in every sample column, one
#' replicate number, no factor - passes [sdrf_validate()], because reserved words are the
#' specification's correct way to say nothing. For grouping results by biology it is the same as
#' having no SDRF at all. mzLib's `SdrfSampleInformativeness.Assess` is the gate for that.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_assessment`: `path`; `verdict` - `"Informative"`, `"Partial"` or
#'   `"Skeleton"`; the three checks `factor_value_varies`, `sample_is_described` and
#'   `biological_replicate_varies`; `row_count` data rows; `caveats`; `column_names`; and
#'   `records`, the evidence - one row per column a check read: `role`, `column_name`, `rows` (data
#'   rows), `filled` (rows with a real answer), `absent` (rows empty or a reserved word),
#'   `distinct_values` (values, compared ignoring case and surrounding space) and `fill_rate`, the
#'   fraction (0-1) `filled / rows`.
#'
#' @section The verdict is how many checks pass:
#'
#' All three is `Informative`, none is `Skeleton`, anything between is `Partial`.
#' `factor_value_varies`: some `factor value[...]` column holds two or more real answers.
#' `sample_is_described`: some `characteristics[...]` column other than organism and biological
#' replicate holds a real answer - organism is left out because a search can fill it without a
#' human. `biological_replicate_varies`: `characteristics[biological replicate]` holds two or more
#' real answers; `FALSE` also when the column is missing. `Partial` is often legitimate - read the
#' `caveats`.
#'
#' @seealso [sdrf_assess_many()], [sdrf_validate()], [sdrf_samples()]
#' @spec sdrf.assess
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' assessment <- sdrf_assess("sdrf_cohort.sdrf.tsv")
#' assessment
#' assessment$records[, c("role", "column_name", "distinct_values", "fill_rate")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_assess <- function(path, timeout = 60) {
  # Built first, so a bad path is refused before the bridge is even looked for.
  args <- sdrf_one_path("assess", path)
  sdrf_parse_assessment(bridge_invoke(args, timeout = timeout))
}

#' Assess many SDRF files in one bridge call
#'
#' The gate for a corpus: which of these documents describe an experimental design?
#'
#' @param paths A character vector of `.sdrf.tsv` paths, in the order to report them. Each may be
#'   given once.
#' @param threads Documents assessed at once. `1` by default; `-1` uses every core. The result is
#'   the same at any value.
#' @param on_error `"fail"`, the default, or `"skip"`, as for [sdrf_validate_many()].
#' @param timeout Seconds to allow for the whole batch, or `NULL`, the default, to wait
#'   indefinitely.
#'
#' @return An `mzlibr_sdrf_assessment_batch`: `file_count` documents given, `read_count`
#'   documents assessed, `failed_count` documents not read, `record_count` rows in `records`,
#'   `verdict_counts` (a list counting documents per verdict: `informative`, `partial`,
#'   `skeleton`), and `caveats`.
#'
#'   `records` has the evidence columns of [sdrf_assess()] - `rows`, `filled` and `absent` in rows,
#'   `distinct_values` in values, `fill_rate` a fraction (0-1) - preceded by `source_index`
#'   (1-based, indexing `files`) and `source_path`. `files` has one row per input: `path`,
#'   `verdict`, the three checks, `row_count`, and `error_kind` / `error_message`, `NA` for a
#'   document that was read; an unread document has `NA` for every other fact.
#'
#' @seealso [sdrf_assess()]
#' @spec sdrf.assess bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' batch <- sdrf_assess_many(
#'   c("sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv")
#' )
#' batch$files[, c("path", "verdict")]
#' unlist(batch$verdict_counts)
#' # The documents fit to group results by biology:
#' batch$files$path[batch$files$verdict != "Skeleton"]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_assess_many <- function(paths, threads = 1, on_error = "fail", timeout = NULL) {
  request <- bulk_request(c("sdrf", "assess"), paths, threads, on_error, fn = "sdrf_assess_many")
  sdrf_parse_batch(
    bridge_invoke(request$args, stdin = request$stdin, timeout = timeout),
    "sdrf assess", "mzlibr_sdrf_assessment_batch"
  )
}

# ---------------------------------------------------------------- samples

#' Each sample's characteristics and factor values, merged over its rows, ages parsed
#'
#' Calls mzLib's `SdrfSampleBlock.BySourceName`: the sample half of the document - `source name`,
#' every `characteristics[...]`, every `factor value[...]` - per sample, merged over the sample's
#' rows. A sample is a `source name`, matched ignoring case and surrounding space, so one sample
#' measured in three fractions is one sample. Where its rows **disagree** - one fraction says
#' `normal` and another `COVID-19` - mzLib names the column and withholds it rather than let the
#' first row win.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_samples`: `path`, `sample_count` samples (distinct source names),
#'   `row_count` data rows, `conflict_count` (sample x column pairs withheld as conflicting),
#'   `problems` (rows mzLib could not place; empty for a well-formed document), `caveats`,
#'   `column_names`, and `records`, a long data.frame with one row per sample x column x position,
#'   samples in document order: `source_name`; `sample_row_count`, the rows carrying the sample;
#'   `column_name`, as the header spells it; `column_kind` (`source_name`, `characteristic` or
#'   `factor_value`); `position` (1 for the first occurrence of a repeated column; `NA` on a
#'   conflicting row); `status` (`"agreed"` or `"conflicting"`); and `value`, the cell verbatim,
#'   `NA` when conflicting.
#'
#'   On the rows for the age characteristic six more columns hold the parsed age, with the meanings
#'   [sdrf_parse_ages()] gives: `age_years`, `age_min_years` and `age_max_years` in years,
#'   `age_precision`, `age_follows_specification` and `age_refusal`. On every other row they are
#'   `NA`.
#'
#' @section Key samples within a document:
#'
#' A source name is unique only within one document - two studies may both have a `"Sample 1"` -
#' so across documents key a sample on the document as well.
#'
#' @seealso [sdrf_samples_many()], [sdrf_parse_ages()], [sdrf_assess()]
#' @spec sdrf.samples
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' samples <- sdrf_samples("sdrf_cohort.sdrf.tsv")
#' samples
#' ages <- samples$records[samples$records$column_name == "characteristics[age]", ]
#' ages[, c("source_name", "value", "age_years", "age_precision", "age_refusal")]
#' samples$records[samples$records$status == "conflicting", c("source_name", "column_name")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_samples <- function(path, timeout = 60) {
  # Built first, so a bad path is refused before the bridge is even looked for.
  args <- sdrf_one_path("samples", path)
  sdrf_parse_samples(bridge_invoke(args, timeout = timeout))
}

#' The samples of many SDRF files, in one bridge call
#'
#' @param paths A character vector of `.sdrf.tsv` paths, in the order to report them. Each may be
#'   given once.
#' @param threads Documents read at once. `1` by default; `-1` uses every core. The result is the
#'   same at any value.
#' @param on_error `"fail"`, the default, or `"skip"`, as for [sdrf_validate_many()].
#' @param timeout Seconds to allow for the whole batch, or `NULL`, the default, to wait
#'   indefinitely.
#'
#' @return An `mzlibr_sdrf_samples_batch`: `file_count` documents given, `read_count` documents
#'   read, `failed_count` documents not read, `sample_count` samples over every document read,
#'   `record_count` rows in `records`, and `caveats`.
#'
#'   `records` has the columns of [sdrf_samples()] - `sample_row_count` in rows; `age_years`,
#'   `age_min_years` and `age_max_years` in years - preceded by `source_index` (1-based, indexing `files`) and `source_path`. **Key a sample on
#'   `source_index` and `source_name` together**: a source name is unique only within a document.
#'   `files` has one row per input: `path`, `sample_count`, `row_count`, `conflict_count`, a
#'   `problems` list column, and `error_kind` / `error_message`, `NA` for a document that was read.
#'
#' @seealso [sdrf_samples()]
#' @spec sdrf.samples bulk
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' batch <- sdrf_samples_many(c("sdrf_cohort.sdrf.tsv", "sdrf_cohort_partner.sdrf.tsv"))
#' batch
#' batch$files[, c("path", "sample_count", "conflict_count")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_samples_many <- function(paths, threads = 1, on_error = "fail", timeout = NULL) {
  request <- bulk_request(c("sdrf", "samples"), paths, threads, on_error, fn = "sdrf_samples_many")
  sdrf_parse_batch(
    bridge_invoke(request$args, stdin = request$stdin, timeout = timeout),
    "sdrf samples", "mzlibr_sdrf_samples_batch", "position"
  )
}

# ---------------------------------------------------------------- parse ages

#' Read SDRF age cells into years, refusing anything that would need a guess
#'
#' Calls mzLib's `SdrfAge.TryParse` on each cell. It reads the specification's grammar (`58Y`,
#' `30Y6M`, `16W`), ranges (`40Y-85Y`, `6-8 weeks`), bounds (`>=90Y`, `<1Y`) and unambiguous words
#' (`3 year`, `4 hour`). It **refuses** a bare number - `63` is 11% of real age cells, and 63 years
#' and 63 days are both plausible in one study - as well as reserved words and free text.
#'
#' One bridge call for any number of cells, so pass them all at once rather than looping.
#'
#' @param cells A character vector of cells, e.g.
#'   the age characteristic's values from [sdrf_value()]. `NA` is sent as an empty cell and comes back
#'   refused as `"empty"`, so the result stays aligned with the input.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_ages`: `cell_count` cells given, `parsed_count` cells read as an age,
#'   `caveats`, `column_names`, and `records`, one row per cell **in input order**: `cell`, as
#'   given; `years`, the single figure to place the sample on an age axis, in years - the age, a
#'   range's midpoint, or a bound; `min_years`, the youngest the cell allows, in years (0 for an
#'   upper bound); `max_years`, the oldest, in years; `precision` (`"Exact"`, `"Range"`,
#'   `"LowerBound"` or `"UpperBound"`); `follows_specification` (`TRUE` for the specification's own
#'   `nYnMnD` / `nW` grammar); and `refusal` - `NA` when the cell was read, else `"empty"`,
#'   `"reserved_word"`, `"no_unit"` or `"unreadable"`.
#'
#'   `years`, `min_years`, `max_years`, `precision` and `follows_specification` are `NA` for a
#'   refused cell. `max_years` is also `NA` for a lower bound such as `>=90Y`, which has no upper
#'   limit; `precision` tells the two apart. A month is 1/12 year, a week 7/365.25, a day 1/365.25.
#'
#' @seealso [sdrf_samples()], which parses the ages of a document's samples for you.
#' @spec sdrf.parse-age
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' ages <- sdrf_parse_ages(c("58Y", "30Y6M", "40Y-85Y", "40Y-40Y", ">=90Y", "<1Y", "6-8 weeks",
#'   "63", "not available", "", "about forty"))
#' ages
#' ages$records[, c("cell", "years", "min_years", "max_years", "precision", "refusal")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_parse_ages <- function(cells, timeout = 60) {
  if (!is.character(cells) && !(is.logical(cells) && all(is.na(cells)))) {
    stop(mzlib_usage_error(paste0(
      "cells must be a character vector of age cells; got ", class(cells)[1L], "."
    )))
  }
  if (length(cells) == 0L) {
    stop(mzlib_usage_error("sdrf_parse_ages() needs at least one cell."))
  }
  cells <- as.character(cells)
  cells[is.na(cells)] <- ""
  if (any(grepl("[\r\n]", cells))) {
    stop(mzlib_usage_error("A cell contains a line break, which separates cells on the wire."))
  }
  # No trailing newline here: system2() writes `input` with writeLines(), which ends the last line
  # itself - so a final blank cell survives the trip, and no phantom blank cell follows it.
  stdin <- paste(cells, collapse = "\n")
  sdrf_parse_ages_result(bridge_invoke(c("sdrf", "parse-age"), stdin = stdin, timeout = timeout))
}

# ---------------------------------------------------------------- print methods

sdrf_print_caveat_count <- function(x) {
  if (length(x$caveats) > 0L) {
    cat("  ", length(x$caveats), " caveat(s) - read x$caveats once\n", sep = "")
  }
}

#' Print an SDRF validation
#'
#' @param x An [sdrf_validate()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_validation <- function(x, ...) {
  cat("<mzlibr_sdrf_validation> ", basename(x$path), ": ",
    if (isTRUE(x$is_valid)) "valid" else "NOT valid", "\n",
    sep = ""
  )
  cat("  ", format(x$error_count), " error(s), ", format(x$warning_count), " warning(s) over ",
    format(x$row_count), " rows\n",
    sep = ""
  )
  records <- x$records
  if (!is.null(records) && nrow(records) > 0L) {
    rules <- table(paste0(records$severity, " ", records$rule))
    for (rule in names(rules)) {
      cat("  ! ", rule, " x", rules[[rule]], "\n", sep = "")
    }
  }
  sdrf_print_caveat_count(x)
  invisible(x)
}

#' Print SDRF drift findings
#'
#' @param x An [sdrf_lint()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_drift <- function(x, ...) {
  cat("<mzlibr_sdrf_drift> ", format(x$document_count), " documents: ",
    paste(x$labels, collapse = ", "), "\n",
    sep = ""
  )
  cat("  ", format(x$finding_count), " finding(s)\n", sep = "")
  records <- x$records
  if (!is.null(records) && nrow(records) > 0L) {
    for (index in unique(records$finding_index)) {
      rows <- records[records$finding_index == index, , drop = FALSE]
      cat("  ! ", rows$kind[[1L]], ": ", paste(rows$value, collapse = " | "), "\n", sep = "")
    }
  }
  sdrf_print_caveat_count(x)
  invisible(x)
}

#' Print an SDRF informativeness assessment
#'
#' @param x An [sdrf_assess()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_assessment <- function(x, ...) {
  cat("<mzlibr_sdrf_assessment> ", basename(x$path), ": ", x$verdict, "\n", sep = "")
  yes_no <- function(value) if (isTRUE(value)) "yes" else "no"
  cat("  factor value varies: ", yes_no(x$factor_value_varies),
    "; sample described: ", yes_no(x$sample_is_described),
    "; biological replicate varies: ", yes_no(x$biological_replicate_varies), "\n",
    sep = ""
  )
  sdrf_print_caveat_count(x)
  invisible(x)
}

#' Print an SDRF document's samples
#'
#' @param x An [sdrf_samples()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_samples <- function(x, ...) {
  cat("<mzlibr_sdrf_samples> ", basename(x$path), "\n", sep = "")
  cat("  ", format(x$sample_count), " samples over ", format(x$row_count), " rows",
    if (isTRUE(x$conflict_count > 0)) {
      paste0("; ", format(x$conflict_count), " column(s) withheld as conflicting")
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  for (problem in x$problems) {
    cat("  ! ", problem, "\n", sep = "")
  }
  sdrf_print_caveat_count(x)
  invisible(x)
}

#' Print parsed SDRF ages
#'
#' @param x An [sdrf_parse_ages()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_ages <- function(x, ...) {
  cat("<mzlibr_sdrf_ages> ", format(x$parsed_count), " of ", format(x$cell_count),
    " cells read as an age\n",
    sep = ""
  )
  records <- x$records
  if (!is.null(records) && "refusal" %in% names(records)) {
    refused <- table(records$refusal[!is.na(records$refusal)])
    for (reason in names(refused)) {
      cat("  refused (", reason, "): ", refused[[reason]], "\n", sep = "")
    }
  }
  invisible(x)
}
