# SDRF analysis: validate, lint, assess, samples and ages (R/sdrf-analysis.R).
#
# Offline except for the live tests at the end. Ported from pyMzLib's test_sdrf_analysis.py; the
# fixtures are pyMzLib's recordings of the real bridge, byte for byte.

sdrf_analysis_fixture <- function(name) {
  path <- fixture_path(name)
  text <- rawToChar(readBin(path, "raw", file.info(path)$size))
  Encoding(text) <- "UTF-8"
  data <- mz$json_parse(text)
  if (is.list(data) && !is.null(data$ok) && "data" %in% names(data)) data$data else data
}

# ---------------------------------------------------------------- validate

test_that("a valid document has warnings and no errors, located by line", {
  result <- mz$sdrf_parse_validation(sdrf_analysis_fixture("sdrf_validate_cohort.json"))
  expect_true(inherits(result, "mzlibr_sdrf_validation"))
  expect_true(isTRUE(result$is_valid))
  expect_equal(result$error_count, 0)
  expect_equal(result$warning_count, 2)
  expect_identical(nrow(result$records), 2L)
  expect_identical(result$records$rule, c("ReservedWordCase", "ReservedWordCase"))
  expect_equal(result$records$line_number, c(10, 11))
  expect_identical(unique(result$records$column_name), "characteristics[age]")
})

test_that("row_index is 1-based, one less than the file line", {
  # The wire's row_index is 0-based into the data rows; the header is line 1.
  result <- mz$sdrf_parse_validation(sdrf_analysis_fixture("sdrf_validate_cohort.json"))
  expect_equal(result$records$row_index, result$records$line_number - 1)
})

test_that("a skeleton is invalid, and a document-wide finding has no line", {
  result <- mz$sdrf_parse_validation(sdrf_analysis_fixture("sdrf_validate_skeleton.json"))
  expect_false(isTRUE(result$is_valid))
  errors <- result$records[result$records$severity == "Error", ]
  expect_identical(errors$rule[[1L]], "RequiredColumn")
  expect_true(is.na(errors$line_number[[1L]]))
  expect_true(is.na(errors$row_index[[1L]]))
  expect_equal(result$error_count, sum(result$records$severity == "Error"))
})

test_that("validate_many keeps one row per input, and an unread one says why", {
  batch <- mz$sdrf_parse_batch(
    sdrf_analysis_fixture("sdrf_validate_bulk.json"), "sdrf validate",
    "mzlibr_sdrf_validation_batch", "row_index"
  )
  expect_true(inherits(batch, "mzlibr_batch"))
  expect_equal(c(batch$file_count, batch$read_count, batch$valid_count), c(3, 2, 1))
  expect_identical(nrow(batch$files), 3L)
  expect_identical(batch$files$is_valid, c(FALSE, TRUE, NA))
  expect_identical(batch$files$error_kind, c(NA, NA, "usage"))
  expect_true(grepl("missing.sdrf.tsv", batch$files$error_message[[3L]], fixed = TRUE))
  # source_index is 1-based, so it indexes files directly.
  expect_true(all(batch$records$source_index %in% c(1, 2)))
  expect_identical(
    unique(basename(batch$files$path[batch$records$source_index])),
    unique(basename(batch$records$source_path))
  )
  expect_equal(batch$record_count, nrow(batch$records))
})

# ---------------------------------------------------------------- lint

test_that("lint is one row per finding x variant, majority first", {
  drift <- mz$sdrf_parse_drift(sdrf_analysis_fixture("sdrf_lint_cohort.json"))
  expect_true(inherits(drift, "mzlibr_sdrf_drift"))
  expect_equal(drift$finding_count, 4)
  expect_identical(drift$labels, c("cohort", "partner"))
  expect_equal(sort(unique(drift$records$finding_index)), c(1, 2, 3, 4))
  expect_equal(min(drift$records$variant_rank), 1)
  first <- drift$records[drift$records$finding_index == 1, ]
  expect_identical(first$kind[[1L]], "AccessionNameConflict")
  expect_identical(first$value, c("Exploris 480", "Orbitrap Exploris 480"))
})

test_that("lint's documents column is a list column of labels", {
  drift <- mz$sdrf_parse_drift(sdrf_analysis_fixture("sdrf_lint_cohort.json"))
  expect_true(is.list(drift$records$documents))
  expect_identical(length(drift$records$documents), nrow(drift$records))
  expect_true(all(unlist(drift$records$documents) %in% c("cohort", "partner")))
  # Filtering the frame keeps the list column in step.
  some <- drift$records[drift$records$variant_rank == 1, ]
  expect_true(is.list(some$documents))
})

test_that("lint's documents are labelled by name, as for sdrf_pool", {
  request <- mz$sdrf_build_pool_request(
    c(cohort = "a.sdrf.tsv", partner = "b.sdrf.tsv"), NULL, NULL, 0
  )
  expect_identical(request$stdin, "a.sdrf.tsv\tcohort\nb.sdrf.tsv\tpartner")
  expect_error(sdrf_lint(character(0)), class = "mzlib_usage_error")
  expect_error(sdrf_lint(c(cohort = "a.sdrf.tsv", "b.sdrf.tsv")), class = "mzlib_usage_error",
    contains = "label")
})

# ---------------------------------------------------------------- assess

test_that("assess gives the verdict, the three checks and the evidence", {
  a <- mz$sdrf_parse_assessment(sdrf_analysis_fixture("sdrf_assess_cohort.json"))
  expect_true(inherits(a, "mzlibr_sdrf_assessment"))
  expect_identical(a$verdict, "Informative")
  expect_true(isTRUE(a$factor_value_varies))
  expect_true(isTRUE(a$sample_is_described))
  expect_true(isTRUE(a$biological_replicate_varies))
  expect_identical(a$records$column_name[[1L]], "factor value[disease]")
  expect_equal(a$records$distinct_values[[1L]], 2)
  expect_true(all(a$records$fill_rate >= 0 & a$records$fill_rate <= 1))
})

test_that("assess_many counts the verdicts, in input order", {
  batch <- mz$sdrf_parse_batch(
    sdrf_analysis_fixture("sdrf_assess_bulk.json"), "sdrf assess", "mzlibr_sdrf_assessment_batch"
  )
  expect_identical(batch$files$verdict, c("Informative", "Skeleton", "Partial"))
  expect_equal(unlist(batch$verdict_counts), c(informative = 1, partial = 1, skeleton = 1))
  expect_equal(batch$record_count, nrow(batch$records))
})

# ---------------------------------------------------------------- samples

test_that("samples merges each sample's rows, withholds a conflict and parses ages", {
  s <- mz$sdrf_parse_samples(sdrf_analysis_fixture("sdrf_samples_cohort.json"))
  expect_true(inherits(s, "mzlibr_sdrf_samples"))
  expect_equal(s$sample_count, 6)
  expect_identical(s$problems, character(0))
  conflicting <- s$records[s$records$status == "conflicting", ]
  expect_identical(conflicting$source_name, "S6")
  expect_identical(conflicting$column_name, "characteristics[disease]")
  expect_true(is.na(conflicting$value))
  expect_true(is.na(conflicting$position))

  ages <- s$records[s$records$column_name == "characteristics[age]", ]
  expect_identical(ages$source_name, c("S1", "S2", "S3", "S4", "S5", "S6"))
  expect_equal(ages$age_years, c(58, 62.5, 90, NA, NA, 40))
  expect_identical(ages$age_refusal, c(NA, NA, NA, "no_unit", "reserved_word", NA))
  # Age columns are NA on every other row.
  others <- s$records[s$records$column_name != "characteristics[age]", ]
  expect_true(all(is.na(others$age_years)))
})

test_that("position is 1-based", {
  s <- mz$sdrf_parse_samples(sdrf_analysis_fixture("sdrf_samples_cohort.json"))
  expect_equal(min(s$records$position, na.rm = TRUE), 1)
})

test_that("samples_many sums the samples and keeps each file's own count", {
  batch <- mz$sdrf_parse_batch(
    sdrf_analysis_fixture("sdrf_samples_bulk.json"), "sdrf samples", "mzlibr_sdrf_samples_batch",
    "position"
  )
  expect_equal(batch$sample_count, 8)
  expect_equal(batch$files$sample_count, c(6, 2))
  expect_true(is.list(batch$files$problems))
})

# ---------------------------------------------------------------- ages

test_that("parse_ages keeps input order and names each refusal", {
  ages <- mz$sdrf_parse_ages_result(sdrf_analysis_fixture("sdrf_parse_age.json"))
  expect_true(inherits(ages, "mzlibr_sdrf_ages"))
  expect_equal(c(ages$parsed_count, ages$cell_count), c(7, 11))
  expect_identical(nrow(ages$records), 11L)
  expect_equal(ages$records$years[1:4], c(58, 30.5, 62.5, 40))
  expect_equal(ages$records$min_years[3], 40)
  expect_equal(ages$records$max_years[3], 85)
  # A lower bound has no upper limit: NA, told apart from a refusal by precision.
  expect_true(is.na(ages$records$max_years[5]))
  expect_identical(ages$records$precision[5], "LowerBound")
  expect_identical(ages$records$refusal[8:11], c("no_unit", "reserved_word", "empty", "unreadable"))
})

test_that("parse_ages refuses what it cannot send", {
  expect_error(sdrf_parse_ages(character(0)), class = "mzlib_usage_error")
  expect_error(sdrf_parse_ages(58), class = "mzlib_usage_error")
  expect_error(sdrf_parse_ages("5\n8Y"), class = "mzlib_usage_error", contains = "line break")
})

# ---------------------------------------------------------------- the _many arguments

test_that("the _many forms validate their arguments before spawning anything", {
  expect_error(sdrf_validate_many(character(0)), class = "mzlib_usage_error")
  expect_error(sdrf_validate_many(1), class = "mzlib_usage_error")
  expect_error(sdrf_assess_many("a.sdrf.tsv", threads = 0), class = "mzlib_usage_error",
    contains = "threads")
  expect_error(sdrf_samples_many("a.sdrf.tsv", on_error = "ignore"), class = "mzlib_usage_error",
    contains = "on_error")
  expect_error(sdrf_validate(""), class = "mzlib_usage_error")
})

test_that("a _many call hands the whole list to one bridge call, on stdin", {
  request <- mz$bulk_request(c("sdrf", "validate"), c("a.sdrf.tsv", "b.sdrf.tsv"), -1, "skip",
    fn = "sdrf_validate_many"
  )
  expect_identical(
    request$args,
    c("sdrf", "validate", "--paths-stdin", "--threads", "-1", "--on-error", "skip")
  )
  expect_identical(request$stdin, "a.sdrf.tsv\nb.sdrf.tsv\n")
})

test_that("printing names the verdict and the findings", {
  out <- capture.output(print(mz$sdrf_parse_validation(sdrf_analysis_fixture("sdrf_validate_skeleton.json"))))
  expect_true(any(grepl("NOT valid", out, fixed = TRUE)))
  out <- capture.output(print(mz$sdrf_parse_drift(sdrf_analysis_fixture("sdrf_lint_cohort.json"))))
  expect_true(any(grepl("ValueCaseVariant", out, fixed = TRUE)))
  out <- capture.output(print(mz$sdrf_parse_ages_result(sdrf_analysis_fixture("sdrf_parse_age.json"))))
  expect_true(any(grepl("7 of 11", out, fixed = TRUE)))
})

# ---------------------------------------------------------------- live

# These verbs need a bridge from pyMzLib 0.2.0 or later; an older one is skipped, not failed.
sdrf_analysis_live <- function(verb) {
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  verbs <- tryCatch(
    {
      options(mzlibr.bridge = live_bridge)
      on.exit(options(mzlibr.bridge = NULL))
      # The raw envelope rather than mzlibr_bridge_version(), so the gate reads the bridge's own
      # verb list whatever this version of the projection exposes.
      unlist(mz$bridge_invoke("version", timeout = 60)$verbs)
    },
    error = function(e) NULL
  )
  skip_if(is.null(verbs) || !verb %in% verbs, paste0("the staged bridge has no '", verb, "'"))
  # The documents themselves: vendored beside the recordings when present, or pyMzLib's own
  # tests/fixtures named by MZLIBR_SDRF_FIXTURES.
  vendored <- tryCatch(dirname(fixture_path("sdrf_cohort.sdrf.tsv")), error = function(e) "")
  dir <- if (nzchar(vendored)) vendored else Sys.getenv("MZLIBR_SDRF_FIXTURES", "")
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "sdrf_cohort.sdrf.tsv")),
    "no SDRF documents to read (set MZLIBR_SDRF_FIXTURES to pyMzLib's tests/fixtures)")
  normalizePath(dir)
}

test_that("LIVE: validate and assess agree with the recordings", {
  dir <- sdrf_analysis_live("sdrf validate")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  cohort <- file.path(dir, "sdrf_cohort.sdrf.tsv")

  live <- sdrf_validate(cohort)
  recorded <- mz$sdrf_parse_validation(sdrf_analysis_fixture("sdrf_validate_cohort.json"))
  expect_identical(live$is_valid, recorded$is_valid)
  expect_identical(live$records$rule, recorded$records$rule)
  expect_equal(live$records$row_index, recorded$records$row_index)

  assessment <- sdrf_assess(cohort)
  expect_identical(assessment$verdict, "Informative")
})

test_that("LIVE: samples, lint and ages agree with the recordings", {
  dir <- sdrf_analysis_live("sdrf samples")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  samples <- sdrf_samples(file.path(dir, "sdrf_cohort.sdrf.tsv"))
  expect_equal(samples$sample_count, 6)

  drift <- sdrf_lint(c(
    cohort = file.path(dir, "sdrf_cohort.sdrf.tsv"),
    partner = file.path(dir, "sdrf_cohort_partner.sdrf.tsv")
  ))
  recorded <- mz$sdrf_parse_drift(sdrf_analysis_fixture("sdrf_lint_cohort.json"))
  expect_identical(drift$records$value, recorded$records$value)
  expect_identical(drift$records$documents, recorded$records$documents)

  ages <- sdrf_parse_ages(c("58Y", NA, ">=90Y", ""))
  # One row per cell, a trailing blank one included and nothing added after it.
  expect_identical(ages$records$refusal, c(NA, "empty", NA, "empty"))

  batch <- sdrf_validate_many(
    file.path(dir, c("sdrf_skeleton.sdrf.tsv", "sdrf_cohort.sdrf.tsv", "missing.sdrf.tsv")),
    on_error = "skip", threads = 2
  )
  expect_identical(batch$files$is_valid, c(FALSE, TRUE, NA))
})
