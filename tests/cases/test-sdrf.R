# SDRF-Proteomics experimental design.
#
# Offline except for the two live tests at the end. Ported from pyMzLib's test_sdrf.py. The
# fixtures are recorded from the real bridge by pyMzLib and shared verbatim, so they carry
# PXD059974's raggedness and PXD000070's eight repetitions of `comment[modification parameters]` -
# both of which are the point rather than an accident.

sdrf_fixture <- function(name) {
  mz$json_parse(paste(readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  ))
}

sdrf_document <- function() mz$sdrf_parse_document(sdrf_fixture("sdrf_read_PXD000070.json"))
sdrf_ragged <- function() mz$sdrf_parse_document(sdrf_fixture("sdrf_read_ragged.json"))
sdrf_pooled <- function() mz$sdrf_parse_pooled(sdrf_fixture("sdrf_pool_two.json"))

# ---------------------------------------------------------------- the shape

test_that("columns keeps repeated names, because SDRF names repeat", {
  doc <- sdrf_document()
  expect_true(inherits(doc, "mzlibr_sdrf"))
  expect_true(sdrf_has_repeated_columns(doc))
  expect_identical(length(sdrf_indexes_of(doc, "comment[modification parameters]")), 8L)
})

test_that("sdrf_all returns every occurrence where sdrf_value returns the first", {
  doc <- sdrf_document()
  every <- sdrf_all(doc, "comment[modification parameters]")[[1L]]
  first <- sdrf_value(doc, "comment[modification parameters]")[[1L]]
  expect_identical(length(every), 8L)
  expect_identical(every[[1L]], first)
  expect_true(startsWith(first, "NT=Carbamidomethyl"))
})

test_that("a CV cell crosses intact rather than being split on its semicolons", {
  # The defect this module exists to close: readers_read_records() joins an SdrfRow's cells
  # with ";", and the SDRF key=value grammar is itself semicolon-delimited.
  doc <- sdrf_document()
  expect_identical(
    sdrf_value(doc, "comment[modification parameters]")[[1L]],
    "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"
  )
})

test_that("an absent column is NA and a reserved word is itself", {
  doc <- sdrf_document()
  expect_identical(sdrf_value(doc, "characteristics[nonesuch]"), rep(NA_character_, length(doc$rows)))
  expect_identical(sdrf_index_of(doc, "characteristics[nonesuch]"), NA_integer_)
  expect_identical(sdrf_indexes_of(doc, "characteristics[nonesuch]"), integer(0))
  expect_identical(sdrf_value(doc, "characteristics[disease]")[[1L]], "not applicable")
})

test_that("positions are 1-based, as everywhere else in R", {
  doc <- sdrf_document()
  expect_identical(sdrf_index_of(doc, doc$columns[[1L]]), 1L)
})

test_that("a short row gives NA rather than shifting later columns", {
  doc <- sdrf_ragged()
  expect_identical(sdrf_ragged_row_count(doc), 17L)
  short <- which(lengths(doc$rows) < length(doc$columns))[[1L]]
  full <- which(lengths(doc$rows) == length(doc$columns))[[1L]]

  last <- doc$columns[[length(doc$columns)]]
  expect_false(is.na(sdrf_value(doc, last)[[full]]))
  expect_true(is.na(sdrf_value(doc, last)[[short]]))
  # A position both rows reach still lines up, which is what "no shifting" means.
  expect_false(is.na(sdrf_value(doc, "characteristics[organism]")[[short]]))
})

test_that("caveats name the raggedness", {
  expect_true(any(grepl("SHORT", sdrf_ragged()$caveats, fixed = TRUE)))
})

test_that("truncated says rows were left behind", {
  # Recorded with limit = 4 against a 24-row merge.
  pooled <- sdrf_pooled()
  expect_true(pooled$truncated)
  expect_identical(pooled$returned_count, 4)
  expect_identical(pooled$row_count, 24)
})

test_that("a complete read is not marked truncated", {
  doc <- sdrf_ragged()
  expect_false(doc$truncated)
  expect_identical(doc$returned_count, doc$row_count)
})

test_that("sdrf_records is offered but loses cells when names repeat", {
  doc <- sdrf_document()
  records <- sdrf_records(doc)
  expect_true(is.data.frame(records))
  expect_identical(nrow(records), length(doc$rows))
  expect_true(ncol(records) < length(doc$columns))
  # Names stay verbatim, brackets and all.
  expect_true("characteristics[disease]" %in% names(records))
  # The last occurrence wins, as in pyMzLib's records.
  expect_identical(
    records[["comment[modification parameters]"]][[1L]],
    tail(sdrf_all(doc, "comment[modification parameters]")[[1L]], 1L)
  )
})

test_that("sdrf_records gives NA where a row is short", {
  doc <- sdrf_ragged()
  records <- sdrf_records(doc)
  short <- which(lengths(doc$rows) < length(doc$columns))[[1L]]
  expect_true(is.na(records[[doc$columns[[length(doc$columns)]]]][[short]]))
})

test_that("a null list on the wire is empty rather than an error", {
  doc <- mz$sdrf_parse_document(list(path = "x", column_names = NA, rows = NA, caveats = NA))
  expect_identical(doc$columns, character(0))
  expect_identical(doc$rows, list())
  expect_identical(doc$caveats, character(0))
})

test_that("printing names the repetition, the truncation and the caveats", {
  output <- paste(capture.output(print(sdrf_pooled())), collapse = "\n")
  expect_true(grepl("mzlibr_pooled_sdrf", output, fixed = TRUE), info = output)
  expect_true(grepl("malaria, colon", output, fixed = TRUE), info = output)
  expect_true(grepl("truncated", output, fixed = TRUE), info = output)
  expect_true(grepl("ANALYSIS table", output, fixed = TRUE), info = output)
})

# ---------------------------------------------------------------- pool

test_that("a pool carries provenance and the labels it was given", {
  pooled <- sdrf_pooled()
  expect_true(inherits(pooled, "mzlibr_sdrf"))
  expect_identical(pooled$document_count, 2)
  expect_identical(pooled$labels, c("malaria", "colon"))
  expect_identical(pooled$path, "")
  expect_true(is.null(pooled$written))
  expect_true(mz$SDRF_SOURCE_DOCUMENT_COLUMN %in% pooled$columns)
  expect_true(all(sdrf_source_documents(pooled) %in% c("malaria", "colon")))
})

test_that("sdrf_source_documents refuses a document that was not pooled", {
  expect_error(sdrf_source_documents(sdrf_document()), class = "mzlib_usage_error")
})

test_that("a pool renders one tab-separated stdin line per document", {
  request <- mz$sdrf_build_pool_request(
    c(malaria = "a.sdrf.tsv", colon = "b.sdrf.tsv"), NULL, NULL, 0
  )
  expect_identical(request$stdin, "a.sdrf.tsv\tmalaria\nb.sdrf.tsv\tcolon")
  expect_identical(request$args, c("sdrf", "pool"))
})

test_that("an unnamed pool sends bare paths", {
  request <- mz$sdrf_build_pool_request(c("a.sdrf.tsv", "b.sdrf.tsv"), NULL, NULL, 0)
  expect_identical(request$stdin, "a.sdrf.tsv\nb.sdrf.tsv")
})

test_that("the fallback-provenance warning is passed through", {
  expect_false(any(grepl("reproducible", sdrf_pooled()$caveats, fixed = TRUE)))
  payload <- sdrf_fixture("sdrf_pool_two.json")
  payload$caveats <- c(payload$caveats, list("... not reproducible on another machine ..."))
  unlabelled <- mz$sdrf_parse_pooled(payload)
  expect_true(any(grepl("reproducible", unlabelled$caveats, fixed = TRUE)))
})

test_that("out writes the whole document and is reported", {
  request <- mz$sdrf_build_pool_request("a.sdrf.tsv", "merged.sdrf.tsv", 4, 0)
  expect_identical(request$args, c("sdrf", "pool", "--limit", "4", "--out", "merged.sdrf.tsv"))

  payload <- sdrf_fixture("sdrf_pool_two.json")
  payload$written <- list(path = "merged.sdrf.tsv", row_count = 24)
  written <- mz$sdrf_parse_pooled(payload)$written
  expect_identical(written, list(path = "merged.sdrf.tsv", row_count = 24))
})

# ---------------------------------------------------------------- argument validation

test_that("an empty selection is refused", {
  expect_error(mz$sdrf_build_pool_request(character(0), NULL, NULL, 0),
    class = "mzlib_usage_error", contains = "At least one"
  )
})

test_that("a list of paths is refused, with the shape that works", {
  expect_error(mz$sdrf_build_pool_request(list("a.sdrf.tsv"), NULL, NULL, 0),
    class = "mzlib_usage_error", contains = "character vector"
  )
})

test_that("a partly named selection is refused rather than silently defaulted", {
  # c(malaria = "a", "b") has names c("malaria", ""): half-labelled input is the trap, and the
  # bridge refuses it too, but failing here costs no subprocess and names the offending path.
  expect_error(mz$sdrf_build_pool_request(c(malaria = "a.sdrf.tsv", "b.sdrf.tsv"), NULL, NULL, 0),
    class = "mzlib_usage_error", contains = "label for 'b.sdrf.tsv'"
  )
})

test_that("a tab in a label is refused because it is the field separator", {
  expect_error(mz$sdrf_build_pool_request(c("mal\taria" = "a.sdrf.tsv"), NULL, NULL, 0),
    class = "mzlib_usage_error", contains = "tab or newline"
  )
})

test_that("a blank path is refused", {
  for (path in list("", "   ", NA_character_, NULL, 7)) {
    expect_error(mz$sdrf_build_read_args(path, NULL, 0), class = "mzlib_usage_error")
  }
  expect_error(mz$sdrf_build_pool_request(c("a.sdrf.tsv", " "), NULL, NULL, 0),
    class = "mzlib_usage_error", contains = "may not be blank"
  )
})

test_that("a bad limit or offset is refused", {
  for (limit in list(-1, 1.5, "3", TRUE, c(1, 2))) {
    expect_error(mz$sdrf_build_read_args("a.sdrf.tsv", limit, 0),
      class = "mzlib_usage_error", contains = "limit must be"
    )
  }
  for (offset in list(-1, 2.5, "0", NA)) {
    expect_error(mz$sdrf_build_read_args("a.sdrf.tsv", NULL, offset),
      class = "mzlib_usage_error", contains = "offset must be"
    )
  }
})

test_that("the window reaches the bridge, and a zero limit is allowed", {
  expect_identical(
    mz$sdrf_build_read_args("a.sdrf.tsv", 2, 1),
    c("sdrf", "read", "--path", "a.sdrf.tsv", "--limit", "2", "--offset", "1")
  )
  # The bridge accepts --limit 0, a header-only answer.
  expect_identical(
    mz$sdrf_build_read_args("a.sdrf.tsv", 0, 0),
    c("sdrf", "read", "--path", "a.sdrf.tsv", "--limit", "0")
  )
})

# ---------------------------------------------------------------- live

sdrf_live_file <- function() {
  root <- Sys.getenv("MZLIB_TEST_FILES")
  path <- file.path(root, "FileReadingTests", "ExternalFileTypes", "PXD000070.sdrf.tsv")
  skip_if(!nzchar(root) || !file.exists(path), "MZLIB_TEST_FILES does not hold PXD000070.sdrf.tsv")
  path
}

test_that("live: the recorded read fixture still matches the bridge", {
  skip_if_no_bridge()
  path <- sdrf_live_file()
  live <- sdrf_read(path)
  recorded <- sdrf_document()
  expect_identical(live$columns, recorded$columns)
  expect_identical(live$rows, recorded$rows)
  expect_identical(live$row_count, recorded$row_count)
  expect_identical(live$caveats, recorded$caveats)
})

test_that("live: pooling one document with itself keeps both copies apart", {
  skip_if_no_bridge()
  path <- sdrf_live_file()
  once <- sdrf_pool(c(first = path))
  twice <- sdrf_pool(c(path, path))
  expect_identical(once$document_count, 1)
  expect_identical(twice$document_count, 2)
  expect_identical(twice$row_count, 2 * once$row_count)
  expect_true(any(grepl("not reproducible", twice$caveats, fixed = TRUE)))
})
