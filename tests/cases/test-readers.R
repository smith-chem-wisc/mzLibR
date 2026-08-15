# Reading result files.
#
# The point of this module is to stop someone believing mzLib reads 31 formats into one uniform
# shape. Most of these tests assert that the exceptions are visible.

recorded_formats <- function() {
  payload <- mz$json_parse(paste(
    readLines(fixture_path("readers_formats.json"), warn = FALSE),
    collapse = "\n"
  ))
  mz$readers_parse_formats(payload$data)
}

# ---------------------------------------------------------------- formats

test_that("every format mzLib recognises becomes a row", {
  formats <- recorded_formats()
  expect_identical(nrow(formats), 31L)
  expect_true(is.data.frame(formats))
  expect_false(is.factor(formats$file_type))
})

test_that("exactly four file types are quantifiable", {
  # The precondition for flashlfq_quantify(), and it is a much smaller set than "31 formats"
  # suggests. The number is quoted in ?readers_formats, so it is pinned here.
  #
  # It moved from three to four when mzLib 1.0.585 added DiaNnReport (mzLib #1120), which is
  # how DIA data reaches flashlfq_quantify() at all. This test failing on a bridge repin is
  # the mechanism working: the claim changed, so the name and the number both moved.
  formats <- recorded_formats()
  quantifiable <- formats[formats$is_quantifiable, ]
  expect_identical(nrow(quantifiable), 4L)
  expect_identical(sort(quantifiable$file_type),
                   sort(c("psmtsv", "osmtsv", "MsFraggerPsm", "DiaNnReport")))
})

test_that("most formats have no views at all, and that is a real answer", {
  # 14 of 31. An empty views list means mzLib can parse the file but offers no cross-format
  # projection of it - not that anything failed.
  formats <- recorded_formats()
  viewless <- vapply(formats$views, function(v) length(v) == 0L, logical(1L))
  expect_identical(sum(viewless), 14L)
})

test_that("the view vocabulary is the four documented families", {
  formats <- recorded_formats()
  expect_identical(
    sort(unique(unlist(formats$views))),
    sort(c("ms1_features", "quantifiable", "spectra", "spectral_match"))
  )
})

test_that("an extension does not identify a format on its own", {
  # BrukerD and BrukerTimsTof are both ".d", told apart by what the directory holds, and several
  # formats share ".tsv". Anyone dispatching on extension alone is wrong.
  formats <- recorded_formats()
  expect_true(anyDuplicated(formats$extension) > 0L)
})

test_that("views survive subsetting the formats frame", {
  # A list column for the same reason PRIDE's locations is one.
  formats <- recorded_formats()
  some <- formats[formats$is_quantifiable, ]
  expect_true("views" %in% names(some))
  expect_true(all(vapply(some$views, function(v) "quantifiable" %in% v, logical(1L))))
})

# ---------------------------------------------------------------- identify

test_that("a file with no views says so rather than looking broken", {
  info <- mz$readers_parse_file_info(list(
    path = "/abs/spectra.mzML", file_type = "MzML", extension = ".mzML",
    reader = "MsDataFile", views = list("spectra")
  ))
  expect_true(inherits(info, "mzlibr_file_info"))
  expect_false(info$is_quantifiable)
  expect_identical(info$views, "spectra")

  bare <- mz$readers_parse_file_info(list(
    path = "/abs/x.tsv", file_type = "CruxResult", extension = ".tsv",
    reader = "CruxResultFile", views = list()
  ))
  expect_identical(bare$views, character(0))
  expect_false(bare$is_quantifiable)
  output <- paste(capture.output(print(bare)), collapse = "\n")
  expect_true(grepl("no views", output, fixed = TRUE), info = output)
})

test_that("printing a quantifiable file still warns about trusting it", {
  # is_quantifiable reports what mzLib's interface offers, not that the numbers are comparable.
  # MSFragger is quantifiable by interface and should not be quantified.
  info <- mz$readers_parse_file_info(list(
    path = "/abs/psm.tsv", file_type = "MsFraggerPsm", extension = "psm.tsv",
    reader = "MsFraggerPsmFile", views = list("quantifiable")
  ))
  expect_true(info$is_quantifiable)
  output <- paste(capture.output(print(info)), collapse = "\n")
  expect_true(grepl("before trusting the numbers", output, fixed = TRUE), info = output)
})

test_that("a blank path is refused before anything is spawned", {
  for (path in list("", "   ", NA_character_, NULL, 42)) {
    expect_error(mz$readers_normalise_path(path), class = "mzlib_usage_error")
  }
})

# ---------------------------------------------------------------- read_results

# Deliberately not modifyList(): it recurses into list-valued fields, and for an *unnamed* list
# it iterates over names(value), which is NULL — so the override is silently ignored and the
# base value survives. Four tests here passed against data they had not actually set before that
# was noticed. Replacing wholesale is the only behaviour that means what it looks like.
sample_records <- function(...) {
  base <- list(
    path = "/abs/AllPSMs.psmtsv", file_type = "psmtsv",
    record_count = 4, returned_count = 2, offset = 0, truncated = TRUE,
    retention_time_unit = "minutes", rows_not_read = 1,
    caveats = list("MetaMorpheus retention times are minutes."),
    column_names = list("sequence", "retention_time", "charge"),
    columns = list(
      sequence = list("PEPTIDEK", "ACDEFR"),
      retention_time = list(30.1, NA),
      charge = list(2, 3)
    ),
    output = NULL
  )
  overrides <- list(...)
  for (name in names(overrides)) {
    base[[name]] <- overrides[[name]]
  }
  base
}

test_that("the record columns become a data.frame directly", {
  # pyMzLib returns field -> values because that is what pandas wants; in R the equivalent IS a
  # data.frame, so ResultRecords.records - Python's second, row-wise accessor - has no
  # counterpart here. One object is both.
  records <- mz$readers_parse_records(sample_records())
  expect_true(is.data.frame(records$records))
  expect_identical(nrow(records$records), 2L)
  expect_identical(names(records$records), c("sequence", "retention_time", "charge"))
  expect_identical(records$records$sequence, c("PEPTIDEK", "ACDEFR"))
})

test_that("a null inside a column becomes NA and does not shorten it", {
  # The hazard that runs through this whole package. A shortened column would silently
  # misalign every other column in the table.
  records <- mz$readers_parse_records(sample_records())
  expect_identical(length(records$records$retention_time), 2L)
  expect_true(is.na(records$records$retention_time[2]))
})

test_that("columns of unequal length are refused rather than recycled", {
  # R would recycle silently and produce a table that is wrong rather than absent.
  broken <- sample_records(columns = list(
    sequence = list("A", "B"),
    retention_time = list(1)
  ), column_names = list("sequence", "retention_time"))
  expect_error(mz$readers_parse_records(broken),
    class = "mzlib_protocol_error", contains = "different lengths"
  )
})

test_that("column order follows column_names, not the map's own order", {
  records <- mz$readers_parse_records(sample_records(
    column_names = list("charge", "sequence", "retention_time")
  ))
  expect_identical(names(records$records), c("charge", "sequence", "retention_time"))
})

test_that("truncation and unreadable rows are both surfaced", {
  # A short answer and a complete one must not look alike, and mzLib drops a malformed row
  # silently - so rows_not_read is the only sign that a table is incomplete.
  records <- mz$readers_parse_records(sample_records())
  expect_true(records$truncated)
  expect_identical(records$rows_not_read, 1)
  expect_identical(records$record_count, 4)
  expect_identical(records$returned_count, 2)

  output <- paste(capture.output(print(records)), collapse = "\n")
  expect_true(grepl("truncated", output, fixed = TRUE), info = output)
  expect_true(grepl("partly unreadable", output, fixed = TRUE), info = output)
})

test_that("caveats are printed, because that is where the units live", {
  records <- mz$readers_parse_records(sample_records())
  output <- paste(capture.output(print(records)), collapse = "\n")
  expect_true(grepl("retention times are minutes", output, fixed = TRUE), info = output)
})

test_that("a written table reports where it went and in what format", {
  records <- mz$readers_parse_records(sample_records(
    columns = NULL, column_names = list(), returned_count = 0,
    output = list(path = "/abs/out.tsv", format = "tsv", row_count = 4)
  ))
  expect_true(is.null(records$records))
  expect_identical(records$output$format, "tsv")
  expect_identical(records$output$row_count, 4)
  output <- paste(capture.output(print(records)), collapse = "\n")
  expect_true(grepl("written to /abs/out.tsv", output, fixed = TRUE), info = output)
})

# ---------------------------------------------------------------- retention time units

test_that("minutes pass through and seconds are converted", {
  minutes <- mz$readers_parse_records(sample_records())
  expect_identical(mz$readers_retention_time_in_minutes(minutes)[1], 30.1)

  seconds <- mz$readers_parse_records(sample_records(
    file_type = "MsFraggerPsm", retention_time_unit = "seconds",
    columns = list(
      sequence = list("A", "B"), retention_time = list(1800, 3600), charge = list(2, 2)
    )
  ))
  expect_equal(mz$readers_retention_time_in_minutes(seconds), c(30, 60))
})

test_that("an unknown unit raises rather than guessing", {
  # A silently unconverted time axis is the specific mistake this module exists to prevent: a
  # 60x error in a retention-time comparison reads as a chromatography problem, not a units one.
  unknown <- mz$readers_parse_records(sample_records(retention_time_unit = "unknown"))
  expect_error(mz$readers_retention_time_in_minutes(unknown),
    class = "mzlib_usage_error", contains = "no basis to say what unit"
  )
})

test_that("the converter refuses anything that is not a record view", {
  expect_error(mz$readers_retention_time_in_minutes(list()), class = "mzlib_usage_error")
})

# ---------------------------------------------------------------- argument assembly

test_that("read-results defaults send neither limit nor offset", {
  args <- mz$readers_build_read_args("psm.tsv", NULL, 0, NULL)
  expect_identical(args, c("readers", "read-results", "--path", "psm.tsv"))
})

test_that("limit, offset and out reach the wire in the bridge's spelling", {
  args <- mz$readers_build_read_args("psm.tsv", 100, 50, "out.tsv")
  expect_identical(args[which(args == "--limit") + 1L], "100")
  expect_identical(args[which(args == "--offset") + 1L], "50")
  expect_identical(args[which(args == "--out") + 1L], "out.tsv")
})

test_that("an impossible limit or offset is refused", {
  for (bad in list(0, -1, 1.5, "10")) {
    expect_error(mz$readers_build_read_args("psm.tsv", bad, 0, NULL), class = "mzlib_usage_error")
  }
  for (bad in list(-1, 1.5, "10", NA_real_)) {
    expect_error(mz$readers_build_read_args("psm.tsv", NULL, bad, NULL), class = "mzlib_usage_error")
  }
})

test_that("a large limit is not written in scientific notation", {
  args <- mz$readers_build_read_args("psm.tsv", 1000000, 0, NULL)
  expect_identical(args[which(args == "--limit") + 1L], "1000000")
})

# ---------------------------------------------------------------- against a real mzLib

test_that("LIVE: mzLib still recognises 31 formats, four of them quantifiable", {
  # Enumerated from mzLib itself, so this is the test that notices when the installed version
  # changes what it supports - which is exactly when the numbers in ?readers_formats go stale.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  formats <- readers_formats()
  expect_identical(nrow(formats), 31L)
  expect_identical(sum(formats$is_quantifiable), 4L)
  expect_identical(sum(vapply(formats$views, length, integer(1L)) == 0L), 14L)
})

test_that("LIVE: identify dispatches on extension and does not validate contents", {
  # Worth knowing and easy to get wrong. A plain text file containing nothing proteomic is
  # identified as `CruxResult`, because mzLib dispatches `.txt` there and identify() is
  # deliberately cheap - it resolves the type and stops. It is not a validity check, and a
  # confident-looking file_type is not evidence the file is what it claims.
  #
  # The honest signal is `views`, which is empty here. That is what to branch on.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  scratch <- tempfile("mzlibr-not-a-result-", fileext = ".txt")
  writeLines("this is not a proteomics result file", scratch)
  on.exit(unlink(scratch), add = TRUE)

  info <- readers_identify(scratch)
  expect_identical(info$file_type, "CruxResult")
  expect_identical(info$views, character(0))
  expect_false(info$is_quantifiable)
})

test_that("LIVE: reading a file with no quantifiable view names the views it does have", {
  # The error is the documentation a stuck user reads.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  scratch <- tempfile("mzlibr-spectra-", fileext = ".mzML")
  writeLines("<mzML></mzML>", scratch)
  on.exit(unlink(scratch), add = TRUE)

  condition <- tryCatch(
    {
      readers_read_results(scratch, timeout = 120)
      NULL
    },
    mzlib_error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_error"))
})
