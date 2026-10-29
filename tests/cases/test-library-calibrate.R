# library_calibrate_rt(), the projection of library calibrate-rt.
#
# Ported from pyMzLib's test_library_calibrate.py and test_library_calibrate_live.py. The payloads are
# pyMzLib's recordings, verbatim, from two DIA-NN reports in mzLib's test data, cut to six columns and
# shipped in extdata/calibration/ (its PROVENANCE.txt). The numbers are mzLib's IrtCalibration, held
# by the bridge's C# tests and re-run here LIVE. These hold the R layer: what it sends, what it
# parses, and what it refuses before starting a process.

calibration_path <- function(name) system.file("extdata", "calibration", name, package = "mzLibR")

calibration_tsv <- function(name) {
  utils::read.delim(calibration_path(name), stringsAsFactors = FALSE)
}

# The article's HeLa call: q-value 0.01 or below are anchors; the rest are queried both ways.
calibration_hela_inputs <- function() {
  rows <- calibration_tsv("diann2_hela_rt.tsv")
  confident <- rows$q_value <= 0.01
  others <- rows[!confident, ]
  list(
    anchors = data.frame(group = rows$run[confident], rt_minutes = rows$retention_time[confident],
      irt = rows$indexed_retention_time[confident], stringsAsFactors = FALSE),
    queries = rbind(
      data.frame(group = others$run, rt_minutes = NA_real_, irt = others$indexed_retention_time, stringsAsFactors = FALSE),
      data.frame(group = others$run, rt_minutes = others$retention_time, irt = NA_real_, stringsAsFactors = FALSE)
    ),
    others = others
  )
}

# The article's and the help page's: one precursor held out of each run's anchors and queried.
calibration_pbmc_inputs <- function() {
  rows <- calibration_tsv("diann_pbmc_rt.tsv")
  held <- rows$precursor_id == "AAAAAAAAGAFAGR2"
  list(
    anchors = data.frame(group = rows$run[!held], rt_minutes = rows$retention_time[!held],
      irt = rows$indexed_retention_time[!held], stringsAsFactors = FALSE),
    queries = data.frame(group = rows$run[held], irt = 11.886096, stringsAsFactors = FALSE),
    observed = stats::setNames(rows$retention_time[held], rows$run[held])
  )
}

calibration_request <- function(anchors, queries = NULL, ...) {
  options <- list(kind = NULL, minimum_anchors = NULL, bandwidth = NULL, knots = NULL, robustness_iterations = NULL,
    on_error = "fail")
  options[names(list(...))] <- list(...)
  mz$library_build_calibrate_request(anchors, queries, options$kind, options$minimum_anchors, options$bandwidth,
    options$knots, options$robustness_iterations, options$on_error)
}

# ---------------------------------------------------------------- what is sent

test_that("calibrate_rt sends the anchors, then a separator, then the queries with the computed time blank", {
  seen <- recorded_bridge_call("library_calibrate_rt_pbmc_linear.json", function() {
    library_calibrate_rt(
      list(list("run1", 10, 20.5), list(group = "run1", rt_minutes = 11, irt = 22)),
      list(list("run1", NA, 30), list(group = "run2", rt_minutes = 12.5), list("run1", NaN, 1))
    )
  })
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("library", "calibrate-rt"))
  expect_identical(seen$stdin, "run1\t10\t20.5\nrun1\t11\t22\n--\nrun1\t\t30\nrun2\t12.5\t\nrun1\t\t1\n")
  expect_identical(seen$timeout, 60)
})

test_that("a data.frame of queries may leave out the time it never gives", {
  request <- calibration_request(data.frame(group = "r", rt_minutes = 1, irt = 2),
    data.frame(group = c("r", "s"), irt = c(5, 6)))
  expect_identical(request$stdin, "r\t1\t2\n--\nr\t\t5\ns\t\t6\n")
  # A factor group, as read.delim() gives one before R 4.0, is its label.
  request <- calibration_request(data.frame(group = factor("r"), rt_minutes = 1, irt = 2))
  expect_identical(request$stdin, "r\t1\t2\n")
})

test_that("no queries sends no separator", {
  expect_identical(calibration_request(list(list("r", 1, 2)))$stdin, "r\t1\t2\n")
  expect_identical(calibration_request(list(list("r", 1, 2)), list())$stdin, "r\t1\t2\n")
})

test_that("a missing anchor time crosses blank, for mzLib to refuse", {
  expect_identical(calibration_request(list(list("r", NA, 2), list("r", 1, NaN)))$stdin, "r\t\t2\nr\t1\t\n")
})

test_that("a time crosses as the exact text mzLib reads back to the same double", {
  # The wire text itself, not as.numeric() of it: R reads 17-digit text through long double, which is
  # only double on Apple silicon, while mzLib (C#) parses it exactly.
  expect_identical(calibration_request(list(list("r", 0.1, 1 / 3)))$stdin, "r\t0.10000000000000001\t0.33333333333333331\n")
})

test_that("options are sent only when given", {
  expect_identical(calibration_request(list(list("r", 1, 2)))$args, c("library", "calibrate-rt"))
  every <- calibration_request(list(list("r", 1, 2)), kind = "lowess", minimum_anchors = 5, bandwidth = 0.5, knots = 7L,
    robustness_iterations = 0, on_error = "skip")
  expect_identical(every$args, c("library", "calibrate-rt", "--kind", "lowess", "--minimum-anchors", "5", "--bandwidth",
    "0.5", "--knots", "7", "--robustness-iterations", "0", "--on-error", "skip"))
  seen <- recorded_bridge_call("library_calibrate_rt_pbmc_linear.json", function() {
    library_calibrate_rt(list(list("r", 1, 2)), timeout = 9)
  })
  expect_identical(seen$timeout, 9)
})

test_that("bad input is refused before the bridge starts", {
  never <- function(args, stdin) stop("the bridge must not be called")
  refused <- function(contains, anchors = list(list("r", 1, 2)), queries = NULL, ...) {
    expect_error(recorded_bridge_call(never, function() library_calibrate_rt(anchors, queries, ...)),
      class = "mzlib_usage_error", contains = contains)
  }
  refused("at least one anchor", list())
  refused("at least one anchor", data.frame(group = character(0), rt_minutes = numeric(0), irt = numeric(0)))
  refused("must be a data.frame with group, rt_minutes and irt columns", "r\t1\t2")
  refused("must be list(group, rt_minutes, irt)", list(list("r", 1)))
  refused("needs a `group` column", data.frame(rt_minutes = 1, irt = 2))
  refused("needs a non-empty group name", list(list("", 1, 2)))
  refused("needs a non-empty group name", list(list(NA_character_, 1, 2)))
  refused("tab or newline", list(list("a\tb", 1, 2)))
  refused("separates anchors from queries", list(list("--", 1, 2)))
  refused("must be a number or NA", list(list("r", "early", 2)))
  refused("must be a number or NA", list(list("r", TRUE, 2)))
  refused("exactly one of rt_minutes and irt", queries = list(list("r", 1, 2)))
  refused("exactly one of rt_minutes and irt", queries = list(list("r", NA, NA)))
  refused("exactly one of rt_minutes and irt", queries = list(list(group = "r")))
  refused("kind must be one of linear, lowess", kind = "cubic")
  refused("minimum_anchors must be a whole number", minimum_anchors = 2.5)
  refused("knots must be a whole number", knots = TRUE)
  refused("bandwidth must be a number", bandwidth = "wide")
  refused("robustness_iterations must be a whole number", robustness_iterations = NA)
  refused("on_error must be 'fail' or 'skip'", on_error = "ignore")
  expect_error(library_calibrate_rt(), class = "mzlib_usage_error", contains = "at least one anchor")
})

test_that("an older bridge is named, not spawned, for calibrate-rt", {
  expect_error(
    recorded_bridge_call(function(args, stdin) stop("reached the verb"), function() {
      library_calibrate_rt(list(list("r", 1, 2)))
    }, version = "readers_formats.json"),
    class = "mzlib_usage_error", contains = "pyMzLib 0.7.0"
  )
})

# ---------------------------------------------------------------- what is parsed

test_that("a real run parses: mzLib's defaults, one group, the anchors' fitted iRT", {
  calibration <- mz$library_parse_calibration(recorded_data("library_calibrate_rt_hela.json"))
  expect_true(inherits(calibration, "mzlibr_rt_calibration"))
  expect_identical(c(calibration$kind, calibration$on_error), c("lowess", "fail"))
  expect_identical(c(calibration$minimum_anchors, calibration$bandwidth, calibration$knots,
    calibration$robustness_iterations), c(20, 0.3, 100, 3))
  expect_identical(c(calibration$group_count, calibration$fitted_count, calibration$failed_count), c(1, 1, 0))
  expect_identical(c(calibration$anchor_count, calibration$query_count, calibration$row_count), c(418, 228, 646))
  expect_identical(calibration$column_names, c("group", "role", "rt_minutes", "irt", "fitted_irt", "is_extrapolated"))
  run <- calibration$groups
  expect_identical(nrow(run), 1L)
  expect_identical(c(run$group, run$kind), c("Fig2HeLa-0-5h_MHRM_R01_T0", "lowess"))
  expect_identical(c(run$anchor_count, run$query_count), c(418, 228))
  expect_true(is.na(run$error_kind) && is.na(run$error_type) && is.na(run$error_message))
  expect_true(abs(run$residual_sd - 5.0035) < 1e-4)
  expect_true(abs(run$first_anchor_minutes - 15.988) < 1e-3)
  expect_true(abs(run$last_anchor_minutes - 40.903) < 1e-3)
  rows <- calibration$records
  expect_identical(nrow(rows), 646L)
  expect_identical(sort(unique(rows$role)), c("anchor", "to_irt", "to_minutes"))
  expect_true(all(is.na(rows$fitted_irt[rows$role != "anchor"])))
  expect_true(is.numeric(rows$fitted_irt) && !anyNA(rows$fitted_irt[rows$role == "anchor"]))
  expect_true(is.logical(rows$is_extrapolated))
  expect_identical(sum(rows$is_extrapolated[rows$role != "anchor"]), 5L)
  expect_identical(calibration$caveats[[length(calibration$caveats)]],
    "5 queries fall outside the anchors' span, where the map is extrapolated (is_extrapolated).")
})

test_that("refused runs parse: each group's error, and no fitted value anywhere", {
  calibration <- mz$library_parse_calibration(recorded_data("library_calibrate_rt_pbmc_skip.json"))
  expect_identical(calibration$on_error, "skip")
  expect_identical(c(calibration$failed_count, calibration$fitted_count), c(3, 0))
  groups <- calibration$groups
  expect_identical(groups$anchor_count, c(13, 11, 10))
  expect_true(all(is.na(groups$kind) & is.na(groups$residual_sd) & is.na(groups$first_anchor_minutes) &
    is.na(groups$last_anchor_minutes)))
  expect_identical(unique(groups$error_kind), "usage")
  expect_identical(unique(groups$error_type), "ArgumentException")
  expect_true(all(startsWith(groups$error_message, "Calibration needs at least 20 anchors")))
  expect_identical(nrow(groups[!is.na(groups$error_kind), ]), 3L)
  expect_true(all(is.na(calibration$records$fitted_irt)) && length(calibration$records$fitted_irt) == calibration$row_count)
  expect_true(all(is.na(calibration$records$is_extrapolated)))
  expect_true(any(startsWith(calibration$caveats, "3 groups were not fitted")))
})

test_that("a line per run parses, with each run's prediction", {
  calibration <- mz$library_parse_calibration(recorded_data("library_calibrate_rt_pbmc_linear.json"))
  expect_identical(c(calibration$kind, calibration$minimum_anchors), c("linear", "10"))
  expect_identical(calibration$groups$kind, rep("linear", 3L))
  queries <- calibration$records[calibration$records$role == "to_minutes", ]
  expect_identical(round(queries$rt_minutes, 2), c(9.22, 9.45, 9.17))
  expect_true(all(queries$irt == 11.886096))
})

test_that("the kinds are mzLib's", {
  expect_identical(mz$LIBRARY_CALIBRATION_KINDS, c("linear", "lowess"))
})

# ---------------------------------------------------------------- end to end, replayed

test_that("the calibrate-rt examples replay end to end, an option choosing the recording", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  hela <- calibration_hela_inputs()
  expect_identical(library_calibrate_rt(hela$anchors, hela$queries)$row_count, 646)
  pbmc <- calibration_pbmc_inputs()
  expect_identical(library_calibrate_rt(pbmc$anchors, pbmc$queries, on_error = "skip")$failed_count, 3)
  lines <- library_calibrate_rt(pbmc$anchors, pbmc$queries, kind = "linear", minimum_anchors = 10)
  expect_identical(lines$fitted_count, 3)
  # Without the options it was recorded with, or with other anchors, no recording answers.
  expect_error(library_calibrate_rt(pbmc$anchors, pbmc$queries, kind = "linear"),
    class = "mzlib_usage_error", contains = "replay bridge")
  expect_error(library_calibrate_rt(pbmc$anchors[-1L, ], pbmc$queries, on_error = "skip"),
    class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

calibration_live <- function() {
  skip_unless_bridge_has("library calibrate-rt")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  invisible(NULL)
}

test_that("LIVE: the calibrate-rt recordings are what the bridge returns today", {
  calibration_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  hela <- calibration_hela_inputs()
  pbmc <- calibration_pbmc_inputs()
  calls <- list(
    library_calibrate_rt_hela.json = function() library_calibrate_rt(hela$anchors, hela$queries),
    library_calibrate_rt_pbmc_skip.json = function() library_calibrate_rt(pbmc$anchors, pbmc$queries, on_error = "skip"),
    library_calibrate_rt_pbmc_linear.json = function() {
      library_calibrate_rt(pbmc$anchors, pbmc$queries, kind = "linear", minimum_anchors = 10)
    }
  )
  for (name in names(calls)) {
    live <- calls[[name]]()
    recorded <- mz$library_parse_calibration(recorded_data(name))
    # Closed-form arithmetic on the anchors, so 1e-9 across platforms.
    expect_identical(records_disagree(live$records, recorded$records), "", info = name)
    expect_identical(records_disagree(live$groups, recorded$groups), "", info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
    for (echo in c("kind", "minimum_anchors", "bandwidth", "knots", "robustness_iterations", "on_error", "group_count",
                   "fitted_count", "failed_count", "anchor_count", "query_count", "row_count")) {
      expect_identical(live[[echo]], recorded[[echo]], info = paste(name, echo))
    }
  }
})

test_that("LIVE: the article's numbers hold", {
  calibration_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  pbmc <- calibration_pbmc_inputs()
  lines <- library_calibrate_rt(pbmc$anchors, pbmc$queries, kind = "linear", minimum_anchors = 10)
  predicted <- lines$records[lines$records$role == "to_minutes", ]
  # Each run's line predicts the held-out precursor within a quarter of a minute.
  expect_true(all(abs(predicted$rt_minutes - pbmc$observed[predicted$group]) < 0.25))
})

test_that("LIVE: mzLib's refusals cross as usage errors", {
  calibration_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  pbmc <- calibration_pbmc_inputs()
  condition <- expect_error(library_calibrate_rt(pbmc$anchors, pbmc$queries), class = "mzlib_usage_error",
    contains = "needs at least 20 anchors; got 13")
  expect_true(grepl("Group '20240316_PBMCs", conditionMessage(condition), fixed = TRUE))
  expect_error(library_calibrate_rt(pbmc$anchors, pbmc$queries, knots = 2, on_error = "skip"),
    class = "mzlib_usage_error", contains = "Knots must be at least 3")
  expect_error(library_calibrate_rt(pbmc$anchors, bandwidth = 1.5), class = "mzlib_usage_error",
    contains = "Bandwidth must be in")
  expect_error(library_calibrate_rt(list(list("r", 1, 1), list("r", 2, NA), list("r", 3, 3)), minimum_anchors = 3),
    class = "mzlib_usage_error", contains = "must be finite")
})
