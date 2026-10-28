# stats_fit_logistic(), stats_fit_mixed() and stats_entrapment_fdp().
#
# Ported from pyMzLib's test_stats_models.py and its _live.py. The payloads are pyMzLib's recordings,
# verbatim; the reference tables are mzLib's, shipped in extdata/stats/. Three layers:
#   1. always: what is sent (nothing for an option left out, so mzLib's default applies), what is
#      parsed, refusals before the bridge, and the recordings against R's glm() and nlme::lme() run
#      here and against the references mzLib's tests freeze;
#   2. replayed end to end through the stand-in bridge;
#   3. LIVE: the bridge today returns what was recorded, and refuses alike.

models_path <- function(name) system.file("extdata", "stats", name, package = "mzLibR")
models_reference <- function(name) utils::read.delim(models_path(name), check.names = FALSE, stringsAsFactors = FALSE)
models_relative <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)

# ---------------------------------------------------------------- fit_logistic

test_that("fit_logistic sends only what was given", {
  request <- mz$stats_build_logistic_request("r.tsv", "d.tsv", c("x", "grp"), NULL, NULL, 1)
  expect_identical(request$args, c("stats", "fit-logistic", "--responses", "r.tsv", "--design", "d.tsv", "--threads", "1"))
  expect_identical(request$stdin, "x\ngrp")
  given <- mz$stats_build_logistic_request("r.tsv", "d.tsv", "x", 25, 1e-8, -1)
  expect_identical(utils::tail(given$args, 6L)[c(1, 2, 3, 4, 5)], c("--threads", "-1", "--max-iterations", "25", "--tolerance"))
  expect_identical(as.numeric(utils::tail(given$args, 1L)), 1e-8)
  expect_true(is.null(formals(stats_fit_logistic)$max_iterations) && is.null(formals(stats_fit_logistic)$tolerance))
})

test_that("bad fit_logistic arguments are refused before the bridge", {
  build <- function(...) {
    given <- list(responses = "r.tsv", design = "d.tsv", coefficients = "x", max_iterations = NULL, tolerance = NULL, threads = 1)
    given[names(list(...))] <- list(...)
    do.call(mz$stats_build_logistic_request, given)
  }
  expect_error(build(coefficients = character(0)), class = "mzlib_usage_error", contains = "at least one coefficient")
  expect_error(build(coefficients = ""), class = "mzlib_usage_error", contains = "non-empty")
  expect_error(build(coefficients = "a\nb"), class = "mzlib_usage_error", contains = "newline")
  expect_error(build(responses = 5), class = "mzlib_usage_error", contains = "responses must be")
  expect_error(build(max_iterations = 2.5), class = "mzlib_usage_error", contains = "max_iterations must be a whole number")
  expect_error(build(tolerance = "tight"), class = "mzlib_usage_error", contains = "tolerance must be a number")
  expect_error(build(threads = "all"), class = "mzlib_usage_error", contains = "threads must be a whole number")
})

test_that("the glm reference recording parses, one block per coefficient", {
  fit <- mz$stats_parse_logistic(recorded_data("stats_fit_logistic_glm.json"))
  expect_true(inherits(fit, "mzlibr_logistic_model_fit"))
  expect_identical(fit$tested_coefficients, c("Intercept", "x", "grp"))
  expect_true(is.na(fit$max_iterations) && is.na(fit$tolerance), info = "left out: mzLib's defaults, echoed as null")
  expect_identical(names(fit$status_counts), mz$STATS_LOGISTIC_STATUSES)
  expect_identical(c(sum(fit$status_counts), fit$feature_count, fit$row_count), c(80, 80, 240))
  expect_identical(fit$status_counts[["separated"]], 4)
  expect_identical(fit$column_names, c("feature", "coefficient", "status", "estimate", "standard_error", "z",
    "p_value", "deviance", "iterations", "observed", "events"))
  x <- fit$records[fit$records$coefficient == "x", ]
  expect_true(all(x$status %in% mz$STATS_LOGISTIC_STATUSES))
  fitted <- x$status == "fitted"
  for (column in c("estimate", "standard_error", "z", "p_value", "deviance")) {
    expect_identical(!is.na(x[[column]]), fitted, info = column)
  }
  expect_true(all(x$events >= 0 & x$events <= x$observed & x$observed <= fit$sample_count))
  expect_identical(c(round(x$estimate[[1L]], 4), round(x$p_value[[1L]], 4), x$events[[1L]], x$observed[[1L]]),
    c(0.4592, 0.2724, 25, 35))
})

test_that("the recorded logistic fit is glm()'s recorded output: 1e-8 in estimates, 1e-6 in errors", {
  fit <- mz$stats_parse_logistic(recorded_data("stats_fit_logistic_glm.json"))
  glm <- models_reference("glm_reference_fit.tsv")
  for (j in 0:2) {
    rows <- fit$records[fit$records$coefficient == c("Intercept", "x", "grp")[j + 1L], ]
    separated <- glm$status == "separated"
    expect_identical(rows$status[separated], rep("separated", sum(separated)))
    expect_true(all(rows$status[!separated] == "fitted"))
    expect_true(max(models_relative(rows$estimate[!separated], glm[[paste0("b", j)]][!separated])) < 1e-8)
    expect_true(max(models_relative(rows$standard_error[!separated], glm[[paste0("se", j)]][!separated])) < 1e-6)
  }
  x <- fit$records[fit$records$coefficient == "x", ]
  fitted <- glm$status == "fitted"
  expect_true(max(models_relative(x$p_value[fitted], glm$p1[fitted])) < 1e-6)
  expect_true(max(models_relative(x$deviance[fitted], glm$deviance[fitted])) < 1e-8)
})

test_that("the recorded logistic fit is glm(family = binomial) run here, feature by feature", {
  fit <- mz$stats_parse_logistic(recorded_data("stats_fit_logistic_glm.json"))
  responses <- models_reference("glm_reference_responses.tsv")
  design <- models_reference("glm_reference_design.tsv")
  y <- as.matrix(responses[, -1L])
  x <- as.matrix(design[match(colnames(y), design[[1L]]), c("Intercept", "x", "grp")])
  rows <- fit$records[fit$records$coefficient == "x", ]
  worst <- 0
  for (i in which(rows$status == "fitted")) {
    keep <- is.finite(y[i, ])
    model <- stats::glm.fit(x[keep, , drop = FALSE], y[i, keep], family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-14, maxit = 100))
    worst <- max(worst, models_relative(rows$estimate[[i]], unname(model$coefficients["x"])))
  }
  expect_true(worst < 1e-8, info = paste("worst relative difference", worst))
})

test_that("the dilution series' detection table separates some oligos, and the caveats say so", {
  fit <- mz$stats_parse_logistic(recorded_data("stats_fit_logistic_malat.json"))
  expect_true(fit$status_counts[["separated"]] > 0)
  expect_true(any(startsWith(fit$caveats, "separated: the design splits")))
  expect_true(any(grepl("were not fitted", fit$caveats, fixed = TRUE) & grepl("separated", fit$caveats, fixed = TRUE)))
})

# ---------------------------------------------------------------- fit_mixed

test_that("fit_mixed sends the estimator only when given", {
  request <- mz$stats_build_mixed_request("r.tsv", "d.tsv", "g.tsv", "age", NULL, 1)
  expect_identical(request$args, c("stats", "fit-mixed", "--responses", "r.tsv", "--design", "d.tsv",
    "--groups", "g.tsv", "--threads", "1"))
  expect_identical(request$stdin, "age")
  expect_identical(utils::tail(mz$stats_build_mixed_request("r.tsv", "d.tsv", "g.tsv", "age", "ml", 1)$args, 2L),
    c("--estimator", "ml"))
  expect_error(mz$stats_build_mixed_request("r.tsv", "d.tsv", "g.tsv", "age", "REML", 1), class = "mzlib_usage_error",
    contains = "estimator must be one of reml, ml")
  expect_error(mz$stats_build_mixed_request("r.tsv", "d.tsv", "g.tsv", "age", "lme4", 1), class = "mzlib_usage_error",
    contains = "estimator must be one of reml, ml")
  expect_error(mz$stats_build_mixed_request("r.tsv", "d.tsv", "g.tsv", "age", NULL, 1.5), class = "mzlib_usage_error",
    contains = "threads must be")
  expect_error(stats_fit_mixed("r.tsv", "d.tsv", NULL, "age"), class = "mzlib_usage_error", contains = "groups must be")
})

test_that("the nlme reference recordings parse, REML by mzLib's default and ML when asked", {
  fit <- mz$stats_parse_mixed(recorded_data("stats_fit_mixed_lmm_reml.json"))
  expect_true(inherits(fit, "mzlibr_mixed_model_fit"))
  expect_identical(fit$estimator, "reml", info = "left out: mzLib's default, as mzLib reports it")
  expect_identical(c(fit$grouping, fit$group_count), c("group", "7"))
  expect_identical(fit$group_names, paste0("d", 1:7))
  expect_identical(names(fit$status_counts), mz$STATS_MIXED_STATUSES)
  expect_identical(c(sum(fit$status_counts), fit$feature_count), c(120, 120))
  expect_identical(fit$column_names, c("feature", "coefficient", "status", "estimate", "standard_error", "df", "t",
    "p_value", "residual_variance", "group_variance", "log_likelihood", "observed", "groups"))
  wp <- fit$records[fit$records$coefficient == "wp", ]
  expect_true(all(wp$df %in% c(4, 5) & wp$groups %in% c(6, 7) & wp$group_variance >= 0))
  ml <- mz$stats_parse_mixed(recorded_data("stats_fit_mixed_lmm_ml.json"))
  expect_identical(ml$estimator, "ml")
  expect_true(any(startsWith(ml$caveats, "ml: residual_variance is the maximum-likelihood estimate")))
})

test_that("the recorded mixed fits are nlme's recorded lme() output, to 1e-5 with identical df", {
  nlme <- models_reference("lmm_reference_fit.tsv")
  for (method in c("REML", "ML")) {
    fit <- mz$stats_parse_mixed(recorded_data(paste0("stats_fit_mixed_lmm_", tolower(method), ".json")))
    reference <- nlme[nlme$method == method, ]
    for (j in 0:2) {
      rows <- fit$records[fit$records$coefficient == c("Intercept", "age", "wp")[j + 1L], ]
      expect_true(all(rows$status == "fitted" & reference$status == "fitted"))
      expect_true(max(models_relative(rows$estimate, reference[[paste0("b", j)]])) < 1e-5, info = method)
      expect_true(max(models_relative(rows$standard_error, reference[[paste0("se", j)]])) < 1e-5, info = method)
      expect_identical(rows$df, as.numeric(reference[[paste0("df", j)]]), info = method)
    }
    age <- fit$records[fit$records$coefficient == "age", ]
    expect_true(max(models_relative(age$log_likelihood, reference$loglik)) < 1e-5, info = method)
    expect_true(all(abs(age$group_variance - reference$tau2) < 1e-5 * reference$sigma2), info = method)
  }
})

test_that("the recorded REML fit is nlme::lme() run here, to 1e-5", {
  skip_if(!requireNamespace("nlme", quietly = TRUE), "nlme is not installed")
  fit <- mz$stats_parse_mixed(recorded_data("stats_fit_mixed_lmm_reml.json"))
  responses <- models_reference("lmm_reference_responses.tsv")
  design <- models_reference("lmm_reference_design.tsv")
  groups <- models_reference("lmm_reference_groups.tsv")
  y <- as.matrix(responses[, -1L])
  samples <- colnames(y)
  frame <- data.frame(design[match(samples, design[[1L]]), c("Intercept", "age", "wp")],
    group = groups[[2L]][match(samples, groups[[1L]])])
  age <- fit$records[fit$records$coefficient == "age", ]
  worst <- 0
  for (i in seq_len(10L)) {
    frame$y <- y[i, ]
    model <- nlme::lme(y ~ 0 + Intercept + age + wp, random = ~ 1 | group, data = frame, method = "REML",
      na.action = stats::na.omit)
    worst <- max(worst, models_relative(age$estimate[[i]], unname(nlme::fixef(model)["age"])))
  }
  expect_true(worst < 1e-5, info = paste("worst relative difference", worst))
})

test_that("a group-level coefficient with no positive df carries t but no p-value", {
  fit <- mz$stats_parse_mixed(recorded_data("stats_fit_mixed_malat.json"))
  no_df <- fit$records[fit$records$status == "fitted" & is.na(fit$records$df), ]
  expect_true(nrow(no_df) > 0L, info = "the dilution series has group-level coefficients for some oligos")
  expect_true(all(!is.na(no_df$t) & is.na(no_df$p_value)))
  expect_true(any(endsWith(fit$caveats, "estimate, standard_error and t are reported.")))
})

# ---------------------------------------------------------------- entrapment_fdp

test_that("entrapment_fdp sends the ratio, and the thresholds only when given", {
  args <- mz$stats_build_entrapment_args("t.tsv", "e.tsv", 1, NULL)
  expect_identical(args, c("stats", "entrapment-fdp", "--targets", "t.tsv", "--entrapments", "e.tsv", "--ratio", "1"))
  given <- mz$stats_build_entrapment_args("t.tsv", "e.tsv", 2, c(0.01, 0.05))
  expect_identical(utils::tail(given, 4L)[1:3], c("--ratio", "2", "--thresholds"))
  expect_identical(as.numeric(strsplit(utils::tail(given, 1L), ",", fixed = TRUE)[[1L]]), c(0.01, 0.05))
  expect_error(stats_entrapment_fdp("t.tsv", "e.tsv"), class = "mzlib_usage_error", contains = "ratio is required")
  expect_error(stats_entrapment_fdp("t.tsv", "e.tsv", ratio = "1"), class = "mzlib_usage_error", contains = "ratio must be a number")
  expect_error(stats_entrapment_fdp("t.tsv", "e.tsv", ratio = TRUE), class = "mzlib_usage_error", contains = "ratio must be a number")
  expect_error(stats_entrapment_fdp("t.tsv", "e.tsv", ratio = 1, thresholds = numeric(0)), class = "mzlib_usage_error",
    contains = "non-empty numeric vector of q-values")
  expect_error(stats_entrapment_fdp("t.tsv", "e.tsv", ratio = 1, thresholds = c("0.01", "5%")), class = "mzlib_usage_error",
    contains = "non-empty numeric vector of q-values")
})

test_that("the entrapment reference recording parses, at mzLib's default thresholds", {
  fdp <- mz$stats_parse_entrapment(recorded_data("stats_entrapment_fdp_reference.json"))
  expect_true(inherits(fdp, "mzlibr_entrapment_fdp_estimates"))
  expect_identical(fdp$ratio, 1)
  expect_false(fdp$thresholds_given)
  expect_true(isTRUE(all.equal(fdp$thresholds, (1:100) / 1000)), info = "mzLib's DefaultQValueThresholds")
  expect_identical(c(fdp$target_candidate_count, fdp$entrapment_candidate_count, fdp$partnered_count, fdp$foreign_count),
    c(7, 4, 3, 0))
  expect_identical(c(fdp$row_count, fdp$threshold_count - fdp$empty_threshold_count), c(100, 100))
  at <- function(q) fdp$records[abs(fdp$records$q_value_threshold - q) < 1e-12, ]
  expect_identical(c(at(0.003)$target_count, at(0.003)$entrapment_count), c(2, 1))
  expect_true(abs(at(0.003)$lower_bound - 1 / 3) < 1e-12 && abs(at(0.003)$paired - 1 / 3) < 1e-12)
  expect_identical(round(c(at(0.01)$lower_bound, at(0.01)$combined, at(0.01)$paired), 5), c(0.42857, 0.85714, 0.85714))
  expect_identical(nrow(at(0.2)), 0L)
})

test_that("paired is NA off a ratio of 1, and the caveat says why", {
  fdp <- mz$stats_parse_entrapment(recorded_data("stats_entrapment_fdp_reference_r2.json"))
  expect_true(fdp$thresholds_given)
  expect_identical(fdp$thresholds, c(0.01, 0.05))
  expect_identical(fdp$records$paired, c(NA, NA))
  expect_true(abs(fdp$records$combined[fdp$records$q_value_threshold == 0.01] - 3 * 1.5 / 7) < 1e-12)
  expect_identical(fdp$caveats[[length(fdp$caveats)]],
    "paired is null at every threshold. Here the ratio is 2 and 0 of 4 entrapments are foreign.")
})

test_that("the recorded estimates are the entrapment project's reference script's, to 1e-5", {
  reference <- models_reference("entrapment_reference_r1.tsv")
  fdp <- mz$stats_parse_entrapment(recorded_data("stats_entrapment_fdp_reference.json"))
  r <- fdp$records[match(round(reference$q_value_threshold, 6), round(fdp$records$q_value_threshold, 6)), ]
  expect_identical(r$target_count, as.numeric(reference$target_count))
  expect_identical(r$entrapment_count, as.numeric(reference$entrapment_count))
  for (column in c("lower_bound", "combined", "paired")) {
    expect_true(max(abs(r[[column]] - reference[[column]])) < 1e-5, info = column)
  }
})

# ---------------------------------------------------------------- end to end, replayed

test_that("the three verbs replay end to end, an estimator or ratio choosing the recording", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  logistic <- stats_fit_logistic("glm_reference_responses.tsv", "glm_reference_design.tsv", c("Intercept", "x", "grp"))
  expect_identical(logistic$feature_count, 80)
  reml <- stats_fit_mixed("lmm_reference_responses.tsv", "lmm_reference_design.tsv", "lmm_reference_groups.tsv",
    c("Intercept", "age", "wp"))
  ml <- stats_fit_mixed("lmm_reference_responses.tsv", "lmm_reference_design.tsv", "lmm_reference_groups.tsv",
    c("Intercept", "age", "wp"), estimator = "ml")
  expect_identical(c(reml$estimator, ml$estimator), c("reml", "ml"))
  one <- stats_entrapment_fdp("entrapment_reference_targets.tsv", "entrapment_reference_entrapments.tsv", ratio = 1)
  two <- stats_entrapment_fdp("entrapment_reference_targets.tsv", "entrapment_reference_entrapments.tsv", ratio = 2,
    thresholds = c(0.01, 0.05))
  expect_identical(c(one$row_count, two$row_count), c(100, 2))
})

# ---------------------------------------------------------------- against a real mzLib

models_live <- function(verb) {
  skip_unless_bridge_has(verb)
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
}

test_that("LIVE: the model recordings are what the bridge returns today", {
  models_live("stats fit-mixed")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  lmm <- function(estimator = NULL) {
    stats_fit_mixed(models_path("lmm_reference_responses.tsv"), models_path("lmm_reference_design.tsv"),
      models_path("lmm_reference_groups.tsv"), c("Intercept", "age", "wp"), estimator = estimator)
  }
  calls <- list(
    stats_fit_logistic_glm.json = function() {
      stats_fit_logistic(models_path("glm_reference_responses.tsv"), models_path("glm_reference_design.tsv"), c("Intercept", "x", "grp"))
    },
    stats_fit_logistic_malat.json = function() {
      stats_fit_logistic(models_path("malat_dilution_detected.tsv"), models_path("malat_dilution_design.tsv"), c("malat_250ng", "malat_125ng"))
    },
    stats_fit_mixed_lmm_reml.json = function() lmm(),
    stats_fit_mixed_lmm_ml.json = function() lmm("ml"),
    stats_fit_mixed_malat.json = function() {
      stats_fit_mixed(models_path("malat_dilution_log2_vs_fluc.tsv"), models_path("malat_dilution_design.tsv"),
        models_path("malat_dilution_day.tsv"), c("malat_250ng", "malat_125ng"))
    },
    stats_entrapment_fdp_reference.json = function() {
      stats_entrapment_fdp(models_path("entrapment_reference_targets.tsv"), models_path("entrapment_reference_entrapments.tsv"), ratio = 1)
    },
    stats_entrapment_fdp_reference_r2.json = function() {
      stats_entrapment_fdp(models_path("entrapment_reference_targets.tsv"), models_path("entrapment_reference_entrapments.tsv"),
        ratio = 2, thresholds = c(0.01, 0.05))
    }
  )
  parsers <- list(stats_fit_logistic = mz$stats_parse_logistic, stats_fit_mixed = mz$stats_parse_mixed,
    stats_entrapment = mz$stats_parse_entrapment)
  for (name in names(calls)) {
    live <- calls[[name]]()
    parse <- parsers[[which(startsWith(name, names(parsers)))]]
    recorded <- parse(recorded_data(name))
    # MixedModel finds its variance ratio by numerical optimisation, whose stopping point differs in
    # the last digits between platforms: a Windows recording replayed on Linux moved estimates by up
    # to 5e-9 relative. 1e-6 still catches any change in what the model computes.
    relative <- if (startsWith(name, "stats_fit_mixed")) 1e-6 else 1e-9
    expect_identical(records_disagree(live$records, recorded$records, relative), "", info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
    if (!is.null(recorded$status_counts)) {
      # Its keys are the bridge's reading of mzLib's code (mzLib#1487); a change must not slip past.
      expect_identical(live$status_counts, recorded$status_counts, info = name)
    }
  }
})

test_that("LIVE: refusals cross as usage errors", {
  models_live("stats fit-logistic")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  dir <- tempfile("models")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  write <- function(name, lines) {
    path <- file.path(dir, name)
    writeLines(lines, path)
    path
  }
  proportions <- write("proportions.tsv", c("feature\tV1\tV2", "A\t0.5\t1"))
  design <- write("design.tsv", c("sample\tintercept", "V1\t1", "V2\t1"))
  expect_error(stats_fit_logistic(proportions, design, "intercept"), class = "mzlib_usage_error",
    contains = c("expected 0, 1 or missing", "feature 'A', sample 'V1'"))
  expect_error(stats_fit_logistic(models_path("glm_reference_responses.tsv"), models_path("glm_reference_design.tsv"), "x",
    max_iterations = 0), class = "mzlib_usage_error", contains = "refuses --max-iterations 0")
  groups <- write("groups.tsv", "sample\tday\tsite")
  expect_error(stats_fit_mixed(models_path("lmm_reference_responses.tsv"), models_path("lmm_reference_design.tsv"), groups, "age"),
    class = "mzlib_usage_error", contains = "needs exactly two")
  expect_error(stats_entrapment_fdp(models_path("entrapment_reference_targets.tsv"),
    models_path("entrapment_reference_entrapments.tsv"), ratio = 0), class = "mzlib_usage_error",
    contains = "ratio must be a positive finite number")
  extra <- write("targets.tsv", c("peptide\tscore\tq_value", "PEPTIDE\t10\t0.01"))
  expect_error(stats_entrapment_fdp(extra, models_path("entrapment_reference_entrapments.tsv"), ratio = 1),
    class = "mzlib_usage_error", contains = "Not a column this file takes: peptide")
})
