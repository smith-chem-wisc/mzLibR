# stats_fit_latent_mixed() and stats_fit_link(): mzLib's latent-scale mixed model and its link.
#
# Ported from pyMzLib's test_stats_latent.py and test_stats_latent_live.py. The payloads are
# pyMzLib's recordings, verbatim, made on mzLib's own test data, which ships in extdata/stats/ with
# lme4's glmer() and glmmTMB's fits of it (lnmm_PROVENANCE.txt). The numbers are mzLib's
# LogitNormalMixedModel; these hold the R layer - what it sends (and nothing for an option left out,
# so mzLib's default applies), what it parses, what it refuses before starting a process - and hold
# the recordings to glmer and glmmTMB, here and LIVE.

latent_path <- function(name) system.file("extdata", "stats", name, package = "mzLibR")

latent_request <- function(link, ...) {
  arguments <- list(design = "d.tsv", groups = "g.tsv", coefficients = "trait", link = link, responses = NULL,
    counts = "c.tsv", depth = NULL, residual_variance = NULL, estimator = NULL, threads = 1)
  arguments[names(list(...))] <- list(...)
  do.call(mz$stats_build_latent_request, arguments)
}
latent_link_sent <- function(link) {
  args <- latent_request(link)$args
  args[[match("--link", args) + 1L]]
}

# ---------------------------------------------------------------- fit_latent_mixed: what is sent

test_that("fit_latent_mixed sends only what was given", {
  seen <- recorded_bridge_call("stats_fit_latent_mixed_scatter.json", function() {
    stats_fit_latent_mixed("d.tsv", "g.tsv", "trait", link = "identity", counts = "c.tsv")
  })
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("stats", "fit-latent-mixed", "--design", "d.tsv", "--groups", "g.tsv",
    "--link", "identity", "--threads", "1", "--counts", "c.tsv"))
  expect_identical(seen$stdin, "trait")
  expect_identical(seen$timeout, 600)
  every <- latent_request(c(-0.4, 0.6, 0.3), coefficients = c("trait", "x"), responses = "r.tsv", depth = "depth.tsv",
    residual_variance = "s.tsv", estimator = "ml", threads = -1)
  expect_identical(every$args[7:length(every$args)], c(
    "--link", "-0.40000000000000002,0.59999999999999998,0.29999999999999999", "--threads", "-1",
    "--responses", "r.tsv", "--counts", "c.tsv", "--depth", "depth.tsv", "--residual-variance", "s.tsv",
    "--estimator", "ml"))
  expect_identical(every$stdin, "trait\nx")
})

test_that("a link sends only the numbers it holds", {
  expect_identical(latent_link_sent("IDENTITY"), "identity")
  # Two numbers send two: the depth slope and overdispersion are mzLib's defaults.
  expect_identical(latent_link_sent(c(0.5, 2)), "0.5,2")
  expect_identical(latent_link_sent(c(1L, 2L, 0.5, 0.25)), "1,2,0.5,0.25")
  fitted <- mz$stats_parse_link_fit(recorded_data("stats_fit_link_simulated.json"))
  expected <- paste(mz$stats_number_text(unlist(unclass(fitted$link))), collapse = ",")
  expect_identical(latent_link_sent(fitted), expected)
  expect_identical(latent_link_sent(fitted$link), expected)
  expect_identical(length(strsplit(expected, ",", fixed = TRUE)[[1L]]), 4L)
})

test_that("fit_latent_mixed refuses bad input without starting a process", {
  never <- function(args, stdin) stop("the bridge must not be called")
  refused <- function(contains, ...) {
    arguments <- utils::modifyList(list(design = "d.tsv", groups = "g.tsv", coefficients = "trait", counts = "c.tsv"),
      list(...))
    expect_error(recorded_bridge_call(never, function() do.call(stats_fit_latent_mixed, arguments)),
      class = "mzlib_usage_error", contains = contains)
  }
  refused("link takes \"identity\"", link = "logit")
  refused("2 to 4 numbers", link = 1)
  refused("2 to 4 numbers", link = c(1, 2, 3, 4, 5))
  refused("2 to 4 numbers", link = list(0, 1))
  refused("link[2] must be a number", link = c(0, NA))
  refused("link is required")
  refused("estimator must be one of reml, ml", link = "identity", estimator = "REML")
  refused("threads must be a whole number", link = "identity", threads = 1.5)
  refused("depth must be a single non-empty file path", link = "identity", depth = 3)
  refused("needs at least one coefficient", link = "identity", coefficients = character(0))
  unconverged <- mz$stats_parse_link_fit(utils::modifyList(recorded_data("stats_fit_link_simulated.json"), list(link = NA)))
  refused("did not converge", link = unconverged)
  expect_error(recorded_bridge_call(never, function() stats_fit_latent_mixed("d.tsv", "g.tsv", "trait", link = "identity")),
    class = "mzlib_usage_error", contains = "needs responses = (measured samples), counts = (counted samples), or both")
})

test_that("an older bridge is named, not spawned, for the latent-scale verbs", {
  for (call in list(
    function() stats_fit_latent_mixed("d.tsv", "g.tsv", "trait", link = "identity", counts = "c.tsv"),
    function() stats_fit_link("cells.tsv")
  )) {
    expect_error(recorded_bridge_call(function(args, stdin) stop("reached the verb"), call, version = "readers_formats.json"),
      class = "mzlib_usage_error", contains = "pyMzLib 0.7.0")
  }
})

# ---------------------------------------------------------------- fit_latent_mixed: what is parsed

test_that("fit_latent_mixed parses the recordings", {
  fit <- mz$stats_parse_latent(recorded_data("stats_fit_latent_mixed_simulated.json"))
  expect_true(inherits(fit, "mzlibr_latent_mixed_model_fit"))
  expect_identical(fit$estimator, "reml", info = "absent, mzLib's default, as mzLib reports it")
  expect_true(inherits(fit$link, "mzlibr_detection_link"))
  expect_identical(unclass(fit$link), list(intercept = 0.4, slope = 0.6, depth_slope = 0, overdispersion = 0))
  expect_identical(c(fit$grouping, fit$group_count), c("animal", "16"))
  expect_identical(c(fit$measured_cells, fit$counted_cells, fit$missing_count), c(96, 96, 0))
  expect_true(is.na(fit$depth_file) && is.na(fit$residual_variance_file))
  expect_false(is.na(fit$responses_file) || is.na(fit$counts_file))
  expect_identical(names(fit$status_counts), mz$STATS_LATENT_MIXED_STATUSES)
  expect_identical(fit$column_names, c("feature", "coefficient", "status", "estimate", "standard_error", "df", "t",
    "p_value", "residual_variance", "residual_variance_given", "group_variance", "log_likelihood", "measured",
    "counted", "groups"))
  expect_identical(names(fit$records), fit$column_names)
  trait <- fit$records[fit$records$coefficient == "trait", ]
  expect_identical(unique(trait$measured), 32)
  expect_identical(unique(trait$counted), 32)
  expect_identical(unique(trait$df), 14)
  expect_identical(unique(trait$residual_variance_given), FALSE)
  # mzLib's own test asks for 3 standard errors of the true 1.0; every one is within 2.
  expect_true(all(abs(trait$estimate - 1) < 2 * trait$standard_error))
  expect_identical(fit$records$df[fit$records$coefficient == "Intercept"], c(48, 48, 48))

  ml <- mz$stats_parse_latent(recorded_data("stats_fit_latent_mixed_link_ml.json"))
  expect_identical(ml$estimator, "ml")
  expect_true(is.na(ml$responses_file))
  expect_identical(ml$measured_cells, 0)
  expect_identical(unclass(ml$link), list(intercept = -0.4, slope = 0.6, depth_slope = 0.3, overdispersion = 0))
  expect_true(all(ml$records$residual_variance_given))
  expect_true(any(startsWith(ml$caveats, "16 of 16 features were fitted at the residual_variance you gave")))
})

test_that("the recordings are lme4's glmer and glmmTMB's fits, as mzLib holds them", {
  # glmer (25-point adaptive quadrature, ML) on mzLib's counted reference input, each feature given
  # mzLib's tiny sigma^2: every estimate within 1e-5 of glmer's standard error, every standard error
  # within 5e-5 relative.
  ml <- mz$stats_parse_latent(recorded_data("stats_fit_latent_mixed_link_ml.json"))
  glmer <- utils::read.delim(latent_path("lnmm_counted_fit.tsv"), stringsAsFactors = FALSE)
  glmer <- glmer[glmer$link == "link" & glmer$method == "ML", ]
  for (j in 0:2) {
    ours <- ml$records[ml$records$coefficient == c("Intercept", "trait", "x")[[j + 1L]], ]
    fitted <- glmer$status == "fitted"
    se <- glmer[[paste0("se", j)]][fitted]
    expect_true(all(abs(ours$estimate[fitted] - glmer[[paste0("b", j)]][fitted]) < 1e-5 * se), info = paste("b", j))
    expect_true(all(abs(ours$standard_error[fitted] - se) / se < 5e-5), info = paste("se", j))
  }
  # glmmTMB (Laplace) with run scatter estimated from counts alone: about 1% in tau^2.
  reml <- mz$stats_parse_latent(recorded_data("stats_fit_latent_mixed_scatter.json"))
  tmb <- utils::read.delim(latent_path("lnmm_scatter_fit.tsv"), stringsAsFactors = FALSE)
  tmb <- tmb[tmb$method == "REML", ]
  ours <- reml$records[reml$records$coefficient == "trait", ]
  expect_true(all(abs(ours$group_variance - tmb$tau2) < 0.03 * tmb$tau2 + 0.003))
  expect_true(all(abs(ours$standard_error - tmb$se1) / tmb$se1 < 0.03))
})

test_that("with the counted runs, the intercept gains df and the trait barely moves", {
  latent <- mz$stats_parse_latent(recorded_data("stats_fit_latent_mixed_simulated.json"))
  measured <- mz$stats_parse_mixed(recorded_data("stats_fit_mixed_lnmm_simulated.json"))
  expect_identical(measured$records$df[measured$records$coefficient == "Intercept"], c(16, 16, 16))
  a <- latent$records[latent$records$coefficient == "trait", ]
  b <- measured$records[measured$records$coefficient == "trait", ]
  expect_identical(a$feature, b$feature)
  expect_identical(sprintf("%+.2f", a$estimate), c("+1.61", "+0.73", "+1.26"))
  expect_identical(sprintf("%+.2f", b$estimate), c("+1.63", "+0.77", "+1.36"))
})

# ---------------------------------------------------------------- fit_link

test_that("fit_link sends each option only when given", {
  seen <- recorded_bridge_call("stats_fit_link_simulated.json", function() stats_fit_link("cells.tsv"))
  expect_identical(seen$args, c("stats", "fit-link", "--cells", "cells.tsv"))
  expect_identical(seen$timeout, 600)
  expect_identical(mz$stats_build_link_args("cells.tsv", TRUE, FALSE, FALSE)[5:10],
    c("--with-depth", "true", "--with-overdispersion", "false", "--detection-only", "false"))
})

test_that("fit_link refuses an option that is not TRUE, FALSE or NULL without starting a process", {
  never <- function(args, stdin) stop("the bridge must not be called")
  for (options in list(list(with_depth = 1), list(detection_only = "false"), list(with_overdispersion = NA))) {
    expect_error(recorded_bridge_call(never, function() do.call(stats_fit_link, c(list("cells.tsv"), options))),
      class = "mzlib_usage_error", contains = "must be TRUE, FALSE or NULL")
  }
  expect_error(stats_fit_link(NULL), class = "mzlib_usage_error", contains = "cells must be")
})

test_that("fit_link parses the recording: mzLib's defaults reported, the link recovered", {
  fit <- mz$stats_parse_link_fit(recorded_data("stats_fit_link_simulated.json"))
  expect_true(inherits(fit, "mzlibr_link_fit"))
  expect_identical(fit$cell_count, 4000)
  expect_identical(c(fit$depth_given, fit$with_depth, fit$with_overdispersion, fit$detection_only),
    c(TRUE, TRUE, FALSE, TRUE), info = "mzLib's defaults, reported")
  expect_true(inherits(fit$link, "mzlibr_detection_link"))
  expect_identical(names(fit$link), c("intercept", "slope", "depth_slope", "overdispersion"))
  expect_true(abs(fit$link$intercept + 0.5) < 0.15 && abs(fit$link$slope - 0.5) < 0.1)
  expect_identical(fit$link$overdispersion, 0)
  expect_true(any(startsWith(fit$caveats, "overdispersion was not fitted")))
  unconverged <- mz$stats_parse_link_fit(utils::modifyList(recorded_data("stats_fit_link_simulated.json"), list(link = NA)))
  expect_true(is.null(unconverged$link) && "link" %in% names(unconverged))
})

# ---------------------------------------------------------------- the replay bridge, which every example runs on

test_that("the replay matches a link by the numbers given, and holds the optional files", {
  data <- recorded_data("stats_fit_latent_mixed_link_ml.json")
  files <- list(design = "lnmm_reference_design.tsv", groups = "lnmm_reference_groups.tsv",
    counts = "lnmm_link_counts.tsv", depth = "lnmm_link_depth.tsv", "residual-variance" = "lnmm_link_residual_variance.tsv")
  expect_identical(mz$replay_mismatch(data, c(files, link = "-0.4,0.6,0.3")), "")
  expect_identical(mz$replay_mismatch(data, c(files, link = "-0.40000000000000002,0.59999999999999998")), "")
  expect_true(grepl("not --link", mz$replay_mismatch(data, c(files, link = "-0.4,0.6,0.2"))))
  expect_true(grepl("not --link", mz$replay_mismatch(data, c(files, link = "identity"))))
  scatter <- recorded_data("stats_fit_latent_mixed_scatter.json")
  expect_identical(mz$replay_mismatch(scatter, list(counts = "lnmm_scatter_counts.tsv", link = "identity")), "")
  expect_identical(mz$replay_mismatch(scatter, list(counts = "lnmm_scatter_counts.tsv", link = "0,1")), "")
  expect_true(grepl("recorded without --responses", mz$replay_mismatch(scatter, list(responses = "r.tsv"))))
  expect_true(grepl("recorded with --counts", mz$replay_mismatch(scatter, list(link = "identity"))))
  expect_true(grepl("recorded without --depth",
    mz$replay_mismatch(scatter, list(counts = "lnmm_scatter_counts.tsv", depth = "d.tsv"))))
  link <- recorded_data("stats_fit_link_simulated.json")
  expect_identical(mz$replay_mismatch(link, list(cells = "link_cells.tsv", "with-depth" = "true")), "")
  expect_true(grepl("with_depth=true", mz$replay_mismatch(link, list(cells = "link_cells.tsv", "with-depth" = "false"))))
})

test_that("the latent-scale examples replay end to end", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  simulated <- stats_fit_latent_mixed(latent_path("lnmm_simulated_design.tsv"), latent_path("lnmm_simulated_groups.tsv"),
    c("Intercept", "trait"), link = c(0.4, 0.6), responses = latent_path("lnmm_simulated_responses.tsv"),
    counts = latent_path("lnmm_simulated_counts.tsv"))
  expect_identical(simulated$feature_count, 3)
  scatter <- stats_fit_latent_mixed(latent_path("lnmm_reference_design.tsv"), latent_path("lnmm_reference_groups.tsv"),
    c("Intercept", "trait", "x"), link = "identity", counts = latent_path("lnmm_scatter_counts.tsv"))
  expect_identical(scatter$feature_count, 20)
  link <- stats_fit_link(latent_path("link_cells.tsv"), with_depth = TRUE)
  expect_identical(link$cell_count, 4000)
  # Leaving with_depth to mzLib is another call, even though mzLib's default happens to be echoed.
  expect_error(stats_fit_link(latent_path("link_cells.tsv")), class = "mzlib_usage_error", contains = "replay bridge")
  expect_error(stats_fit_latent_mixed(latent_path("lnmm_reference_design.tsv"), latent_path("lnmm_reference_groups.tsv"),
    c("Intercept", "trait", "x"), link = c(0, 0.9), counts = latent_path("lnmm_scatter_counts.tsv")),
    class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

latent_live <- function() {
  skip_unless_bridge_has("stats fit-latent-mixed")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  invisible(NULL)
}

# LogitNormalMixedModel finds its variances by numerical search, so a recording's last digits may
# move between platforms: 1e-6, as pyMzLib's live tests hold them. A variance near its boundary sits
# where the likelihood is flat, so the two variance columns are compared on the scale of the
# feature's larger variance, as mzLib's own reference test compares tau^2 on max(tau^2, sigma^2).
latent_disagree <- function(live, recorded) {
  scale <- pmax(abs(ifelse(is.na(recorded$residual_variance), 0, recorded$residual_variance)),
    abs(ifelse(is.na(recorded$group_variance), 0, recorded$group_variance)))
  rest <- setdiff(names(recorded), c("residual_variance", "group_variance"))
  why <- records_disagree(live[, rest, drop = FALSE], recorded[, rest, drop = FALSE], relative = 1e-6)
  if (nzchar(why) || !identical(names(live), names(recorded))) {
    return(if (nzchar(why)) why else "columns")
  }
  for (name in c("residual_variance", "group_variance")) {
    a <- live[[name]]
    b <- recorded[[name]]
    known <- !is.na(b)
    if (!identical(is.na(a), is.na(b)) || any(abs(a[known] - b[known]) > 1e-6 * pmax(scale[known], 1e-300))) {
      return(name)
    }
  }
  ""
}

test_that("LIVE: the latent-mixed recordings are what the bridge returns today", {
  latent_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  design <- latent_path("lnmm_reference_design.tsv")
  groups <- latent_path("lnmm_reference_groups.tsv")
  calls <- list(
    stats_fit_latent_mixed_scatter.json = function() {
      stats_fit_latent_mixed(design, groups, c("Intercept", "trait", "x"), link = "identity",
        counts = latent_path("lnmm_scatter_counts.tsv"), threads = -1)
    },
    stats_fit_latent_mixed_link_ml.json = function() {
      stats_fit_latent_mixed(design, groups, c("Intercept", "trait", "x"), link = c(-0.4, 0.6, 0.3),
        counts = latent_path("lnmm_link_counts.tsv"), depth = latent_path("lnmm_link_depth.tsv"),
        residual_variance = latent_path("lnmm_link_residual_variance.tsv"), estimator = "ml", threads = -1)
    },
    stats_fit_latent_mixed_simulated.json = function() {
      stats_fit_latent_mixed(latent_path("lnmm_simulated_design.tsv"), latent_path("lnmm_simulated_groups.tsv"),
        c("Intercept", "trait"), link = c(0.4, 0.6), responses = latent_path("lnmm_simulated_responses.tsv"),
        counts = latent_path("lnmm_simulated_counts.tsv"), threads = -1)
    }
  )
  for (name in names(calls)) {
    live <- calls[[name]]()
    recorded <- mz$stats_parse_latent(recorded_data(name))
    expect_identical(latent_disagree(live$records, recorded$records), "", info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
    # Its keys are the bridge's reading of mzLib's code (mzLib#1487); a change must not slip past.
    expect_identical(live$status_counts, recorded$status_counts, info = name)
    expect_identical(live$link, recorded$link, info = name)
    expect_identical(c(live$measured_cells, live$counted_cells, live$missing_count),
      c(recorded$measured_cells, recorded$counted_cells, recorded$missing_count), info = name)
  }
})

test_that("LIVE: the link recording is what the bridge returns today", {
  latent_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  live <- stats_fit_link(latent_path("link_cells.tsv"), with_depth = TRUE)
  recorded <- mz$stats_parse_link_fit(recorded_data("stats_fit_link_simulated.json"))
  # FitLink's estimate is Nelder-Mead's stopping point, with no polishing step: 1e-4, as pyMzLib holds it.
  for (name in names(recorded$link)) {
    expect_true(abs(live$link[[name]] - recorded$link[[name]]) <= 1e-4 * max(abs(recorded$link[[name]]), 1e-300),
      info = name)
  }
  expect_identical(live$caveats, recorded$caveats)
  expect_identical(c(live$with_depth, live$with_overdispersion, live$detection_only),
    c(recorded$with_depth, recorded$with_overdispersion, recorded$detection_only))
})

test_that("LIVE: a fitted link goes straight into fit_latent_mixed", {
  latent_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  # Without a depth term: these samples have no depth table, and mzLib refuses a depth slope without one.
  link <- stats_fit_link(latent_path("link_cells.tsv"), with_depth = FALSE)
  expect_identical(link$link$depth_slope, 0)
  fit <- stats_fit_latent_mixed(latent_path("lnmm_simulated_design.tsv"), latent_path("lnmm_simulated_groups.tsv"),
    "trait", link = link, responses = latent_path("lnmm_simulated_responses.tsv"),
    counts = latent_path("lnmm_simulated_counts.tsv"))
  # mzLib echoes the link it used: the four numbers sent. Compared to 1e-12, not exactly: R reads
  # 17-digit wire text through long double, which is only double on Apple silicon.
  expect_true(isTRUE(all.equal(unclass(fit$link), unclass(link$link), tolerance = 1e-12)))
})

test_that("LIVE: refusals cross as usage errors", {
  latent_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  design <- latent_path("lnmm_simulated_design.tsv")
  groups <- latent_path("lnmm_simulated_groups.tsv")
  counts <- tempfile(fileext = ".tsv")
  writeLines(c("feature\ts1", "A\t6/5"), counts)
  expect_error(stats_fit_latent_mixed(design, groups, "trait", link = "identity",
    responses = latent_path("lnmm_simulated_responses.tsv"), counts = counts),
    class = "mzlib_usage_error", contains = "must list exactly the responses' features")
  expect_error(stats_fit_latent_mixed(design, groups, "trait", link = c(0, 0), counts = latent_path("lnmm_simulated_counts.tsv")),
    class = "mzlib_usage_error", contains = "non-zero slope")
  expect_error(stats_fit_latent_mixed(design, groups, "trait", link = c(0, 1, 0.2), counts = latent_path("lnmm_simulated_counts.tsv")),
    class = "mzlib_usage_error", contains = "supply one depth per sample")
  cells <- tempfile(fileext = ".tsv")
  writeLines(c("latent_mean\tlatent_sd\tsuccesses\ttrials\tfeature", "0\t1\t1\t2\tA"), cells)
  expect_error(stats_fit_link(cells), class = "mzlib_usage_error", contains = "Not a column this file takes: feature")
})
