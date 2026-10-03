# Differential abundance: stats_fit(), stats_adjust(), stats_meta().
#
# Ported from pyMzLib's test_stats.py and test_stats_live.py. The payloads were recorded from the
# real bridge by pyMzLib and are shared verbatim. The reference tables in inst/extdata/stats/ are
# pyMzLib's, verbatim: mzLib's own limma 3.68.5 and metafor 5.2.1 reference data and output (see
# limma_metafor_PROVENANCE.txt), and the MALAT1 dilution series derived from mzLib's RNA test table.
#
# Three layers of agreement, each to 1e-8 relative:
#   1. always: the recorded bridge output against the recorded limma / metafor output, and BH
#      against stats::p.adjust(), which is base R;
#   2. when limma or metafor is installed: the recorded bridge output against limma or metafor run
#      HERE, on the same inputs - the check an R user would make;
#   3. LIVE, with a 1.0.593 bridge: a fresh fit against both.

stats_fixture <- function(name) {
  parsed <- mz$json_parse(paste(readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}
# The reference tables ship with the package, so the examples can read them too.
stats_reference_path <- function(name) {
  path <- system.file("extdata", "stats", name, package = "mzLibR")
  if (!nzchar(path)) stop("reference table not installed: ", name)
  path
}
stats_reference <- function(name) {
  utils::read.delim(stats_reference_path(name), check.names = FALSE, stringsAsFactors = FALSE)
}

# The largest relative difference, |a - b| / |b|, over paired values; both must be finite.
stats_worst <- function(ours, theirs) {
  stopifnot(length(ours) == length(theirs), all(is.finite(ours)), all(is.finite(theirs)))
  max(abs(ours - theirs) / pmax(abs(theirs), .Machine$double.xmin))
}

stats_limma_fit <- function() mz$stats_parse_fit(stats_fixture("stats_fit_limma.json"))
stats_malat_fit <- function() mz$stats_parse_fit(stats_fixture("stats_fit_malat.json"))
stats_metafor <- function() mz$stats_parse_meta(stats_fixture("stats_meta_metafor.json"))

# What limma and metafor say, computed here when they are installed.
limma_here <- function() {
  skip_if(!requireNamespace("limma", quietly = TRUE), "limma is not installed")
  responses <- stats_reference("limma_reference_responses.tsv")
  design <- stats_reference("limma_reference_design.tsv")
  y <- as.matrix(responses[, -1L])
  rownames(y) <- responses[[1L]]
  x <- as.matrix(design[match(colnames(y), design$sample), -1L])
  fit <- limma::eBayes(limma::lmFit(y, x), legacy = TRUE)
  list(
    t = unname(fit$t[, "age_decades"]),
    p = unname(fit$p.value[, "age_decades"]),
    bh = stats::p.adjust(unname(fit$p.value[, "age_decades"]), method = "BH"),
    df_total = unname(fit$df.total),
    s2_post = unname(fit$s2.post),
    s2_prior = rep(fit$s2.prior, length.out = nrow(y)),
    estimate = unname(fit$coefficients[, "age_decades"]),
    df_prior = fit$df.prior
  )
}

metafor_here <- function() {
  skip_if(!requireNamespace("metafor", quietly = TRUE), "metafor is not installed")
  inputs <- stats_reference("metafor_dl_inputs.tsv")
  cases <- unique(inputs$case)
  rows <- lapply(cases, function(case) {
    one <- inputs[inputs$case == case, ]
    m <- metafor::rma(yi = one$yi, sei = one$sei, method = "DL")
    c(estimate = as.numeric(m$b), se = m$se, tau2 = m$tau2, q = m$QE, i2 = m$I2 / 100,
      ci_lb = m$ci.lb, ci_ub = m$ci.ub, pval = m$pval)
  })
  out <- as.data.frame(do.call(rbind, rows))
  out$case <- cases
  out
}

# ---------------------------------------------------------------- stats_fit: the payload

test_that("the statuses are mzLib's three", {
  expect_identical(mz$STATS_FIT_STATUSES, c("fitted", "too_few_observations", "rank_deficient"))
})

test_that("the limma reference fit projects every feature, one block per coefficient", {
  fit <- stats_limma_fit()
  expect_true(inherits(fit, "mzlibr_moderated_fit"))
  expect_identical(c(fit$feature_count, fit$sample_count, fit$row_count), c(400, 12, 400))
  expect_identical(fit$tested_coefficients, "age_decades")
  expect_identical(fit$coefficient_names, c("Intercept", "age_decades", "sex"))
  expect_true(fit$residual_df_differ)
  expect_identical(nrow(fit$records), 400L)
  expect_identical(names(fit$records), fit$column_names)
  expect_identical(unique(fit$records$coefficient), "age_decades")
  expect_true(is.numeric(fit$status_counts) && !is.null(names(fit$status_counts)))
  expect_true(is.list(fit$prior) && all(c("df", "df_infinite", "trended", "spline_basis_count", "scale") %in% names(fit$prior)))
})

test_that("the recorded fit agrees with limma's recorded eBayes(legacy = TRUE) to 1e-8", {
  fit <- stats_limma_fit()
  limma <- stats_reference("limma_ebayes_notrend.tsv")
  expect_identical(nrow(limma), nrow(fit$records))
  expect_true(stats_worst(fit$records$t, limma$t_age) < 1e-8)
  expect_true(stats_worst(fit$records$p_value, limma$p_age) < 1e-8)
  expect_true(stats_worst(fit$records$bh_adjusted, limma$bh_age) < 1e-8)
  expect_true(stats_worst(fit$records$df_total, limma$df_total) < 1e-8)
  expect_true(stats_worst(fit$records$posterior_variance, limma$s2_post) < 1e-8)
  expect_true(stats_worst(fit$records$prior_variance, limma$s2_prior) < 1e-8)
})

test_that("bh_adjusted is stats::p.adjust(method = \"BH\") of p_value, to 1e-8", {
  fit <- stats_limma_fit()
  expect_true(stats_worst(fit$records$bh_adjusted, stats::p.adjust(fit$records$p_value, method = "BH")) < 1e-8)
})

test_that("the recorded fit agrees with limma run here, to 1e-8", {
  here <- limma_here()
  fit <- stats_limma_fit()
  expect_true(stats_worst(fit$records$t, here$t) < 1e-8)
  expect_true(stats_worst(fit$records$p_value, here$p) < 1e-8)
  expect_true(stats_worst(fit$records$bh_adjusted, here$bh) < 1e-8)
  expect_true(stats_worst(fit$records$df_total, here$df_total) < 1e-8)
  expect_true(stats_worst(fit$records$posterior_variance, here$s2_post) < 1e-8)
  expect_true(stats_worst(fit$records$prior_variance, here$s2_prior) < 1e-8)
  expect_true(stats_worst(fit$records$estimate, here$estimate) < 1e-8)
  expect_true(stats_worst(fit$prior$df, here$df_prior) < 1e-8)
})

test_that("the MALAT1 dilution reads as about -1 per halving, unfitted features kept", {
  fit <- stats_malat_fit()
  expect_identical(c(fit$feature_count, fit$sample_count), c(492, 27))
  expect_identical(fit$status_counts[["fitted"]], 212)
  expect_identical(fit$status_counts[["too_few_observations"]], 214)
  expect_identical(fit$status_counts[["rank_deficient"]], 66)
  expect_identical(fit$row_count, 2 * 492)
  expect_true(fit$residual_df_differ)
  expect_identical(round(fit$prior$df, 1), 19.6)
  rows <- fit$records[fit$records$coefficient == "malat_250ng" & fit$records$status == "fitted" &
    startsWith(fit$records$feature, "MALAT1:"), ]
  expect_identical(nrow(rows), 116L)
  expect_identical(round(stats::median(rows$estimate), 2), -0.89)
  unfitted <- fit$records[fit$records$status != "fitted", ]
  expect_true(all(is.na(unfitted$t) & is.na(unfitted$p_value) & is.na(unfitted$bh_adjusted)))
})

# ---------------------------------------------------------------- stats_fit: arguments

test_that("the fit request names the files, threads and coefficients on stdin", {
  request <- mz$stats_build_fit_request("r.tsv", "d.tsv", c("a", "b"), FALSE, NULL, 1)
  expect_identical(request$args, c("stats", "fit", "--responses", "r.tsv", "--design", "d.tsv", "--threads", "1"))
  expect_identical(request$stdin, "a\nb")
  trended <- mz$stats_build_fit_request("r.tsv", "d.tsv", "a", TRUE, 3, -1)
  expect_identical(utils::tail(trended$args, 5L), c("--threads", "-1", "--trend", "--spline-basis", "3"))
})

test_that("bad fit arguments are refused before the bridge", {
  build <- function(...) {
    given <- list(responses = "r.tsv", design = "d.tsv", coefficients = "a", trend = FALSE, spline_basis = NULL, threads = 1)
    given[names(list(...))] <- list(...)
    do.call(mz$stats_build_fit_request, given)
  }
  expect_error(build(coefficients = character(0)), class = "mzlib_usage_error", contains = "coefficients must be")
  expect_error(build(coefficients = 1), class = "mzlib_usage_error", contains = "coefficients must be")
  expect_error(build(coefficients = c("a", "")), class = "mzlib_usage_error", contains = "non-empty")
  expect_error(build(coefficients = "a\nb"), class = "mzlib_usage_error", contains = "newline")
  expect_error(build(trend = "yes"), class = "mzlib_usage_error", contains = "TRUE or FALSE")
  expect_error(build(threads = 1.5), class = "mzlib_usage_error", contains = "threads must be a whole number")
  expect_error(build(spline_basis = 2.5), class = "mzlib_usage_error", contains = "spline_basis must be")
  expect_error(build(responses = ""), class = "mzlib_usage_error")
})

# ---------------------------------------------------------------- stats_adjust

test_that("adjust keeps every position, and NA is not counted in m", {
  adjusted <- mz$stats_parse_adjusted(stats_fixture("stats_adjust.json"))
  expect_true(inherits(adjusted, "mzlibr_adjusted"))
  expect_identical(c(adjusted$line_count, adjusted$tested_count), c(8, 6))
  expect_identical(round(adjusted$records$bh_adjusted, 4), c(0.0012, 0.012, 0.038, NA, 0.0465, 0.24, NA, 0.74))
  tested <- !is.na(adjusted$records$p_value)
  expect_true(stats_worst(adjusted$records$bh_adjusted[tested],
                          stats::p.adjust(adjusted$records$p_value[tested], method = "BH")) < 1e-8)
})

test_that("NA and NaN cross as empty lines, and numbers survive the trip exactly", {
  expect_identical(mz$stats_adjust_stdin(c(0.0002, NA, NaN, 0.74)),
                   "0.00020000000000000001\n\n\n0.73999999999999999")
  # 17 significant digits. R's own parser is not correctly rounded everywhere (macOS arm64 misreads
  # 1e-300), so it is the check here only at ordinary magnitudes; .NET reads every one exactly.
  for (x in c(1 / 3, 0.1 + 0.2, 4.9534166563722163e-08)) {
    expect_identical(as.numeric(strsplit(mz$stats_adjust_stdin(x), "\n")[[1L]]), x)
  }
  expect_identical(mz$stats_adjust_stdin(NA), "")
  expect_error(stats_adjust(character(0)), class = "mzlib_usage_error", contains = "numeric vector")
  expect_error(stats_adjust(numeric(0)), class = "mzlib_usage_error", contains = "at least one")
  expect_error(stats_adjust("0.01"), class = "mzlib_usage_error", contains = "numeric vector")
})

# ---------------------------------------------------------------- stats_meta

test_that("the recorded pooling agrees with metafor's recorded DerSimonian-Laird to 1e-8", {
  pooled <- stats_metafor()
  expect_true(inherits(pooled, "mzlibr_meta_analysis"))
  reference <- stats_reference("metafor_dl_results.tsv")
  r <- pooled$records[match(reference$case, pooled$records$feature), ]
  expect_identical(r$feature, reference$case)
  expect_true(stats_worst(r$estimate, reference$estimate) < 1e-8)
  expect_true(stats_worst(r$standard_error, reference$se) < 1e-8)
  expect_true(stats_worst(r$confidence_low, reference$ci_lb) < 1e-8)
  expect_true(stats_worst(r$confidence_high, reference$ci_ub) < 1e-8)
  expect_true(stats_worst(r$p_value, reference$pval) < 1e-8)
  expect_true(stats_worst(r$q, reference$q) < 1e-8)
  positive <- reference$tau2 > 0
  expect_true(stats_worst(r$tau2[positive], reference$tau2[positive]) < 1e-8)
  expect_identical(r$tau2[!positive], reference$tau2[!positive])
})

test_that("the recorded pooling agrees with metafor run here, to 1e-8", {
  here <- metafor_here()
  r <- stats_metafor()$records
  r <- r[match(here$case, r$feature), ]
  for (pair in list(c("estimate", "estimate"), c("standard_error", "se"), c("confidence_low", "ci_lb"),
                    c("confidence_high", "ci_ub"), c("p_value", "pval"), c("q", "q"))) {
    expect_true(stats_worst(r[[pair[1L]]], here[[pair[2L]]]) < 1e-8, info = pair[1L])
  }
  positive <- here$tau2 > 0
  expect_true(stats_worst(r$tau2[positive], here$tau2[positive]) < 1e-8)
  expect_true(stats_worst(r$i_squared[positive], here$i2[positive]) < 1e-8)
})

test_that("the pooled rows report direction and the cases' studies", {
  r <- stats_metafor()$records
  expect_identical(r$feature, c("homogeneous", "heterogeneous", "two_studies", "random_eight"))
  expect_identical(r$studies, c(5, 6, 2, 8))
  expect_identical(r$direction_agree, c(5, 5, 1, 6))
  expect_identical(r$direction_disagree, c(0, 1, 1, 2))
})

test_that("meta studies travel as tab-separated lines, exactly", {
  studies <- data.frame(feature = c("a", "a"), estimate = c(0.3, -1 / 3), standard_error = c(0.1, 0.2))
  lines <- strsplit(mz$stats_meta_stdin(studies), "\n", fixed = TRUE)[[1L]]
  expect_identical(length(lines), 2L)
  fields <- strsplit(lines[[2L]], "\t", fixed = TRUE)[[1L]]
  expect_identical(fields[[1L]], "a")
  expect_identical(as.numeric(fields[[2L]]), -1 / 3)
})

test_that("bad meta arguments are refused before the bridge", {
  ok <- data.frame(feature = "a", estimate = 0.1, standard_error = 0.1)
  expect_error(stats_meta(list(feature = "a")), class = "mzlib_usage_error", contains = "data.frame")
  expect_error(stats_meta(ok[, 1:2]), class = "mzlib_usage_error", contains = "standard_error")
  expect_error(stats_meta(ok[0, ]), class = "mzlib_usage_error", contains = "at least one")
  expect_error(stats_meta(transform(ok, feature = "")), class = "mzlib_usage_error", contains = "non-empty feature")
  expect_error(stats_meta(transform(ok, feature = "a\tb")), class = "mzlib_usage_error", contains = "tab or newline")
  expect_error(stats_meta(transform(ok, estimate = NA_real_)), class = "mzlib_usage_error", contains = "no NA")
  expect_error(stats_meta(ok, confidence = "high"), class = "mzlib_usage_error", contains = "confidence")
})

test_that("all three replay end to end, verb check included, and the fit picks its recording by file", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(stats_fit("limma_reference_responses.tsv", "limma_reference_design.tsv", "age_decades")$feature_count, 400)
  expect_identical(stats_fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
                             c("malat_250ng", "malat_125ng"))$feature_count, 492)
  expect_error(stats_fit("other.tsv", "d.tsv", "a"), class = "mzlib_usage_error", contains = "replay bridge")
  expect_identical(stats_adjust(c(0.0002, 0.004, 0.019, NA, 0.031, 0.2, NaN, 0.74))$tested_count, 6)
  inputs <- stats_reference("metafor_dl_inputs.tsv")
  expect_identical(stats_meta(data.frame(feature = inputs$case, estimate = inputs$yi, standard_error = inputs$sei))$feature_count, 4)
})

# ---------------------------------------------------------------- against a real mzLib

test_that("LIVE: a fresh fit on limma's reference data agrees with limma's recorded output", {
  skip_unless_bridge_has("stats fit")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  fit <- stats_fit(stats_reference_path("limma_reference_responses.tsv"),
                   stats_reference_path("limma_reference_design.tsv"), "age_decades", threads = -1)
  limma <- stats_reference("limma_ebayes_notrend.tsv")
  expect_true(stats_worst(fit$records$t, limma$t_age) < 1e-8)
  expect_true(stats_worst(fit$records$bh_adjusted, limma$bh_age) < 1e-8)
  expect_identical(fit$records$t, stats_limma_fit()$records$t, info = "the same at any thread count")
})

test_that("LIVE: a fresh fit agrees with limma run here", {
  skip_unless_bridge_has("stats fit")
  here <- limma_here()
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  fit <- stats_fit(stats_reference_path("limma_reference_responses.tsv"),
                   stats_reference_path("limma_reference_design.tsv"), "age_decades")
  expect_true(stats_worst(fit$records$t, here$t) < 1e-8)
  expect_true(stats_worst(fit$records$p_value, here$p) < 1e-8)
})

test_that("LIVE: adjust and meta agree with p.adjust and metafor's recorded output", {
  skip_unless_bridge_has("stats meta")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  p <- c(0.0002, 0.004, 0.019, NA, 0.031, 0.2, NaN, 0.74)
  adjusted <- stats_adjust(p)
  tested <- !is.na(p)
  expect_true(stats_worst(adjusted$records$bh_adjusted[tested], stats::p.adjust(p[tested], method = "BH")) < 1e-8)
  inputs <- stats_reference("metafor_dl_inputs.tsv")
  pooled <- stats_meta(data.frame(feature = inputs$case, estimate = inputs$yi, standard_error = inputs$sei))
  reference <- stats_reference("metafor_dl_results.tsv")
  r <- pooled$records[match(reference$case, pooled$records$feature), ]
  expect_true(stats_worst(r$estimate, reference$estimate) < 1e-8)
  expect_true(stats_worst(r$standard_error, reference$se) < 1e-8)
  expect_error(stats_adjust(1.5), class = "mzlib_usage_error")
  expect_error(stats_meta(data.frame(feature = "a", estimate = 0.1, standard_error = 0)), class = "mzlib_usage_error")
})
