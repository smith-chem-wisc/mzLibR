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
#   3. LIVE, with a bridge that has the stats verbs: a fresh fit against both, and against the
#      recordings to 1e-9.

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
  # A trend with no spline_basis asks for none: mzLib then sizes the spline as limma does
  # (pyMzLib #82). Nothing in mzLibR may fill in a number of its own.
  expect_true(is.null(formals(stats_fit)$spline_basis))
  unsized <- mz$stats_build_fit_request("r.tsv", "d.tsv", "a", TRUE, NULL, 1)
  expect_identical(utils::tail(unsized$args, 1L), "--trend")
  expect_false("--spline-basis" %in% unsized$args)
})

test_that("bad fit arguments are refused before the bridge", {
  build <- function(...) {
    given <- list(responses = "r.tsv", design = "d.tsv", coefficients = "a", trend = FALSE, spline_basis = NULL, threads = 1)
    given[names(list(...))] <- list(...)
    do.call(mz$stats_build_fit_request, given)
  }
  expect_error(build(coefficients = character(0)), class = "mzlib_usage_error", contains = "or a contrasts file")
  expect_error(build(coefficients = 1), class = "mzlib_usage_error", contains = "coefficients must be")
  expect_error(build(coefficients = c("a", "")), class = "mzlib_usage_error", contains = "non-empty")
  expect_error(build(coefficients = "a\nb"), class = "mzlib_usage_error", contains = "newline")
  expect_error(build(trend = "yes"), class = "mzlib_usage_error", contains = "TRUE or FALSE")
  expect_error(build(threads = 1.5), class = "mzlib_usage_error", contains = "threads must be a whole number")
  expect_error(build(spline_basis = 2.5), class = "mzlib_usage_error", contains = "spline_basis must be")
  expect_error(build(responses = ""), class = "mzlib_usage_error")
})

# ---------------------------------------------------------------- stats_fit: contrasts and the estimator

stats_malat_contrasts <- function() mz$stats_parse_fit(stats_fixture("stats_fit_malat_contrasts.json"))

test_that("contrasts, confidence and estimator are sent only when given, and contrasts alone send no stdin", {
  request <- mz$stats_build_fit_request("r.tsv", "d.tsv", character(0), FALSE, NULL, 1,
    contrasts = "contrasts_1.tsv", confidence = 0.9, estimator = "marginal_likelihood"
  )
  args <- request$args
  expect_identical(args[match("--contrasts", args) + 1L], "contrasts_1.tsv")
  expect_identical(as.numeric(args[match("--confidence", args) + 1L]), 0.9)
  expect_identical(args[match("--estimator", args) + 1L], "marginal_likelihood")
  expect_true(is.null(request$stdin), info = "contrasts alone: no coefficient on stdin, not even a blank line")
  plain <- mz$stats_build_fit_request("r.tsv", "d.tsv", "malat_250ng", FALSE, NULL, 1)
  expect_false(any(c("--contrasts", "--confidence", "--estimator") %in% plain$args),
    info = "left out, they are not sent, so mzLib's defaults apply"
  )
  for (name in c("contrasts", "confidence", "estimator")) {
    expect_true(is.null(formals(stats_fit)[[name]]), info = paste(name, "defaults to NULL: mzLib's default"))
  }
})

test_that("the recorded contrasts parse into their own table, the coefficient table empty", {
  fit <- stats_malat_contrasts()
  expect_identical(fit$tested_coefficients, character(0))
  expect_identical(c(fit$row_count, nrow(fit$records)), c(0, 0))
  expect_identical(names(fit$records), fit$column_names)
  expect_true(endsWith(fit$contrasts_file, "contrasts_1.tsv"))
  pairs <- fit$contrasts
  expect_true(inherits(pairs, "mzlibr_moderated_contrasts"))
  expect_identical(pairs$tested_contrasts, c("c1", "c2", "c3"))
  expect_identical(pairs$confidence_level, 0.95)
  expect_identical(fit$coefficient_names, c("intercept", "malat_load=250ng", "malat_load=500ng"))
  expect_identical(pairs$weights, list(c1 = c(0, 1, 0), c2 = c(0, 0, 1), c3 = c(0, -1, 1)))
  expect_identical(pairs$row_count, 3 * fit$feature_count)
  expect_identical(nrow(pairs$records), as.integer(pairs$row_count))
  expect_identical(names(pairs$records)[1:7], c(
    "feature", "contrast", "status", "estimate", "standard_error", "confidence_low", "confidence_high"
  ))
  rows <- pairs$records[pairs$records$contrast == "c3", ]
  expect_identical(nrow(rows), as.integer(fit$feature_count))
  expect_identical(rows$feature, pairs$records$feature[seq_len(fit$feature_count)])
  fitted <- rows$status == "fitted"
  expect_identical(!is.na(rows$confidence_low), fitted)
  expect_identical(!is.na(rows$bh_adjusted), fitted)
  expect_true(all(rows$confidence_low[fitted] < rows$estimate[fitted] & rows$estimate[fitted] < rows$confidence_high[fitted]))
  expect_true(any(startsWith(fit$caveats, "Contrasts are tested with the same variance prior")))
})

test_that("a fit without contrasts has none, and reports mzLib's default estimator", {
  fit <- stats_malat_fit()
  expect_true(is.null(fit$contrasts))
  expect_true(is.na(fit$contrasts_file))
  expect_identical(fit$prior$estimator, "moments_legacy")
  expect_true(startsWith(fit$caveats[[1L]], "moments_legacy:"))
})

test_that("the marginal-likelihood recording reports its estimator and a lower prior df", {
  fit <- mz$stats_parse_fit(stats_fixture("stats_fit_malat_marginal.json"))
  expect_identical(fit$prior$estimator, "marginal_likelihood")
  expect_true(fit$prior$estimator %in% mz$STATS_PRIOR_ESTIMATORS)
  expect_true(startsWith(fit$caveats[[1L]], "marginal_likelihood:"))
  expect_identical(round(fit$prior$df, 1), 11)
  expect_identical(round(stats_malat_fit()$prior$df, 1), 19.6)
  expect_identical(mz$STATS_PRIOR_ESTIMATORS, c("moments_legacy", "marginal_likelihood"))
})

test_that("bad contrast options are refused before the bridge", {
  build <- function(...) mz$stats_build_fit_request("r.tsv", "d.tsv", "x", FALSE, NULL, 1, ...)
  expect_error(build(estimator = "reml"), class = "mzlib_usage_error",
    contains = "estimator must be one of moments_legacy, marginal_likelihood")
  expect_error(build(estimator = 1), class = "mzlib_usage_error", contains = "estimator must be one of")
  expect_error(build(confidence = "0.9"), class = "mzlib_usage_error", contains = "confidence must be a number")
  expect_error(build(confidence = TRUE), class = "mzlib_usage_error", contains = "confidence must be a number")
  expect_error(build(contrasts = 3), class = "mzlib_usage_error", contains = "contrasts must be")
  expect_error(stats_fit("r.tsv", "d.tsv"), class = "mzlib_usage_error", contains = "or a contrasts file")
})

# What a bridge from before pyMzLib 0.6.0 answers: no contrasts, no estimator.
stats_older <- function(name) {
  data <- stats_fixture(name)
  data$contrasts_file <- NULL
  data$contrasts <- NULL
  data$prior$estimator <- NULL
  data
}

test_that("an older bridge that ignored contrasts or the estimator is named, not trusted", {
  older <- function(...) stats_older("stats_fit_malat.json")
  request <- mz$stats_build_fit_request("r.tsv", "d.tsv", "malat_250ng", FALSE, NULL, 1, contrasts = "c.tsv")
  expect_error(mz$stats_fit_answer(request, "c.tsv", NULL, 60, invoke = older), class = "mzlib_usage_error",
    contains = c("stats_fit(contrasts = ...) needs the bridge from pyMzLib 0.6.0", "ignored --contrasts"))
  expect_error(mz$stats_fit_answer(request, NULL, "moments_legacy", 60, invoke = older), class = "mzlib_usage_error",
    contains = "stats_fit(estimator = ...) needs the bridge from pyMzLib 0.6.0")
  # Without the new options, an older bridge's answer is read as it always was.
  plain <- mz$stats_parse_fit(mz$stats_fit_answer(request, NULL, NULL, 60, invoke = older))
  expect_true(is.na(plain$prior$estimator))
  expect_true(is.null(plain$contrasts))
})

test_that("an older bridge that wanted a coefficient is told apart from a real refusal", {
  refuse <- function(...) stop(mz$mzlib_usage_error("No coefficient to test was given on stdin. Give one name per line."))
  contrasts_only <- mz$stats_build_fit_request("r.tsv", "d.tsv", character(0), FALSE, NULL, 1, contrasts = "c.tsv")
  expect_error(mz$stats_fit_answer(contrasts_only, "c.tsv", NULL, 60, invoke = refuse), class = "mzlib_usage_error",
    contains = "needs the bridge from pyMzLib 0.6.0")
  other <- function(...) stop(mz$mzlib_usage_error("The design is not of full column rank."))
  expect_error(mz$stats_fit_answer(contrasts_only, "c.tsv", NULL, 60, invoke = other), class = "mzlib_usage_error",
    contains = "not of full column rank")
})

test_that("the contrasts chain replays end to end: diff_design(out = ), then stats_fit(contrasts = )", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
  design <- diff_design(samples, factors = list(factor_spec("malat_load")), out = file.path(tempdir(), "malat_pairs"))
  stratum <- design$strata[[1L]]
  pairs <- stats_fit("malat_dilution_log2_vs_fluc.tsv", stratum$design_file, contrasts = stratum$contrasts_file)
  expect_identical(pairs$contrasts$tested_contrasts, c("c1", "c2", "c3"))
  marginal <- stats_fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv", c("malat_250ng", "malat_125ng"),
    estimator = "marginal_likelihood")
  expect_identical(marginal$prior$estimator, "marginal_likelihood")
  # Without the estimator the default recording answers, never the marginal one.
  default <- stats_fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv", c("malat_250ng", "malat_125ng"))
  expect_identical(default$prior$estimator, "moments_legacy")
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

  # Bit for bit across thread counts: every core and one give the same answer.
  single <- stats_fit(stats_reference_path("limma_reference_responses.tsv"),
                      stats_reference_path("limma_reference_design.tsv"), "age_decades", threads = 1)
  expect_identical(fit$records$t, single$records$t, info = "the same at any thread count")
  expect_identical(fit$records$bh_adjusted, single$records$bh_adjusted, info = "the same at any thread count")

  # Against the recording, to 1e-9 relative, the tolerance pyMzLib's live check uses: still far
  # below the 1e-8 promised against limma, so drift between the recorded bridge and this one is
  # caught. Not bit for bit, because the recording was made at mzLib 1.0.593, and 1.0.595's prior
  # arithmetic differs from it by about 1e-14 relative (prior df 4.531798758719062 against
  # 4.5317987587192796), with identical estimates.
  recorded <- stats_limma_fit()
  expect_true(stats_worst(fit$records$t, recorded$records$t) <= 1e-9, info = "t against the recording")
  expect_true(stats_worst(fit$records$bh_adjusted, recorded$records$bh_adjusted) <= 1e-9,
    info = "bh_adjusted against the recording"
  )
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

test_that("LIVE: a small trend fit sizes its spline as limma does", {
  # mzLib 1.0.594 (#1398) follows limma's fitFDist rule, 1 + (n >= 3) + (n >= 6) + (n >= 30) basis
  # functions for n features, when it is not asked for a number. The bridge used to ask for 4 every
  # time (pyMzLib #82). Twelve features of limma's reference table must therefore get 3.
  skip_unless_bridge_has("stats fit")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  responses <- stats_reference("limma_reference_responses.tsv")[1:12, ]
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  utils::write.table(responses, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  fit <- stats_fit(path, stats_reference_path("limma_reference_design.tsv"), "age_decades", trend = TRUE)
  expect_true(fit$prior$trended)
  expect_identical(fit$prior$spline_basis_count, 3)
  asked <- stats_fit(path, stats_reference_path("limma_reference_design.tsv"), "age_decades",
                     trend = TRUE, spline_basis = 2)
  expect_identical(asked$prior$spline_basis_count, 2)
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

# ---------------------------------------------------------------- LIVE: contrasts and the prior estimator

stats_live_bridge <- function(verb) {
  skip_unless_bridge_has(verb)
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
}

# diff design on the dilution series with no reference level: every pair of loads a contrast.
stats_pairs_design <- function() {
  samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
  out <- tempfile("malat_pairs")
  diff_design(samples, factors = list(factor_spec("malat_load")), out = out)$strata[[1L]]
}

test_that("LIVE: the contrasts recording is what the chain returns today", {
  stats_live_bridge("diff design")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  stratum <- stats_pairs_design()
  live <- stats_fit(stats_reference_path("malat_dilution_log2_vs_fluc.tsv"), stratum$design_file,
    contrasts = stratum$contrasts_file)
  recorded <- stats_malat_contrasts()
  expect_identical(live$contrasts$tested_contrasts, recorded$contrasts$tested_contrasts)
  expect_identical(live$contrasts$weights, recorded$contrasts$weights)
  expect_identical(live$contrasts$confidence_level, recorded$contrasts$confidence_level)
  expect_identical(records_disagree(live$contrasts$records, recorded$contrasts$records), "")
  expect_identical(live$caveats, recorded$caveats)
})

test_that("LIVE: the recorded fits are what the bridge returns today, the repinned caveats included", {
  stats_live_bridge("stats fit")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  data <- system.file("extdata", "stats", package = "mzLibR")
  cases <- list(
    stats_fit_malat.json = list(estimator = NULL, relative = 1e-9),
    stats_fit_malat_marginal.json = list(estimator = "marginal_likelihood", relative = 1e-6)
  )
  for (name in names(cases)) {
    live <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"), file.path(data, "malat_dilution_design.tsv"),
      c("malat_250ng", "malat_125ng"), estimator = cases[[name]]$estimator)
    recorded <- mz$stats_parse_fit(stats_fixture(name))
    # The marginal-likelihood prior is a minimum found by Brent's search over log d0, flat near the
    # minimum, so the last bits of Math.Log move it by about sqrt(eps) across platforms: 1e-6.
    expect_identical(records_disagree(live$records, recorded$records, cases[[name]]$relative), "", info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
    expect_identical(live$prior$estimator, recorded$prior$estimator, info = name)
  }
})

test_that("LIVE: a contrast of one coefficient is that coefficient's test", {
  stats_live_bridge("diff design")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  stratum <- stats_pairs_design()
  fit <- stats_fit(stats_reference_path("malat_dilution_log2_vs_fluc.tsv"), stratum$design_file,
    "malat_load=250ng", contrasts = stratum$contrasts_file)
  expect_identical(fit$contrasts$weights$c1, c(0, 1, 0))
  coefficient <- fit$records
  contrast <- fit$contrasts$records[fit$contrasts$records$contrast == "c1", ]
  for (column in c("estimate", "standard_error", "t", "p_value", "bh_adjusted")) {
    expect_identical(coefficient[[column]], contrast[[column]], info = column)
  }
})

test_that("LIVE: a lower confidence gives narrower intervals, and mzLib refuses one outside 0 to 1", {
  stats_live_bridge("diff design")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  stratum <- stats_pairs_design()
  responses <- stats_reference_path("malat_dilution_log2_vs_fluc.tsv")
  wide <- stats_fit(responses, stratum$design_file, contrasts = stratum$contrasts_file)
  narrow <- stats_fit(responses, stratum$design_file, contrasts = stratum$contrasts_file, confidence = 0.8)
  expect_identical(narrow$contrasts$confidence_level, 0.8)
  w <- wide$contrasts$records
  n <- narrow$contrasts$records
  fitted <- w$status == "fitted"
  expect_true(all(w$confidence_low[fitted] < n$confidence_low[fitted] & n$confidence_low[fitted] < n$confidence_high[fitted] &
    n$confidence_high[fitted] < w$confidence_high[fitted]))
  expect_error(stats_fit(responses, stratum$design_file, contrasts = stratum$contrasts_file, confidence = 1.5),
    class = "mzlib_usage_error", contains = "strictly between 0 and 1")
  expect_error(stats_fit(responses, stratum$design_file, "malat_load=250ng", confidence = 0.8),
    class = "mzlib_usage_error", contains = "give --contrasts as well")
})

test_that("LIVE: the default estimator is the one mzLib reports, and gives the legacy numbers", {
  stats_live_bridge("stats fit")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  args <- list(stats_reference_path("malat_dilution_log2_vs_fluc.tsv"), stats_reference_path("malat_dilution_design.tsv"),
    "malat_250ng")
  default <- do.call(stats_fit, args)
  legacy <- do.call(stats_fit, c(args, list(estimator = "moments_legacy")))
  expect_identical(default$prior$estimator, "moments_legacy")
  expect_identical(default$records$t, legacy$records$t)
})
