# stats_combine_p() and stats_correlate().
#
# Ported from pyMzLib's test_stats_combine_correlate.py and its _live.py. The payloads are pyMzLib's
# recordings, verbatim; the reference tables are mzLib's, shipped in extdata/stats/. Three layers:
#   1. always: what is sent, what is parsed, refusals before the bridge, and the recordings against
#      R's own answers - combination_results.tsv, and stats::cor() / cor.test() run here;
#   2. replayed end to end through the stand-in bridge;
#   3. LIVE: the bridge today returns what was recorded, still agrees with R, and refuses alike.

combine_reference <- function(name) {
  utils::read.delim(system.file("extdata", "stats", name, package = "mzLibR"), check.names = FALSE,
    stringsAsFactors = FALSE)
}
combine_path <- function(name) system.file("extdata", "stats", name, package = "mzLibR")

# mzLib's signed cohorts (pyMzLib's COHORTS): three features, each three studies.
COMBINE_COHORTS <- data.frame(
  feature = rep(c("up_everywhere", "flips_sign", "no_effect"), each = 3L),
  p_value = c(0.003, 0.02, 0.3, 0.001, 0.004, 0.4, 0.6, 0.45, 0.85),
  weight = c(4, 5, 3, 4, 5, 3, 4, 5, 3),
  effect = c(0.8, 0.5, 0.3, 1.2, -0.9, 0.2, -0.1, 0.1, -0.05),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------- combine_p: what is sent

test_that("combine_p writes one line per test, a blank field for NA", {
  tests <- data.frame(feature = c("A", "A", "B", "B"), p_value = c(0.01, NA, 0.2, NaN), weight = c(4, 5, 3, 2.5),
    effect = c(0.5, -1, NA, 0.1))
  lines <- strsplit(mz$stats_combine_stdin(tests), "\n", fixed = TRUE)[[1L]]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  expect_identical(length(lines), 4L)
  expect_identical(vapply(fields, `[`, "", 1L), c("A", "A", "B", "B"))
  expect_identical(as.numeric(fields[[1L]][2:4]), c(0.01, 4, 0.5))
  expect_true(grepl("^A\t\t5\t-1$", lines[[2L]]), info = "an NA p-value crosses as a blank field")
  expect_true(endsWith(lines[[3L]], "\t"), info = "an NA effect is a blank fourth field")
})

test_that("combine_p sends only the columns tests has; an effect with no weight keeps field 4", {
  expect_identical(mz$stats_combine_stdin(data.frame(feature = c("A", "B"), p_value = c(0.5, 0.25))), "A\t0.5\nB\t0.25")
  expect_identical(mz$stats_combine_stdin(data.frame(feature = "A", p_value = 0.5, effect = -2)), "A\t0.5\t\t-2")
  expect_identical(mz$stats_combine_stdin(data.frame(feature = "A", p_value = 0.5, weight = NA)), "A\t0.5\t")
})

test_that("bad combine_p arguments are refused before the bridge", {
  ok <- data.frame(feature = "A", p_value = 0.1)
  expect_error(stats_combine_p(ok, method = "Fisher"), class = "mzlib_usage_error", contains = "method must be one of fisher, stouffer")
  expect_error(stats_combine_p(ok), class = "mzlib_usage_error", contains = "no default")
  expect_error(stats_combine_p(ok[0, ], method = "fisher"), class = "mzlib_usage_error", contains = "at least one test")
  expect_error(stats_combine_p("A\t0.1", method = "fisher"), class = "mzlib_usage_error", contains = "must be a data.frame")
  expect_error(stats_combine_p(list(feature = "A", p_value = 0.1), method = "fisher"), class = "mzlib_usage_error",
    contains = "must be a data.frame")
  expect_error(stats_combine_p(data.frame(feature = "A"), method = "fisher"), class = "mzlib_usage_error",
    contains = "no p_value column")
  expect_error(stats_combine_p(data.frame(feature = "", p_value = 0.1), method = "fisher"), class = "mzlib_usage_error",
    contains = "non-empty feature")
  expect_error(stats_combine_p(data.frame(feature = "A\tB", p_value = 0.1), method = "fisher"), class = "mzlib_usage_error",
    contains = "tab or newline")
  expect_error(stats_combine_p(data.frame(feature = "A", p_value = "0.1"), method = "fisher"), class = "mzlib_usage_error",
    contains = "tests$p_value must be numbers")
  expect_error(stats_combine_p(data.frame(feature = "A", p_value = 0.1, weight = "heavy"), method = "stouffer"),
    class = "mzlib_usage_error", contains = "tests$weight must be numbers")
  expect_error(stats_combine_p(data.frame(feature = "A", p_value = 0.1, effect = "up"), method = "stouffer"),
    class = "mzlib_usage_error", contains = "tests$effect must be numbers")
})

# ---------------------------------------------------------------- combine_p: what is parsed

test_that("the Fisher recording parses, one row per feature in order of first appearance", {
  combined <- mz$stats_parse_combined(recorded_data("stats_combine_p_fisher.json"))
  expect_true(inherits(combined, "mzlibr_combined_p_values"))
  expect_identical(c(combined$method, combined$weighted, combined$signed), c("fisher", "FALSE", "FALSE"))
  r <- combined$records
  expect_identical(r$feature, c("few", "tiny", "one", "many", "mixed", "tail"))
  expect_identical(r$studies, c(3, 4, 1, 12, 5, 3))
  expect_identical(r$omitted, rep(0, 6))
  expect_identical(c(combined$line_count, combined$feature_count, combined$row_count), c(28, 6, 6))
  expect_true(abs(r$p_value[r$feature == "one"] - 0.037) < 1e-12)
  expect_identical(combined$column_names, c("feature", "studies", "omitted", "statistic", "p_value"))
  expect_identical(length(combined$caveats), 3L)
})

test_that("the signed recording cancels a feature whose studies disagree", {
  signed <- mz$stats_parse_combined(recorded_data("stats_combine_p_signed.json"))
  expect_true(signed$signed && signed$weighted)
  expect_identical(signed$method, "stouffer")
  p <- stats::setNames(signed$records$p_value, signed$records$feature)
  expect_true(p[["flips_sign"]] > 0.4 && 0.4 > p[["up_everywhere"]])
  expect_true(any(startsWith(signed$caveats, "Each p-value was read as two-sided")))
  unsigned <- mz$stats_parse_combined(recorded_data("stats_combine_p_unsigned.json"))
  expect_true(unsigned$records$p_value[unsigned$records$feature == "flips_sign"] < 0.01,
    info = "without effects, opposite directions reinforce")
})

test_that("the recorded combinations are R's pchisq() and qnorm() answers, to 1e-8", {
  inputs <- combine_reference("combination_inputs.tsv")
  results <- combine_reference("combination_results.tsv")
  pairs <- list(c("stats_combine_p_fisher.json", "fisher_p"), c("stats_combine_p_weighted.json", "weighted_p"))
  for (pair in pairs) {
    r <- mz$stats_parse_combined(recorded_data(pair[[1L]]))$records
    theirs <- results[[pair[[2L]]]][match(r$feature, results$case)]
    expect_true(max(abs(r$p_value - theirs) / theirs) < 1e-8, info = pair[[1L]])
  }
  # ... and the same rules written in base R here.
  fisher <- vapply(split(inputs$p, factor(inputs$case, levels = unique(inputs$case))), function(p) {
    stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
  }, numeric(1L))
  r <- mz$stats_parse_combined(recorded_data("stats_combine_p_fisher.json"))$records
  expect_true(max(abs(r$p_value - fisher[r$feature]) / fisher[r$feature]) < 1e-8)
  stouffer <- vapply(split(inputs, factor(inputs$case, levels = unique(inputs$case))), function(s) {
    stats::pnorm(sum(s$weight * stats::qnorm(s$p, lower.tail = FALSE)) / sqrt(sum(s$weight^2)), lower.tail = FALSE)
  }, numeric(1L))
  w <- mz$stats_parse_combined(recorded_data("stats_combine_p_weighted.json"))$records
  expect_true(max(abs(w$p_value - stouffer[w$feature]) / stouffer[w$feature]) < 1e-8)
})

# ---------------------------------------------------------------- correlate

test_that("correlate names the files, and given only when set", {
  args <- mz$stats_build_correlate_args("r.tsv", "t.tsv", NULL)
  expect_identical(args, c("stats", "correlate", "--responses", "r.tsv", "--trait", "t.tsv"))
  expect_identical(utils::tail(mz$stats_build_correlate_args("r.tsv", "t.tsv", "g.tsv"), 2L), c("--given", "g.tsv"))
  expect_error(stats_correlate(5, "t.tsv"), class = "mzlib_usage_error", contains = "responses must be")
  expect_error(stats_correlate("r.tsv", "t.tsv", given = ""), class = "mzlib_usage_error", contains = "given must be")
})

test_that("the MALAT1 correlation parses; rho and p are NA exactly when not estimable", {
  result <- mz$stats_parse_correlations(recorded_data("stats_correlate_malat.json"))
  expect_true(inherits(result, "mzlibr_correlations"))
  expect_identical(result$trait_name, "malat_ng")
  expect_false(result$partial)
  expect_true(is.na(result$given_file) && is.na(result$given_incomplete_count))
  expect_identical(result$covariate_names, character(0))
  expect_identical(sort(names(result$method_counts)), sort(mz$STATS_SPEARMAN_METHODS))
  expect_identical(c(sum(result$method_counts), result$feature_count, result$row_count), c(492, 492, 492))
  expect_identical(result$method_counts[["exact"]], 4)
  expect_identical(result$column_names, c("feature", "rho", "p_value", "n", "has_ties", "method"))
  r <- result$records
  estimable <- r$method != "not_estimable"
  expect_identical(!is.na(r$rho), estimable)
  expect_identical(!is.na(r$p_value), estimable)
  exact <- r$method == "exact"
  expect_true(all(r$n[exact] <= 9 & !r$has_ties[exact]))
  med <- function(prefix) stats::median(r$rho[estimable & startsWith(r$feature, prefix)])
  expect_identical(sprintf("%+.2f", c(med("MALAT1:"), med("FLuc:"))), c("+0.43", "+0.00"))
})

test_that("the partial recording reports its covariate", {
  result <- mz$stats_parse_correlations(recorded_data("stats_correlate_limma_age_given_sex.json"))
  expect_true(result$partial)
  expect_identical(result$covariate_names, "sex")
  expect_true(endsWith(result$given_file, "limma_reference_sex.tsv"))
  expect_identical(result$given_incomplete_count, 0)
  expect_identical(result$method_counts[c("exact", "asymptotic", "not_estimable")], c(exact = 0, asymptotic = 400, not_estimable = 0))
})

test_that("the recorded Spearman correlations are stats::cor.test()'s, run here, to 1e-8", {
  # The limma reference responses (missing values in every feature) against age_decades: each row's
  # rho is cor(method = "spearman") on the complete samples, and its p-value cor.test()'s - exact
  # for the "exact" rows, exact = FALSE (t on n - 2 df) for the "asymptotic" ones.
  r <- mz$stats_parse_correlations(recorded_data("stats_correlate_limma_age.json"))$records
  responses <- combine_reference("limma_reference_responses.tsv")
  age <- combine_reference("limma_reference_age.tsv")
  y <- as.matrix(responses[, -1L])
  x <- age$age_decades[match(colnames(y), age$sample)]
  worst <- 0
  for (i in seq_len(nrow(y))[r$method != "not_estimable"]) {
    keep <- is.finite(y[i, ]) & is.finite(x)
    test <- suppressWarnings(stats::cor.test(y[i, keep], x[keep], method = "spearman", exact = r$method[[i]] == "exact"))
    worst <- max(worst, abs(r$rho[[i]] - unname(test$estimate)) / max(abs(test$estimate), 1e-300))
    if (abs(r$rho[[i]]) == 1 && r$method[[i]] == "asymptotic") {
      # A perfect rank correlation has an infinite t: mzLib reports p = 0, where cor.test() prints
      # a p-value below 1e-60. Both say the same thing; neither is a number to compare.
      expect_true(r$p_value[[i]] == 0 && test$p.value < 1e-60, info = r$feature[[i]])
      next
    }
    worst <- max(worst, abs(r$p_value[[i]] - test$p.value) / test$p.value)
  }
  expect_true(worst < 1e-8, info = paste("worst relative difference", worst))
})

# ---------------------------------------------------------------- end to end, replayed

test_that("combine_p and correlate replay end to end, each call answered by its own recording", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  cases <- combine_reference("combination_inputs.tsv")
  fisher <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p), method = "fisher")
  expect_identical(fisher$method, "fisher")
  weighted <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p, weight = cases$weight), method = "stouffer")
  expect_true(weighted$weighted && !weighted$signed)
  signed <- stats_combine_p(COMBINE_COHORTS, method = "stouffer")
  expect_true(signed$signed)
  unsigned <- stats_combine_p(COMBINE_COHORTS[, 1:3], method = "stouffer")
  expect_false(unsigned$signed)
  expect_identical(stats_correlate("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_load_ng.tsv")$trait_name, "malat_ng")
  plain <- stats_correlate("limma_reference_responses.tsv", "limma_reference_age.tsv")
  partial <- stats_correlate("limma_reference_responses.tsv", "limma_reference_age.tsv", given = "limma_reference_sex.tsv")
  expect_identical(c(plain$partial, partial$partial), c(FALSE, TRUE))
})

# ---------------------------------------------------------------- against a real mzLib

combine_live <- function(verb) {
  skip_unless_bridge_has(verb)
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
}

test_that("LIVE: the combine-p and correlate recordings are what the bridge returns today", {
  combine_live("stats combine-p")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  cases <- combine_reference("combination_inputs.tsv")
  calls <- list(
    stats_combine_p_fisher.json = function() stats_combine_p(data.frame(feature = cases$case, p_value = cases$p), method = "fisher"),
    stats_combine_p_weighted.json = function() {
      stats_combine_p(data.frame(feature = cases$case, p_value = cases$p, weight = cases$weight), method = "stouffer")
    },
    stats_combine_p_signed.json = function() stats_combine_p(COMBINE_COHORTS, method = "stouffer"),
    stats_combine_p_unsigned.json = function() stats_combine_p(COMBINE_COHORTS[, 1:3], method = "stouffer"),
    stats_correlate_malat.json = function() stats_correlate(combine_path("malat_dilution_log2_vs_fluc.tsv"), combine_path("malat_dilution_load_ng.tsv")),
    stats_correlate_limma_age.json = function() stats_correlate(combine_path("limma_reference_responses.tsv"), combine_path("limma_reference_age.tsv")),
    stats_correlate_limma_age_given_sex.json = function() {
      stats_correlate(combine_path("limma_reference_responses.tsv"), combine_path("limma_reference_age.tsv"),
        given = combine_path("limma_reference_sex.tsv"))
    }
  )
  for (name in names(calls)) {
    live <- calls[[name]]()
    parse <- if (startsWith(name, "stats_combine")) mz$stats_parse_combined else mz$stats_parse_correlations
    recorded <- parse(recorded_data(name))
    expect_identical(records_disagree(live$records, recorded$records), "", info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
    if (!is.null(recorded$method_counts)) {
      # Its keys come from an mzLib enum; a repin that grows it must not slip past the recording.
      expect_identical(live$method_counts, recorded$method_counts, info = name)
    }
  }
})

test_that("LIVE: combine-p still agrees with R, and refusals cross as usage errors", {
  combine_live("stats combine-p")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  cases <- combine_reference("combination_inputs.tsv")
  results <- combine_reference("combination_results.tsv")
  plain <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p), method = "stouffer")
  theirs <- results$stouffer_p[match(plain$records$feature, results$case)]
  expect_true(max(abs(plain$records$p_value - theirs) / theirs) < 1e-8)
  expect_error(stats_combine_p(data.frame(feature = c("A", "A"), p_value = c(0.1, 0.2), weight = c(2, NA)), method = "stouffer"),
    class = "mzlib_usage_error", contains = "every line or on none")
  expect_error(stats_combine_p(data.frame(feature = "A", p_value = 0.1, weight = 2), method = "fisher"),
    class = "mzlib_usage_error", contains = "takes no weights")
  expect_error(stats_combine_p(data.frame(feature = "A", p_value = 1.5), method = "fisher"),
    class = "mzlib_usage_error", contains = "outside [0, 1]")
  trait <- tempfile(fileext = ".tsv")
  on.exit(unlink(trait), add = TRUE)
  writeLines("sample\tage\tsex", trait)
  expect_error(stats_correlate(combine_path("limma_reference_responses.tsv"), trait), class = "mzlib_usage_error",
    contains = "has 2 value columns")
})
