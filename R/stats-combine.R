# One answer from several studies, and features against a per-sample trait, computed by mzLib.
#
#   stats_combine_p()   PValueCombination.Fisher / .Stouffer: one p-value per feature from several
#                       independent tests, optionally signed (PValueCombination.OneSided)
#   stats_correlate()   SpearmanCorrelation.Correlate / .PartialCorrelate: every feature against a
#                       per-sample trait, optionally given covariates
#
# Ported from pyMzLib's `stats.py` (combine_p, correlate; mzLib #1367, #1446). Every number is
# mzLib's: this file writes the tests as lines on stdin, or names the files, and projects each
# result into a data.frame, `records`. The tests hold the recordings to R's own answers
# (combination_results.tsv; stats::cor.test() run here).

# The first pyMzLib whose bridge has `stats combine-p` and `stats correlate` (mzLib 1.0.595).
STATS_COMBINE_CORRELATE_SINCE <- "0.6.0"

# The rules stats_combine_p() takes as `method`: mzLib's PValueCombinationMethod. There is no
# default; you choose.
STATS_COMBINATION_METHODS <- c("fisher", "stouffer")

# Every value of stats_correlate()'s `method` column, as mzLib's SpearmanPValueMethod is written on
# the wire. Only "not_estimable" rows have NA for `rho` and `p_value`.
STATS_SPEARMAN_METHODS <- c("exact", "asymptotic", "not_estimable")

# ---------------------------------------------------------------- argument shaping

# A numeric column of `tests` as stdin fields: NA and NaN are blank (missing), a number exact.
stats_cells <- function(values, what) {
  if (!is.numeric(values) && !(is.logical(values) && all(is.na(values)))) {
    stop(mzlib_usage_error(paste0(what, " must be numbers, NA for none; got ", class(values)[1L], ".")))
  }
  values <- as.numeric(values)
  ifelse(is.na(values), "", stats_number_text(values))
}

# One line per test: feature, p_value, then weight and effect when `tests` has those columns. With
# an effect but no weight column the weight field is blank, so the effect stays in field 4.
stats_combine_stdin <- function(tests) {
  if (!is.data.frame(tests)) {
    stop(mzlib_usage_error(paste0(
      "tests must be a data.frame with feature and p_value columns, and optionally weight and ",
      "effect, one row per test; got ", class(tests)[1L], "."
    )))
  }
  missing <- setdiff(c("feature", "p_value"), names(tests))
  if (length(missing) > 0L) {
    stop(mzlib_usage_error(paste0(
      "tests has no ", paste(missing, collapse = " or "), " column; it needs feature and p_value, ",
      "and optionally weight and effect."
    )))
  }
  if (nrow(tests) == 0L) {
    stop(mzlib_usage_error("stats_combine_p() needs at least one test."))
  }
  feature <- as.character(tests$feature)
  for (i in seq_along(feature)) {
    if (is.na(feature[[i]]) || !nzchar(trimws(feature[[i]]))) {
      stop(mzlib_usage_error(paste0("tests row ", i, " needs a non-empty feature name.")))
    }
    if (grepl("[\t\r\n]", feature[[i]])) {
      stop(mzlib_usage_error(paste0("tests row ", i, "'s feature name contains a tab or newline: '", feature[[i]], "'.")))
    }
  }
  fields <- list(feature, stats_cells(tests$p_value, "tests$p_value"))
  has_effect <- "effect" %in% names(tests)
  if ("weight" %in% names(tests) || has_effect) {
    fields[[3L]] <- if ("weight" %in% names(tests)) stats_cells(tests$weight, "tests$weight") else rep("", nrow(tests))
  }
  if (has_effect) {
    fields[[4L]] <- stats_cells(tests$effect, "tests$effect")
  }
  # No trailing newline: system2(input =) ends the last line itself.
  paste(do.call(paste, c(fields, sep = "\t")), collapse = "\n")
}

stats_build_correlate_args <- function(responses, trait, given) {
  args <- c(
    "stats", "correlate",
    "--responses", readers_normalise_path(responses, "responses"),
    "--trait", readers_normalise_path(trait, "trait")
  )
  if (!is.null(given)) {
    args <- c(args, "--given", readers_normalise_path(given, "given"))
  }
  args
}

# ---------------------------------------------------------------- parsing

stats_parse_combined <- function(data) {
  structure(
    list(
      method = wire_field(data, "method", "character", NA_character_),
      weighted = isTRUE(data[["weighted"]]),
      signed = isTRUE(data[["signed"]]),
      line_count = stats_count(data, "line_count"),
      feature_count = stats_count(data, "feature_count"),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_combined_p_values"
  )
}

stats_parse_correlations <- function(data) {
  incomplete <- data[["given_incomplete_count"]]
  structure(
    list(
      responses_file = wire_field(data, "responses_file", "character", NA_character_),
      trait_file = wire_field(data, "trait_file", "character", NA_character_),
      given_file = wire_field(data, "given_file", "character", NA_character_),
      trait_name = wire_field(data, "trait_name", "character", NA_character_),
      covariate_names = wire_strings(data[["covariate_names"]]),
      partial = isTRUE(data[["partial"]]),
      feature_count = stats_count(data, "feature_count"),
      sample_count = stats_count(data, "sample_count"),
      sample_names = wire_strings(data[["sample_names"]]),
      missing_count = stats_count(data, "missing_count"),
      trait_missing_count = stats_count(data, "trait_missing_count"),
      given_incomplete_count = if (wire_null(incomplete)) NA_real_ else as.numeric(incomplete),
      method_counts = proteins_parse_counts(data[["method_counts"]]),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_correlations"
  )
}

# ---------------------------------------------------------------- the public surface

#' Combine several independent p-values per feature into one, by Fisher's or Stouffer's method
#'
#' Each feature is combined on its own by mzLib's `PValueCombination`: Fisher (1925),
#' X = -2 sum ln p against chi-square on 2k degrees of freedom, or Stouffer (1949),
#' Z = sum w z / sqrt(sum w^2) against the standard normal. mzLib holds both to R's `pchisq()` and
#' `qnorm()` to 1e-8 relative. Use it when you have one p-value per study - per dataset, cohort or
#' laboratory - for the same feature, and want one answer from all of them.
#'
#' **Direction.** Both rules combine *one-sided* p-values. A two-sided p-value from [stats_fit()]
#' does not say which way the effect went, so two studies that disagree look like two that agree.
#' Give every test its `effect` (its estimate; only the sign is used) and mzLib first turns each p
#' into the one-sided p-value for effect > 0 (`PValueCombination.OneSided`: p/2 for a positive
#' effect, 1 - p/2 for a negative one), so opposite effects cancel. The combined `p_value` then
#' tests effect > 0; for Stouffer, a `p_value` near 1 means a consistently negative effect.
#'
#' Neither rule estimates an effect size: [stats_meta()] does that, from estimates and standard
#' errors. Neither adjusts across features: send `records$p_value` to [stats_adjust()].
#'
#' @param tests A data.frame with one row per test: `feature` and `p_value`, and optionally `weight`
#'   and `effect` columns. A feature's tests need not be adjacent. `p_value` is a fraction from 0 to
#'   1, or `NA` for a test that could not be run (it is omitted, and counted in the `omitted`
#'   column). `weight` (Stouffer only) is a positive number, for example the square root of the
#'   study's sample size; give it for every test or for none (`NA` is no weight). `effect` is the
#'   test's estimate in any unit; give it for every test or for none, and `NA` for a test with no
#'   direction (omitted). An `effect` column without a `weight` column sends no weights.
#' @param method `"fisher"` or `"stouffer"`. Required: mzLib has no preferred rule. Fisher is driven
#'   by the strongest single result and takes no weights; Stouffer weighs every study, and takes
#'   weights.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_combined_p_values`:
#'
#'   - `records`, a data.frame with one row per feature in order of first appearance: `feature`;
#'     `studies`, the tests combined (tests with a finite p-value, after the signed conversion when
#'     effects were given); `omitted`, the feature's tests that were not combined (tests): no
#'     p-value, or (signed) an effect of 0 or `NA`; `statistic`, Fisher's X or Stouffer's Z (fisher:
#'     chi-square statistic on 2 x studies degrees of freedom; stouffer: standard normal deviates),
#'     `NA` when it is not finite; and
#'     `p_value`, a fraction (0 to 1), the upper tail of the statistic, so ONE direction (effect > 0
#'     when signed), `NA` when there is nothing to combine, or a p of 0 and a p of 1 together
#'     (Stouffer). An `NA` statistic with a p-value is an infinite one: `p_value` 0 means plus
#'     infinity (a p-value of 0 among the tests), 1 means minus infinity (Stouffer, a p-value of 1).
#'   - `method`, as you asked; `weighted` (Stouffer with weights); `signed` (effects were given,
#'     so every p-value was read as two-sided and converted before combining); `line_count` tests
#'     (lines) sent over every feature; `feature_count` features; `row_count` rows; `column_names`.
#'   - `caveats`: what the combination assumes, plus any that apply to this input (omitted tests,
#'     infinite statistics).
#'
#' @section The trap:
#'
#' Two-sided p-values given without their effects are combined as if one-sided, so two studies
#' whose effects point opposite ways reinforce each other instead of cancelling. Give `effect` for
#' that. mzLib computes no two-sided combined p-value.
#'
#' @seealso [stats_meta()] to pool effect sizes, [stats_adjust()] to adjust across features
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own reference cases (constructed p-values, not proteomics data), by both rules:
#' cases <- read.delim(system.file("extdata", "stats", "combination_inputs.tsv", package = "mzLibR"))
#' fisher <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p), method = "fisher")
#' stouffer <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p,
#'   weight = cases$weight), method = "stouffer")
#' data.frame(feature = fisher$records$feature, studies = fisher$records$studies,
#'   fisher = signif(fisher$records$p_value, 3),
#'   weighted_stouffer = signif(stouffer$records$p_value, 3))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.combine-p
#' @export
stats_combine_p <- function(tests, method, timeout = 60) {
  if (missing(method) || !is.character(method) || length(method) != 1L || !method %in% STATS_COMBINATION_METHODS) {
    stop(mzlib_usage_error(paste0(
      "method must be one of ", paste(STATS_COMBINATION_METHODS, collapse = ", "),
      " (mzLib's PValueCombinationMethod); mzLib prefers neither, so there is no default",
      if (!missing(method)) paste0("; got ", paste(deparse(method), collapse = " ")) else "", "."
    )))
  }
  stdin <- stats_combine_stdin(tests)
  bridge_require_verb("stats combine-p", STATS_COMBINE_CORRELATE_SINCE)
  stats_parse_combined(bridge_invoke(c("stats", "combine-p", "--method", method),
    stdin = stdin, timeout = timeout
  ))
}

#' Spearman-correlate every feature with a per-sample trait, optionally given covariates
#'
#' For each feature (row) of the responses table, mzLib's `SpearmanCorrelation.Correlate` ranks
#' the feature's values and the trait over the samples where both are known, and tests rho = 0,
#' two-sided. With `given`, `SpearmanCorrelation.PartialCorrelate` first removes from both sets of
#' ranks what the covariates' ranks explain, as `ppcor::pcor.test(method = "spearman")` does. mzLib
#' holds both to R's `cor.test()` and to ppcor 1.1 to 1e-8 relative.
#'
#' Spearman uses ranks only. It needs no log transform, one extreme value cannot drive it, and it
#' finds any relationship that only ever rises (or only ever falls), not just a straight line.
#'
#' @param responses [stats_fit()]'s responses table: a TSV with a header row, the feature id first,
#'   then one column per sample, named by its header. Blank, `NA` or `NaN` is missing. **0 is an
#'   observation**, ranked like any other value.
#' @param trait A TSV with a header row and one row per sample: the sample name (exactly the
#'   responses' sample columns, in any order), then ONE value column, whose header names the trait -
#'   age, dose, loading. A blank or `NA` cell leaves that sample out.
#' @param given A TSV of the same shape with one or more covariate columns (batch, run order,
#'   loading) to partial out. Leave it `NULL` for a plain correlation.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_correlations`:
#'
#'   - `records`, a data.frame with one row per feature in the responses file's order: `feature`;
#'     `rho`, Spearman's rho (unitless correlation, -1 to 1), ties at their average rank, the partial rho with
#'     `given`; `p_value`, a fraction (0 to 1), two-sided against rho = 0 and **not adjusted across
#'     features** - both `NA` exactly when `method` is `"not_estimable"`; `n`, the complete samples
#'     used (samples: response, trait and every covariate finite); `has_ties`, whether any variable
#'     was tied among them (`FALSE` too when there were too few to rank); and `method`, `"exact"`
#'     (n <= 9, no ties, plain), `"asymptotic"` (t on n - 2 df, or n - 2 - k with k covariates kept)
#'     or `"not_estimable"`.
#'   - `responses_file`, `trait_file` and `given_file` (`NA` without `given`), as absolute paths;
#'     `trait_name`, the trait column's header; `covariate_names` (empty without `given`);
#'     `partial`, `TRUE` with `given`; `feature_count` features and `sample_count` samples read;
#'     `sample_names`; `missing_count` response cells (cells) that were blank, `NA` or `NaN`;
#'     `trait_missing_count`, samples with no trait value (samples); `given_incomplete_count`,
#'     samples missing a covariate (samples), `NA` without `given`; `method_counts`, features per
#'     method; `row_count` rows; `column_names`.
#'   - `caveats`: how the p-values were computed and how they differ from R's defaults, plus a count
#'     of not_estimable features when there are any.
#'
#' @section Against R's cor.test():
#'
#' For 10 <= n < 1290 without ties, R's `cor.test()` default is an Edgeworth series (AS 89), which
#' mzLib does not use, so p-values there differ slightly from R's default output; they match
#' `cor.test(exact = FALSE)`. Send `records$p_value` to [stats_adjust()] for Benjamini-Hochberg.
#'
#' @seealso [stats_fit()] for a linear model on the same responses, [stats_adjust()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # The MALAT1 dilution series stats_fit() analyses, correlated with each run's MALAT1 load:
#' data <- system.file("extdata", "stats", package = "mzLibR")
#' result <- stats_correlate(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
#'   file.path(data, "malat_dilution_load_ng.tsv"))
#' result
#' rho <- result$records$rho
#' c(MALAT1 = median(rho[startsWith(result$records$feature, "MALAT1:")], na.rm = TRUE),
#'   FLuc = median(rho[startsWith(result$records$feature, "FLuc:")], na.rm = TRUE))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.correlate
#' @export
stats_correlate <- function(responses, trait, given = NULL, timeout = 600) {
  args <- stats_build_correlate_args(responses, trait, given)
  bridge_require_verb("stats correlate", STATS_COMBINE_CORRELATE_SINCE)
  stats_parse_correlations(bridge_invoke(args, timeout = timeout))
}

# ---------------------------------------------------------------- printing

#' Print combined p-values
#'
#' @param x A [stats_combine_p()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_combined_p_values <- function(x, ...) {
  cat("<mzlibr_combined_p_values> ", format(x$feature_count), " features from ", format(x$line_count),
    " tests, ", x$method, if (x$weighted) ", weighted" else "", if (x$signed) ", signed" else "", "\n",
    sep = ""
  )
  invisible(x)
}

#' Print Spearman correlations
#'
#' @param x A [stats_correlate()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_correlations <- function(x, ...) {
  cat("<mzlibr_correlations> ", format(x$feature_count), " features x ", format(x$sample_count),
    " samples against ", x$trait_name,
    if (x$partial) paste0(", given ", paste(x$covariate_names, collapse = ", ")) else "", "\n",
    sep = ""
  )
  if (length(x$method_counts) > 0L) {
    cat("  ", paste0(names(x$method_counts), " ", x$method_counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
