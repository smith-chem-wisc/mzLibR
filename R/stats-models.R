# Per-feature models beyond limma's, and an entrapment check of a search's FDR, computed by mzLib.
#
#   stats_fit_logistic()    LogisticRegression.Fit: glm(family = binomial) per feature of a 0/1 table
#   stats_fit_mixed()       MixedModel.Fit: nlme::lme(random = ~ 1 | group) per feature
#   stats_entrapment_fdp()  EntrapmentFdp.Sweep: Wen et al. 2025's FDP estimates per q-value
#
# Ported from pyMzLib's `stats.py` (fit_logistic, fit_mixed, entrapment_fdp; mzLib #1367, #1385).
# Every number is mzLib's: this file names the files, sends the coefficients on stdin and projects
# each result into a data.frame, `records`. The tests hold the recordings to R's own glm() and
# nlme::lme(), run here, and to the references mzLib's tests freeze.

# The first pyMzLib whose bridge has `stats fit-logistic`, `stats fit-mixed` and
# `stats entrapment-fdp` (mzLib 1.0.595).
STATS_MODELS_SINCE <- "0.6.0"

# Every value of stats_fit_logistic()'s `status` column: the outcomes mzLib's LogisticRegression.Fit
# can give a feature, read from its code (mzLib does not list them per model yet). Only "fitted"
# rows carry statistics.
STATS_LOGISTIC_STATUSES <- c("fitted", "too_few_observations", "rank_deficient", "separated", "not_converged")

# Every value of stats_fit_mixed()'s `status` column: the outcomes mzLib's MixedModel.Fit can give
# a feature, read from its code. Only "fitted" rows carry statistics.
STATS_MIXED_STATUSES <- c(
  "fitted", "too_few_observations", "rank_deficient", "too_few_groups", "variance_ratio_at_limit",
  "not_converged"
)

# The variance estimators stats_fit_mixed() takes as `estimator`: mzLib's VarianceEstimator Reml
# and MaximumLikelihood. Leave `estimator` NULL for mzLib's default.
STATS_MIXED_ESTIMATORS <- c("reml", "ml")

# ---------------------------------------------------------------- argument shaping

stats_build_logistic_request <- function(responses, design, coefficients, max_iterations, tolerance, threads) {
  coefficients <- stats_coefficient_names(coefficients, caller = "stats_fit_logistic")
  args <- c(
    "stats", "fit-logistic",
    "--responses", readers_normalise_path(responses, "responses"),
    "--design", readers_normalise_path(design, "design"),
    "--threads", stats_whole(threads, "threads")
  )
  # Each is sent only when given, so mzLib's own default applies otherwise.
  if (!is.null(max_iterations)) {
    args <- c(args, "--max-iterations", stats_whole(max_iterations, "max_iterations"))
  }
  if (!is.null(tolerance)) {
    args <- c(args, "--tolerance", stats_number(tolerance, "tolerance"))
  }
  list(args = args, stdin = paste(coefficients, collapse = "\n"))
}

stats_build_mixed_request <- function(responses, design, groups, coefficients, estimator, threads) {
  coefficients <- stats_coefficient_names(coefficients, caller = "stats_fit_mixed")
  args <- c(
    "stats", "fit-mixed",
    "--responses", readers_normalise_path(responses, "responses"),
    "--design", readers_normalise_path(design, "design"),
    "--groups", readers_normalise_path(groups, "groups"),
    "--threads", stats_whole(threads, "threads")
  )
  if (!is.null(estimator)) {
    if (!is.character(estimator) || length(estimator) != 1L || !estimator %in% STATS_MIXED_ESTIMATORS) {
      stop(mzlib_usage_error(paste0(
        "estimator must be one of ", paste(STATS_MIXED_ESTIMATORS, collapse = ", "),
        " (mzLib's VarianceEstimator Reml or MaximumLikelihood), or NULL for mzLib's default; got ",
        paste(deparse(estimator), collapse = " "), "."
      )))
    }
    args <- c(args, "--estimator", estimator)
  }
  list(args = args, stdin = paste(coefficients, collapse = "\n"))
}

stats_build_entrapment_args <- function(targets, entrapments, ratio, thresholds) {
  args <- c(
    "stats", "entrapment-fdp",
    "--targets", readers_normalise_path(targets, "targets"),
    "--entrapments", readers_normalise_path(entrapments, "entrapments"),
    "--ratio", stats_number(ratio, "ratio")
  )
  if (!is.null(thresholds)) {
    if (!is.numeric(thresholds) || length(thresholds) == 0L || anyNA(thresholds)) {
      stop(mzlib_usage_error(
        "thresholds takes a non-empty numeric vector of q-values, with no NA, or NULL for mzLib's default."
      ))
    }
    args <- c(args, "--thresholds", paste(stats_number_text(thresholds), collapse = ","))
  }
  args
}

# ---------------------------------------------------------------- parsing

stats_parse_logistic <- function(data) {
  limit <- data[["max_iterations"]]
  tolerance <- data[["tolerance"]]
  structure(
    list(
      responses_file = wire_field(data, "responses_file", "character", NA_character_),
      design_file = wire_field(data, "design_file", "character", NA_character_),
      feature_count = stats_count(data, "feature_count"),
      sample_count = stats_count(data, "sample_count"),
      sample_names = wire_strings(data[["sample_names"]]),
      coefficient_names = wire_strings(data[["coefficient_names"]]),
      tested_coefficients = wire_strings(data[["tested_coefficients"]]),
      max_iterations = if (wire_null(limit)) NA_real_ else as.numeric(limit),
      tolerance = if (wire_null(tolerance)) NA_real_ else as.numeric(tolerance),
      threads = stats_count(data, "threads"),
      missing_count = stats_count(data, "missing_count"),
      status_counts = proteins_parse_counts(data[["status_counts"]]),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_logistic_model_fit"
  )
}

stats_parse_mixed <- function(data) {
  structure(
    list(
      responses_file = wire_field(data, "responses_file", "character", NA_character_),
      design_file = wire_field(data, "design_file", "character", NA_character_),
      groups_file = wire_field(data, "groups_file", "character", NA_character_),
      grouping = wire_field(data, "grouping", "character", NA_character_),
      group_names = wire_strings(data[["group_names"]]),
      group_count = stats_count(data, "group_count"),
      feature_count = stats_count(data, "feature_count"),
      sample_count = stats_count(data, "sample_count"),
      sample_names = wire_strings(data[["sample_names"]]),
      coefficient_names = wire_strings(data[["coefficient_names"]]),
      tested_coefficients = wire_strings(data[["tested_coefficients"]]),
      estimator = wire_field(data, "estimator", "character", NA_character_),
      threads = stats_count(data, "threads"),
      missing_count = stats_count(data, "missing_count"),
      zero_count = stats_count(data, "zero_count"),
      status_counts = proteins_parse_counts(data[["status_counts"]]),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_mixed_model_fit"
  )
}

stats_parse_entrapment <- function(data) {
  thresholds <- data[["thresholds"]]
  structure(
    list(
      targets_file = wire_field(data, "targets_file", "character", NA_character_),
      entrapments_file = wire_field(data, "entrapments_file", "character", NA_character_),
      ratio = stats_count(data, "ratio"),
      thresholds = if (is.list(thresholds)) vapply(thresholds, as.numeric, numeric(1L)) else numeric(0),
      thresholds_given = isTRUE(data[["thresholds_given"]]),
      target_candidate_count = stats_count(data, "target_candidate_count"),
      entrapment_candidate_count = stats_count(data, "entrapment_candidate_count"),
      foreign_count = stats_count(data, "foreign_count"),
      partnered_count = stats_count(data, "partnered_count"),
      threshold_count = stats_count(data, "threshold_count"),
      empty_threshold_count = stats_count(data, "empty_threshold_count"),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_entrapment_fdp_estimates"
  )
}

# ---------------------------------------------------------------- the public surface

#' Fit one logistic regression per feature of a 0/1 table, such as a detection table
#'
#' mzLib's `LogisticRegression.Fit` models logit P(outcome = 1) = X beta for every feature (row)
#' separately, by maximum likelihood, and tests each coefficient with a two-sided Wald z. It is R's
#' `glm(y ~ 0 + X, family = binomial)` per feature, which mzLib reproduces to 1e-8 relative in the
#' coefficients (1e-6 in standard errors and p-values).
#'
#' The usual question in proteomics: **is a feature's chance of being measured related to the
#' condition?** Write 1 where a protein was quantified and 0 where it was not; a coefficient far
#' from 0 says its missing values are not missing at random.
#'
#' A feature is reported - never dropped - when it cannot be fitted. `"separated"` is the one to
#' read carefully: the design splits its 0s from its 1s (or every outcome is the same), so no finite
#' estimate exists. R's `glm()` prints a huge coefficient and a warning there; mzLib reports the
#' status and `NA`. Look at `events` and `observed` for those rows.
#'
#' @param responses A TSV with a header row: the feature id first, then one column per sample, named
#'   by its header - [stats_fit()]'s layout - with every cell 0, 1, or blank / `NA` / `NaN` for
#'   missing (left out of that feature's fit). Any other value is refused.
#' @param design [stats_fit()]'s design: a TSV with a header row and one row per sample, the sample
#'   name (matching the responses header, in any order) then one numeric column per coefficient.
#'   Include an intercept column yourself if you want one.
#' @param coefficients Design column names to report, each as one block of rows. Every design column
#'   is fitted whatever you name.
#' @param max_iterations The IRLS iteration limit per feature (iterations). Leave it `NULL` for
#'   mzLib's default (100 at mzLib `d10017ce`); a feature still changing at the limit is
#'   `"not_converged"`.
#' @param tolerance Stop when the relative change in deviance is below this (R `glm()`'s rule).
#'   Leave it `NULL` for mzLib's default (1e-12 at mzLib `d10017ce`), which is tighter than
#'   `glm()`'s 1e-8.
#' @param threads Features fitted at once (threads), or `-1` for every core. The answer is identical
#'   at any value.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_logistic_model_fit`:
#'
#'   - `records`, a data.frame with one row per (reported coefficient, feature) - one block of
#'     `feature_count` rows per coefficient you named, in that order, each in the responses file's
#'     feature order: `feature`; `coefficient`; `status`; `estimate`, the maximum-likelihood
#'     coefficient in log-odds per coefficient unit, so `exp(estimate)` is the **odds ratio**;
#'     `standard_error` (log-odds per coefficient unit), from the inverse Fisher information at the
#'     estimate; `z` (standard errors), the Wald statistic; `p_value`, a fraction (0 to 1),
#'     two-sided and **not adjusted across features**; `deviance` (deviance, -2 log L), the
#'     feature's residual deviance; `iterations`, IRLS iterations used (0 when no fit was
#'     attempted); `observed`, samples with a 0 or 1; and `events`, observed samples whose outcome
#'     is 1. Statistics are `NA` where `status` is not `"fitted"`.
#'   - `responses_file` and `design_file`, as absolute paths; `feature_count` features and
#'     `sample_count` samples read; `sample_names`; `coefficient_names` (every design column) and
#'     `tested_coefficients`; `max_iterations` (iterations) and `tolerance` (relative change in
#'     deviance) as you passed them, `NA` when left to mzLib; `threads`; `missing_count` response
#'     cells (cells); `status_counts`, features per status; `row_count` rows; `column_names`.
#'   - `caveats`: what the numbers are and are not, with a count of unfitted features by status.
#'
#' @section What it is not:
#'
#' Wald tests are the least reliable of the likelihood tests near separation and in small samples;
#' no likelihood-ratio test is offered, and Firth's penalized likelihood (logistf), which gives
#' finite estimates for separated features, is not implemented. Send one coefficient's `p_value`
#' column to [stats_adjust()] for Benjamini-Hochberg.
#'
#' @seealso [stats_fit()] for the intensities themselves, [stats_adjust()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own reference input (80 simulated binary features, 40 samples), and R's glm() output
#' # for it, which this package ships:
#' data <- system.file("extdata", "stats", package = "mzLibR")
#' fit <- stats_fit_logistic(file.path(data, "glm_reference_responses.tsv"),
#'   file.path(data, "glm_reference_design.tsv"), c("Intercept", "x", "grp"))
#' fit
#' x <- fit$records[fit$records$coefficient == "x", ]
#' glm <- read.delim(file.path(data, "glm_reference_fit.tsv"))
#' fitted <- glm$status == "fitted"
#' max(abs(x$estimate[fitted] - glm$b1[fitted]) / abs(glm$b1[fitted]))
#' x[1, c("estimate", "p_value", "events", "observed")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.fit-logistic
#' @export
stats_fit_logistic <- function(responses, design, coefficients, max_iterations = NULL, tolerance = NULL,
                               threads = 1, timeout = 600) {
  request <- stats_build_logistic_request(responses, design, coefficients, max_iterations, tolerance, threads)
  bridge_require_verb("stats fit-logistic", STATS_MODELS_SINCE)
  stats_parse_logistic(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Fit one linear mixed model per feature, with a random intercept for each group of samples
#'
#' mzLib's `MixedModel.Fit` models y = X beta + u_group + e for every feature separately: the fixed
#' effects in your design, plus one baseline per group (a dataset, a donor, an acquisition day)
#' drawn from N(0, tau^2). Each coefficient is tested by a two-sided Wald t on nlme's "containment"
#' degrees of freedom. It is `nlme::lme(y ~ ..., random = ~ 1 | group)` per feature, which mzLib
#' reproduces to 1e-5 relative (nlme's own precision), with identical df.
#'
#' Use it when samples come in groups that share a baseline - runs acquired on the same day,
#' replicates of one donor, studies in a pooled analysis - and you want each coefficient's effect
#' *within* groups. [stats_fit()] treats every sample as independent, so a day-to-day shift there
#' either inflates the residual variance or, if the days are unbalanced, leaks into the effect.
#'
#' @param responses [stats_fit()]'s responses table: a TSV with a header row, the feature id first,
#'   then one column per sample, named by its header. Log-transform intensities first. Blank, `NA`
#'   or `NaN` is missing; **0 is an observation**.
#' @param design [stats_fit()]'s design: one row per sample, the sample name then one numeric column
#'   per fixed-effect coefficient. Do not add a column per group: the random intercept replaces it.
#' @param groups A TSV with a header row and one row per sample: the sample name (exactly the
#'   responses' samples, in any order), then ONE column with the sample's group label, whose header
#'   names the grouping (`day`, `donor`, `dataset`).
#' @param coefficients Design column names to report, each as one block of rows. Every design column
#'   is fitted whatever you name.
#' @param estimator `"reml"` or `"ml"`: how the variances are estimated. Leave it `NULL` for mzLib's
#'   default (`"reml"` at mzLib `d10017ce`); `estimator` in the result says which was used. `"ml"`
#'   gives log-likelihoods you can compare between designs.
#' @param threads Features fitted at once (threads), or `-1` for every core. The answer is identical
#'   at any value.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_mixed_model_fit`:
#'
#'   - `records`, a data.frame with one row per (reported coefficient, feature) - one block of
#'     `feature_count` rows per coefficient you named, in that order: `feature`; `coefficient`;
#'     `status`; `estimate`, the fixed effect within groups in response units per coefficient unit
#'     (a **log2 fold change** for a 0/1 column on log2 intensities); `standard_error` (response
#'     units per coefficient unit); `df` (degrees of freedom), nlme's containment df, a whole
#'     number, `NA` when not fitted or when no positive df exists (a group-level coefficient with
#'     too few groups); `t` (standard errors), the Wald statistic; `p_value`, a fraction (0 to 1),
#'     two-sided on `df`, **not adjusted across features**, `NA` also when `df` is; `residual_variance`
#'     (squared response units), sigma^2 within groups; `group_variance` (squared response units),
#'     tau^2 between the groups' baselines, where 0 is the boundary (no group effect);
#'     `log_likelihood` (log-likelihood), restricted under `"reml"`; `observed`, samples with a
#'     finite response; and `groups`, distinct groups among them. Statistics are `NA` where
#'     `status` is not `"fitted"`.
#'   - `responses_file`, `design_file` and `groups_file`, as absolute paths; `grouping`, the groups
#'     file's label column header; `group_names`, in order of first appearance; `group_count`
#'     groups; `feature_count` features and `sample_count` samples read; `sample_names`;
#'     `coefficient_names` and `tested_coefficients`; `estimator`, the estimator mzLib used;
#'     `threads`; `missing_count` and `zero_count` response cells (cells); `status_counts`,
#'     features per status; `row_count` rows; `column_names`.
#'   - `caveats`: what the model is and is not, plus counts for this input: unfitted features, rows
#'     with no df, features at the boundary, zero cells.
#'
#' @section What it is not:
#'
#' A random intercept per group and nothing more: no random slopes, nested or crossed groupings, or
#' residual correlation. Satterthwaite and Kenward-Roger df (lmerTest, pbkrtest) are not
#' implemented, so p-values differ from lmerTest's. A feature the fixed effects fit exactly is
#' reported fitted, with a `residual_variance` at rounding level and a meaningless `t` and
#' `p_value`: check `residual_variance` before trusting a very small `p_value`.
#'
#' @seealso [stats_fit()] when samples are independent, [stats_adjust()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own reference input (120 simulated features, 29 samples in 7 groups of unequal size;
#' # `wp` is constant within each group, `age` varies within them), and nlme's output for it:
#' data <- system.file("extdata", "stats", package = "mzLibR")
#' fit <- stats_fit_mixed(file.path(data, "lmm_reference_responses.tsv"),
#'   file.path(data, "lmm_reference_design.tsv"), file.path(data, "lmm_reference_groups.tsv"),
#'   c("Intercept", "age", "wp"))
#' fit
#' # `wp` is a group-level covariate, so it is tested on the groups, not the samples:
#' table(fit$records$df[fit$records$coefficient == "wp"])
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.fit-mixed
#' @export
stats_fit_mixed <- function(responses, design, groups, coefficients, estimator = NULL, threads = 1,
                            timeout = 600) {
  request <- stats_build_mixed_request(responses, design, groups, coefficients, estimator, threads)
  bridge_require_verb("stats fit-mixed", STATS_MODELS_SINCE)
  stats_parse_mixed(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Estimate a search's false discovery proportion from its entrapment hits
#'
#' An entrapment experiment adds peptides that cannot be in the sample (from another species, or
#' made from the target peptides) to the search database. Every entrapment the search reports is
#' then a known false discovery, which lets you check the search's FDR claim instead of trusting it.
#' mzLib's `EntrapmentFdp.Sweep` computes the three valid estimators of Wen et al. 2025 (the method
#' FDRBench implements) at each q-value threshold:
#'
#' - `lower_bound` = N_E / (N_T + N_E). Above the FDR you claimed, the search **fails** to control
#'   it. Below, it proves nothing.
#' - `combined` = N_E (1 + 1/r) / (N_T + N_E). An upper bound: at or below the claimed FDR, it is
#'   **evidence of control**.
#' - `paired`: a tighter upper bound when each target has exactly one entrapment made from it
#'   (`ratio = 1`), using whether each discovered entrapment's own target was discovered too.
#'
#' The "sample" estimator N_E / (r N_T) is invalid in both directions and is deliberately absent.
#'
#' @param targets A TSV whose header names exactly `score` and `q_value`, in any order: one row per
#'   original-target candidate (precursor, peptide or protein group, whichever you are estimating
#'   for), decoys removed, each with its best q-value. Higher scores are better.
#' @param entrapments A TSV whose header names `score`, `q_value`, `partner_score` and
#'   `partner_q_value`, and optionally `foreign`: one row per entrapment candidate. The partner
#'   columns copy the score and q-value of the one target this entrapment was made from, when the
#'   search identified it; leave both blank otherwise. `foreign` is 1 for an entrapment with no
#'   partner by construction (another species); blank or absent is mzLib's default, not foreign.
#'   Other columns are refused, so select these first.
#' @param ratio r, the entrapment-to-target database size ratio (unitless). Required: mzLib has no
#'   default. `paired` is reported only when it is exactly 1, so for 1:1 pairing pass 1, not a ratio
#'   computed from database sizes after exclusions.
#' @param thresholds The q-value thresholds to estimate at, as fractions (0.01 for 1%). Leave it
#'   `NULL` for mzLib's default sweep (0.001 to 0.100 in steps of 0.001 at mzLib `d10017ce`),
#'   echoed in `thresholds`.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_entrapment_fdp_estimates`:
#'
#'   - `records`, a data.frame with one row per q-value threshold with a discovery:
#'     `q_value_threshold` (a q-value, fraction), where a candidate is discovered at or below it;
#'     `target_count`, N_T, targets discovered (discoveries); `entrapment_count`, N_E, entrapments
#'     discovered (discoveries); and three estimates of the false discovery proportion:
#'     `lower_bound` (fraction of discoveries), a **lower bound** (above your claimed FDR, the
#'     search fails to control it); `combined` (fraction of discoveries), an **upper bound** (at or
#'     below your claimed FDR, evidence of control; can exceed 1); and `paired` (fraction of
#'     discoveries), a tighter upper bound for 1:1 pairing (can exceed 1), `NA` where it does not
#'     apply: a ratio not exactly 1, or every entrapment foreign.
#'   - `targets_file` and `entrapments_file`, as absolute paths; `ratio`, r as you gave it (entrapment-to-target database size ratio,
#'     unitless);
#'     `thresholds`, the q-value thresholds estimated at, in order (yours or mzLib's default);
#'     `thresholds_given`; `target_candidate_count` and `entrapment_candidate_count`, the rows of the
#'     two files (candidates); `foreign_count` and `partnered_count` entrapments (candidates);
#'     `threshold_count` thresholds; `empty_threshold_count`, thresholds with no discoveries
#'     (thresholds), which have no row; `row_count` rows; `column_names`.
#'   - `caveats`: what each estimate can and cannot show, plus any that apply to this input.
#'
#' @section The input traps:
#'
#' Give one row per unit estimated, decoys removed, each with its best q-value: mzLib does not
#' deduplicate, and it does not check a partner's copied score against the targets file. The k-fold
#' paired estimator for k entrapments per target is not implemented.
#'
#' @seealso [stats_adjust()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own hand-built search (seven targets, four entrapments, one for each case the paired
#' # estimate distinguishes), checked against the entrapment project's reference script:
#' data <- system.file("extdata", "stats", package = "mzLibR")
#' fdp <- stats_entrapment_fdp(file.path(data, "entrapment_reference_targets.tsv"),
#'   file.path(data, "entrapment_reference_entrapments.tsv"), ratio = 1)
#' fdp
#' fdp$records[fdp$records$q_value_threshold == 0.01, ]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.entrapment-fdp
#' @export
stats_entrapment_fdp <- function(targets, entrapments, ratio, thresholds = NULL, timeout = 600) {
  if (missing(ratio)) {
    stop(mzlib_usage_error("ratio is required: mzLib has no default for the entrapment-to-target database size ratio."))
  }
  args <- stats_build_entrapment_args(targets, entrapments, ratio, thresholds)
  bridge_require_verb("stats entrapment-fdp", STATS_MODELS_SINCE)
  stats_parse_entrapment(bridge_invoke(args, timeout = timeout))
}

# ---------------------------------------------------------------- printing

#' Print a logistic fit
#'
#' @param x A [stats_fit_logistic()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_logistic_model_fit <- function(x, ...) {
  cat("<mzlibr_logistic_model_fit> ", format(x$feature_count), " features x ", format(x$sample_count),
    " samples, reporting ", paste(x$tested_coefficients, collapse = ", "), "\n",
    sep = ""
  )
  if (length(x$status_counts) > 0L) {
    cat("  ", paste0(names(x$status_counts), " ", x$status_counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Print a mixed-model fit
#'
#' @param x A [stats_fit_mixed()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_mixed_model_fit <- function(x, ...) {
  cat("<mzlibr_mixed_model_fit> ", format(x$feature_count), " features x ", format(x$sample_count),
    " samples in ", format(x$group_count), " groups (", x$grouping, "), ", x$estimator, ", reporting ",
    paste(x$tested_coefficients, collapse = ", "), "\n",
    sep = ""
  )
  if (length(x$status_counts) > 0L) {
    cat("  ", paste0(names(x$status_counts), " ", x$status_counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Print entrapment FDP estimates
#'
#' @param x A [stats_entrapment_fdp()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_entrapment_fdp_estimates <- function(x, ...) {
  cat("<mzlibr_entrapment_fdp_estimates> ", format(x$target_candidate_count), " targets, ",
    format(x$entrapment_candidate_count), " entrapments, ratio ", format(x$ratio), ", ",
    format(x$row_count), " of ", format(x$threshold_count), " thresholds with a discovery\n",
    sep = ""
  )
  invisible(x)
}
