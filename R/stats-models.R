# Per-feature models beyond limma's, and an entrapment check of a search's FDR, computed by mzLib.
#
#   stats_fit_logistic()    LogisticRegression.Fit: glm(family = binomial) per feature of a 0/1 table
#   stats_fit_mixed()       MixedModel.Fit: nlme::lme(random = ~ 1 | group) per feature
#   stats_fit_latent_mixed()  LogitNormalMixedModel.Fit: the same model on a latent scale, from
#                           measured and counted (k of n) samples, lme4::glmer and lmer in one
#   stats_fit_link()        LogitNormalMixedModel.FitLink: the detection link it reads counts through
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

# ---------------------------------------------------------------- the latent scale (pyMzLib 0.7.0)

# Ported from pyMzLib's `stats.py` (fit_latent_mixed, fit_link, LatentMixedModelFit, LinkFit,
# DetectionLink; mzLib #1469). A verb of its own, not an option of stats_fit_mixed(): mzLib has two
# entry points with different optimisers, statuses and ML standard errors, so the caller picks the
# model by picking the function, and stats_fit_mixed()'s answer does not move.

# The first pyMzLib whose bridge has `stats fit-latent-mixed` and `stats fit-link` (mzLib 1.0.596).
STATS_LATENT_SINCE <- "0.7.0"

# Every value of stats_fit_latent_mixed()'s `status` column: the outcomes mzLib's
# LogitNormalMixedModel.Fit can give a feature, read from its code (mzLib does not list them per model
# yet, mzLib#1487). It never reports "variance_ratio_at_limit". Only "fitted" rows carry statistics.
STATS_LATENT_MIXED_STATUSES <- c(
  "fitted", "too_few_observations", "rank_deficient", "separated", "not_converged", "too_few_groups"
)

# The fields of mzLib's DetectionLink, in its positional order.
STATS_LINK_FIELDS <- c("intercept", "slope", "depth_slope", "overdispersion")

# `--link`: "identity", or the link's numbers comma-separated - only those given, so a depth slope or
# overdispersion left out is mzLib's own default.
stats_link_argument <- function(link) {
  if (inherits(link, "mzlibr_link_fit")) {
    if (is.null(link$link)) {
      stop(mzlib_usage_error("link is a stats_fit_link() result whose search did not converge (link is NULL)."))
    }
    link <- link$link
  }
  if (inherits(link, "mzlibr_detection_link")) {
    return(paste(stats_number_text(vapply(STATS_LINK_FIELDS, function(f) as.numeric(link[[f]]), numeric(1L))),
      collapse = ","))
  }
  if (is.character(link) && length(link) == 1L && !is.na(link)) {
    if (identical(tolower(trimws(link)), "identity")) {
      return("identity")
    }
    stop(mzlib_usage_error(paste0(
      "link takes \"identity\", a stats_fit_link() result (or its $link), or c(intercept, slope[, ",
      "depth_slope[, overdispersion]]); got ", paste(deparse(link), collapse = " "), "."
    )))
  }
  if (!is.numeric(link) || length(link) < 2L || length(link) > 4L) {
    stop(mzlib_usage_error(paste0(
      "link takes \"identity\", a stats_fit_link() result (or its $link), or 2 to 4 numbers ",
      "c(intercept, slope[, depth_slope[, overdispersion]]); got ", paste(deparse(link), collapse = " "), "."
    )))
  }
  for (i in seq_along(link)) {
    if (is.na(link[[i]])) {
      stop(mzlib_usage_error(paste0("link[", i, "] must be a number; got NA.")))
    }
  }
  paste(stats_number_text(unname(link)), collapse = ",")
}

# `--<name> true|false` when given; nothing when NULL, so mzLib's declared default applies.
stats_bool_option <- function(value, name) {
  if (is.null(value)) {
    return(character(0))
  }
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(mzlib_usage_error(paste0(
      gsub("-", "_", name, fixed = TRUE), " must be TRUE, FALSE or NULL; got ", paste(deparse(value), collapse = " "), "."
    )))
  }
  c(paste0("--", name), if (value) "true" else "false")
}

stats_build_latent_request <- function(design, groups, coefficients, link, responses, counts, depth,
                                       residual_variance, estimator, threads) {
  coefficients <- stats_coefficient_names(coefficients, caller = "stats_fit_latent_mixed")
  if (is.null(responses) && is.null(counts)) {
    stop(mzlib_usage_error(
      "stats_fit_latent_mixed() needs responses = (measured samples), counts = (counted samples), or both."
    ))
  }
  if (is.null(link)) {
    stop(mzlib_usage_error(paste0(
      "link is required: mzLib's DetectionLink has no default. Give \"identity\", c(intercept, slope[, ",
      "depth_slope[, overdispersion]]), or a stats_fit_link() result."
    )))
  }
  args <- c(
    "stats", "fit-latent-mixed",
    "--design", readers_normalise_path(design, "design"),
    "--groups", readers_normalise_path(groups, "groups"),
    "--link", stats_link_argument(link),
    "--threads", stats_whole(threads, "threads")
  )
  # Each optional file is sent only when given.
  files <- list(responses = responses, counts = counts, depth = depth, "residual-variance" = residual_variance)
  for (name in names(files)) {
    if (!is.null(files[[name]])) {
      args <- c(args, paste0("--", name), readers_normalise_path(files[[name]], gsub("-", "_", name, fixed = TRUE)))
    }
  }
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

stats_build_link_args <- function(cells, with_depth, with_overdispersion, detection_only) {
  c(
    "stats", "fit-link", "--cells", readers_normalise_path(cells, "cells"),
    stats_bool_option(with_depth, "with-depth"),
    stats_bool_option(with_overdispersion, "with-overdispersion"),
    stats_bool_option(detection_only, "detection-only")
  )
}

# A wire DetectionLink as an `mzlibr_detection_link`, or NULL for a null link.
stats_parse_link <- function(value) {
  if (!is.list(value) || is.null(names(value))) {
    return(NULL)
  }
  link <- lapply(STATS_LINK_FIELDS, function(f) stats_count(value, f))
  names(link) <- STATS_LINK_FIELDS
  structure(link, class = "mzlibr_detection_link")
}

stats_parse_latent <- function(data) {
  structure(
    list(
      responses_file = wire_field(data, "responses_file", "character", NA_character_),
      counts_file = wire_field(data, "counts_file", "character", NA_character_),
      design_file = wire_field(data, "design_file", "character", NA_character_),
      groups_file = wire_field(data, "groups_file", "character", NA_character_),
      depth_file = wire_field(data, "depth_file", "character", NA_character_),
      residual_variance_file = wire_field(data, "residual_variance_file", "character", NA_character_),
      grouping = wire_field(data, "grouping", "character", NA_character_),
      group_names = wire_strings(data[["group_names"]]),
      group_count = stats_count(data, "group_count"),
      feature_count = stats_count(data, "feature_count"),
      sample_count = stats_count(data, "sample_count"),
      sample_names = wire_strings(data[["sample_names"]]),
      coefficient_names = wire_strings(data[["coefficient_names"]]),
      tested_coefficients = wire_strings(data[["tested_coefficients"]]),
      estimator = wire_field(data, "estimator", "character", NA_character_),
      link = stats_parse_link(data[["link"]]),
      threads = stats_count(data, "threads"),
      measured_cells = stats_count(data, "measured_cells"),
      counted_cells = stats_count(data, "counted_cells"),
      missing_count = stats_count(data, "missing_count"),
      status_counts = proteins_parse_counts(data[["status_counts"]]),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_latent_mixed_model_fit"
  )
}

stats_parse_link_fit <- function(data) {
  structure(
    list(
      cells_file = wire_field(data, "cells_file", "character", NA_character_),
      cell_count = stats_count(data, "cell_count"),
      depth_given = isTRUE(data[["depth_given"]]),
      with_depth = isTRUE(data[["with_depth"]]),
      with_overdispersion = isTRUE(data[["with_overdispersion"]]),
      detection_only = isTRUE(data[["detection_only"]]),
      link = stats_parse_link(data[["link"]]),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_link_fit"
  )
}

#' Fit a mixed model per feature on a latent scale, from samples measured on it, counted through it, or both
#'
#' Answers **does a value follow a trait, when some runs measured it and others only counted
#' successes, k of n?** A PTM site's occupancy is the usual case: where the modified and unmodified
#' peptides were both quantified you have a number; where they were not, you may know only that 3 of
#' 8 covering peptidoforms carried the modification. [stats_fit_mixed()] has to drop the counted
#' runs. This uses both kinds.
#'
#' mzLib's `LogitNormalMixedModel.Fit` (mzLib #1469) is [stats_fit_mixed()]'s model - fixed effects
#' from your design, one baseline per group drawn from N(0, tau^2), run-to-run scatter sigma^2 - on a
#' **latent** value, such as the log-odds of a site's occupancy. A sample can be **measured** (its
#' latent value is in `responses`) or **counted**: k of its n trials succeeded (`counts`), each with
#' a chance g set by the latent value through the `link`, logit g = intercept + slope x latent +
#' depth_slope x depth. It is `lme4::glmer` (a binomial GLMM, with an observation-level random effect
#' for sigma^2) and `lmer` in one model, which mzLib checks against both and against `glmmTMB`.
#'
#' **Which cells are counted, and the link, are your choices**: mzLibR builds neither.
#' [stats_fit_link()] fits a link from cells.
#'
#' @param design [stats_fit()]'s design: a TSV with a header row and one row per sample, the sample
#'   name then one numeric column per fixed-effect coefficient. Do not add a column per group: the
#'   random intercept replaces it. mzLib refuses a design that is not of full column rank.
#' @param groups [stats_fit_mixed()]'s groups table: one row per sample, the sample name then ONE
#'   column with its group label (an animal, a donor, a dataset), whose header names the grouping.
#' @param coefficients Design column names to report, each as one block of rows, sent one per stdin
#'   line. Every design column is fitted whatever you name.
#' @param link mzLib's `DetectionLink`, in natural log-odds: `"identity"` (`DetectionLink.Identity`,
#'   intercept 0 and slope 1: a count samples the same probability a measurement reads), a
#'   [stats_fit_link()] result or its `$link`, or 2 to 4 numbers `c(intercept, slope[, depth_slope[,
#'   overdispersion]])`, sent as given: a value you leave out is mzLib's default (0 at the pin).
#'   Required: mzLib has no default link.
#' @param responses The measured samples, in latent units (logit occupancy, for example, not a
#'   fraction): [stats_fit_mixed()]'s responses table, a TSV with the feature id first and one column
#'   per sample. Blank, `NA` or `NaN` is not measured. Give `responses`, `counts` or both. `NULL`
#'   (the default): nothing was measured.
#' @param counts The counted samples, in trials: a table laid out as `responses`, each cell
#'   `successes/trials` (`3/8`), or blank where the sample was not counted. With `responses` too,
#'   the two must list the same features and samples, and no cell may be in both. `NULL` (the
#'   default): nothing was counted, and the model is [stats_fit_mixed()]'s.
#' @param depth A TSV with a header row and one row per sample: the sample name, then its depth in
#'   your depth units (log2 target PSMs, say), for the link's depth term. Needed when the link's
#'   `depth_slope` is not 0. `NULL` (the default): none is sent.
#' @param residual_variance A TSV with a header row and one row per feature: the feature id, then
#'   sigma^2 to use for it (squared latent units) instead of estimating it, for a feature with too
#'   few measured samples. Blank, or a feature not listed, is estimated. `NULL` (the default): every
#'   sigma^2 is estimated.
#' @param estimator `"reml"` or `"ml"`. Leave it `NULL` (the default) for mzLib's default (`"reml"`
#'   at the pin); `estimator` in the result says which was used.
#' @param threads Features fitted at once (threads), or `-1` for every core. The answer is identical
#'   at any value.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_latent_mixed_model_fit`:
#'
#'   - `records`, a data.frame with one row per (reported coefficient, feature) - one block of
#'     `feature_count` rows per coefficient you named, in that order, each in the feature order of
#'     the responses file (of the counts file, without one): `feature`; `coefficient`; `status`, one
#'     of `"fitted"`, `"too_few_observations"`, `"rank_deficient"`, `"separated"` (no finite
#'     estimate exists: a status, never a p-value), `"not_converged"` or `"too_few_groups"`;
#'     `estimate`, the fixed effect within groups on the latent scale (latent units per coefficient
#'     unit: **natural log-odds** for a logit latent, so `exp(estimate)` is an odds ratio);
#'     `standard_error` (latent units per coefficient unit), from the inverse Hessian at the
#'     estimated variances; `df` (degrees of freedom), nlme's containment df over the samples with a
#'     response of either kind, `NA` when not fitted or when no positive df exists; `t` (standard
#'     errors), the Wald statistic; `p_value`, a fraction (0 to 1), two-sided on `df`, **not
#'     adjusted across features**, `NA` also when `df` is; `residual_variance` (squared latent
#'     units), sigma^2, run to run within a group, estimated or the one you gave;
#'     `residual_variance_given`, `TRUE` when it is the one you gave; `group_variance` (squared
#'     latent units), tau^2, where 0 is the boundary; `log_likelihood` (log-likelihood, natural log),
#'     Laplace-approximated restricted under `"reml"`, full under `"ml"`; `measured`, samples with a
#'     finite response (samples); `counted`, samples with a count (samples); and `groups`, distinct
#'     groups among both (groups). Statistics are `NA` where `status` is not `"fitted"`.
#'   - `responses_file`, `counts_file`, `depth_file` and `residual_variance_file` (`NA` when not
#'     given), `design_file` and `groups_file`, as absolute paths; `grouping`, the groups file's
#'     label column header; `group_names`, in order of first appearance; `group_count` groups;
#'     `feature_count` features and `sample_count` samples read; `sample_names`;
#'     `coefficient_names` and `tested_coefficients`; `estimator`, the estimator mzLib used;
#'     `link`, the link every count was read through, as mzLib reports it: an
#'     `mzlibr_detection_link` with `intercept` (natural log-odds), `slope` (natural log-odds per
#'     latent unit), `depth_slope` (natural log-odds per depth unit) and `overdispersion` (fraction,
#'     the intra-cell correlation, 0 to below 1); `threads` (threads); `measured_cells`,
#'     `counted_cells` and `missing_count` (cells); `status_counts`, features per status
#'     (features); `row_count` rows; `column_names`.
#'   - `caveats`: what the model is and is not, the link and the estimator used, plus counts for this
#'     input: unfitted features, rows with no df, features at the boundary, residual variances given.
#'
#' @section What it is not:
#'
#' A random intercept per group and nothing more: no random slopes, nested or crossed groupings.
#' Responses and estimates are on the latent scale, never fractions or intensities. Standard errors
#' are at the estimated variances and do not include their uncertainty, so they are not `glmmTMB`'s
#' `vcov()`. Under `"reml"` the fixed effects are integrated out by a Laplace approximation, as
#' `glmmTMB`'s `REML = TRUE` does; with nothing counted it is exactly [stats_fit_mixed()]'s REML,
#' and the two agree to about 1e-5. Send one coefficient's `p_value` column to [stats_adjust()] for
#' Benjamini-Hochberg.
#'
#' @seealso [stats_fit_link()] for the link, [stats_fit_mixed()] when every sample is measured
#' @spec stats.fit-latent-mixed
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own reference input: 16 simulated features, 23 samples in 8 animals, every sample
#' # counted, read through the link (-0.4, 0.6, 0.3) with a depth term; and lme4's glmer() fit of it
#' # (25-point adaptive quadrature, maximum likelihood), which this package ships.
#' data <- system.file("extdata", "stats", package = "mzLibR")
#' fit <- stats_fit_latent_mixed(file.path(data, "lnmm_reference_design.tsv"),
#'   file.path(data, "lnmm_reference_groups.tsv"), c("Intercept", "trait", "x"),
#'   link = c(-0.4, 0.6, 0.3), counts = file.path(data, "lnmm_link_counts.tsv"),
#'   depth = file.path(data, "lnmm_link_depth.tsv"),
#'   residual_variance = file.path(data, "lnmm_link_residual_variance.tsv"), estimator = "ml")
#' fit
#' trait <- fit$records[fit$records$coefficient == "trait", ]
#' glmer <- read.delim(file.path(data, "lnmm_counted_fit.tsv"))
#' glmer <- glmer[glmer$link == "link" & glmer$method == "ML", ]
#' # Every trait estimate is within 1e-5 standard errors of glmer's:
#' max(abs(trait$estimate - glmer$b1) / glmer$se1) < 1e-5
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
stats_fit_latent_mixed <- function(design, groups, coefficients, link, responses = NULL, counts = NULL,
                                   depth = NULL, residual_variance = NULL, estimator = NULL, threads = 1,
                                   timeout = 600) {
  request <- stats_build_latent_request(design, groups, coefficients, if (missing(link)) NULL else link,
    responses, counts, depth, residual_variance, estimator, threads)
  bridge_require_verb("stats fit-latent-mixed", STATS_LATENT_SINCE)
  stats_parse_latent(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Fit the detection link: how a cell's latent value sets the chance that a trial succeeds
#'
#' Answers **how does a count's chance of success follow the latent value?** mzLib's
#' `LogitNormalMixedModel.FitLink` (mzLib #1469) finds the `DetectionLink` - logit g = intercept +
#' slope x latent + depth_slope x depth - that makes your cells' counts most likely, each cell's
#' latent value taken as a normal with the mean and SD you give, by maximum likelihood. By default
#' it fits on **detection**: whether any trial succeeded, which is what makes a covered zero. Pass
#' the result to [stats_fit_latent_mixed()] as `link`.
#'
#' A slope below 1 says the chance of a success rises more slowly than the latent: on PTM data, for
#' example, whether a modified form is seen depends on how much of it there is, not only on its
#' share. **The cells, and their latent means and SDs, are your choice**: mzLibR does not build them.
#' mzLib's own example is a covered cell's feature mean and SD over its measured cells.
#'
#' @param cells A TSV whose header names `latent_mean`, `latent_sd`, `successes` and `trials`, and
#'   optionally `depth`, in any order and no other column: one row per covered cell. The latent
#'   mean and SD are in latent units, successes of trials is the cell's count, and depth is in your
#'   depth units. Without a depth column every depth is mzLib's default (0 at the pin).
#' @param with_depth Fit the depth slope (mzLib's `withDepth`). `NULL` (the default) sends nothing:
#'   mzLib's declared default (`FALSE` at the pin), which `with_depth` in the result reports.
#' @param with_overdispersion Fit rho, the beta-binomial correlation of the trials in a cell (mzLib's
#'   `withOverdispersion`). `NULL` (the default): mzLib's declared default (`FALSE` at the pin).
#' @param detection_only `TRUE` fits on whether any trial succeeded, `FALSE` on the full counts
#'   (mzLib's `detectionOnly`). `NULL` (the default): mzLib's declared default (`TRUE` at the pin).
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_link_fit`: `cells_file`, the cells table's absolute path; `cell_count`, rows
#'   read (cells); `depth_given`, whether the cells file has a `depth` column; `with_depth`,
#'   `with_overdispersion` and `detection_only`, as passed: what you gave, else mzLib's declared
#'   default; `link`, the fitted `mzlibr_detection_link`, ready for [stats_fit_latent_mixed()]'s
#'   `link`, or `NULL` when mzLib's search did not converge - `intercept` (natural log-odds, the
#'   per-trial log-odds at latent 0 and depth 0), `slope` (natural log-odds per latent unit),
#'   `depth_slope` (natural log-odds per depth unit; exactly 0 when `with_depth` is `FALSE`) and
#'   `overdispersion` (fraction, the intra-cell correlation, 0 to below 1; exactly 0 when
#'   `with_overdispersion` is `FALSE`); and `caveats`, what the fit is and is not.
#'
#' @section What it is not:
#'
#' A point estimate only: mzLib returns no standard errors. The cells and their latent normals are
#' modelling choices, and the link fitted from them is only as good as they are. With
#' `with_depth = TRUE` and no depth column, every depth is 0 and the data say nothing about
#' `depth_slope`.
#'
#' @seealso [stats_fit_latent_mixed()]
#' @spec stats.fit-link
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own test cells: 4,000 cells simulated through the link (-0.5, 0.5, 0.3), which this
#' # package ships. The link is recovered from detection alone:
#' cells <- system.file("extdata", "stats", "link_cells.tsv", package = "mzLibR")
#' fit <- stats_fit_link(cells, with_depth = TRUE)
#' fit
#' round(unlist(fit$link), 2)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
stats_fit_link <- function(cells, with_depth = NULL, with_overdispersion = NULL, detection_only = NULL,
                           timeout = 600) {
  args <- stats_build_link_args(cells, with_depth, with_overdispersion, detection_only)
  bridge_require_verb("stats fit-link", STATS_LATENT_SINCE)
  stats_parse_link_fit(bridge_invoke(args, timeout = timeout))
}

# ---------------------------------------------------------------- printing

# One line for a detection link, as the print methods show it.
stats_link_text <- function(x) {
  paste0("logit g = ", format(x$intercept, digits = 4), " + ", format(x$slope, digits = 4), " x latent + ",
    format(x$depth_slope, digits = 4), " x depth, overdispersion ", format(x$overdispersion, digits = 4))
}

#' Print a latent-scale mixed-model fit
#'
#' @param x A [stats_fit_latent_mixed()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_latent_mixed_model_fit <- function(x, ...) {
  cat("<mzlibr_latent_mixed_model_fit> ", format(x$feature_count), " features x ", format(x$sample_count),
    " samples in ", format(x$group_count), " groups (", x$grouping, "), ", x$estimator, ", reporting ",
    paste(x$tested_coefficients, collapse = ", "), "\n",
    sep = ""
  )
  cat("  ", format(x$measured_cells), " measured and ", format(x$counted_cells), " counted cells, through ",
    stats_link_text(x$link), "\n",
    sep = ""
  )
  if (length(x$status_counts) > 0L) {
    cat("  ", paste0(names(x$status_counts), " ", x$status_counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Print a fitted detection link
#'
#' @param x A [stats_fit_link()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_link_fit <- function(x, ...) {
  cat("<mzlibr_link_fit> ", format(x$cell_count), " cells, fitted on ",
    if (x$detection_only) "detection (any success or none)" else "the full counts", "\n",
    sep = ""
  )
  cat("  ", if (is.null(x$link)) "no link: mzLib's search did not converge" else stats_link_text(x$link), "\n", sep = "")
  invisible(x)
}

#' Print a detection link
#'
#' @param x A detection link: `$link` of a [stats_fit_link()] or [stats_fit_latent_mixed()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_detection_link <- function(x, ...) {
  cat("<mzlibr_detection_link> ", stats_link_text(x), "\n", sep = "")
  invisible(x)
}

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
