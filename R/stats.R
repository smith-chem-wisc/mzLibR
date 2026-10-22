# Differential abundance, computed by mzLib: limma's moderated t-test, Benjamini-Hochberg, and
# random-effects meta-analysis.
#
#   stats_fit()     LinearModel.Fit + EmpiricalBayes.Moderate: limma's lmFit then
#                   eBayes(legacy = TRUE), one moderated test per named coefficient; with a
#                   contrasts file, EmpiricalBayes.ModerateContrasts as well (limma's contrasts.fit
#                   and topTable(confint = ...)); and either of mzLib's variance-prior estimators
#   stats_adjust()  MultipleTesting.BenjaminiHochberg over p-values from anywhere
#   stats_meta()    RandomEffectsMeta.Pool: DerSimonian-Laird random-effects pooling per feature
#
# Ported from pyMzLib's `stats.py` (mzLib #1341, #1357, #1368, #1452). Every number is mzLib's: this
# file shapes arguments into lines on stdin and projects each result into a data.frame, `records`.
# The tests hold the recorded results to limma, metafor and stats::p.adjust() when those are
# installed. The other stats verbs are in stats-combine.R and stats-models.R.
#
# R has limma and metafor already. These functions exist so that an R analysis can state that its
# statistics are the ones MetaMorpheus and the other bindings report, from the same code, and so a
# pipeline that runs in Python, Rust and R gives the same moderated t in all three.

# The pyMzLib release whose bridge first dispatches the stats verbs (their specs' since.pymzlib).
STATS_SINCE <- "0.3.0"

# Every value of the `status` column, as mzLib's FeatureFitStatus is written on the wire. Only
# "fitted" rows carry statistics or enter the Benjamini-Hochberg family.
STATS_FIT_STATUSES <- c("fitted", "too_few_observations", "rank_deficient")

# The first pyMzLib whose `stats fit` takes --contrasts, --confidence and --estimator (mzLib
# 1.0.595). An older bridge lists the verb but ignores the options, so stats_fit() checks the
# answer for what they add.
STATS_FIT_OPTIONS_SINCE <- "0.6.0"

# The variance-prior estimators stats_fit() takes as `estimator`: mzLib's VariancePriorEstimator
# MomentsLegacy (limma's eBayes(legacy = TRUE)) and MarginalLikelihood (each feature's own residual
# df, for unequal df). Leave `estimator` NULL for mzLib's default; `prior$estimator` says which.
STATS_PRIOR_ESTIMATORS <- c("moments_legacy", "marginal_likelihood")

# ---------------------------------------------------------------- argument shaping

# A number as text that .NET parses back to exactly the same double: 17 significant digits always
# suffice, and .NET's parser rounds correctly. Not "the shortest that round-trips in R": R's own
# parser is not correctly rounded on every platform, so it cannot be the judge.
stats_number_text <- function(x) {
  sprintf("%.17g", x)
}

stats_one_line <- function(value, what) {
  if (grepl("[\r\n]", value)) {
    stop(mzlib_usage_error(paste0(what, " contains a newline, which separates entries: '", value, "'.")))
  }
  value
}

stats_whole <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value != round(value)) {
    stop(mzlib_usage_error(paste0(name, " must be a whole number; got ", paste(deparse(value), collapse = " "), ".")))
  }
  formatC(value, format = "d")
}

# One design column name per stdin line, for every per-coefficient fit. Empty is allowed only
# where the caller says so (stats_fit() with contrasts).
stats_coefficient_names <- function(coefficients, allow_empty = FALSE, caller = "stats_fit") {
  if (is.null(coefficients)) {
    coefficients <- character(0)
  }
  if (!is.character(coefficients)) {
    stop(mzlib_usage_error(
      "coefficients must be a character vector of design column names, e.g. \"age_decades\"."
    ))
  }
  if (length(coefficients) == 0L && !allow_empty) {
    stop(mzlib_usage_error(paste0(
      caller, "() needs at least one coefficient to report, as a character vector."
    )))
  }
  for (name in coefficients) {
    if (is.na(name) || !nzchar(trimws(name))) {
      stop(mzlib_usage_error("Every coefficient must be a non-empty name; got a blank or NA."))
    }
    stats_one_line(name, "A coefficient name")
  }
  unname(coefficients)
}

# A single finite or infinite number, as .NET reads it back exactly.
stats_number <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    stop(mzlib_usage_error(paste0(name, " must be a number; got ", paste(deparse(value), collapse = " "), ".")))
  }
  stats_number_text(value)
}

stats_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(mzlib_usage_error(paste0(name, " must be TRUE or FALSE; got ", paste(deparse(value), collapse = " "), ".")))
  }
  value
}

stats_build_fit_request <- function(responses, design, coefficients, trend, spline_basis, threads,
                                    contrasts = NULL, confidence = NULL, estimator = NULL) {
  if (!is.character(coefficients) && !is.null(coefficients)) {
    stop(mzlib_usage_error(
      "coefficients must be a character vector of design column names to test, e.g. \"age_decades\"."
    ))
  }
  if (length(coefficients) == 0L && is.null(contrasts)) {
    stop(mzlib_usage_error(
      "stats_fit() needs at least one coefficient to test, as a character vector, or a contrasts file."
    ))
  }
  coefficients <- stats_coefficient_names(coefficients, allow_empty = TRUE)
  if (!is.logical(trend) || length(trend) != 1L || is.na(trend)) {
    stop(mzlib_usage_error("trend must be TRUE or FALSE."))
  }
  args <- c(
    "stats", "fit",
    "--responses", readers_normalise_path(responses, "responses"),
    "--design", readers_normalise_path(design, "design"),
    "--threads", stats_whole(threads, "threads")
  )
  if (isTRUE(trend)) {
    args <- c(args, "--trend")
  }
  if (!is.null(spline_basis)) {
    args <- c(args, "--spline-basis", stats_whole(spline_basis, "spline_basis"))
  }
  # Each is sent only when given, so mzLib's own default applies otherwise.
  if (!is.null(contrasts)) {
    args <- c(args, "--contrasts", readers_normalise_path(contrasts, "contrasts"))
  }
  if (!is.null(confidence)) {
    args <- c(args, "--confidence", stats_number(confidence, "confidence"))
  }
  if (!is.null(estimator)) {
    if (!is.character(estimator) || length(estimator) != 1L || !estimator %in% STATS_PRIOR_ESTIMATORS) {
      stop(mzlib_usage_error(paste0(
        "estimator must be one of ", paste(STATS_PRIOR_ESTIMATORS, collapse = ", "),
        " (mzLib's VariancePriorEstimator MomentsLegacy or MarginalLikelihood), or NULL for mzLib's ",
        "default; got ", paste(deparse(estimator), collapse = " "), "."
      )))
    }
    args <- c(args, "--estimator", estimator)
  }
  # No trailing newline: system2(input =) ends the last line itself (see stats_adjust_stdin()).
  # Contrasts alone send nothing at all, not one blank line.
  list(args = args, stdin = if (length(coefficients) > 0L) paste(coefficients, collapse = "\n") else NULL)
}

# The error for a `stats fit` option the bridge in use predates, and so ignored.
stats_older_fit_bridge <- function(option) {
  exe <- tryCatch(mzlibr_bridge_path(), error = function(e) "(unknown)")
  mzlib_usage_error(paste0(
    "stats_fit(", option, " = ...) needs the bridge from pyMzLib ", STATS_FIT_OPTIONS_SINCE,
    " or later, and the bridge at '", exe, "' predates it, so it ignored --", option, ". Point ",
    MZLIB_BRIDGE_ENV_VAR, " (or options(", MZLIB_BRIDGE_OPTION, ")) at a newer bridge, or install ",
    "one with mzlibr_install_bridge()."
  ))
}

# Run the fit and hold the answer to the options it was asked: a bridge from before pyMzLib 0.6.0
# lists `stats fit` but ignores --contrasts and --estimator, and would otherwise return a fit that
# silently did not do what was asked. `invoke` is the seam the tests replace.
stats_fit_answer <- function(request, contrasts, estimator, timeout, invoke = bridge_invoke) {
  data <- tryCatch(
    invoke(request$args, stdin = request$stdin, timeout = timeout),
    mzlib_usage_error = function(e) {
      # Such a bridge ignores --contrasts and so, with no coefficient, says it was given none: the
      # option it did not read is the real cause.
      if (!is.null(contrasts) && is.null(request$stdin) &&
        grepl("No coefficient to test", conditionMessage(e), fixed = TRUE)) {
        stop(stats_older_fit_bridge("contrasts"))
      }
      stop(e)
    }
  )
  if (!is.null(contrasts) && wire_null(data[["contrasts"]])) {
    stop(stats_older_fit_bridge("contrasts"))
  }
  prior <- data[["prior"]]
  if (!is.null(estimator) && (!is.list(prior) || wire_null(prior[["estimator"]]))) {
    stop(stats_older_fit_bridge("estimator"))
  }
  data
}

# One line per p-value; NA or NaN crosses as an empty line, which the bridge keeps as untested.
stats_adjust_stdin <- function(p_values) {
  if (!is.numeric(p_values) && !(is.logical(p_values) && all(is.na(p_values)))) {
    stop(mzlib_usage_error(paste0(
      "p_values must be a numeric vector, NA for an untested feature; got ", class(p_values)[1L], "."
    )))
  }
  if (length(p_values) == 0L) {
    stop(mzlib_usage_error("stats_adjust() needs at least one p-value."))
  }
  values <- as.numeric(p_values)
  lines <- ifelse(is.na(values), "", stats_number_text(values))
  # No trailing newline here: system2(input =) writes the text with writeLines(), which ends the
  # last line itself. One more would be read as one more untested p-value - blank lines are kept.
  paste(lines, collapse = "\n")
}

stats_meta_stdin <- function(studies) {
  if (!is.data.frame(studies)) {
    stop(mzlib_usage_error(
      "studies must be a data.frame with feature, estimate and standard_error columns, one row per study."
    ))
  }
  missing <- setdiff(c("feature", "estimate", "standard_error"), names(studies))
  if (length(missing) > 0L) {
    stop(mzlib_usage_error(paste0(
      "studies has no ", paste(missing, collapse = ", "), " column; it needs feature, estimate and standard_error."
    )))
  }
  if (nrow(studies) == 0L) {
    stop(mzlib_usage_error("stats_meta() needs at least one study."))
  }
  feature <- as.character(studies$feature)
  for (i in seq_along(feature)) {
    if (is.na(feature[[i]]) || !nzchar(trimws(feature[[i]]))) {
      stop(mzlib_usage_error(paste0("studies row ", i, " needs a non-empty feature name.")))
    }
    if (grepl("[\t\r\n]", feature[[i]])) {
      stop(mzlib_usage_error(paste0("studies row ", i, "'s feature name contains a tab or newline: '", feature[[i]], "'.")))
    }
  }
  for (column in c("estimate", "standard_error")) {
    if (!is.numeric(studies[[column]]) || anyNA(studies[[column]])) {
      stop(mzlib_usage_error(paste0("studies$", column, " must be numbers, with no NA.")))
    }
  }
  paste(feature, stats_number_text(studies$estimate), stats_number_text(studies$standard_error),
    sep = "\t", collapse = "\n"
  )
}

# ---------------------------------------------------------------- parsing

stats_table <- function(data) {
  table <- readers_parse_records_table(data)
  columns <- data[["columns"]]
  empty <- !is.list(columns) || all(vapply(columns, function(v) !is.list(v) || length(v) == 0L, logical(1L)))
  if (is.null(table) || empty) {
    names <- wire_strings(data[["column_names"]])
    table <- as.data.frame(structure(rep(list(logical(0)), length(names)), names = names),
      stringsAsFactors = FALSE, optional = TRUE
    )
  }
  table
}

stats_count <- function(data, name) {
  as.numeric(wire_field(data, name, "numeric", NA_real_))
}

# The contrasts block of a fit, or NULL when no contrasts file was given.
stats_parse_contrasts <- function(block) {
  if (!is.list(block) || is.null(names(block))) {
    return(NULL)
  }
  weights <- block[["weights"]]
  weights <- if (is.list(weights) && length(weights) > 0L) {
    lapply(weights, function(w) vapply(w, function(v) if (wire_null(v)) NA_real_ else as.numeric(v), numeric(1L)))
  } else {
    structure(list(), names = character(0))
  }
  structure(
    list(
      confidence_level = stats_count(block, "confidence_level"),
      tested_contrasts = wire_strings(block[["tested_contrasts"]]),
      weights = weights,
      row_count = stats_count(block, "row_count"),
      column_names = wire_strings(block[["column_names"]]),
      records = stats_table(block)
    ),
    class = "mzlibr_moderated_contrasts"
  )
}

stats_parse_fit <- function(data) {
  prior <- data[["prior"]]
  prior <- if (is.list(prior)) lapply(prior, function(v) if (wire_null(v)) NA else v) else NULL
  # A bridge before pyMzLib 0.6.0 does not report the estimator: NA says "not reported".
  if (is.list(prior) && is.null(prior[["estimator"]])) {
    prior <- c(list(estimator = NA_character_), prior)
  }
  structure(
    list(
      responses_file = wire_field(data, "responses_file", "character", NA_character_),
      design_file = wire_field(data, "design_file", "character", NA_character_),
      contrasts_file = wire_field(data, "contrasts_file", "character", NA_character_),
      feature_count = stats_count(data, "feature_count"),
      sample_count = stats_count(data, "sample_count"),
      sample_names = wire_strings(data[["sample_names"]]),
      coefficient_names = wire_strings(data[["coefficient_names"]]),
      tested_coefficients = wire_strings(data[["tested_coefficients"]]),
      trend = isTRUE(data[["trend"]]),
      threads = stats_count(data, "threads"),
      missing_count = stats_count(data, "missing_count"),
      zero_count = stats_count(data, "zero_count"),
      status_counts = proteins_parse_counts(data[["status_counts"]]),
      residual_df_differ = isTRUE(data[["residual_df_differ"]]),
      prior = prior,
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      contrasts = stats_parse_contrasts(data[["contrasts"]]),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_moderated_fit"
  )
}

stats_parse_adjusted <- function(data) {
  structure(
    list(
      line_count = stats_count(data, "line_count"),
      tested_count = stats_count(data, "tested_count"),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_adjusted"
  )
}

stats_parse_meta <- function(data) {
  structure(
    list(
      study_count = stats_count(data, "study_count"),
      feature_count = stats_count(data, "feature_count"),
      confidence = stats_count(data, "confidence"),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = stats_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_meta_analysis"
  )
}

# ---------------------------------------------------------------- the public surface

#' Fit one linear model per feature and test coefficients, and contrasts, with limma's moderated t
#'
#' Runs mzLib's `LinearModel.Fit` once - every design column is fitted - and then
#' `EmpiricalBayes.Moderate` for each coefficient you name, which also applies Benjamini-Hochberg
#' across that coefficient's fitted features. This is limma's `lmFit()` then
#' `eBayes(legacy = TRUE)` (Smyth 2004), which mzLib holds to limma to 1e-8 relative on limma's own
#' reference output. With `contrasts`, `EmpiricalBayes.ModerateContrasts` also tests each row of a
#' contrasts file - "250 ng against 125 ng", "the mean of two groups against a third" - with its
#' confidence interval: limma's `contrasts.fit()` then `eBayes()` then `topTable(confint = ...)`.
#'
#' A feature is fitted on the samples where it was observed, and is reported - never dropped - when
#' it cannot be: `"too_few_observations"` when it has no more observed samples than the design has
#' coefficients, `"rank_deficient"` when its observed samples cannot separate the coefficients (a
#' feature seen in only one group cannot have a group effect).
#'
#' @param responses A TSV with a header row. First column the feature id; every other column one
#'   sample, named by its header. Use the values as you want them modelled - log2 intensities,
#'   usually. Blank, `NA` or `NaN` is missing and left out of that feature's fit, never imputed;
#'   **0 is an observation**, so blank out any 0 that means "not measured".
#' @param design A TSV with a header row and one row per sample: first column the sample name
#'   (matching the responses header exactly, in any order), then one numeric column per
#'   coefficient. Include an intercept column yourself if you want one. [diff_design()] writes one
#'   (`design_<i>.tsv`) from a sample sheet.
#' @param coefficients Design column names to test, each with its own moderated test and its own
#'   Benjamini-Hochberg family. May be empty (`character(0)`) when you give `contrasts`.
#' @param contrasts A TSV of contrasts, as [diff_design()] writes it beside the design
#'   (`contrasts_<i>.tsv`): a header of the contrast id column, then exactly the design's coefficient
#'   columns (in any order); then one row per contrast, its id and one weight per coefficient, e.g.
#'   `-1` and `1` for the second level minus the first. Each contrast gets its own moderated test,
#'   Benjamini-Hochberg family and confidence interval, in `contrasts`. `NULL` sends none. Needs the
#'   pyMzLib 0.6.0 bridge.
#' @param confidence The two-sided level of the contrasts' intervals, as a fraction (0 to 1),
#'   strictly between 0 and 1; only with `contrasts`. Leave it `NULL` for mzLib's default (0.95 at
#'   mzLib `d10017ce`); `contrasts$confidence_level` says which level was used.
#' @param estimator How the variance prior is estimated, one of `"moments_legacy"` (limma's
#'   `eBayes(legacy = TRUE)`) or `"marginal_likelihood"` (mzLib's own estimator, which uses each
#'   feature's own residual degrees of freedom, so it handles the unequal df that missing values
#'   cause exactly). Leave it `NULL` for mzLib's default (`"moments_legacy"` at mzLib `d10017ce`);
#'   `prior$estimator` says which was used.
#' @param trend Let the prior variance follow average response (limma's `trend = TRUE`), for data
#'   whose low-abundance features are noisier.
#' @param spline_basis Basis functions of the trend spline, intercept included. Only with
#'   `trend = TRUE`. Leave it `NULL` to let mzLib choose, as limma does, for the number of
#'   features; mzLib may also use fewer than you ask for. Either way `prior$spline_basis_count`
#'   reports what it used.
#' @param threads Features fitted at once (threads), or `-1` for every core. The answer is identical
#'   at any value.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_moderated_fit`:
#'
#'   - `records`, a data.frame with one row per (tested coefficient, feature) - one block of
#'     `feature_count` rows per coefficient, in the order you named them, and no rows when you named
#'     none - with `feature`,
#'     `coefficient`, `status`, `estimate` (response units per coefficient unit: a **log2 fold
#'     change** for a 0/1 group column on log2 intensities), `standard_error` (response units per
#'     coefficient unit), `t` (standard errors), `df_total` (degrees of freedom), `p_value` and
#'     `bh_adjusted` (fractions, 0 to 1), `posterior_variance` and `prior_variance` (squared
#'     response units), `sigma` (response units), `df_residual`, `observed` (samples) and
#'     `average_response` (response units). Statistics are `NA` where `status` is not `"fitted"`.
#'   - `responses_file`, `design_file` and `contrasts_file` (`NA` when you gave no `contrasts`), as
#'     absolute paths; `feature_count` features and `sample_count` samples read; `sample_names`,
#'     `coefficient_names` (every design column, in the order each contrast's weights follow) and
#'     `tested_coefficients` (empty when you named none); `missing_count` and `zero_count` response
#'     cells; `status_counts`, features per status; `threads`; `row_count` rows.
#'   - `residual_df_differ`: `TRUE` when the fitted features do not all have the same residual
#'     degrees of freedom. Then default limma (3.61 and later) uses a different prior estimator and
#'     reports different moderated statistics; under `"moments_legacy"` these are `legacy = TRUE`,
#'     and `estimator = "marginal_likelihood"` is mzLib's estimator for this case.
#'   - `prior`, the fitted variance prior, shared by every coefficient and contrast: `estimator`,
#'     the estimator used (`NA` only from a bridge before pyMzLib 0.6.0, which did not report it);
#'     `df` (degrees of freedom, `NA` when infinite; under `"marginal_likelihood"` mzLib searches d0
#'     from about 0.05 upward, so a d0 there is that floor, not an optimum), `df_infinite`,
#'     `trended`, `spline_basis_count` (basis functions) and `scale` (squared response units), s0^2
#'     when every fitted feature shares one. `scale` is `NA` when it differs by feature, which is
#'     then in the `prior_variance` column: with a trend, or with fewer than 3 features. Then mzLib
#'     fits no prior (mzLib #1398, as limma's `squeezeVar` does): `df` is 0 and each feature keeps
#'     its own variance. Bridges before pyMzLib 0.6.0 report the first feature's own variance as
#'     `scale` in that case.
#'   - `contrasts`, `NULL` when you gave no `contrasts`; otherwise an `mzlibr_moderated_contrasts`:
#'     `confidence_level`, the two-sided level of every interval as a fraction (0 to 1);
#'     `tested_contrasts`, the contrast ids in the file's order; `weights`, a list from contrast id
#'     to its weights (contrast weight per coefficient, unitless), one per `coefficient_names`
#'     entry, whatever order the file's columns were in; `row_count` rows, `feature_count` x the
#'     contrasts; `column_names`; and `records`, one row per (contrast, feature), a block per
#'     contrast in file order, with `feature`, `contrast` (its id), `status`, `estimate` (response
#'     units per contrast unit: a **log2 fold change** for one level against another on log2
#'     intensities), `standard_error`, `confidence_low` and `confidence_high` (response units per
#'     contrast unit; the two-sided `confidence_level` interval, not adjusted across features), and
#'     the coefficient table's `t`, `df_total`, `p_value`, `bh_adjusted` (over this contrast's
#'     fitted features), `posterior_variance`, `prior_variance`, `sigma`, `df_residual`, `observed`
#'     and `average_response`, in the same units. Statistics are `NA` where `status` is not
#'     `"fitted"`.
#'   - `caveats`, including any that apply to this input: the estimator, unfitted features, zeros,
#'     `residual_df_differ`, contrasts.
#'
#' @section What it is not:
#'
#' By default it is limma's *legacy* estimator; `estimator = "marginal_likelihood"` is mzLib's own,
#' close to default limma's but not equal. `robust = TRUE`, the B-statistic and observation weights
#' are not implemented. With missing values a contrast's standard error is exact per feature, where
#' limma's `contrasts.fit()` approximates it, so on such data the contrasts differ from limma's by
#' design. `bh_adjusted` controls the false discovery rate among the features tested; it is not a
#' target-decoy q-value.
#'
#' @section Bridges before pyMzLib 0.6.0:
#'
#' Such a bridge has `stats fit` but ignores `--contrasts`, `--confidence` and `--estimator`. Rather
#' than return a fit that silently did not do what you asked, `stats_fit()` raises an
#' `mzlib_usage_error` naming the 0.6.0 bridge when you gave `contrasts` or `estimator` and the answer
#' lacks what they add.
#'
#' @seealso [diff_design()], which writes the design and contrasts files; [stats_adjust()],
#'   [stats_meta()], [stats_fit_mixed()] for samples in groups
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # A dilution series from mzLib's own RNA test data: MALAT1 at 500, 250 and 125 ng against a
#' # constant FLuc spike, log2 intensities normalised to FLuc. Halving MALAT1 should read as about -1.
#' data <- system.file("extdata", "stats", package = "mzLibR")
#' fit <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
#'   file.path(data, "malat_dilution_design.tsv"), c("malat_250ng", "malat_125ng"))
#' fit
#' fit$status_counts
#' c(estimator = fit$prior$estimator, prior_df = round(fit$prior$df, 1))
#' malat <- fit$records[fit$records$coefficient == "malat_250ng" & fit$records$status == "fitted" &
#'   startsWith(fit$records$feature, "MALAT1:"), ]
#' c(oligos = nrow(malat), median_log2_fold_change = median(malat$estimate))
#'
#' # Every pair of loads, with intervals: diff_design() with no reference level writes one contrast
#' # per pair beside the design, and `contrasts` tests them all.
#' samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
#' design <- diff_design(samples, factors = list(factor_spec("malat_load")),
#'   out = file.path(tempdir(), "malat_pairs"))
#' stratum <- design$strata[[1]]
#' pairs <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"), stratum$design_file,
#'   contrasts = stratum$contrasts_file)
#' pairs$contrasts
#' pairs$contrasts$weights$c3
#' pairs$coefficient_names
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.fit
#' @export
stats_fit <- function(responses, design, coefficients = character(0), contrasts = NULL,
                      confidence = NULL, estimator = NULL, trend = FALSE, spline_basis = NULL,
                      threads = 1, timeout = 600) {
  request <- stats_build_fit_request(responses, design, coefficients, trend, spline_basis, threads,
    contrasts = contrasts, confidence = confidence, estimator = estimator
  )
  bridge_require_verb("stats fit", STATS_SINCE)
  stats_parse_fit(stats_fit_answer(request, contrasts, estimator, timeout))
}

#' Benjamini-Hochberg adjust p-values, keeping every position
#'
#' For the i-th smallest of m p-values the adjusted value is the minimum over j >= i of
#' p_(j) m / j, capped at 1 (Benjamini and Hochberg 1995), computed by mzLib's
#' `MultipleTesting.BenjaminiHochberg`. `NA` (or `NaN`) marks a feature that was not tested: it
#' stays `NA` and is **not counted in m**, which is the honest family. Do not write 1 for an
#' untested feature; that enlarges m.
#'
#' [stats_fit()] already reports `bh_adjusted`; this is for p-values from anywhere else.
#'
#' @param p_values A numeric vector of p-values, each a fraction from 0 to 1, or `NA` for untested.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_adjusted`: `records`, a data.frame with row i for `p_values[i]` - `p_value`
#'   as read and `bh_adjusted`, both fractions (0 to 1), `NA` for an untested entry; `line_count`
#'   lines (values) given; `tested_count`, m, the p-values the adjustment is over; `row_count` rows;
#'   and `caveats`.
#'
#' @seealso [stats_fit()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' adjusted <- stats_adjust(c(0.0002, 0.004, 0.019, NA, 0.031, 0.2, NaN, 0.74))
#' adjusted$tested_count
#' adjusted$records
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.adjust
#' @export
stats_adjust <- function(p_values, timeout = 60) {
  stdin <- stats_adjust_stdin(p_values)
  bridge_require_verb("stats adjust", STATS_SINCE)
  stats_parse_adjusted(bridge_invoke(c("stats", "adjust"), stdin = stdin, timeout = timeout))
}

#' Pool one effect size per study into a random-effects estimate per feature
#'
#' Each feature is pooled on its own by mzLib's `RandomEffectsMeta.Pool`, with the
#' DerSimonian-Laird (1986) moment estimator of the between-study variance, as metafor's
#' `rma(method = "DL")` does (mzLib matches metafor to 1e-8 relative). Every row also says how many
#' studies agree in direction and how far the estimate moves when any one is dropped - the two
#' checks a reader of a pooled result asks for first.
#'
#' @param studies A data.frame with one row per study and columns `feature`, `estimate` and
#'   `standard_error`. A feature's studies need not be adjacent. The standard error must be
#'   positive.
#' @param confidence Coverage of the reported interval, as a fraction (0 to 1).
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_meta_analysis`: `records`, a data.frame with one row per feature in order of
#'   first appearance - `feature`; `studies` pooled; `estimate`, `standard_error`, `confidence_low`
#'   and `confidence_high` (units of the study estimates); `p_value` (a fraction, 0 to 1, normal
#'   reference, not adjusted across features); `tau2` (squared units of the study estimates); `q`
#'   (a chi-square statistic) and `i_squared` (a fraction, 0 to 1), both `NA` for one study;
#'   `direction_agree` and `direction_disagree` (studies); and `leave_one_out_max_delta` (units of
#'   the study estimates, `NA` for one study). Also `study_count` studies, `feature_count`
#'   features, `confidence` (a fraction, 0 to 1), `row_count` rows and `caveats`.
#'
#' @section With few studies:
#'
#' DerSimonian-Laird with a normal reference: the interval and p-value do not widen for few studies
#' (no Hartung-Knapp adjustment), so treat a pooled estimate from two or three studies with care.
#' Send `records$p_value` to [stats_adjust()] to adjust across features.
#'
#' @seealso [stats_fit()], [stats_adjust()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # metafor's own DerSimonian-Laird reference cases, from mzLib's test data:
#' inputs <- read.delim(system.file("extdata", "stats", "metafor_dl_inputs.tsv", package = "mzLibR"))
#' studies <- data.frame(feature = inputs$case, estimate = inputs$yi, standard_error = inputs$sei)
#' pooled <- stats_meta(studies)
#' pooled$records[, c("feature", "studies", "estimate", "tau2", "direction_agree")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec stats.meta
#' @export
stats_meta <- function(studies, confidence = 0.95, timeout = 60) {
  if (!is.numeric(confidence) || length(confidence) != 1L || is.na(confidence)) {
    stop(mzlib_usage_error("confidence must be a number strictly between 0 and 1."))
  }
  stdin <- stats_meta_stdin(studies)
  bridge_require_verb("stats meta", STATS_SINCE)
  stats_parse_meta(bridge_invoke(c("stats", "meta", "--confidence", stats_number_text(confidence)),
    stdin = stdin, timeout = timeout
  ))
}

# ---------------------------------------------------------------- printing

#' Print a moderated fit
#'
#' @param x A [stats_fit()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_moderated_fit <- function(x, ...) {
  tested <- if (length(x$tested_coefficients) > 0L) paste(x$tested_coefficients, collapse = ", ") else "no coefficient"
  cat("<mzlibr_moderated_fit> ", format(x$feature_count), " features x ", format(x$sample_count),
    " samples, testing ", tested,
    if (!is.null(x$contrasts)) paste0(" and ", length(x$contrasts$tested_contrasts), " contrast(s)") else "",
    "\n",
    sep = ""
  )
  if (length(x$status_counts) > 0L) {
    cat("  ", paste0(names(x$status_counts), " ", x$status_counts, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$prior)) {
    estimator <- x$prior$estimator
    cat("  prior df ", if (isTRUE(x$prior$df_infinite)) "infinite" else format(signif(x$prior$df, 4)),
      if (isTRUE(x$prior$trended)) ", trended" else "",
      if (!is.null(estimator) && !is.na(estimator)) paste0(" (", estimator, ")") else "", "\n",
      sep = ""
    )
  }
  if (isTRUE(x$residual_df_differ)) {
    cat("  ! residual df differ between features",
      if (identical(x$prior$estimator, "marginal_likelihood")) {
        ": the marginal-likelihood prior uses each feature's own df\n"
      } else {
        ": these are eBayes(legacy = TRUE) statistics\n"
      },
      sep = ""
    )
  }
  invisible(x)
}

#' Print the contrasts of a moderated fit
#'
#' @param x The `contrasts` of a [stats_fit()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_moderated_contrasts <- function(x, ...) {
  cat("<mzlibr_moderated_contrasts> ", length(x$tested_contrasts), " contrast(s) x ",
    format(if (length(x$tested_contrasts) > 0L) x$row_count / length(x$tested_contrasts) else 0),
    " features, ", format(100 * x$confidence_level), "% intervals\n",
    sep = ""
  )
  for (id in x$tested_contrasts) {
    cat("  ", id, ": ", paste(format(x$weights[[id]]), collapse = " "), "\n", sep = "")
  }
  invisible(x)
}

#' Print Benjamini-Hochberg adjusted p-values
#'
#' @param x A [stats_adjust()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_adjusted <- function(x, ...) {
  cat("<mzlibr_adjusted> ", format(x$line_count), " values, ", format(x$tested_count),
    " tested (m)\n",
    sep = ""
  )
  invisible(x)
}

#' Print a meta-analysis
#'
#' @param x A [stats_meta()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_meta_analysis <- function(x, ...) {
  cat("<mzlibr_meta_analysis> ", format(x$feature_count), " features pooled from ",
    format(x$study_count), " studies, ", format(100 * x$confidence), "% intervals\n",
    sep = ""
  )
  invisible(x)
}
