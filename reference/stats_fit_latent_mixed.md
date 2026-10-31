# Fit a mixed model per feature on a latent scale, from samples measured on it, counted through it, or both

Fit one random-intercept mixed model per feature on a latent scale (such
as the log-odds of a site's occupancy), where each sample is either
measured on that scale or counted through it (k of n trials succeed),
and test each named coefficient with a Wald t.

Answers **does a value follow a trait, when some runs measured it and
others only counted successes, k of n?** A PTM site's occupancy is the
usual case: where the modified and unmodified peptides were both
quantified you have a number; where they were not, you may know only
that 3 of 8 covering peptidoforms carried the modification.
[`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)
has to drop the counted runs. This uses both kinds.

## Usage

``` r
stats_fit_latent_mixed(design, groups, coefficients, link, responses = NULL,
  counts = NULL, depth = NULL, residual_variance = NULL, estimator = NULL, threads = 1,
  timeout = 600)
```

## Arguments

- design:

  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  design: a TSV with a header row and one row per sample, the sample
  name then one numeric column per fixed-effect coefficient. Do not add
  a column per group: the random intercept replaces it. mzLib refuses a
  design that is not of full column rank.

- groups:

  [`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)'s
  groups table: one row per sample, the sample name then ONE column with
  its group label (an animal, a donor, a dataset), whose header names
  the grouping.

- coefficients:

  Design column names to report, each as one block of rows, sent one per
  stdin line. Every design column is fitted whatever you name.

- link:

  mzLib's `DetectionLink`, in natural log-odds: `"identity"`
  (`DetectionLink.Identity`, intercept 0 and slope 1: a count samples
  the same probability a measurement reads), a
  [`stats_fit_link`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md)
  result or its `$link`, or 2 to 4 numbers
  `c(intercept, slope[, depth_slope[, overdispersion]])`, sent as given:
  a value you leave out is mzLib's default (0 at the pin). Required:
  mzLib has no default link.

- responses:

  The measured samples, in latent units (logit occupancy, for example,
  not a fraction):
  [`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)'s
  responses table, a TSV with the feature id first and one column per
  sample. Blank, `NA` or `NaN` is not measured. Give `responses`,
  `counts` or both. `NULL` (the default): nothing was measured.

- counts:

  The counted samples, in trials: a table laid out as `responses`, each
  cell `successes/trials` (`3/8`), or blank where the sample was not
  counted. With `responses` too, the two must list the same features and
  samples, and no cell may be in both. `NULL` (the default): nothing was
  counted, and the model is
  [`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)'s.

- depth:

  A TSV with a header row and one row per sample: the sample name, then
  its depth in your depth units (log2 target PSMs, say), for the link's
  depth term. Needed when the link's `depth_slope` is not 0. `NULL` (the
  default): none is sent.

- residual_variance:

  A TSV with a header row and one row per feature: the feature id, then
  sigma^2 to use for it (squared latent units) instead of estimating it,
  for a feature with too few measured samples. Blank, or a feature not
  listed, is estimated. `NULL` (the default): every sigma^2 is
  estimated.

- estimator:

  `"reml"` or `"ml"`. Leave it `NULL` (the default) for mzLib's default
  (`"reml"` at the pin); `estimator` in the result says which was used.

- threads:

  Features fitted at once (threads), or `-1` for every core. The answer
  is identical at any value.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

mzLib's `LogitNormalMixedModel.Fit` (mzLib \#1469) is
[`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)'s
model - fixed effects from your design, one baseline per group drawn
from N(0, tau^2), run-to-run scatter sigma^2 - on a **latent** value,
such as the log-odds of a site's occupancy. A sample can be **measured**
(its latent value is in `responses`) or **counted**: k of its n trials
succeeded (`counts`), each with a chance g set by the latent value
through the `link`, logit g = intercept + slope x latent + depth_slope x
depth. It is `lme4::glmer` (a binomial GLMM, with an observation-level
random effect for sigma^2) and `lmer` in one model, which mzLib checks
against both and against `glmmTMB`.

**Which cells are counted, and the link, are your choices**: mzLibR
builds neither.
[`stats_fit_link`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md)
fits a link from cells.

## Value

An `mzlibr_latent_mixed_model_fit`:

\- `records`, a data.frame with one row per (reported coefficient,
feature) - one block of `feature_count` rows per coefficient you named,
in that order, each in the feature order of the responses file (of the
counts file, without one): `feature`; `coefficient`; `status`, one of
`"fitted"`, `"too_few_observations"`, `"rank_deficient"`, `"separated"`
(no finite estimate exists: a status, never a p-value),
`"not_converged"` or `"too_few_groups"`; `estimate`, the fixed effect
within groups on the latent scale (latent units per coefficient unit:
**natural log-odds** for a logit latent, so `exp(estimate)` is an odds
ratio); `standard_error` (latent units per coefficient unit), from the
inverse Hessian at the estimated variances; `df` (degrees of freedom),
nlme's containment df over the samples with a response of either kind,
`NA` when not fitted or when no positive df exists; `t` (standard
errors), the Wald statistic; `p_value`, a fraction (0 to 1), two-sided
on `df`, **not adjusted across features**, `NA` also when `df` is;
`residual_variance` (squared latent units), sigma^2, run to run within a
group, estimated or the one you gave; `residual_variance_given`, `TRUE`
when it is the one you gave; `group_variance` (squared latent units),
tau^2, where 0 is the boundary; `log_likelihood` (log-likelihood,
natural log), Laplace-approximated restricted under `"reml"`, full under
`"ml"`; `measured`, samples with a finite response (samples); `counted`,
samples with a count (samples); and `groups`, distinct groups among both
(groups). Statistics are `NA` where `status` is not `"fitted"`. -
`responses_file`, `counts_file`, `depth_file` and
`residual_variance_file` (`NA` when not given), `design_file` and
`groups_file`, as absolute paths; `grouping`, the groups file's label
column header; `group_names`, in order of first appearance;
`group_count` groups; `feature_count` features and `sample_count`
samples read; `sample_names`; `coefficient_names` and
`tested_coefficients`; `estimator`, the estimator mzLib used; `link`,
the link every count was read through, as mzLib reports it: an
`mzlibr_detection_link` with `intercept` (natural log-odds), `slope`
(natural log-odds per latent unit), `depth_slope` (natural log-odds per
depth unit) and `overdispersion` (fraction, the intra-cell correlation,
0 to below 1); `threads` (threads); `measured_cells`, `counted_cells`
and `missing_count` (cells); `status_counts`, features per status
(features); `row_count` rows; `column_names`. - `caveats`: what the
model is and is not, the link and the estimator used, plus counts for
this input: unfitted features, rows with no df, features at the
boundary, residual variances given.

## What it is not

A random intercept per group and nothing more: no random slopes, nested
or crossed groupings. Responses and estimates are on the latent scale,
never fractions or intensities. Standard errors are at the estimated
variances and do not include their uncertainty, so they are not
`glmmTMB`'s [`vcov()`](https://rdrr.io/r/stats/vcov.html). Under
`"reml"` the fixed effects are integrated out by a Laplace
approximation, as `glmmTMB`'s `REML = TRUE` does; with nothing counted
it is exactly
[`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)'s
REML, and the two agree to about 1e-5. Send one coefficient's `p_value`
column to
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
for Benjamini-Hochberg.

## Wraps

Wire verb `stats fit-latent-mixed`. Generated from the bridge's verb
spec `stats.fit-latent-mixed.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`LogitNormalMixedModel.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`LogitNormalMixedModelFit`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`DetectionLink`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`BinomialCount`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`LogitNormalMixedModel.SeparationLimit`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`VarianceEstimator`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/MixedModel.cs)
  in `mzLib/StatisticalModels/MixedModel.cs` at mzLib `cc67b7a2`

- [`FeatureFitStatus`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LinearModel.cs)
  in `mzLib/StatisticalModels/LinearModel.cs` at mzLib `cc67b7a2`

## Parameters: units, ranges and defaults

- `responses`:

  path; in **latent units (e.g. natural log-odds of occupancy)**;
  default absent. The measured samples: stats fit-mixed's responses
  table, a TSV with a header row, the feature id first (unique) and then
  one column per sample, named by its header. Values are on the LATENT
  scale (logit occupancy, for example), not fractions or intensities. A
  blank cell, NA or NaN (any case) is not measured. At least one of
  responses and counts is required. Absent, no sample is measured: the
  bridge passes an all-missing matrix of the counts' shape.

- `counts`:

  path; in **trials (successes of trials)**; default absent; range
  `each cell k/n, two whole numbers`. The counted samples (mzLib's
  BinomialCount): a TSV laid out as responses, each cell
  successes/trials such as 3/8, or blank, NA or NaN for a sample not
  counted. When both tables are given they must list exactly the same
  feature ids and sample columns, in any order (the responses' order is
  kept), and a cell may be measured or counted, not both (mzLib's rule).
  mzLib refuses a count with no trials, or more successes than trials.
  Absent, no sample is counted, and the model is stats fit-mixed's.

- `design`:

  path; required. stats fit's design table: one row per sample, the
  sample name (matching the fit's sample columns exactly, rows matched
  by name) then one column per fixed-effect coefficient, every cell a
  finite number. Include an intercept column if you want one; do not add
  a column per group, which the random intercept replaces. mzLib refuses
  a design that is not of full column rank over all samples.

- `groups`:

  path; required. stats fit-mixed's groups table: one row per sample,
  the sample name then EXACTLY ONE column, the sample's group label (an
  animal, a donor, a dataset), whose header names the grouping. The
  sample names must be exactly the fit's sample columns. Labels are
  trimmed and compared exactly; mzLib refuses an empty one.

- `link`:

  string; in **natural log-odds**; required; range
  `identity, or 2 to 4 comma-separated numbers: intercept,slope[,depth_slope[,overdispersion]]`.
  mzLib's DetectionLink: how the latent value sets the chance g that one
  trial of a counted cell succeeds, logit g = intercept + slope x
  latent + depth_slope x depth, with trials in a cell correlated by
  overdispersion (rho, a beta-binomial; 0 is the binomial). REQUIRED:
  mzLib's link parameter has no default. identity is
  DetectionLink.Identity, (0, 1): a counted cell then samples the same
  probability a measured one reads. Numbers are passed as the record's
  positional arguments, so a depth_slope or overdispersion left out is
  the record's own default (0 at cc67b7a2). stats fit-link fits one from
  cells. mzLib refuses a non-finite intercept or depth slope, a slope
  that is 0 or not finite, and an overdispersion outside \[0, 1).

- `depth`:

  path; in **depth units (the caller's, e.g. log2 target PSMs)**;
  default absent. A TSV with a header row and one row per sample: the
  sample name (exactly the fit's samples, in any order), then one
  number, the sample's depth for the link's depth term. Needed when the
  link's depth_slope is not 0 (mzLib refuses otherwise). Absent, the
  bridge passes none.

- `residual_variance` (wire `--residual-variance`):

  path; in **squared latent units**; default absent; range
  `each > 0 and finite, or blank to estimate`. A TSV with a header row
  and one row per feature: the feature id, then sigma^2 to use for it
  instead of estimating it (e.g. borrowed when a feature has too few
  measured samples). A blank, NA or NaN cell, or a feature not listed,
  is NaN, which mzLib reads as estimate it. An id that is not one of the
  fit's features is refused. Absent, the bridge passes none and every
  sigma^2 is estimated.

- `coefficients` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line; each a design column name, none repeated`. The
  coefficients to report, one design column name per line. The model
  always fits every design column; the table holds one block of rows per
  coefficient named here, in this order.

- `estimator`:

  string; default absent; range `reml | ml`. reml is
  VarianceEstimator.Reml (the fixed effects integrated out by a Laplace
  approximation, as glmmTMB's REML; exact when nothing is counted), ml
  is VarianceEstimator.MaximumLikelihood. Absent, the bridge sends none
  and LogitNormalMixedModel.Fit's default applies (Reml at cc67b7a2);
  the envelope's estimator says which was used.

- `threads`:

  int; in **threads**; default `1`; range `>= 1, or -1 for every core`.
  Features fitted at once. The bridge default is 1 (BULK.md); mzLib's
  own default is -1. Results are identical at any value (mzLib).

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `responses_file`:

  string; `NA` when no `--responses`: nothing was measured. Absolute
  path of the responses table.

- `counts_file`:

  string; `NA` when no `--counts`: nothing was counted. Absolute path of
  the counts table.

- `design_file`:

  string; never `NA`. Absolute path of the design table.

- `groups_file`:

  string; never `NA`. Absolute path of the groups table.

- `depth_file`:

  string; `NA` when no `--depth` was given. Absolute path of the depth
  table.

- `residual_variance_file`:

  string; `NA` when no `--residual-variance`: every sigma^2 was
  estimated. Absolute path of the residual-variance table.

- `grouping`:

  string; never `NA`. The groups file's label column header.

- `group_names`:

  string\[\]; never `NA`. The distinct labels, in order of first
  appearance over the samples (the order mzLib numbers them in).

- `group_count`:

  int; in **groups**; never `NA`. Distinct labels over all samples. A
  feature's own count is its groups column.

- `feature_count`:

  int; in **features**; never `NA`. Feature rows read.

- `sample_count`:

  int; in **samples**; never `NA`. Sample columns read.

- `sample_names`:

  string\[\]; never `NA`. In the responses header's order (the counts
  header's, with counts alone): the order the design, group and depth
  rows were put in.

- `coefficient_names`:

  string\[\]; never `NA`. Every design column, in order.

- `tested_coefficients`:

  string\[\]; never `NA`. The stdin coefficients, in order; the table
  holds one block of feature_count rows per entry.

- `estimator`:

  string; never `NA`. reml or ml: the estimator mzLib reports it fitted
  with (LogitNormalMixedModelFit.Estimator), the default when
  `--estimator` was absent.

- `link`:

  object; never `NA`. The link every count cell was read through, as
  mzLib reports it (LogitNormalMixedModelFit.Link); fields under
  result.tables.link.

- `threads`:

  int; in **threads**; never `NA`. The thread count used.

- `measured_cells`:

  int; in **cells**; never `NA`. Feature x sample cells with a finite
  response.

- `counted_cells`:

  int; in **cells**; never `NA`. Feature x sample cells with a count.

- `missing_count`:

  int; in **cells**; never `NA`. Cells with neither a response nor a
  count, left out of their feature's fit.

- `status_counts`:

  map\<string,int\>; in **features**; never `NA`. Features per status,
  always these six keys: fitted, too_few_observations, rank_deficient,
  separated, not_converged, too_few_groups. They are the outcomes
  LogitNormalMixedModel.Fit's body can return at cc67b7a2; a stop-gap
  until mzLib says per model (mzLib#1487).

- `row_count`:

  int; in **rows**; never `NA`. feature_count x the number of tested
  coefficients.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (tested coefficient, feature):
  coefficient blocks in stdin order, features in file order.

- `caveats`:

  string\[\]; never `NA`. Seven static caveats (one names the link used,
  one the estimator), plus one each, with counts, when any feature was
  not fitted, when a fitted row has no df, when a fitted feature's
  group_variance is 0, when any residual variance was given, and when
  nothing was counted.

## Columns of `records`

- `feature`:

  string; never `NA`. The feature id from the first column.

- `coefficient`:

  string; never `NA`. The coefficient this row reports.

- `status`:

  string; never `NA`. fitted \| too_few_observations (no more samples
  with a response of either kind than coefficients) \| rank_deficient
  (the design over those samples is not of full rank) \| separated (no
  estimate exists: a coefficient beyond +-SeparationLimit, 30 at
  cc67b7a2, on the latent scale, a Hessian that is not positive
  definite, or nothing measured and every group's counts one-sided, all
  0 or all full) \| not_converged (the variance search failed, or the
  likelihood could not be evaluated at its optimum) \| too_few_groups
  (fewer than two groups, or none with two samples). Only fitted rows
  carry statistics.

- `estimate`:

  float; in **latent units per coefficient unit (natural log-odds for a
  logit latent)**; `NA` when status is not fitted. The fixed-effect
  coefficient at the estimated variances: the effect within groups, on
  the latent scale.

- `standard_error`:

  float; in **latent units per coefficient unit**; `NA` when not fitted.
  sqrt of the diagonal of H^-1, H = -d2 l / d beta2 at the estimate (l
  the log-likelihood with the random effects integrated out). With
  nothing counted, (X' V^-1 X)^-1 as stats fit-mixed reports under reml.

- `df`:

  float; in **degrees of freedom**; `NA` when not fitted, or the
  containment df is not positive (a group-level coefficient with too few
  groups). Denominator df of the t test, by nlme's containment rule over
  the samples with a response of either kind, a whole number.

- `t`:

  float; in **standard errors**; `NA` when not fitted. Wald statistic,
  estimate / standard_error. Reported even when df is null.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not fitted, or df is null.
  Two-sided p-value of the t test on df. Not adjusted across features.

- `residual_variance`:

  float; in **squared latent units**; `NA` when not fitted. sigma^2, the
  run-to-run variance of the latent within a group: estimated, or the
  value given (residual_variance_given). The same in every coefficient
  block.

- `residual_variance_given`:

  bool; never `NA`. True when residual_variance is the value the
  `--residual-variance` file gave, not an estimate. Reported whatever
  the status.

- `group_variance`:

  float; in **squared latent units**; `NA` when not fitted. tau^2, the
  between-group variance of the random intercept. 0 is the boundary
  estimate. The same in every coefficient block.

- `log_likelihood`:

  float; in **log-likelihood (natural log)**; `NA` when not fitted. The
  maximized log-likelihood: under reml the Laplace-approximated
  restricted likelihood (exactly the restricted likelihood when nothing
  is counted), under ml the full one. The same in every coefficient
  block.

- `measured`:

  int; in **samples**; never `NA`. Samples with a finite response for
  this feature.

- `counted`:

  int; in **samples**; never `NA`. Samples with a count for this
  feature.

- `groups`:

  int; in **groups**; never `NA`. Distinct groups among this feature's
  samples with a response of either kind.

## Fields of `link`

an `mzlibr_detection_link`, a named list of these four numbers.

- `link$intercept` (wire `intercept`):

  float; in **natural log-odds**; never `NA`. a: the per-trial log-odds
  at latent 0 and depth 0.

- `link$slope` (wire `slope`):

  float; in **natural log-odds per latent unit**; never `NA`. b: how
  fast the per-trial log-odds follow the latent; 1 for the identity
  link, below 1 when the chance of a success rises more slowly than the
  latent.

- `link$depth_slope` (wire `depth_slope`):

  float; in **natural log-odds per depth unit**; never `NA`. c: the
  depth term's coefficient; 0 when the link has none.

- `link$overdispersion` (wire `overdispersion`):

  float; in **fraction (intra-cell correlation, 0 to \< 1)**; never
  `NA`. rho, the correlation of trials within one cell (beta-binomial);
  0 is the binomial.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  neither `--responses` nor `--counts`; design, groups or link missing;
  a file given without a value or not found; a table empty, header-only
  or with too few columns; a row whose field count differs from its
  header; an empty or repeated feature id, sample or coefficient name; a
  response cell that is not a number (other than blank, NA, NaN); a
  count cell that is not two whole numbers k/n; a design cell that is
  not a finite number

- `mzlib_usage_error` (usage):

  the counts' features or samples are not exactly the responses' (each
  difference named); the design's, the groups' or the depth file's
  samples are not exactly the fit's; the groups or depth file with other
  than two columns; a depth that is not a number; a residual-variance id
  that is not a feature, or a value that is not a number; no coefficient
  on stdin, one named twice, or one not in the design; `--link` not
  identity or 2 to 4 numbers; `--estimator` other than reml or ml;
  `--threads` 0 or below -1

- `mzlib_usage_error` (usage):

  mzLib refuses the input (its ArgumentException, reclassified): a cell
  both measured and counted, or a count with no trials or more successes
  than trials, named by 0-based position and then by feature and sample
  name; a link with a slope of 0, a non-finite value or overdispersion
  outside \[0, 1); a link with a depth term and no `--depth`; a depth
  that is not finite; a given residual variance that is not positive and
  finite; a design not of full column rank; an empty group label

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- A random intercept per group on a latent scale: a measured sample
  reads the latent; a counted sample has k of n trials succeed with
  logit g = intercept + slope x latent + depth_slope x depth, binomially
  or beta-binomially. No random slopes, nested or crossed groupings.

- Responses and estimates are on the latent scale (natural log-odds for
  a logit latent), not fractions or intensities. Which cells are
  counted, and the link, are the caller's modelling choices; the bridge
  builds neither (stats fit-link fits a link).

- Counted samples are integrated by adaptive Gauss-Hermite quadrature
  (40 nodes over the group intercept, 12 over each count's residual),
  measured ones exactly. mzLib agrees with lme4's glmer (nAGQ 25, ML) to
  about 1e-6, with lmer to about 1e-5 when nothing is counted, and with
  glmmTMB to its Laplace error (about 1% in tau^2 and standard errors,
  5% in sigma^2) when sigma^2 is estimated from counts alone
  (ReferenceComparisonLogitNormalTests).

- Fixed effects are tested by Wald t on nlme's containment df, over
  samples with a response of either kind. standard_error is at the
  estimated variances and does not include their uncertainty, so it is
  not glmmTMB's vcov().

- Under reml the fixed effects are integrated out by a Laplace
  approximation (glmmTMB's REML = TRUE); it is exact when nothing is
  counted. Under REML glmmTMB reports beta as the joint mode of (beta,
  u), so compare its estimates with this verb's ml ones only.

- With nothing counted the model is stats fit-mixed's, fitted by a
  different optimiser: under reml the numbers agree with fit-mixed's to
  about 1e-5; under ml estimates and variances agree, but standard
  errors are at the ML variances, where fit-mixed's scale by RSS /
  (observed - coefficients) as nlme's summary does.

- separated is a status, never a p-value: the counts are
  (quasi-)completely separated by the design and the estimate runs to
  infinity.

- p_value is not adjusted across features. Send one coefficient's
  p_value column to stats adjust for Benjamini-Hochberg.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.fit_latent_mixed`

- Rust (mzLibRust): `mzlib::stats::fit_latent_mixed_with` with
  `FitLatentMixedOptions`

- R (mzLibR): `stats_fit_latent_mixed`

## Since

Wire protocol 1; pyMzLib 0.7.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284):
  generalized linear mixed models, the class of model this verb fits: a
  latent linear predictor with a random intercept, observed through a
  binomial link (Breslow and Clayton, J Am Stat Assoc 88:9, 1993)

- [doi:10.1093/biomet/81.3.624](https://doi.org/10.1093/biomet/81.3.624):
  adaptive Gauss-Hermite quadrature, centred at the integrand's mode and
  scaled by its curvature, which mzLib uses to integrate the counted
  samples (Liu and Pierce, Biometrika 81:624, 1994)

- [doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01):
  lme4's lmer and glmer (nAGQ 25), the references mzLib's fits are held
  to (Bates, Maechler, Bolker and Walker, J Stat Softw 67:1, 2015)

- [doi:10.32614/RJ-2017-066](https://doi.org/10.32614/RJ-2017-066):
  glmmTMB, whose Laplace-approximated REML mzLib's REML follows and
  whose fits with run scatter are a reference (Brooks et al., R Journal
  9:378, 2017)

- [doi:10.1093/biomet/58.3.545](https://doi.org/10.1093/biomet/58.3.545):
  restricted maximum likelihood, the default estimator (Patterson and
  Thompson, Biometrika 58:545, 1971)

- [doi:10.1111/j.2517-6161.1948.tb00014.x](https://doi.org/10.1111/j.2517-6161.1948.tb00014.x):
  the beta-binomial distribution, which a link with overdispersion gives
  a count cell (Skellam, J R Stat Soc B 10:257, 1948)

## See also

[`stats_fit_link`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md)
for the link,
[`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)
when every sample is measured

## Examples

``` r

# mzLib's own reference input: 16 simulated features, 23 samples in 8 animals, every sample
# counted, read through the link (-0.4, 0.6, 0.3) with a depth term; and lme4's glmer() fit of it
# (25-point adaptive quadrature, maximum likelihood), which this package ships.
data <- system.file("extdata", "stats", package = "mzLibR")
fit <- stats_fit_latent_mixed(file.path(data, "lnmm_reference_design.tsv"),
  file.path(data, "lnmm_reference_groups.tsv"), c("Intercept", "trait", "x"),
  link = c(-0.4, 0.6, 0.3), counts = file.path(data, "lnmm_link_counts.tsv"),
  depth = file.path(data, "lnmm_link_depth.tsv"),
  residual_variance = file.path(data, "lnmm_link_residual_variance.tsv"), estimator = "ml")
fit
#> <mzlibr_latent_mixed_model_fit> 16 features x 23 samples in 8 groups (group), ml, reporting Intercept, trait, x
#>   0 measured and 338 counted cells, through logit g = -0.4 + 0.6 x latent + 0.3 x depth, overdispersion 0
#>   fitted 16, too_few_observations 0, rank_deficient 0, separated 0, not_converged 0, too_few_groups 0
trait <- fit$records[fit$records$coefficient == "trait", ]
glmer <- read.delim(file.path(data, "lnmm_counted_fit.tsv"))
glmer <- glmer[glmer$link == "link" & glmer$method == "ML", ]
# Every trait estimate is within 1e-5 standard errors of glmer's:
max(abs(trait$estimate - glmer$b1) / glmer$se1) < 1e-5
#> [1] TRUE
```
