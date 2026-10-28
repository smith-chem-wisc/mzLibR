# Fit one linear model per feature and test coefficients, and contrasts, with limma's moderated t

Fit one linear model per feature of a feature-by-sample table and test
each named coefficient, and each contrast of a contrasts file with its
confidence interval, with limma's empirical-Bayes moderated t,
Benjamini-Hochberg adjusted, in one call.

Runs mzLib's `LinearModel.Fit` once - every design column is fitted -
and then `EmpiricalBayes.Moderate` for each coefficient you name, which
also applies Benjamini-Hochberg across that coefficient's fitted
features. This is limma's `lmFit()` then `eBayes(legacy = TRUE)` (Smyth
2004), which mzLib holds to limma to 1e-8 relative on limma's own
reference output. With `contrasts`, `EmpiricalBayes.ModerateContrasts`
also tests each row of a contrasts file - "250 ng against 125 ng", "the
mean of two groups against a third" - with its confidence interval:
limma's `contrasts.fit()` then `eBayes()` then
`topTable(confint = ...)`.

## Usage

``` r
stats_fit(responses, design, coefficients = character(0), contrasts = NULL,
  confidence = NULL, estimator = NULL, trend = FALSE, spline_basis = NULL, threads = 1,
  timeout = 600)
```

## Arguments

- responses:

  A TSV with a header row. First column the feature id; every other
  column one sample, named by its header. Use the values as you want
  them modelled - log2 intensities, usually. Blank, `NA` or `NaN` is
  missing and left out of that feature's fit, never imputed; **0 is an
  observation**, so blank out any 0 that means "not measured".

- design:

  A TSV with a header row and one row per sample: first column the
  sample name (matching the responses header exactly, in any order),
  then one numeric column per coefficient. Include an intercept column
  yourself if you want one.
  [`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
  writes one (`design_<i>.tsv`) from a sample sheet.

- coefficients:

  Design column names to test, each with its own moderated test and its
  own Benjamini-Hochberg family. May be empty (`character(0)`) when you
  give `contrasts`.

- contrasts:

  A TSV of contrasts, as
  [`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
  writes it beside the design (`contrasts_<i>.tsv`): a header of the
  contrast id column, then exactly the design's coefficient columns (in
  any order); then one row per contrast, its id and one weight per
  coefficient, e.g. `-1` and `1` for the second level minus the first.
  Each contrast gets its own moderated test, Benjamini-Hochberg family
  and confidence interval, in `contrasts`. `NULL` sends none. Needs the
  pyMzLib 0.6.0 bridge.

- confidence:

  The two-sided level of the contrasts' intervals, as a fraction (0 to
  1), strictly between 0 and 1; only with `contrasts`. Leave it `NULL`
  for mzLib's default (0.95 at mzLib `d10017ce`);
  `contrasts$confidence_level` says which level was used.

- estimator:

  How the variance prior is estimated, one of `"moments_legacy"`
  (limma's `eBayes(legacy = TRUE)`) or `"marginal_likelihood"` (mzLib's
  own estimator, which uses each feature's own residual degrees of
  freedom, so it handles the unequal df that missing values cause
  exactly). Leave it `NULL` for mzLib's default (`"moments_legacy"` at
  mzLib `d10017ce`); `prior$estimator` says which was used.

- trend:

  Let the prior variance follow average response (limma's
  `trend = TRUE`), for data whose low-abundance features are noisier.

- spline_basis:

  Basis functions of the trend spline, intercept included. Only with
  `trend = TRUE`. Leave it `NULL` to let mzLib choose, as limma does,
  for the number of features; mzLib may also use fewer than you ask for.
  Either way `prior$spline_basis_count` reports what it used.

- threads:

  Features fitted at once (threads), or `-1` for every core. The answer
  is identical at any value.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

A feature is fitted on the samples where it was observed, and is
reported - never dropped - when it cannot be: `"too_few_observations"`
when it has no more observed samples than the design has coefficients,
`"rank_deficient"` when its observed samples cannot separate the
coefficients (a feature seen in only one group cannot have a group
effect).

## Value

An `mzlibr_moderated_fit`:

\- `records`, a data.frame with one row per (tested coefficient,
feature) - one block of `feature_count` rows per coefficient, in the
order you named them, and no rows when you named none - with `feature`,
`coefficient`, `status`, `estimate` (response units per coefficient
unit: a **log2 fold change** for a 0/1 group column on log2
intensities), `standard_error` (response units per coefficient unit),
`t` (standard errors), `df_total` (degrees of freedom), `p_value` and
`bh_adjusted` (fractions, 0 to 1), `posterior_variance` and
`prior_variance` (squared response units), `sigma` (response units),
`df_residual`, `observed` (samples) and `average_response` (response
units). Statistics are `NA` where `status` is not `"fitted"`. -
`responses_file`, `design_file` and `contrasts_file` (`NA` when you gave
no `contrasts`), as absolute paths; `feature_count` features and
`sample_count` samples read; `sample_names`, `coefficient_names` (every
design column, in the order each contrast's weights follow) and
`tested_coefficients` (empty when you named none); `missing_count` and
`zero_count` response cells; `status_counts`, features per status;
`threads`; `row_count` rows. - `residual_df_differ`: `TRUE` when the
fitted features do not all have the same residual degrees of freedom.
Then default limma (3.61 and later) uses a different prior estimator and
reports different moderated statistics; under `"moments_legacy"` these
are `legacy = TRUE`, and `estimator = "marginal_likelihood"` is mzLib's
estimator for this case. - `prior`, the fitted variance prior, shared by
every coefficient and contrast: `estimator`, the estimator used (`NA`
only from a bridge before pyMzLib 0.6.0, which did not report it); `df`
(degrees of freedom, `NA` when infinite; under `"marginal_likelihood"`
mzLib searches d0 from about 0.05 upward, so a d0 there is that floor,
not an optimum), `df_infinite`, `trended`, `spline_basis_count` (basis
functions) and `scale` (squared response units), s0^2 when every fitted
feature shares one. `scale` is `NA` when it differs by feature, which is
then in the `prior_variance` column: with a trend, or with fewer than 3
features. Then mzLib fits no prior (mzLib \#1398, as limma's
`squeezeVar` does): `df` is 0 and each feature keeps its own variance.
Bridges before pyMzLib 0.6.0 report the first feature's own variance as
`scale` in that case. - `contrasts`, `NULL` when you gave no
`contrasts`; otherwise an `mzlibr_moderated_contrasts`:
`confidence_level`, the two-sided level of every interval as a fraction
(0 to 1); `tested_contrasts`, the contrast ids in the file's order;
`weights`, a list from contrast id to its weights (contrast weight per
coefficient, unitless), one per `coefficient_names` entry, whatever
order the file's columns were in; `row_count` rows, `feature_count` x
the contrasts; `column_names`; and `records`, one row per (contrast,
feature), a block per contrast in file order, with `feature`, `contrast`
(its id), `status`, `estimate` (response units per contrast unit: a
**log2 fold change** for one level against another on log2 intensities),
`standard_error`, `confidence_low` and `confidence_high` (response units
per contrast unit; the two-sided `confidence_level` interval, not
adjusted across features), and the coefficient table's `t`, `df_total`,
`p_value`, `bh_adjusted` (over this contrast's fitted features),
`posterior_variance`, `prior_variance`, `sigma`, `df_residual`,
`observed` and `average_response`, in the same units. Statistics are
`NA` where `status` is not `"fitted"`. - `caveats`, including any that
apply to this input: the estimator, unfitted features, zeros,
`residual_df_differ`, contrasts.

## What it is not

By default it is limma's \*legacy\* estimator;
`estimator = "marginal_likelihood"` is mzLib's own, close to default
limma's but not equal. `robust = TRUE`, the B-statistic and observation
weights are not implemented. With missing values a contrast's standard
error is exact per feature, where limma's `contrasts.fit()` approximates
it, so on such data the contrasts differ from limma's by design.
`bh_adjusted` controls the false discovery rate among the features
tested; it is not a target-decoy q-value.

## Bridges before pyMzLib 0.6.0

Such a bridge has `stats fit` but ignores `--contrasts`, `--confidence`
and `--estimator`. Rather than return a fit that silently did not do
what you asked, `stats_fit()` raises an `mzlib_usage_error` naming the
0.6.0 bridge when you gave `contrasts` or `estimator` and the answer
lacks what they add.

## Wraps

Wire verb `stats fit`. Generated from the bridge's verb spec
`stats.fit.yaml` (bridge commit `60b1d25968a1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`LinearModel.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/LinearModel.cs)
  in `mzLib/StatisticalModels/LinearModel.cs` at mzLib `0a808fec`

- [`EmpiricalBayes.Moderate`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `d10017ce`

- [`EmpiricalBayes.ModerateContrasts`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `d10017ce`

- [`ContrastWeights`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `d10017ce`

- [`ModeratedContrast`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `d10017ce`

- [`VariancePriorEstimator`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `d10017ce`

- [`EmpiricalBayes.FitPrior`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `0a808fec`

- [`EmpiricalBayes.DefaultSplineBasisCount`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `d10017ce`

- [`MultipleTesting.BenjaminiHochberg`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/MultipleTesting.cs)
  in `mzLib/StatisticalModels/MultipleTesting.cs` at mzLib `0a808fec`

## Parameters: units, ranges and defaults

- `responses`:

  path; required. A TSV with a header row. The first column is the
  feature id (unique, any header); every other column is one sample,
  named by its header. Cells are used as given, so log-transform
  intensities first. A blank cell, NA or NaN (any case) is missing and
  is omitted from that feature's fit, never imputed. A 0 is an
  OBSERVATION. Blank lines are skipped.

- `design`:

  path; required. A TSV with a header row and one row per sample. The
  first column is the sample name and must match the responses header
  exactly (rows are matched by name, so their order does not matter);
  every other column is one coefficient, named by its header, and every
  cell must be a finite number. The design is used as written: include
  an intercept column if you want one. diff design `--out` writes one
  (design_i.tsv).

- `coefficients` (wire `--stdin`):

  string\[\]; default absent; range
  `>= 1 line, or none with --contrasts; each a design column name, none repeated`.
  The coefficients to test, one design column name per line. Each gets
  its own moderated test and its own Benjamini-Hochberg family. The
  model always fits every design column. With `--contrasts` it may be
  empty, and the coefficient table then has no rows.

- `contrasts`:

  path; default absent. A TSV of contrasts to test, as diff design
  `--out` writes it (contrasts_i.tsv, beside its design_i.tsv). A header
  row: the contrast id column (any header; diff design writes contrast),
  then exactly the design's coefficient columns, matched by name, so in
  any order. Then one row per contrast: a non-empty id (diff design's
  c1, c2, ...), then one finite weight per coefficient; the contrast is
  the weighted sum of the coefficients, e.g. -1 and 1 for the second
  minus the first. A design column the file lacks, or a column the
  design lacks, is refused. mzLib refuses a contrast whose weights are
  all 0, and two contrasts with one id. Each contrast is moderated with
  the same prior as the coefficients (EmpiricalBayes.ModerateContrasts)
  and has its own Benjamini-Hochberg family; its rows are the contrasts
  table. Absent, there is no contrasts table. Needs the pyMzLib 0.6.0
  bridge; an older bridge ignores the option.

- `confidence`:

  float; in **fraction (0 to 1)**; default absent; range
  `strictly between 0 and 1 (mzLib's rule); only with --contrasts`. The
  two-sided level of each contrast's confidence interval. Absent, the
  bridge sends none and ModerateContrasts' own default applies
  (confidenceLevel = 0.95 at d10017ce); contrasts.confidence_level
  reports the level used either way. A value mzLib refuses is a usage
  error with mzLib's message. Needs the pyMzLib 0.6.0 bridge.

- `estimator`:

  string; default absent; range `moments_legacy | marginal_likelihood`.
  How the variance prior (d0, s0^2) is estimated, mzLib's
  VariancePriorEstimator by a wire name. moments_legacy is
  MomentsLegacy: the method of moments on log variances (Smyth 2004),
  limma's eBayes(legacy = TRUE), exact when every feature has the same
  residual df. marginal_likelihood is MarginalLikelihood: maximum
  marginal likelihood with each feature's own residual df, which handles
  the unequal df that omitting missing values causes exactly (mzLib
  \#1368); it is mzLib's own estimator, not limma's default
  fitFDistUnequalDF1, so its numbers are close to default limma's but
  not equal. The one estimator serves the coefficients and the
  contrasts. Absent, the bridge sends none and mzLib's default applies
  (MomentsLegacy at d10017ce, kept for limma legacy parity);
  prior.estimator reports the estimator used either way. Needs the
  pyMzLib 0.6.0 bridge.

- `trend`:

  flag; default `FALSE`. Let the prior variance trend with each
  feature's average response (limma's trend = TRUE), as a natural cubic
  spline.

- `spline_basis` (wire `--spline-basis`):

  int; in **basis functions**; default absent; range
  `>= 1; only with --trend`. Basis functions of the trend spline,
  intercept included. Absent, the bridge sends none and mzLib chooses by
  limma's fitFDist rule (EmpiricalBayes.DefaultSplineBasisCount, mzLib
  1.0.594 \#1398): 1 + (n \>= 3) + (n \>= 6) + (n \>= 30) for n fitted
  features, so 4 only from 30 features. A given value is capped at the
  fitted features minus one. mzLib caps either at the distinct average
  responses and may merge coinciding knots; prior.spline_basis_count
  reports what was used.

- `threads`:

  int; in **threads**; default `1`; range `>= 1, or -1 for every core`.
  Features fitted at once. The bridge default is 1 (BULK.md); mzLib's
  own default is -1. Results are identical at any value.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `responses_file`:

  string; never `NA`. Absolute path of the responses table.

- `design_file`:

  string; never `NA`. Absolute path of the design table.

- `contrasts_file`:

  string; `NA` when no `--contrasts` was given. Absolute path of the
  contrasts table. Absent from bridges before pyMzLib 0.6.0 (optional
  forever).

- `feature_count`:

  int; in **features**; never `NA`. Feature rows read.

- `sample_count`:

  int; in **samples**; never `NA`. Sample columns read.

- `sample_names`:

  string\[\]; never `NA`. In responses header order, which is the order
  the design rows were put in.

- `coefficient_names`:

  string\[\]; never `NA`. Every design column, in order: the order of
  each contrast's weights.

- `tested_coefficients`:

  string\[\]; never `NA`. The stdin coefficients, in order; the table
  holds one block of feature_count rows per entry. Empty when only
  `--contrasts` was given.

- `trend`:

  bool; never `NA`. Whether `--trend` was given.

- `threads`:

  int; in **threads**; never `NA`. The thread count used.

- `missing_count`:

  int; in **cells**; never `NA`. Response cells that were blank, NA or
  NaN.

- `zero_count`:

  int; in **cells**; never `NA`. Response cells that are exactly 0,
  fitted as observations.

- `status_counts`:

  map\<string,int\>; in **features**; never `NA`. Features per status:
  fitted, too_few_observations, rank_deficient.

- `residual_df_differ`:

  bool; never `NA`. True when the fitted features do not all have the
  same residual df (missing values omitted). Then default limma (\>=
  3.61, eBayes legacy = FALSE) would use a different prior estimator and
  report different moderated statistics; under moments_legacy these
  follow legacy = TRUE, and marginal_likelihood is mzLib's estimator for
  this case. A caveat says so.

- `prior`:

  object; never `NA`. The fitted variance prior, shared by every tested
  coefficient and contrast; fields under result.tables.prior.

- `row_count`:

  int; in **rows**; never `NA`. feature_count x the number of tested
  coefficients.

- `column_names`:

  string\[\]; never `NA`. The coefficient table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (tested coefficient, feature):
  coefficient blocks in stdin order, features in file order. The same
  with or without `--contrasts`.

- `contrasts`:

  object; `NA` when no `--contrasts` was given. The tested contrasts, as
  a second table beside the coefficient table: fields under
  result.tables.contrasts, its rows' columns under
  result.tables.contrasts.columns. A separate table because a contrast
  row is keyed by a contrast id and carries a confidence interval, which
  a coefficient row has neither of; one table holding both would leave
  half of every row's keys null, and the coefficient table stays as it
  was without `--contrasts`. Absent from bridges before pyMzLib 0.6.0
  (optional forever).

- `caveats`:

  string\[\]; never `NA`. Three static caveats (the first names the
  prior estimator), plus one each when residual_df_differ, when any
  feature was not fitted (with counts), when any cell is 0, and when the
  prior df is infinite; and two more with `--contrasts`.

## Columns of `records`

- `feature`:

  string; never `NA`. The feature id from the responses file's first
  column.

- `coefficient`:

  string; never `NA`. The coefficient this row tests.

- `status`:

  string; never `NA`. fitted \| too_few_observations (fewer observed
  samples than coefficients plus one) \| rank_deficient (the design
  restricted to this feature's observed samples is not of full rank).
  Only fitted rows carry statistics or enter the BH family.

- `estimate`:

  float; in **response units per coefficient unit (e.g. log2 fold change
  for a 0/1 group column on log2 intensities)**; `NA` when status is not
  fitted. The least-squares coefficient; moderation does not change it.

- `standard_error`:

  float; in **response units per coefficient unit**; `NA` when not
  fitted. Moderated standard error: sqrt(posterior_variance) x the
  coefficient's unscaled standard deviation.

- `t`:

  float; in **standard errors**; `NA` when not fitted. Moderated
  t-statistic.

- `df_total`:

  float; in **degrees of freedom**; `NA` when not fitted. Degrees of
  freedom of the moderated t: df_residual + prior df, capped at the
  pooled residual df.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not fitted. Two-sided
  p-value.

- `bh_adjusted`:

  float; in **fraction (0 to 1)**; `NA` when not fitted (not in the
  family). Benjamini-Hochberg adjusted p-value over this coefficient's
  fitted features. Not a target-decoy q-value.

- `posterior_variance`:

  float; in **squared response units**; `NA` when not fitted. The
  moderated residual variance (d0 s0^2 + d s^2) / (d0 + d).

- `prior_variance`:

  float; in **squared response units**; `NA` when not fitted (the
  feature did not enter the prior). s0^2 for this feature: one value for
  every feature without `--trend`, varying with average_response with
  it. With fewer than 3 features it is the feature's own variance,
  because no prior is fitted.

- `sigma`:

  float; in **response units**; `NA` when not fitted. Residual standard
  deviation, sqrt(RSS / df_residual), before moderation.

- `df_residual`:

  int; `NA` when not fitted. Observed samples minus coefficients.

- `observed`:

  int; in **samples**; never `NA`. Samples with a finite response for
  this feature.

- `average_response`:

  float; in **response units**; `NA` when the feature has no finite
  response at all. Mean of the feature's observed responses: the
  covariate `--trend` fits on.

## Fields of `prior`

- `prior$estimator` (wire `estimator`):

  string; never `NA`. moments_legacy or marginal_likelihood: the
  estimator mzLib reports it fitted the prior with
  (VariancePrior.Estimator), its default when `--estimator` was absent.
  Absent from bridges before pyMzLib 0.6.0.

- `prior$df` (wire `df`):

  float; in **degrees of freedom**; `NA` when infinite: df_infinite is
  true. Prior degrees of freedom d0. 0 with fewer than 3 features, where
  mzLib fits no prior (mzLib#1398, as limma's squeezeVar), so moderation
  leaves each feature's variance and t as they are. Under
  marginal_likelihood mzLib searches d0 from about 0.05 (e^-3) upward,
  so a d0 there is the search's floor, not an optimum.

- `prior$df_infinite` (wire `df_infinite`):

  bool; never `NA`. True when the residual variances are no more
  dispersed than sampling alone explains, so every feature takes the
  prior variance.

- `prior$trended` (wire `trended`):

  bool; never `NA`. Whether s0^2 varies with average_response.

- `prior$spline_basis_count` (wire `spline_basis_count`):

  int; in **basis functions**; never `NA`. Basis functions the trend
  used, intercept included; 1 without a trend.

- `prior$scale` (wire `scale`):

  float; in **squared response units**; `NA` when s0^2 differs by
  feature, so it is in the prior_variance column: trended, or no prior
  was fitted (fewer than 3 features; mzLib#1398, as limma's squeezeVar).
  s0^2, the prior variance, when every fitted feature shares one. With
  fewer than 3 features mzLib fits no prior: df is 0 and each feature
  keeps its own variance. Bridges before pyMzLib 0.6.0 reported the
  first feature's own variance here in that case.

## Fields of `contrasts`

- `contrasts$confidence_level` (wire `confidence_level`):

  float; in **fraction (0 to 1)**; never `NA`. The two-sided level of
  every confidence_low and confidence_high: `--confidence`, or mzLib's
  default when it was absent (ModeratedContrast.ConfidenceLevel).

- `contrasts$tested_contrasts` (wire `tested_contrasts`):

  string\[\]; never `NA`. The contrast ids, in file order; the table
  holds one block of feature_count rows per entry.

- `contrasts$weights` (wire `weights`):

  map\<string,float\[\]\>; in **contrast weight per coefficient
  (unitless)**; never `NA`. Each contrast id's weights, one per
  coefficient_names entry, in that order (the file's columns matched by
  name): what mzLib tested.

- `contrasts$row_count` (wire `row_count`):

  int; in **rows**; never `NA`. feature_count x the number of contrasts.

- `contrasts$column_names` (wire `column_names`):

  string\[\]; never `NA`. The contrasts table's columns, in order
  (result.tables.contrasts.columns).

- `contrasts$records` (wire `columns`):

  table; never `NA`. One row per (contrast, feature): contrast blocks in
  file order, features in file order.

## Fields of `contrasts.columns`

- `feature`:

  string; never `NA`. The feature id from the responses file's first
  column.

- `contrast`:

  string; never `NA`. The contrast id this row tests, from the contrasts
  file's first column. diff design numbers its contrasts within each
  stratum; its label column says what each id compares.

- `status`:

  string; never `NA`. The feature's fit status, as in the coefficient
  table. Only fitted rows carry statistics or enter the contrast's BH
  family.

- `estimate`:

  float; in **response units per contrast unit (e.g. log2 fold change
  for a level-against-level contrast on log2 intensities)**; `NA` when
  not fitted. The contrast's least-squares estimate, the weighted sum of
  the coefficients; moderation does not change it.

- `standard_error`:

  float; in **response units per contrast unit**; `NA` when not fitted.
  Moderated standard error: sqrt(posterior_variance) x sqrt(c'(X'X)^-1
  c), with the feature's own (X'X)^-1 over its observed samples.

- `confidence_low`:

  float; in **response units per contrast unit**; `NA` when not fitted.
  Lower bound of the two-sided confidence_level interval: estimate -
  t-quantile((1 + level) / 2, df_total) x standard_error. Per feature,
  not adjusted across features.

- `confidence_high`:

  float; in **response units per contrast unit**; `NA` when not fitted.
  Upper bound: estimate + the same half-width.

- `t`:

  float; in **standard errors**; `NA` when not fitted. Moderated
  t-statistic.

- `df_total`:

  float; in **degrees of freedom**; `NA` when not fitted. df_residual +
  prior df, capped at the pooled residual df, as for a coefficient.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not fitted. Two-sided
  p-value.

- `bh_adjusted`:

  float; in **fraction (0 to 1)**; `NA` when not fitted (not in the
  family). Benjamini-Hochberg adjusted p-value over this contrast's
  fitted features. Not a target-decoy q-value.

- `posterior_variance`:

  float; in **squared response units**; `NA` when not fitted. The
  moderated residual variance, as in the coefficient table.

- `prior_variance`:

  float; in **squared response units**; `NA` when not fitted (the
  feature did not enter the prior). s0^2 for this feature, as in the
  coefficient table.

- `sigma`:

  float; in **response units**; `NA` when not fitted. Residual standard
  deviation before moderation.

- `df_residual`:

  int; `NA` when not fitted. Observed samples minus coefficients.

- `observed`:

  int; in **samples**; never `NA`. Samples with a finite response for
  this feature.

- `average_response`:

  float; in **response units**; `NA` when the feature has no finite
  response at all. Mean of the feature's observed responses.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  responses, design or contrasts missing, blank or not found (or
  `--contrasts` without a value); a table empty, header-only or with too
  few columns; a row whose field count differs from its header; an empty
  or repeated feature id, sample, coefficient name or contrasts column;
  a response cell that is not a number (other than blank, NA, NaN); a
  design cell or contrast weight that is not a finite number

- `mzlib_usage_error` (usage):

  the design's samples are not exactly the responses' sample columns, or
  the contrasts file's columns are not exactly the design's coefficients
  (both differences are named); an empty contrast id; no coefficient on
  stdin without `--contrasts`, one named twice, or one not in the
  design; `--spline-basis` without `--trend`, or below 1; `--confidence`
  without `--contrasts`, without a value or not a number; `--estimator`
  other than moments_legacy or marginal_likelihood, or without a value;
  `--threads` 0 or below -1

- `mzlib_usage_error` (usage):

  mzLib refuses the input (its ArgumentException, reclassified, with
  mzLib's message): the design is not of full column rank; fewer than 2
  features could be fitted, so no variance prior exists; a contrast
  whose weights are all 0; two contrasts with one id; a confidence level
  not strictly between 0 and 1

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Under the default estimator, moments_legacy, this is limma's LEGACY
  empirical-Bayes estimator (eBayes(legacy = TRUE)). Since limma 3.61,
  eBayes defaults to a different estimator whenever residual df differ
  between features, which omitting missing values causes;
  residual_df_differ says when that applies. For that case mzLib offers
  marginal_likelihood, which handles unequal df exactly
  (VariancePriorEstimator.MarginalLikelihood, 'this library's answer to
  that gap' in mzLib's EmpiricalBayes docs). It is mzLib's own
  estimator, not limma's fitFDistUnequalDF1, so its numbers differ from
  default limma's too; mzLib keeps moments_legacy the default for limma
  legacy parity.

- Not implemented: robust = TRUE, the B-statistic, observation weights.

- Contrasts (`--contrasts`) use the same prior as the coefficients, and
  each has its own Benjamini-Hochberg family. confidence_low and
  confidence_high are per feature, not adjusted across features. With
  missing values a contrast's standard error is exact per feature, where
  limma's contrasts.fit approximates it, so on such data the contrasts
  differ from limma's by design; with complete data they agree, and a
  contrast of one coefficient gives exactly that coefficient's row.

- A 0 is an observation. Log-transform intensities first, and blank any
  cell that means 'not measured'.

- bh_adjusted adjusts for the false discovery rate among the features
  tested. It is not a target-decoy q-value for identifications.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.fit`

- Rust (mzLibRust): `mzlib::stats::fit_with` with `FitOptions`

- R (mzLibR): `stats_fit`

## Since

Wire protocol 1; pyMzLib 0.3.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## References

- [doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027):
  the empirical-Bayes moderated t this verb computes, for coefficients
  and contrasts (Smyth, Stat Appl Genet Mol Biol 3:3, 2004)

- [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007): limma,
  whose eBayes(legacy = TRUE), and contrasts.fit with topTable(confint),
  mzLib reproduces to 1e-8 relative (Ritchie et al., Nucleic Acids Res
  43:e47, 2015)

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the Benjamini-Hochberg adjustment in bh_adjusted (J R Stat Soc B
  57:289, 1995)

## See also

[`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md),
which writes the design and contrasts files;
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md),
[`stats_meta`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md),
[`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)
for samples in groups

## Examples

``` r

# A dilution series from mzLib's own RNA test data: MALAT1 at 500, 250 and 125 ng against a
# constant FLuc spike, log2 intensities normalised to FLuc. Halving MALAT1 should read as about -1.
data <- system.file("extdata", "stats", package = "mzLibR")
fit <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
  file.path(data, "malat_dilution_design.tsv"), c("malat_250ng", "malat_125ng"))
fit
#> <mzlibr_moderated_fit> 492 features x 27 samples, testing malat_250ng, malat_125ng
#>   fitted 212, too_few_observations 214, rank_deficient 66
#>   prior df 19.56 (moments_legacy)
#>   ! residual df differ between features: these are eBayes(legacy = TRUE) statistics
fit$status_counts
#>               fitted too_few_observations       rank_deficient 
#>                  212                  214                   66 
c(estimator = fit$prior$estimator, prior_df = round(fit$prior$df, 1))
#>        estimator         prior_df 
#> "moments_legacy"           "19.6" 
malat <- fit$records[fit$records$coefficient == "malat_250ng" & fit$records$status == "fitted" &
  startsWith(fit$records$feature, "MALAT1:"), ]
c(oligos = nrow(malat), median_log2_fold_change = median(malat$estimate))
#>                  oligos median_log2_fold_change 
#>             116.0000000              -0.8914901 

# Every pair of loads, with intervals: diff_design() with no reference level writes one contrast
# per pair beside the design, and `contrasts` tests them all.
samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
design <- diff_design(samples, factors = list(factor_spec("malat_load")),
  out = file.path(tempdir(), "malat_pairs"))
stratum <- design$strata[[1]]
pairs <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"), stratum$design_file,
  contrasts = stratum$contrasts_file)
pairs$contrasts
#> <mzlibr_moderated_contrasts> 3 contrast(s) x 492 features, 95% intervals
#>   c1: 0 1 0
#>   c2: 0 0 1
#>   c3:  0 -1  1
pairs$contrasts$weights$c3
#> [1]  0 -1  1
pairs$coefficient_names
#> [1] "intercept"        "malat_load=250ng" "malat_load=500ng"
```
