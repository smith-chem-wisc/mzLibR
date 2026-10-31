# Fit one linear mixed model per feature, with a random intercept for each group of samples

Fit one linear mixed model per feature of a feature-by-sample table,
with a random intercept per group of samples (dataset, donor, day), and
test each named coefficient with a Wald t.

mzLib's `MixedModel.Fit` models y = X beta + u_group + e for every
feature separately: the fixed effects in your design, plus one baseline
per group (a dataset, a donor, an acquisition day) drawn from N(0,
tau^2). Each coefficient is tested by a two-sided Wald t on nlme's
"containment" degrees of freedom. It is
`nlme::lme(y ~ ..., random = ~ 1 | group)` per feature, which mzLib
reproduces to 1e-5 relative (nlme's own precision), with identical df.

## Usage

``` r
stats_fit_mixed(responses, design, groups, coefficients, estimator = NULL, threads = 1,
  timeout = 600)
```

## Arguments

- responses:

  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  responses table: a TSV with a header row, the feature id first, then
  one column per sample, named by its header. Log-transform intensities
  first. Blank, `NA` or `NaN` is missing; **0 is an observation**.

- design:

  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  design: one row per sample, the sample name then one numeric column
  per fixed-effect coefficient. Do not add a column per group: the
  random intercept replaces it.

- groups:

  A TSV with a header row and one row per sample: the sample name
  (exactly the responses' samples, in any order), then ONE column with
  the sample's group label, whose header names the grouping (`day`,
  `donor`, `dataset`).

- coefficients:

  Design column names to report, each as one block of rows. Every design
  column is fitted whatever you name.

- estimator:

  `"reml"` or `"ml"`: how the variances are estimated. Leave it `NULL`
  for mzLib's default (`"reml"` at mzLib `d10017ce`); `estimator` in the
  result says which was used. `"ml"` gives log-likelihoods you can
  compare between designs.

- threads:

  Features fitted at once (threads), or `-1` for every core. The answer
  is identical at any value.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

Use it when samples come in groups that share a baseline - runs acquired
on the same day, replicates of one donor, studies in a pooled analysis -
and you want each coefficient's effect \*within\* groups.
[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
treats every sample as independent, so a day-to-day shift there either
inflates the residual variance or, if the days are unbalanced, leaks
into the effect.

## Value

An `mzlibr_mixed_model_fit`:

\- `records`, a data.frame with one row per (reported coefficient,
feature) - one block of `feature_count` rows per coefficient you named,
in that order: `feature`; `coefficient`; `status`; `estimate`, the fixed
effect within groups in response units per coefficient unit (a **log2
fold change** for a 0/1 column on log2 intensities); `standard_error`
(response units per coefficient unit); `df` (degrees of freedom), nlme's
containment df, a whole number, `NA` when not fitted or when no positive
df exists (a group-level coefficient with too few groups); `t` (standard
errors), the Wald statistic; `p_value`, a fraction (0 to 1), two-sided
on `df`, **not adjusted across features**, `NA` also when `df` is;
`residual_variance` (squared response units), sigma^2 within groups;
`group_variance` (squared response units), tau^2 between the groups'
baselines, where 0 is the boundary (no group effect); `log_likelihood`
(log-likelihood), restricted under `"reml"`; `observed`, samples with a
finite response; and `groups`, distinct groups among them. Statistics
are `NA` where `status` is not `"fitted"`. - `responses_file`,
`design_file` and `groups_file`, as absolute paths; `grouping`, the
groups file's label column header; `group_names`, in order of first
appearance; `group_count` groups; `feature_count` features and
`sample_count` samples read; `sample_names`; `coefficient_names` and
`tested_coefficients`; `estimator`, the estimator mzLib used; `threads`;
`missing_count` and `zero_count` response cells (cells);
`status_counts`, features per status; `row_count` rows;
`column_names`. - `caveats`: what the model is and is not, plus counts
for this input: unfitted features, rows with no df, features at the
boundary, zero cells.

## What it is not

A random intercept per group and nothing more: no random slopes, nested
or crossed groupings, or residual correlation. Satterthwaite and
Kenward-Roger df (lmerTest, pbkrtest) are not implemented, so p-values
differ from lmerTest's. A feature the fixed effects fit exactly is
reported fitted, with a `residual_variance` at rounding level and a
meaningless `t` and `p_value`: check `residual_variance` before trusting
a very small `p_value`.

## Wraps

Wire verb `stats fit-mixed`. Generated from the bridge's verb spec
`stats.fit-mixed.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`MixedModel.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/MixedModel.cs)
  in `mzLib/StatisticalModels/MixedModel.cs` at mzLib `d10017ce`

- [`MixedModelFit`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/MixedModel.cs)
  in `mzLib/StatisticalModels/MixedModel.cs` at mzLib `d10017ce`

- [`VarianceEstimator`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/MixedModel.cs)
  in `mzLib/StatisticalModels/MixedModel.cs` at mzLib `d10017ce`

- [`FeatureFitStatus`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/LinearModel.cs)
  in `mzLib/StatisticalModels/LinearModel.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `responses`:

  path; required. stats fit's responses table. A TSV with a header row:
  the first column is the feature id (unique, any header) and every
  other column is one sample, named by its header. Cells are used as
  given, so log-transform intensities first. A blank cell, NA or NaN
  (any case) is missing and is omitted from that feature's fit, never
  imputed. A 0 is an OBSERVATION. Blank lines are skipped.

- `design`:

  path; required. stats fit's design table: a TSV with a header row and
  one row per sample, the sample name (matching the responses header
  exactly, rows matched by name) then one column per fixed-effect
  coefficient, every cell a finite number. Include an intercept column
  if you want one; do not add a column per group, which the random
  intercept replaces. mzLib refuses a design that is not of full column
  rank over all samples.

- `groups`:

  path; required. A TSV with a header row and one row per sample: the
  sample name, then EXACTLY ONE column, the sample's group label (a
  dataset accession, a donor, an acquisition day), whose header names
  the grouping. The sample names must be exactly the responses' sample
  columns, in any order. Labels are text, trimmed and compared exactly;
  NA is a label like any other. An empty label is refused by mzLib
  (MixedModel.Fit), naming the sample.

- `coefficients` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line; each a design column name, none repeated`. The
  coefficients to report, one design column name per line. The model
  always fits every design column; the table holds one block of rows per
  coefficient named here, in this order.

- `estimator`:

  string; default absent; range `reml | ml`. How the variances are
  estimated: reml is mzLib's VarianceEstimator.Reml (restricted maximum
  likelihood, not biased downward by the number of coefficients), ml is
  VarianceEstimator.MaximumLikelihood (biased variances in small
  samples, but log-likelihoods comparable between designs). Absent, the
  bridge sends none and MixedModel.Fit's default applies (Reml at
  d10017ce); the envelope's estimator says which was used.

- `threads`:

  int; in **threads**; default `1`; range `>= 1, or -1 for every core`.
  Features fitted at once. The bridge default is 1 (BULK.md); mzLib's
  own default is -1. Results are identical at any value (mzLib).

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `responses_file`:

  string; never `NA`. Absolute path of the responses table.

- `design_file`:

  string; never `NA`. Absolute path of the design table.

- `groups_file`:

  string; never `NA`. Absolute path of the groups table.

- `grouping`:

  string; never `NA`. The groups file's label column header.

- `group_names`:

  string\[\]; never `NA`. The distinct labels, in order of first
  appearance over the responses' samples (the order mzLib numbers them
  in).

- `group_count`:

  int; in **groups**; never `NA`. Distinct labels over all samples. A
  feature's own count, over its observed samples, is its groups column.

- `feature_count`:

  int; in **features**; never `NA`. Feature rows read.

- `sample_count`:

  int; in **samples**; never `NA`. Sample columns read.

- `sample_names`:

  string\[\]; never `NA`. In responses header order, which is the order
  the design and group rows were put in.

- `coefficient_names`:

  string\[\]; never `NA`. Every design column, in order.

- `tested_coefficients`:

  string\[\]; never `NA`. The stdin coefficients, in order; the table
  holds one block of feature_count rows per entry.

- `estimator`:

  string; never `NA`. reml or ml: the estimator mzLib reports it fitted
  with (MixedModelFit.Estimator), the default when `--estimator` was
  absent.

- `threads`:

  int; in **threads**; never `NA`. The thread count used.

- `missing_count`:

  int; in **cells**; never `NA`. Response cells that were blank, NA or
  NaN.

- `zero_count`:

  int; in **cells**; never `NA`. Response cells that are exactly 0,
  fitted as observations.

- `status_counts`:

  map\<string,int\>; in **features**; never `NA`. Features per status,
  always these six keys: fitted, too_few_observations, rank_deficient,
  too_few_groups, variance_ratio_at_limit, not_converged. They are the
  outcomes MixedModel.Fit's body can return at d10017ce; a stop-gap
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

  string\[\]; never `NA`. Seven static caveats (one of them says which
  estimator), plus one each, with counts, when any feature was not
  fitted, when a fitted row has no df, when a fitted feature's
  group_variance is 0, and when any cell is 0.

## Columns of `records`

- `feature`:

  string; never `NA`. The feature id from the responses file's first
  column.

- `coefficient`:

  string; never `NA`. The coefficient this row reports.

- `status`:

  string; never `NA`. fitted \| too_few_observations (no more observed
  samples than coefficients) \| rank_deficient (the design restricted to
  the observed samples is not of full rank) \| too_few_groups (fewer
  than two observed groups, or no group observed twice, so group and
  residual variance cannot be told apart) \| variance_ratio_at_limit
  (group_variance / residual_variance ran to the top of mzLib's search,
  e^12, typically near-identical replicates within each group) \|
  not_converged (the likelihood at the chosen variance ratio could not
  be evaluated: a non-finite value, or a residual sum of squares of 0 or
  below). Only fitted rows carry statistics.

- `estimate`:

  float; in **response units per coefficient unit (e.g. log2 fold change
  for a 0/1 column on log2 intensities)**; `NA` when status is not
  fitted. The generalized-least-squares fixed-effect coefficient at the
  estimated variances: the effect within groups, each group keeping its
  own baseline.

- `standard_error`:

  float; in **response units per coefficient unit**; `NA` when not
  fitted. sqrt of the diagonal of s^2 (X' V^-1 X)^-1 with s^2 = RSS /
  (observed - coefficients), under both estimators, as nlme's summary
  reports.

- `df`:

  float; in **degrees of freedom**; `NA` when not fitted, or the
  containment df is not positive (a group-level coefficient with too few
  observed groups to test it). Denominator df of the t test, by nlme's
  containment rule for one grouping level, a whole number: observed -
  observed groups - within-group coefficients for a coefficient that
  varies within some group and for the intercept; observed groups -
  group-level coefficients, minus 1 more with an intercept, for one
  constant within every group.

- `t`:

  float; in **standard errors**; `NA` when not fitted. Wald statistic,
  estimate / standard_error. Reported even when df is null.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not fitted, or df is null.
  Two-sided p-value of the t test on df. Not adjusted across features.

- `residual_variance`:

  float; in **squared response units**; `NA` when not fitted. sigma^2,
  the within-group variance: RSS / (observed - coefficients) under reml,
  RSS / observed under ml. The same in every coefficient block.

- `group_variance`:

  float; in **squared response units**; `NA` when not fitted. tau^2, the
  between-group variance of the random intercept. 0 is the boundary
  estimate (no between-group variance), which nlme cannot reach and
  reports as a tiny positive value. The same in every coefficient block.

- `log_likelihood`:

  float; in **log-likelihood (natural log)**; `NA` when not fitted. The
  maximized log-likelihood, restricted under reml. The same in every
  coefficient block.

- `observed`:

  int; in **samples**; never `NA`. Samples with a finite response for
  this feature.

- `groups`:

  int; in **groups**; never `NA`. Distinct groups among this feature's
  observed samples.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  responses, design or groups missing, blank or not found; a table
  empty, header-only or with too few columns; the groups file with other
  than two columns or an empty label header; a row whose field count
  differs from its header; an empty or repeated feature id, sample or
  coefficient name; a response cell that is not a number (other than
  blank, NA, NaN); a design cell that is not a finite number

- `mzlib_usage_error` (usage):

  the design's or the groups' samples are not exactly the responses'
  sample columns (both differences are named); no coefficient on stdin,
  one named twice, or one not in the design; `--estimator` other than
  reml or ml, or without a value; `--threads` 0 or below -1

- `mzlib_usage_error` (usage):

  mzLib refuses the input (its ArgumentException, reclassified): a
  design not of full column rank; a sample whose group label is empty,
  named by 0-based position and then by its name

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- A random intercept per group and nothing more, as nlme::lme(y ~ ...,
  random = ~ 1 \| group): no random slopes, nested or crossed groupings,
  or residual correlation.

- Fixed effects are tested by Wald t on nlme's containment degrees of
  freedom. Satterthwaite and Kenward-Roger (lmerTest, pbkrtest) are not
  implemented, so p-values differ from lmerTest's. mzLib matches nlme
  3.1.169 to 1e-5 relative (nlme's own optimizer precision) in
  estimates, standard errors, variances, p-values and log-likelihoods,
  with identical df.

- group_variance 0 is the boundary estimate. nlme works on the log
  scale, cannot reach 0, and reports a tiny positive value there
  instead.

- Under reml, log_likelihood is the restricted log-likelihood: compare
  it only between models with the same fixed effects, and use
  `--estimator` ml to compare designs. Under ml, residual_variance is
  RSS / observed, but standard errors use RSS / (observed -
  coefficients), as nlme's summary reports ML fits.

- p_value is not adjusted across features. Send one coefficient's
  p_value column to stats adjust for Benjamini-Hochberg.

- A 0 is an observation. Log-transform intensities first, and blank any
  cell that means 'not measured'.

- A feature the fixed effects fit exactly (a constant one, say) is
  reported fitted, with a residual_variance at rounding level (about
  1e-16) and a meaningless t and p_value: MixedModel's not_converged
  guard catches only a residual sum of squares of exactly 0 (d10017ce).
  Check residual_variance before trusting a very small p_value. Pinned
  by
  StatisticsTests.FitMixed_AnExactFit_IsReportedFitted_WithRoundingLevelVariance.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.fit_mixed`

- Rust (mzLibRust): `mzlib::stats::fit_mixed_with` with
  `FitMixedOptions`

- R (mzLibR): `stats_fit_mixed`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.2307/2529876](https://doi.org/10.2307/2529876): the
  random-effects linear model y = X beta + Z u + e this verb fits with
  one random intercept (Laird and Ware, Biometrics 38:963, 1982)

- [doi:10.1093/biomet/58.3.545](https://doi.org/10.1093/biomet/58.3.545):
  restricted maximum likelihood, the default estimator (Patterson and
  Thompson, Biometrika 58:545, 1971)

## See also

[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
when samples are independent,
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)

## Examples

``` r

# mzLib's own reference input (120 simulated features, 29 samples in 7 groups of unequal size;
# `wp` is constant within each group, `age` varies within them), and nlme's output for it:
data <- system.file("extdata", "stats", package = "mzLibR")
fit <- stats_fit_mixed(file.path(data, "lmm_reference_responses.tsv"),
  file.path(data, "lmm_reference_design.tsv"), file.path(data, "lmm_reference_groups.tsv"),
  c("Intercept", "age", "wp"))
fit
#> <mzlibr_mixed_model_fit> 120 features x 29 samples in 7 groups (group), reml, reporting Intercept, age, wp
#>   fitted 120, too_few_observations 0, rank_deficient 0, too_few_groups 0, variance_ratio_at_limit 0, not_converged 0
# `wp` is a group-level covariate, so it is tested on the groups, not the samples:
table(fit$records$df[fit$records$coefficient == "wp"])
#> 
#>   4   5 
#>   1 119 
```
