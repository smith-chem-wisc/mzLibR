# Fit one logistic regression per feature of a 0/1 table, such as a detection table

Fit one logistic regression per feature of a 0/1 feature-by-sample
table, such as a detection table, and test each named coefficient with a
Wald z.

mzLib's `LogisticRegression.Fit` models logit P(outcome = 1) = X beta
for every feature (row) separately, by maximum likelihood, and tests
each coefficient with a two-sided Wald z. It is R's
`glm(y ~ 0 + X, family = binomial)` per feature, which mzLib reproduces
to 1e-8 relative in the coefficients (1e-6 in standard errors and
p-values).

## Usage

``` r
stats_fit_logistic(responses, design, coefficients, max_iterations = NULL,
  tolerance = NULL, threads = 1, timeout = 600)
```

## Arguments

- responses:

  A TSV with a header row: the feature id first, then one column per
  sample, named by its header -
  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  layout - with every cell 0, 1, or blank / `NA` / `NaN` for missing
  (left out of that feature's fit). Any other value is refused.

- design:

  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  design: a TSV with a header row and one row per sample, the sample
  name (matching the responses header, in any order) then one numeric
  column per coefficient. Include an intercept column yourself if you
  want one.

- coefficients:

  Design column names to report, each as one block of rows. Every design
  column is fitted whatever you name.

- max_iterations:

  The IRLS iteration limit per feature (iterations). Leave it `NULL` for
  mzLib's default (100 at mzLib `d10017ce`); a feature still changing at
  the limit is `"not_converged"`.

- tolerance:

  Stop when the relative change in deviance is below this (R
  [`glm()`](https://rdrr.io/r/stats/glm.html)'s rule). Leave it `NULL`
  for mzLib's default (1e-12 at mzLib `d10017ce`), which is tighter than
  [`glm()`](https://rdrr.io/r/stats/glm.html)'s 1e-8.

- threads:

  Features fitted at once (threads), or `-1` for every core. The answer
  is identical at any value.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

The usual question in proteomics: **is a feature's chance of being
measured related to the condition?** Write 1 where a protein was
quantified and 0 where it was not; a coefficient far from 0 says its
missing values are not missing at random.

A feature is reported - never dropped - when it cannot be fitted.
`"separated"` is the one to read carefully: the design splits its 0s
from its 1s (or every outcome is the same), so no finite estimate
exists. R's [`glm()`](https://rdrr.io/r/stats/glm.html) prints a huge
coefficient and a warning there; mzLib reports the status and `NA`. Look
at `events` and `observed` for those rows.

## Value

An `mzlibr_logistic_model_fit`:

\- `records`, a data.frame with one row per (reported coefficient,
feature) - one block of `feature_count` rows per coefficient you named,
in that order, each in the responses file's feature order: `feature`;
`coefficient`; `status`; `estimate`, the maximum-likelihood coefficient
in log-odds per coefficient unit, so `exp(estimate)` is the **odds
ratio**; `standard_error` (log-odds per coefficient unit), from the
inverse Fisher information at the estimate; `z` (standard errors), the
Wald statistic; `p_value`, a fraction (0 to 1), two-sided and **not
adjusted across features**; `deviance` (deviance, -2 log L), the
feature's residual deviance; `iterations`, IRLS iterations used (0 when
no fit was attempted); `observed`, samples with a 0 or 1; and `events`,
observed samples whose outcome is 1. Statistics are `NA` where `status`
is not `"fitted"`. - `responses_file` and `design_file`, as absolute
paths; `feature_count` features and `sample_count` samples read;
`sample_names`; `coefficient_names` (every design column) and
`tested_coefficients`; `max_iterations` (iterations) and `tolerance`
(relative change in deviance) as you passed them, `NA` when left to
mzLib; `threads`; `missing_count` response cells (cells);
`status_counts`, features per status; `row_count` rows;
`column_names`. - `caveats`: what the numbers are and are not, with a
count of unfitted features by status.

## What it is not

Wald tests are the least reliable of the likelihood tests near
separation and in small samples; no likelihood-ratio test is offered,
and Firth's penalized likelihood (logistf), which gives finite estimates
for separated features, is not implemented. Send one coefficient's
`p_value` column to
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
for Benjamini-Hochberg.

## Wraps

Wire verb `stats fit-logistic`. Generated from the bridge's verb spec
`stats.fit-logistic.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`LogisticRegression.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/LogisticRegression.cs)
  in `mzLib/StatisticalModels/LogisticRegression.cs` at mzLib `d10017ce`

- [`LogisticModelFit`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/LogisticRegression.cs)
  in `mzLib/StatisticalModels/LogisticRegression.cs` at mzLib `d10017ce`

- [`LogisticRegression.SeparationThreshold`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/LogisticRegression.cs)
  in `mzLib/StatisticalModels/LogisticRegression.cs` at mzLib `d10017ce`

- [`FeatureFitStatus`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/LinearModel.cs)
  in `mzLib/StatisticalModels/LinearModel.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `responses`:

  path; required. stats fit's responses table, with outcomes for cells.
  A TSV with a header row: the first column is the feature id (unique,
  any header) and every other column is one sample, named by its header.
  Every cell is 0, 1, or blank, NA or NaN (any case) for missing, which
  is omitted from that feature's fit. mzLib refuses any other value (a
  proportion is not a binary outcome); the refusal names the feature and
  sample. A 0 is an outcome, not a missing value. Blank lines are
  skipped.

- `design`:

  path; required. stats fit's design table. A TSV with a header row and
  one row per sample: the sample name, matching the responses header
  exactly (rows are matched by name, in any order), then one column per
  coefficient, every cell a finite number. Used as written: include an
  intercept column if you want one. mzLib refuses a design that is not
  of full column rank over all samples.

- `coefficients` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line; each a design column name, none repeated`. The
  coefficients to report, one design column name per line. The model
  always fits every design column; the table holds one block of rows per
  coefficient named here, in this order.

- `max_iterations` (wire `--max-iterations`):

  int; in **iterations**; default absent; range
  `>= 1 (mzLib's rule; anything else is mzLib's refusal, a usage error)`.
  The iteratively reweighted least squares iteration limit per feature
  (LogisticRegression.Fit maxIterations). A feature still changing at
  the limit is not_converged. Absent, the bridge sends none and mzLib's
  default applies (100 at d10017ce), and the envelope's max_iterations
  is null.

- `tolerance`:

  float; in **relative change in deviance (unitless)**; default absent;
  range
  `> 0 (mzLib's rule; anything else is mzLib's refusal, a usage error)`.
  Convergence rule (LogisticRegression.Fit tolerance): iteration stops
  when \|dev - dev_prev\| / (\|dev\| + 0.1) is below it, R glm's rule.
  Absent, the bridge sends none and mzLib's default applies (1e-12 at
  d10017ce, tighter than glm's default epsilon of 1e-8), and the
  envelope's tolerance is null.

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

- `feature_count`:

  int; in **features**; never `NA`. Feature rows read.

- `sample_count`:

  int; in **samples**; never `NA`. Sample columns read.

- `sample_names`:

  string\[\]; never `NA`. In responses header order, which is the order
  the design rows were put in.

- `coefficient_names`:

  string\[\]; never `NA`. Every design column, in order.

- `tested_coefficients`:

  string\[\]; never `NA`. The stdin coefficients, in order; the table
  holds one block of feature_count rows per entry.

- `max_iterations`:

  int; in **iterations**; `NA` when not given: mzLib's default was used
  (100 at d10017ce). The iteration limit sent to mzLib.

- `tolerance`:

  float; in **relative change in deviance (unitless)**; `NA` when not
  given: mzLib's default was used (1e-12 at d10017ce). The convergence
  tolerance sent to mzLib.

- `threads`:

  int; in **threads**; never `NA`. The thread count used.

- `missing_count`:

  int; in **cells**; never `NA`. Response cells that were blank, NA or
  NaN.

- `status_counts`:

  map\<string,int\>; in **features**; never `NA`. Features per status,
  always these five keys: fitted, too_few_observations, rank_deficient,
  separated, not_converged. They are the outcomes
  LogisticRegression.Fit's body can return at d10017ce; a stop-gap until
  mzLib says per model (mzLib#1487).

- `row_count`:

  int; in **rows**; never `NA`. feature_count x the number of tested
  coefficients.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (tested coefficient, feature):
  coefficient blocks in stdin order, features in file order.

- `caveats`:

  string\[\]; never `NA`. Four static caveats, plus one with the counts
  when any feature was not fitted.

## Columns of `records`

- `feature`:

  string; never `NA`. The feature id from the responses file's first
  column.

- `coefficient`:

  string; never `NA`. The coefficient this row reports.

- `status`:

  string; never `NA`. fitted \| too_few_observations (no more observed
  samples than coefficients) \| rank_deficient (the design restricted to
  the observed samples is not of full rank) \| separated (the observed
  outcomes are all 0, all 1, or split by the design, so no finite
  maximum-likelihood estimate exists; mzLib calls a fit separated when a
  fitted probability comes within
  LogisticRegression.SeparationThreshold, 1e-10 at d10017ce, of 0 or 1)
  \| not_converged (still changing at max-iterations). Only fitted rows
  carry statistics.

- `estimate`:

  float; in **log-odds per coefficient unit (the log odds ratio for a
  0/1 column)**; `NA` when status is not fitted. The maximum-likelihood
  coefficient on the logit scale. exp(estimate) is the odds ratio.

- `standard_error`:

  float; in **log-odds per coefficient unit**; `NA` when not fitted.
  From the inverse Fisher information at the estimate.

- `z`:

  float; in **standard errors**; `NA` when not fitted. Wald statistic,
  estimate / standard_error.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not fitted. Two-sided Wald
  p-value against a coefficient of 0, from the standard normal. Not
  adjusted across features.

- `deviance`:

  float; in **deviance (-2 x log-likelihood, natural log)**; `NA` when
  not fitted. The feature's residual deviance; the same in every
  coefficient block.

- `iterations`:

  int; in **iterations**; never `NA`. IRLS iterations used. 0 when no
  fit was attempted (too_few_observations, rank_deficient, or every
  outcome the same); max-iterations when not_converged; a separated
  feature that IRLS ran on reports the iterations it ran.

- `observed`:

  int; in **samples**; never `NA`. Samples with a 0 or 1 for this
  feature.

- `events`:

  int; in **samples**; never `NA`. Observed samples whose outcome is 1.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  responses or design missing, blank or not found; a table empty,
  header-only or with too few columns; a row whose field count differs
  from its header; an empty or repeated feature id, sample or
  coefficient name; a response cell that is not a number (other than
  blank, NA, NaN); a design cell that is not a finite number; the
  design's samples are not exactly the responses' sample columns (both
  differences are named)

- `mzlib_usage_error` (usage):

  no coefficient on stdin, one named twice, or one not in the design;
  `--max-iterations` or `--tolerance` without a value or not a number;
  `--threads` 0 or below -1

- `mzlib_usage_error` (usage):

  mzLib refuses the input (its ArgumentException, reclassified): a
  response that is not 0, 1 or missing, named by 0-based feature and
  sample and then by the names in the files; a design not of full column
  rank; `--max-iterations` below 1 or `--tolerance` not above 0

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Maximum likelihood by iteratively reweighted least squares, started
  and stopped as R's glm does; standard errors from the inverse Fisher
  information at the estimate; Wald z tests. mzLib matches R 4.6.1's
  glm(family = binomial, epsilon = 1e-14) to 1e-8 relative in
  coefficients and deviance and to 1e-6 relative in standard errors and
  p-values (summary.glm takes its weights from the step before the
  estimate). Wald tests are the least reliable of the likelihood tests
  near separation and in small samples; no likelihood-ratio test is
  offered.

- separated features carry no statistics. R's glm reports large,
  meaningless coefficients and a warning for them instead. A separated
  feature can be the strongest effect in the table (detected in every
  run of one condition and in none of another): read its events and
  observed counts. Firth's penalized likelihood (logistf), which gives
  finite estimates there, is not implemented.

- p_value is not adjusted across features. Send one coefficient's
  p_value column to stats adjust for Benjamini-Hochberg.

- A 0 is an outcome (for a detection table: measured as absent), not a
  missing value. Missing cells are omitted per feature, never imputed.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.fit_logistic`

- Rust (mzLibRust): `mzlib::stats::fit_logistic_with` with
  `FitLogisticOptions`

- R (mzLibR): `stats_fit_logistic`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.2307/2344614](https://doi.org/10.2307/2344614): generalized
  linear models and their fitting by iteratively reweighted least
  squares, the estimator mzLib implements (Nelder and Wedderburn, J R
  Stat Soc A 135:370, 1972)

- [doi:10.1093/biomet/71.1.1](https://doi.org/10.1093/biomet/71.1.1):
  separation: when the outcomes are split by the design the
  maximum-likelihood estimate does not exist, which is the separated
  status (Albert and Anderson, Biometrika 71:1, 1984)

## See also

[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
for the intensities themselves,
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)

## Examples

``` r

# mzLib's own reference input (80 simulated binary features, 40 samples), and R's glm() output
# for it, which this package ships:
data <- system.file("extdata", "stats", package = "mzLibR")
fit <- stats_fit_logistic(file.path(data, "glm_reference_responses.tsv"),
  file.path(data, "glm_reference_design.tsv"), c("Intercept", "x", "grp"))
fit
#> <mzlibr_logistic_model_fit> 80 features x 40 samples, reporting Intercept, x, grp
#>   fitted 76, too_few_observations 0, rank_deficient 0, separated 4, not_converged 0
x <- fit$records[fit$records$coefficient == "x", ]
glm <- read.delim(file.path(data, "glm_reference_fit.tsv"))
fitted <- glm$status == "fitted"
max(abs(x$estimate[fitted] - glm$b1[fitted]) / abs(glm$b1[fitted]))
#> [1] 4.203008e-12
x[1, c("estimate", "p_value", "events", "observed")]
#>     estimate   p_value events observed
#> 81 0.4592209 0.2723941     25       35
```
