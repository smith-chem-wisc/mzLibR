# Fit one linear model per feature and test coefficients with limma's moderated t

Fit one linear model per feature of a feature-by-sample table and test
each named coefficient with limma's empirical-Bayes moderated t,
Benjamini-Hochberg adjusted, in one call.

Runs mzLib's `LinearModel.Fit` once - every design column is fitted -
and then `EmpiricalBayes.Moderate` for each coefficient you name, which
also applies Benjamini-Hochberg across that coefficient's fitted
features. This is limma's `lmFit()` then `eBayes(legacy = TRUE)` (Smyth
2004), which mzLib holds to limma to 1e-8 relative on limma's own
reference output.

## Usage

``` r
stats_fit(responses, design, coefficients, trend = FALSE, spline_basis = NULL,
  threads = 1, timeout = 600)
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

- coefficients:

  Design column names to test, each with its own moderated test and its
  own Benjamini-Hochberg family.

- trend:

  Let the prior variance follow average response (limma's
  `trend = TRUE`), for data whose low-abundance features are noisier.

- spline_basis:

  Basis functions of the trend spline, intercept included (the bridge's
  default is 4). Only with `trend = TRUE`; mzLib may use fewer, and
  `prior$spline_basis_count` reports what it used.

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
order you named them - with `feature`, `coefficient`, `status`,
`estimate` (response units per coefficient unit: a **log2 fold change**
for a 0/1 group column on log2 intensities), `standard_error` (response
units per coefficient unit), `t` (standard errors), `df_total` (degrees
of freedom), `p_value` and `bh_adjusted` (fractions, 0 to 1),
`posterior_variance` and `prior_variance` (squared response units),
`sigma` (response units), `df_residual`, `observed` (samples) and
`average_response` (response units). Statistics are `NA` where `status`
is not `"fitted"`. - `feature_count` features and `sample_count` samples
read; `sample_names`, `coefficient_names` (every design column) and
`tested_coefficients`; `missing_count` and `zero_count` response cells;
`status_counts`, features per status; `threads`; `row_count` rows. -
`residual_df_differ`: `TRUE` when the fitted features do not all have
the same residual degrees of freedom. Then default limma (3.61 and
later) uses a different prior estimator and reports different moderated
statistics; these are `legacy = TRUE`. - `prior`, the fitted variance
prior, shared by every coefficient: `df` (degrees of freedom, `NA` when
infinite), `df_infinite`, `trended`, `spline_basis_count` (basis
functions) and `scale` (squared response units, `NA` when trended). -
`caveats`, including any that apply to this input.

## What it is not

It is limma's \*legacy\* estimator. `robust = TRUE`, contrasts, the
B-statistic and observation weights are not implemented; write a
contrast as its own design column instead. `bh_adjusted` controls the
false discovery rate among the features tested; it is not a target-decoy
q-value.

## Wraps

Wire verb `stats fit`. Generated from the bridge's verb spec
`stats.fit.yaml` (bridge commit `c0cc92372aad`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`LinearModel.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/LinearModel.cs)
  in `mzLib/StatisticalModels/LinearModel.cs` at mzLib `0a808fec`

- [`EmpiricalBayes.Moderate`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `0a808fec`

- [`EmpiricalBayes.FitPrior`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/EmpiricalBayes.cs)
  in `mzLib/StatisticalModels/EmpiricalBayes.cs` at mzLib `0a808fec`

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
  an intercept column if you want one.

- `coefficients` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line; each a design column name, none repeated`. The
  coefficients to test, one design column name per line. Each gets its
  own moderated test and its own Benjamini-Hochberg family. The model
  always fits every design column.

- `trend`:

  flag; default `FALSE`. Let the prior variance trend with each
  feature's average response (limma's trend = TRUE), as a natural cubic
  spline.

- `spline_basis` (wire `--spline-basis`):

  int; in **basis functions**; default `4`; range
  `>= 1; only with --trend`. Basis functions of the trend spline,
  intercept included. mzLib caps it at the distinct average responses
  and the fitted features minus one; prior.spline_basis_count reports
  what was used.

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
  report different moderated statistics; these follow legacy = TRUE. A
  caveat says so.

- `prior`:

  object; never `NA`. The fitted variance prior, shared by every tested
  coefficient; fields under result.tables.prior.

- `row_count`:

  int; in **rows**; never `NA`. feature_count x the number of tested
  coefficients.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (tested coefficient, feature):
  coefficient blocks in stdin order, features in file order.

- `caveats`:

  string\[\]; never `NA`. Three static caveats, plus one each when
  residual_df_differ, when any feature was not fitted (with counts),
  when any cell is 0, and when the prior df is infinite.

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
  it.

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

- `prior$df` (wire `df`):

  float; in **degrees of freedom**; `NA` when infinite: df_infinite is
  true. Prior degrees of freedom d0.

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

  float; in **squared response units**; `NA` when trended: s0^2 is per
  feature, in the prior_variance column. s0^2, the prior variance.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  responses or design missing, blank or not found; a table empty,
  header-only or with too few columns; a row whose field count differs
  from its header; an empty or repeated feature id, sample or
  coefficient name; a response cell that is not a number (other than
  blank, NA, NaN); a design cell that is not a finite number

- `mzlib_usage_error` (usage):

  the design's samples are not exactly the responses' sample columns
  (both differences are named); no coefficient on stdin, one named
  twice, or one not in the design; `--spline-basis` without `--trend`,
  or below 1; `--threads` 0 or below -1

- `mzlib_usage_error` (usage):

  mzLib refuses the input: the design is not of full column rank, or
  fewer than 2 features could be fitted, so no variance prior exists
  (mzLib's ArgumentException, reclassified)

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- This is limma's LEGACY empirical-Bayes estimator (eBayes(legacy =
  TRUE)). Since limma 3.61, eBayes defaults to a different estimator
  whenever residual df differ between features, which omitting missing
  values causes; residual_df_differ says when that applies.

- Not implemented: robust = TRUE, contrasts, the B-statistic,
  observation weights. Write a contrast as its own design column
  instead.

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

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Rust and R spellings are proposals until those bindings ship the verb.

- stats_fit_malat.json is the largest recorded fixture (about 258 KB):
  492 features x 2 coefficients. It is real data (mzLib's RNA oligo
  table, \#1388), derived by a rule held in
  StatisticsTests.MalatDilutionFixture_IsLog2OfMzLibsOligoTable.

## References

- [doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027):
  the empirical-Bayes moderated t this verb computes (Smyth, Stat Appl
  Genet Mol Biol 3:3, 2004)

- [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007): limma,
  whose eBayes(legacy = TRUE) mzLib reproduces to 1e-8 relative (Ritchie
  et al., Nucleic Acids Res 43:e47, 2015)

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the Benjamini-Hochberg adjustment in bh_adjusted (J R Stat Soc B
  57:289, 1995)

## See also

[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md),
[`stats_meta`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)

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
#>   prior df 19.56
#>   ! residual df differ between features: these are eBayes(legacy = TRUE) statistics
fit$status_counts
#>               fitted too_few_observations       rank_deficient 
#>                  212                  214                   66 
malat <- fit$records[fit$records$coefficient == "malat_250ng" & fit$records$status == "fitted" &
  startsWith(fit$records$feature, "MALAT1:"), ]
c(oligos = nrow(malat), median_log2_fold_change = median(malat$estimate))
#>                  oligos median_log2_fold_change 
#>             116.0000000              -0.8914901 
```
