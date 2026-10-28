# Spearman-correlate every feature with a per-sample trait, optionally given covariates

Spearman's rank correlation of every feature of a feature-by-sample
table with one per-sample trait, with a two-sided p-value, optionally
partial given per-sample covariates.

For each feature (row) of the responses table, mzLib's
`SpearmanCorrelation.Correlate` ranks the feature's values and the trait
over the samples where both are known, and tests rho = 0, two-sided.
With `given`, `SpearmanCorrelation.PartialCorrelate` first removes from
both sets of ranks what the covariates' ranks explain, as
`ppcor::pcor.test(method = "spearman")` does. mzLib holds both to R's
[`cor.test()`](https://rdrr.io/r/stats/cor.test.html) and to ppcor 1.1
to 1e-8 relative.

## Usage

``` r
stats_correlate(responses, trait, given = NULL, timeout = 600)
```

## Arguments

- responses:

  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  responses table: a TSV with a header row, the feature id first, then
  one column per sample, named by its header. Blank, `NA` or `NaN` is
  missing. **0 is an observation**, ranked like any other value.

- trait:

  A TSV with a header row and one row per sample: the sample name
  (exactly the responses' sample columns, in any order), then ONE value
  column, whose header names the trait - age, dose, loading. A blank or
  `NA` cell leaves that sample out.

- given:

  A TSV of the same shape with one or more covariate columns (batch, run
  order, loading) to partial out. Leave it `NULL` for a plain
  correlation.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

Spearman uses ranks only. It needs no log transform, one extreme value
cannot drive it, and it finds any relationship that only ever rises (or
only ever falls), not just a straight line.

## Value

An `mzlibr_correlations`:

\- `records`, a data.frame with one row per feature in the responses
file's order: `feature`; `rho`, Spearman's rho (unitless correlation, -1
to 1), ties at their average rank, the partial rho with `given`;
`p_value`, a fraction (0 to 1), two-sided against rho = 0 and **not
adjusted across features** - both `NA` exactly when `method` is
`"not_estimable"`; `n`, the complete samples used (samples: response,
trait and every covariate finite); `has_ties`, whether any variable was
tied among them (`FALSE` too when there were too few to rank); and
`method`, `"exact"` (n \<= 9, no ties, plain), `"asymptotic"` (t on n -
2 df, or n - 2 - k with k covariates kept) or `"not_estimable"`. -
`responses_file`, `trait_file` and `given_file` (`NA` without `given`),
as absolute paths; `trait_name`, the trait column's header;
`covariate_names` (empty without `given`); `partial`, `TRUE` with
`given`; `feature_count` features and `sample_count` samples read;
`sample_names`; `missing_count` response cells (cells) that were blank,
`NA` or `NaN`; `trait_missing_count`, samples with no trait value
(samples); `given_incomplete_count`, samples missing a covariate
(samples), `NA` without `given`; `method_counts`, features per method;
`row_count` rows; `column_names`. - `caveats`: how the p-values were
computed and how they differ from R's defaults, plus a count of
not_estimable features when there are any.

## Against R's cor.test()

For 10 \<= n \< 1290 without ties, R's
[`cor.test()`](https://rdrr.io/r/stats/cor.test.html) default is an
Edgeworth series (AS 89), which mzLib does not use, so p-values there
differ slightly from R's default output; they match
`cor.test(exact = FALSE)`. Send `records$p_value` to
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
for Benjamini-Hochberg.

## Wraps

Wire verb `stats correlate`. Generated from the bridge's verb spec
`stats.correlate.yaml` (bridge commit `60b1d25968a1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SpearmanCorrelation.Correlate`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/SpearmanCorrelation.cs)
  in `mzLib/StatisticalModels/SpearmanCorrelation.cs` at mzLib
  `d10017ce`

- [`SpearmanCorrelation.PartialCorrelate`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/SpearmanCorrelation.cs)
  in `mzLib/StatisticalModels/SpearmanCorrelation.cs` at mzLib
  `d10017ce`

- [`SpearmanCorrelation.ExactMaxN`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/SpearmanCorrelation.cs)
  in `mzLib/StatisticalModels/SpearmanCorrelation.cs` at mzLib
  `d10017ce`

- [`SpearmanResult`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/SpearmanCorrelation.cs)
  in `mzLib/StatisticalModels/SpearmanCorrelation.cs` at mzLib
  `d10017ce`

- [`SpearmanPValueMethod`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/SpearmanCorrelation.cs)
  in `mzLib/StatisticalModels/SpearmanCorrelation.cs` at mzLib
  `d10017ce`

## Parameters: units, ranges and defaults

- `responses`:

  path; required. stats fit's responses table. A TSV with a header row:
  the first column is the feature id (unique, any header) and every
  other column is one sample, named by its header. A blank cell, NA or
  NaN (any case) is missing. Only ranks are used, so a monotone
  transform such as log2 changes nothing; but a 0 is an observation,
  ranked like any other value, so blank any 0 that means not measured.
  Blank lines are skipped.

- `trait`:

  path; required. A TSV with a header row and one row per sample: the
  sample name, then EXACTLY ONE value column, whose header names the
  trait (for example age or dose). The sample names must be exactly the
  responses' sample columns, in any order: rows are matched by name. A
  blank, NA or NaN cell is a sample with no trait value, left out of
  every feature's correlation.

- `given`:

  path; default absent. Optional. A TSV with the trait file's shape but
  one or more covariate columns, each named by its header: per-sample
  nuisance variables (batch, run order, loading) to partial out. With
  it, every feature is correlated by PartialCorrelate instead of
  Correlate. The sample names must be exactly the responses' sample
  columns; a blank, NA or NaN cell leaves that sample out of every
  feature's correlation. Absent, the bridge sends none and the
  correlation is plain.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `responses_file`:

  string; never `NA`. Absolute path of the responses table.

- `trait_file`:

  string; never `NA`. Absolute path of the trait table.

- `given_file`:

  string; `NA` when no `--given`: the correlation is plain. Absolute
  path of the covariates table.

- `trait_name`:

  string; never `NA`. The trait file's value column header.

- `covariate_names`:

  string\[\]; never `NA`. The `--given` file's column headers, in order;
  empty without `--given`.

- `partial`:

  bool; never `NA`. True with `--given`: rows are partial correlations
  (PartialCorrelate).

- `feature_count`:

  int; in **features**; never `NA`. Feature rows read; the table's row
  count.

- `sample_count`:

  int; in **samples**; never `NA`. Sample columns read.

- `sample_names`:

  string\[\]; never `NA`. In responses header order, which is the order
  the trait and covariate rows were put in.

- `missing_count`:

  int; in **cells**; never `NA`. Response cells that were blank, NA or
  NaN.

- `trait_missing_count`:

  int; in **samples**; never `NA`. Samples with no finite trait value,
  so in no feature's correlation.

- `given_incomplete_count`:

  int; in **samples**; `NA` when no `--given`. Samples missing at least
  one covariate, so in no feature's correlation.

- `method_counts`:

  map\<string,int\>; in **features**; never `NA`. Features per method:
  exact, asymptotic, not_estimable.

- `row_count`:

  int; in **rows**; never `NA`. Equal to feature_count.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per feature, in responses file order.

- `caveats`:

  string\[\]; never `NA`. Four static caveats (five with `--given`),
  plus one when any feature is not_estimable, with the count.

## Columns of `records`

- `feature`:

  string; never `NA`. The feature id from the responses file's first
  column.

- `rho`:

  float; in **unitless correlation (-1 to 1)**; `NA` when not_estimable
  (see method). Spearman's rho: Pearson's correlation of the ranks, ties
  given their average rank. With `--given`, the partial rho: the
  correlation of the feature's and the trait's ranks after each is
  regressed on the covariates' ranks, all ranked over the complete
  samples.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not_estimable. Two-sided
  p-value against rho = 0, computed as method says. 0 when rho is
  exactly -1 or 1 with the asymptotic test. Not adjusted across
  features.

- `n`:

  int; in **samples**; never `NA`. Complete samples the correlation
  used: the response, the trait and (with `--given`) every covariate
  finite.

- `has_ties`:

  bool; never `NA`. Whether any variable (the response, the trait or,
  with `--given`, a covariate) has tied values among the complete
  samples. Also false when there were too few complete samples to rank
  at all (n below 3, or below 4 with `--given`), where mzLib does not
  look.

- `method`:

  string; never `NA`. exact (the full permutation distribution of S =
  sum d^2, for n \<= 9 with no ties and no `--given`, as R's cor.test
  exact = TRUE) \| asymptotic (t = rho sqrt(df / (1 - rho^2)) on df =
  n - 2, or n - 2 - k for k covariates kept, as cor.test exact = FALSE
  and ppcor::pcor.test) \| not_estimable (fewer than 3 complete samples,
  or k + 3 with `--given`; a constant variable; or covariates that carry
  the response or the trait entirely).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  responses or trait missing, blank or not found; `--given` given
  without a value or not found; a table empty, header-only or with too
  few columns; the trait file with other than exactly one value column;
  a row whose field count differs from its header; an empty or repeated
  feature id, sample or column name; a cell that is not a number (other
  than blank, NA, NaN)

- `mzlib_usage_error` (usage):

  the trait's or the covariates' samples are not exactly the responses'
  sample columns (both differences are named)

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Plain: p_value is exact, from the full permutation distribution, for n
  \<= 9 without ties, as R's cor.test(method = 'spearman'); otherwise it
  is asymptotic, t with n - 2 degrees of freedom, as cor.test(exact =
  FALSE).

- Plain: for 10 \<= n \< 1290 without ties, R's cor.test default is an
  Edgeworth series (AS 89), which mzLib does not use, so p-values there
  differ slightly from R's default output.

- With `--given`: every variable is ranked over the samples complete in
  all of them, rho is the correlation of the feature's and the trait's
  ranks after each is regressed on the covariates' ranks, and p_value is
  from t with n - 2 - k degrees of freedom, as ppcor::pcor.test(method =
  'spearman'). mzLib drops a covariate that the others carry entirely,
  and k counts the covariates kept; mzLib does not report k, so neither
  does this verb.

- With `--given`: give only nuisance covariates. A covariate that
  follows the trait takes the trait's own signal with it.

- A sample enters a feature's correlation only when its response, its
  trait and every covariate are finite: missing values are omitted,
  never imputed. n is the number of samples used.

- p_value is not adjusted across features. Send the p_value column to
  stats adjust for Benjamini-Hochberg.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.correlate`

- Rust (mzLibRust): `mzlib::stats::correlate_with` with
  `CorrelateOptions`

- R (mzLibR): `stats_correlate`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## References

- [doi:10.2307/1412159](https://doi.org/10.2307/1412159): the rank
  correlation coefficient rho (Spearman, Am J Psychol 15:72, 1904)

- [doi:10.5351/CSAM.2015.22.6.665](https://doi.org/10.5351/CSAM.2015.22.6.665):
  ppcor, whose pcor.test(method = 'spearman') the partial form
  (`--given`) reproduces (Kim, Commun Stat Appl Methods 22:665, 2015)

- [doi:10.2307/2347111](https://doi.org/10.2307/2347111): the Edgeworth
  series (AS 89) R's cor.test uses by default for 10 \<= n \< 1290,
  which mzLib does NOT use (Best and Roberts, Appl Stat 24:377, 1975)

## See also

[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
for a linear model on the same responses,
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)

## Examples

``` r

# The MALAT1 dilution series stats_fit() analyses, correlated with each run's MALAT1 load:
data <- system.file("extdata", "stats", package = "mzLibR")
result <- stats_correlate(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
  file.path(data, "malat_dilution_load_ng.tsv"))
result
#> <mzlibr_correlations> 492 features x 27 samples against malat_ng
#>   exact 4, asymptotic 308, not_estimable 180
rho <- result$records$rho
c(MALAT1 = median(rho[startsWith(result$records$feature, "MALAT1:")], na.rm = TRUE),
  FLuc = median(rho[startsWith(result$records$feature, "FLuc:")], na.rm = TRUE))
#>    MALAT1      FLuc 
#> 0.4330127 0.0000000 
```
