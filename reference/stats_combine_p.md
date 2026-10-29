# Combine several independent p-values per feature into one, by Fisher's or Stouffer's method

Combine several independent p-values for each feature into one, by
Fisher's or Stouffer's method, optionally weighted and optionally signed
by each test's effect.

Each feature is combined on its own by mzLib's `PValueCombination`:
Fisher (1925), X = -2 sum ln p against chi-square on 2k degrees of
freedom, or Stouffer (1949), Z = sum w z / sqrt(sum w^2) against the
standard normal. mzLib holds both to R's
[`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) and
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) to 1e-8 relative. Use
it when you have one p-value per study - per dataset, cohort or
laboratory - for the same feature, and want one answer from all of them.

## Usage

``` r
stats_combine_p(tests, method, timeout = 60)
```

## Arguments

- tests:

  A data.frame with one row per test: `feature` and `p_value`, and
  optionally `weight` and `effect` columns. A feature's tests need not
  be adjacent. `p_value` is a fraction from 0 to 1, or `NA` for a test
  that could not be run (it is omitted, and counted in the `omitted`
  column). `weight` (Stouffer only) is a positive number, for example
  the square root of the study's sample size; give it for every test or
  for none (`NA` is no weight). `effect` is the test's estimate in any
  unit; give it for every test or for none, and `NA` for a test with no
  direction (omitted). An `effect` column without a `weight` column
  sends no weights.

- method:

  `"fisher"` or `"stouffer"`. Required: mzLib has no preferred rule.
  Fisher is driven by the strongest single result and takes no weights;
  Stouffer weighs every study, and takes weights.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**Direction.** Both rules combine \*one-sided\* p-values. A two-sided
p-value from
[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
does not say which way the effect went, so two studies that disagree
look like two that agree. Give every test its `effect` (its estimate;
only the sign is used) and mzLib first turns each p into the one-sided
p-value for effect \> 0 (`PValueCombination.OneSided`: p/2 for a
positive effect, 1 - p/2 for a negative one), so opposite effects
cancel. The combined `p_value` then tests effect \> 0; for Stouffer, a
`p_value` near 1 means a consistently negative effect.

Neither rule estimates an effect size:
[`stats_meta`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
does that, from estimates and standard errors. Neither adjusts across
features: send `records$p_value` to
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md).

## Value

An `mzlibr_combined_p_values`:

\- `records`, a data.frame with one row per feature in order of first
appearance: `feature`; `studies`, the tests combined (tests with a
finite p-value, after the signed conversion when effects were given);
`omitted`, the feature's tests that were not combined (tests): no
p-value, or (signed) an effect of 0 or `NA`; `statistic`, Fisher's X or
Stouffer's Z (fisher: chi-square statistic on 2 x studies degrees of
freedom; stouffer: standard normal deviates), `NA` when it is not
finite; and `p_value`, a fraction (0 to 1), the upper tail of the
statistic, so ONE direction (effect \> 0 when signed), `NA` when there
is nothing to combine, or a p of 0 and a p of 1 together (Stouffer). An
`NA` statistic with a p-value is an infinite one: `p_value` 0 means plus
infinity (a p-value of 0 among the tests), 1 means minus infinity
(Stouffer, a p-value of 1). - `method`, as you asked; `weighted`
(Stouffer with weights); `signed` (effects were given, so every p-value
was read as two-sided and converted before combining); `line_count`
tests (lines) sent over every feature; `feature_count` features;
`row_count` rows; `column_names`. - `caveats`: what the combination
assumes, plus any that apply to this input (omitted tests, infinite
statistics).

## The trap

Two-sided p-values given without their effects are combined as if
one-sided, so two studies whose effects point opposite ways reinforce
each other instead of cancelling. Give `effect` for that. mzLib computes
no two-sided combined p-value.

## Wraps

Wire verb `stats combine-p`. Generated from the bridge's verb spec
`stats.combine-p.yaml` (bridge commit `aea14604d5f7`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PValueCombination.Fisher`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/PValueCombination.cs)
  in `mzLib/StatisticalModels/PValueCombination.cs` at mzLib `d10017ce`

- [`PValueCombination.Stouffer`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/PValueCombination.cs)
  in `mzLib/StatisticalModels/PValueCombination.cs` at mzLib `d10017ce`

- [`PValueCombination.OneSided`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/PValueCombination.cs)
  in `mzLib/StatisticalModels/PValueCombination.cs` at mzLib `d10017ce`

- [`CombinedPValue`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/PValueCombination.cs)
  in `mzLib/StatisticalModels/PValueCombination.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `method`:

  string; required; range `fisher | stouffer`. The combination rule.
  REQUIRED: mzLib has the two as separate methods and prefers neither,
  so the bridge picks none. fisher is PValueCombination.Fisher, Fisher
  (1925): X = -2 sum ln p, referred to chi-square with 2k degrees of
  freedom; it is driven by the smallest p-values and takes no weights.
  stouffer is PValueCombination.Stouffer, Stouffer et al. (1949): Z =
  sum w z / sqrt(sum w^2) with z = the standard normal quantile of 1 -
  p, referred to N(0, 1); it weighs every test.

- `tests` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line of 2 to 4 tab-separated fields`. One test per line,
  tab-separated: feature, p_value, then optionally weight and effect. No
  header; blank lines are skipped; a UTF-8 BOM is dropped. A feature's
  lines need not be adjacent: features are reported in order of first
  appearance, and each is combined by its own mzLib call. p_value is a
  fraction in \[0, 1\]; a blank, NA or NaN p_value is a test that could
  not be run, which mzLib omits. weight (stouffer only; a positive
  number, for example the square root of the study's sample size) must
  be on every line or on none; a blank or NA third field is no weight.
  effect (the test's estimate, any unit: only its sign is used) must be
  a fourth field on every line or on none. With effects, each p_value is
  read as TWO-SIDED and converted by PValueCombination.OneSided to the
  one-sided p-value for effect \> 0 (p/2 for a positive effect, 1 - p/2
  for a negative one) before combining; an effect that is 0, blank, NA
  or NaN has no direction, so mzLib returns NaN and the test is omitted.
  Write feature, p, a blank weight and the effect to give an effect
  without a weight.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `method`:

  string; never `NA`. fisher \| stouffer, as given.

- `weighted`:

  bool; never `NA`. Whether the lines carried weights (stouffer only).

- `signed`:

  bool; never `NA`. Whether the lines carried effects, so each p-value
  was read as two-sided and converted to the one-sided p-value for
  effect \> 0.

- `line_count`:

  int; in **lines**; never `NA`. Test lines read from stdin, over every
  feature.

- `feature_count`:

  int; in **features**; never `NA`. Distinct features; the table's row
  count.

- `row_count`:

  int; in **rows**; never `NA`. Equal to feature_count.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per feature, in order of first appearance.

- `caveats`:

  string\[\]; never `NA`. Three static caveats, plus: two when signed
  (one only for fisher), and one each when any test was omitted (with
  the count), when a feature has nothing left to combine, when a feature
  holds both a p-value of 0 and of 1 (stouffer), and when a statistic is
  infinite.

## Columns of `records`

- `feature`:

  string; never `NA`. As given on stdin.

- `studies`:

  int; in **tests**; never `NA`. Tests combined: mzLib's
  CombinedPValue.Studies, the finite p-values after any signed
  conversion.

- `omitted`:

  int; in **tests**; never `NA`. This feature's lines that were not
  combined: a blank, NA or NaN p_value, or (signed) an effect that is 0,
  blank, NA or NaN. studies + omitted is the feature's line count.

- `statistic`:

  float; in **fisher: chi-square statistic on 2 x studies degrees of
  freedom; stouffer: standard normal deviates**; `NA` when not finite.
  With studies 0, nothing was combined (p_value null too). Otherwise
  p_value says which: 0 means plus infinity (a p-value of 0 among the
  tests), 1 means minus infinity (stouffer, a p-value of 1), and null
  means undefined (stouffer, a p-value of 0 and one of 1 together).
  Fisher's X = -2 sum ln p, or Stouffer's Z = sum w z / sqrt(sum w^2).

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when studies is 0 (nothing to
  combine), or stouffer met a p-value of 0 and one of 1 (z of plus and
  minus infinity). Combined p-value: the UPPER tail of the statistic's
  reference distribution, so it tests one direction (effect \> 0 when
  signed). Fisher's is the regularized upper incomplete gamma Q(k, X/2),
  so a tiny value keeps its digits. Not adjusted across features.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--method` missing, given without a value, or not fisher or stouffer;
  stdin empty; a line with fewer than 2 or more than 4 tab-separated
  fields; an empty feature name; a p_value, weight or effect that is not
  a number (other than blank, NA, NaN)

- `mzlib_usage_error` (usage):

  a weight on some lines and not others; an effect field on some lines
  and not others; any weight with `--method` fisher

- `mzlib_usage_error` (usage):

  mzLib refuses a feature's input, named in the message: a p-value
  outside \[0, 1\], or a weight that is not finite and positive (mzLib's
  ArgumentException, reclassified)

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Both rules assume the tests are independent and test the same
  one-sided hypothesis. Fisher is driven by the smallest p-values;
  Stouffer weighs every test, and with weights (for example the square
  root of each study's sample size) lets larger studies count more.

- p_value is the upper tail of the statistic, so it tests one direction.
  Two-sided p-values given without their effects are combined as if
  one-sided, and opposite effects then reinforce each other instead of
  cancelling: give each test's effect for that.

- Neither rule estimates an effect size; pool the estimates with stats
  meta for that. p_value is not adjusted across features; send the
  p_value column to stats adjust.

- Signed: p_value tests effect \> 0, and for stouffer a p_value near 1
  means a consistently negative effect. mzLib computes no two-sided
  combined p-value. For a negative effect whose two-sided p-value is
  below about 2e-16, 1 - p/2 rounds to 1, so that test's z is minus
  infinity and decides a stouffer feature (p_value 1).

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.combine_p`

- Rust (mzLibRust): `mzlib::stats::combine_p`

- R (mzLibR): `stats_combine_p`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.1111/j.1420-9101.2005.00917.x](https://doi.org/10.1111/j.1420-9101.2005.00917.x):
  the weighted Z (Stouffer) method and why it outperforms Fisher's for a
  shared effect (Whitlock, J Evol Biol 18:1368, 2005)

## See also

[`stats_meta`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
to pool effect sizes,
[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
to adjust across features

## Examples

``` r

# mzLib's own reference cases (constructed p-values, not proteomics data), by both rules:
cases <- read.delim(system.file("extdata", "stats", "combination_inputs.tsv", package = "mzLibR"))
fisher <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p), method = "fisher")
stouffer <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p,
  weight = cases$weight), method = "stouffer")
data.frame(feature = fisher$records$feature, studies = fisher$records$studies,
  fisher = signif(fisher$records$p_value, 3),
  weighted_stouffer = signif(stouffer$records$p_value, 3))
#>   feature studies   fisher weighted_stouffer
#> 1     few       3 8.71e-02          0.150000
#> 2    tiny       4 3.14e-11          0.002280
#> 3     one       1 3.70e-02          0.037000
#> 4    many      12 1.13e-01          0.032300
#> 5   mixed       5 1.91e-02          0.713000
#> 6    tail       3 2.14e-18          0.000836
```
