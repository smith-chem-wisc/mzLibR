# Benjamini-Hochberg adjust p-values, keeping every position

Adjust a list of p-values for the false discovery rate with
Benjamini-Hochberg, keeping every input position and leaving untested
entries out of the family.

For the i-th smallest of m p-values the adjusted value is the minimum
over j \>= i of p\_(j) m / j, capped at 1 (Benjamini and Hochberg 1995),
computed by mzLib's `MultipleTesting.BenjaminiHochberg`. `NA` (or `NaN`)
marks a feature that was not tested: it stays `NA` and is **not counted
in m**, which is the honest family. Do not write 1 for an untested
feature; that enlarges m.

## Usage

``` r
stats_adjust(p_values, timeout = 60)
```

## Arguments

- p_values:

  A numeric vector of p-values, each a fraction from 0 to 1, or `NA` for
  untested.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
already reports `bh_adjusted`; this is for p-values from anywhere else.

## Value

An `mzlibr_adjusted`: `records`, a data.frame with row i for
`p_values[i]` - `p_value` as read and `bh_adjusted`, both fractions (0
to 1), `NA` for an untested entry; `line_count` lines (values) given;
`tested_count`, m, the p-values the adjustment is over; `row_count`
rows; and `caveats`.

## Wraps

Wire verb `stats adjust`. Generated from the bridge's verb spec
`stats.adjust.yaml` (bridge commit `e76157831b15`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`MultipleTesting.BenjaminiHochberg`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/StatisticalModels/MultipleTesting.cs)
  in `mzLib/StatisticalModels/MultipleTesting.cs` at mzLib `0a808fec`

## Parameters: units, ranges and defaults

- `p_values` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line; each in [0, 1], blank, NA or NaN`. One p-value per line.
  BLANK LINES ARE KEPT (as in sdrf parse-age) so output row i is line i;
  one trailing newline does not add a line. A blank, NA or NaN line is a
  feature that was not tested: it stays null and is not counted in m. A
  UTF-8 BOM is dropped.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `line_count`:

  int; in **lines**; never `NA`. Lines read from stdin; the length of
  every column.

- `tested_count`:

  int; in **p-values**; never `NA`. m: the finite p-values, the family
  the adjustment is over.

- `row_count`:

  int; in **rows**; never `NA`. Equal to line_count.

- `column_names`:

  string\[\]; never `NA`. The table's columns.

- `records` (wire `columns`):

  table; never `NA`. One row per input line, in order.

- `caveats`:

  string\[\]; never `NA`. Two static caveats.

## Columns of `records`

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when the line was blank, NA or
  NaN: not tested. The p-value as read.

- `bh_adjusted`:

  float; in **fraction (0 to 1)**; `NA` when not tested. min over j \>=
  i of p\_(j) m / j, capped at 1, for the i-th smallest p-value.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  stdin empty; a line that is not a number (other than blank, NA, NaN);
  a p-value outside \[0, 1\] (mzLib's refusal, reclassified)

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Benjamini-Hochberg controls the false discovery rate among the
  features tested. It is not a target-decoy q-value for identifications.

- m is the number of finite p-values. Leave a line blank, or write NA,
  for a feature that was not tested rather than writing 1, which would
  enlarge the family.

- stats fit already reports bh_adjusted per coefficient; use this verb
  for p-values from anywhere else.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.adjust`

- Rust (mzLibRust): `mzlib::stats::adjust`

- R (mzLibR): `stats_adjust`

## Since

Wire protocol 1; pyMzLib 0.3.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Rust and R spellings are proposals until those bindings ship the verb.

## References

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the step-up adjustment (Benjamini and Hochberg, J R Stat Soc B 57:289,
  1995)

## See also

[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)

## Examples

``` r

adjusted <- stats_adjust(c(0.0002, 0.004, 0.019, NA, 0.031, 0.2, NaN, 0.74))
adjusted$tested_count
#> [1] 6
adjusted$records
#>   p_value bh_adjusted
#> 1  0.0002      0.0012
#> 2  0.0040      0.0120
#> 3  0.0190      0.0380
#> 4      NA          NA
#> 5  0.0310      0.0465
#> 6  0.2000      0.2400
#> 7      NA          NA
#> 8  0.7400      0.7400
```
