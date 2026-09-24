# Read SDRF age cells into years, refusing anything that would need a guess

Read characteristics\[age\] cells into years with an honest precision,
refusing any cell that would need a guess.

Calls mzLib's `SdrfAge.TryParse` on each cell. It reads the
specification's grammar (`58Y`, `30Y6M`, `16W`), ranges (`40Y-85Y`,
`6-8 weeks`), bounds (`>=90Y`, `<1Y`) and unambiguous words (`3 year`,
`4 hour`). It **refuses** a bare number - `63` is 11% of real age cells,
and 63 years and 63 days are both plausible in one study - as well as
reserved words and free text.

## Usage

``` r
sdrf_parse_ages(cells, timeout = 60)
```

## Arguments

- cells:

  A character vector of cells, e.g. the age characteristic's values from
  [`sdrf_value`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md).
  `NA` is sent as an empty cell and comes back refused as `"empty"`, so
  the result stays aligned with the input.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

One bridge call for any number of cells, so pass them all at once rather
than looping.

## Value

An `mzlibr_sdrf_ages`: `cell_count` cells given, `parsed_count` cells
read as an age, `caveats`, `column_names`, and `records`, one row per
cell **in input order**: `cell`, as given; `years`, the single figure to
place the sample on an age axis, in years - the age, a range's midpoint,
or a bound; `min_years`, the youngest the cell allows, in years (0 for
an upper bound); `max_years`, the oldest, in years; `precision`
(`"Exact"`, `"Range"`, `"LowerBound"` or `"UpperBound"`);
`follows_specification` (`TRUE` for the specification's own `nYnMnD` /
`nW` grammar); and `refusal` - `NA` when the cell was read, else
`"empty"`, `"reserved_word"`, `"no_unit"` or `"unreadable"`.

`years`, `min_years`, `max_years`, `precision` and
`follows_specification` are `NA` for a refused cell. `max_years` is also
`NA` for a lower bound such as `>=90Y`, which has no upper limit;
`precision` tells the two apart. A month is 1/12 year, a week 7/365.25,
a day 1/365.25.

## Wraps

Wire verb `sdrf parse-age`. Generated from the bridge's verb spec
`sdrf.parse-age.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfAge.TryParse`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfAge.cs)
  in `mzLib/Readers/Sdrf/SdrfAge.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `cells` (wire `--stdin`):

  string\[\]; required; range `>= 1 line`. One cell per line. BLANK
  LINES ARE KEPT (unlike every other stdin verb) so output row i is line
  i; one trailing newline does not add a cell. A UTF-8 BOM is dropped.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `cell_count`:

  int; in **cells**; never `NA`. Cells read from stdin; the length of
  every column.

- `parsed_count`:

  int; in **cells**; never `NA`. Cells read as an age.

- `column_names`:

  string\[\]; never `NA`. The table's columns.

- `records` (wire `columns`):

  table; never `NA`. One row per input line, in order.

- `caveats`:

  string\[\]; never `NA`. Four static caveats.

## Columns of `records`

- `cell`:

  string; never `NA`. The line exactly as given (mzLib's SdrfAge.Cell is
  the trimmed form).

- `years`:

  float; in **years**; `NA` when the cell was refused; refusal says why.
  The age, a range's MIDPOINT, or a bound.

- `min_years`:

  float; in **years**; `NA` when refused. Youngest the cell allows; 0
  for UpperBound.

- `max_years`:

  float; in **years**; `NA` when refused, OR precision is LowerBound:
  UNBOUNDED (mzLib's +infinity, which JSON cannot carry). Oldest the
  cell allows.

- `precision`:

  string; `NA` when refused. Exact (incl. degenerate 40Y-40Y since
  \#1333) \| Range \| LowerBound \| UpperBound.

- `follows_specification`:

  bool; `NA` when refused. true for the spec grammar, false for
  unambiguous words ("3 year", "6-8 weeks", "4 hour").

- `refusal`:

  string; `NA` when the cell was read. empty \| reserved_word (any case)
  \| no_unit (a bare number) \| unreadable.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  stdin empty

- `mzlib_bridge_error` (correctness):

  never expected: TryParse does not throw

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Units are years: month 1/12, week 7/365.25, day 1/365.25, hour
  1/(365.25 x 24).

- A bare number (63, 11% of real corpus age cells) is refused, never
  assumed to be years.

- years for a wide Range is a poor stand-in for any one sample; filter
  on precision.

- Every cell costs one line, not one process: send all cells in one
  call.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.parse_ages`

- Rust (mzLibRust): `mzlib::sdrf::parse_ages`

- R (mzLibR): `sdrf_parse_ages`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The refusal reason is a bridge-side classification of TryParse's
  false; see sdrf.samples.yaml. An upstream TryParse(cell, out age, out
  reason) would retire it.

- Binding name is parse_ages (plural) for the verb parse-age (one cell
  per line). Kept because the function takes a list; revisit if a
  binding wants the wire name.

- Rust and R spellings are intended, not yet ported. R: NA_character\_
  in the input should map to an empty line, as Python maps None.

- since.pymzlib is null until released.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics age grammar nYnMnD (Dai et al., Nat Commun 12, 5854,
  2021)

## See also

[`sdrf_samples`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md),
which parses the ages of a document's samples for you.

## Examples

``` r

ages <- sdrf_parse_ages(c("58Y", "30Y6M", "40Y-85Y", "40Y-40Y", ">=90Y", "<1Y", "6-8 weeks",
  "63", "not available", "", "about forty"))
ages
#> <mzlibr_sdrf_ages> 7 of 11 cells read as an age
#>   refused (empty): 1
#>   refused (no_unit): 1
#>   refused (reserved_word): 1
#>   refused (unreadable): 1
ages$records[, c("cell", "years", "min_years", "max_years", "precision", "refusal")]
#>             cell      years  min_years  max_years  precision       refusal
#> 1            58Y 58.0000000 58.0000000 58.0000000      Exact          <NA>
#> 2          30Y6M 30.5000000 30.5000000 30.5000000      Exact          <NA>
#> 3        40Y-85Y 62.5000000 40.0000000 85.0000000      Range          <NA>
#> 4        40Y-40Y 40.0000000 40.0000000 40.0000000      Exact          <NA>
#> 5          >=90Y 90.0000000 90.0000000         NA LowerBound          <NA>
#> 6            <1Y  1.0000000  0.0000000  1.0000000 UpperBound          <NA>
#> 7      6-8 weeks  0.1341547  0.1149897  0.1533196      Range          <NA>
#> 8             63         NA         NA         NA       <NA>       no_unit
#> 9  not available         NA         NA         NA       <NA> reserved_word
#> 10                       NA         NA         NA       <NA>         empty
#> 11   about forty         NA         NA         NA       <NA>    unreadable
```
