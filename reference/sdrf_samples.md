# Each sample's characteristics and factor values, merged over its rows, ages parsed

One row per sample, column and position: what each source name's rows
agree on, the columns they disagree about named and withheld, and
characteristics\[age\] read into years.

Calls mzLib's `SdrfSampleBlock.BySourceName`: the sample half of the
document - `source name`, every `characteristics[...]`, every
`factor value[...]` - per sample, merged over the sample's rows. A
sample is a `source name`, matched ignoring case and surrounding space,
so one sample measured in three fractions is one sample. Where its rows
**disagree** - one fraction says `normal` and another `COVID-19` - mzLib
names the column and withholds it rather than let the first row win.

## Usage

``` r
sdrf_samples(path, timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_sdrf_samples`: `path`, `sample_count` samples (distinct
source names), `row_count` data rows, `conflict_count` (sample x column
pairs withheld as conflicting), `problems` (rows mzLib could not place;
empty for a well-formed document), `caveats`, `column_names`, and
`records`, a long data.frame with one row per sample x column x
position, samples in document order: `source_name`; `sample_row_count`,
the rows carrying the sample; `column_name`, as the header spells it;
`column_kind` (`source_name`, `characteristic` or `factor_value`);
`position` (1 for the first occurrence of a repeated column; `NA` on a
conflicting row); `status` (`"agreed"` or `"conflicting"`); and `value`,
the cell verbatim, `NA` when conflicting.

On the rows for the age characteristic six more columns hold the parsed
age, with the meanings
[`sdrf_parse_ages`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md)
gives: `age_years`, `age_min_years` and `age_max_years` in years,
`age_precision`, `age_follows_specification` and `age_refusal`. On every
other row they are `NA`.

## Key samples within a document

A source name is unique only within one document - two studies may both
have a `"Sample 1"` - so across documents key a sample on the document
as well.

## Wraps

Wire verb `sdrf samples`. Generated from the bridge's verb spec
`sdrf.samples.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfSampleBlock.BySourceName`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfSampleBlock.cs)
  in `mzLib/Readers/Sdrf/SdrfSampleBlock.cs` at mzLib `23c2490e`

- [`SdrfAge.TryParse`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfAge.cs)
  in `mzLib/Readers/Sdrf/SdrfAge.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; default absent. One .sdrf.tsv file. Exactly one of path and
  paths-stdin.

- wire `--paths-stdin`:

  flag; default absent. Many paths on stdin, as for sdrf validate. *Not
  an argument here: sdrf_samples_many() is the many-documents form.*

- wire `--threads`:

  int; in **documents**; default `1`; range `>= 1, or -1`. Documents
  read at once; output identical at any value. *Not an argument here:
  only sdrf_samples_many() takes threads.*

- wire `--on-error`:

  string; default `fail`; range `fail | skip`. As for sdrf validate.
  *Not an argument here: only sdrf_samples_many() takes on_error.*

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The absolute path read.

- `sample_count`:

  int; in **samples**; never `NA`. Distinct source names.

- `row_count`:

  int; in **rows**; never `NA`. Data rows in the document.

- `conflict_count`:

  int; never `NA`. Sample x column pairs withheld as conflicting.

- `problems`:

  string\[\]; never `NA`. mzLib's BySourceName problems: no source name
  column, or a row with a blank source name. Empty for a well-formed
  document.

- `column_names`:

  string\[\]; never `NA`. The table's columns.

- `records` (wire `columns`):

  table; never `NA`. Samples in order of first appearance in the
  document; within a sample: source name, characteristics, factor values
  (header order), then conflicting columns.

- `caveats`:

  string\[\]; never `NA`. Four static caveats.

## Columns of `records`

- `source_name`:

  string; never `NA`. The sample as its first row spelled it, trimmed.
  Matched case-insensitively (mzLib's key).

- `sample_row_count`:

  int; in **rows**; never `NA`. Rows carrying this sample.

- `column_name`:

  string; never `NA`. The header's own spelling; characteristics\[age\]
  and characteristics\[Age\] are two columns.

- `column_kind`:

  string; never `NA`. source_name \| characteristic \| factor_value.

- `position`:

  int; `NA` when status is conflicting: the column was withheld, so it
  has no positions. 0-based occurrence of a repeated column.

- `status`:

  string; never `NA`. agreed \| conflicting (the sample's rows disagree;
  mzLib withheld the column).

- `value`:

  string; `NA` when status is conflicting (withheld, not absent). The
  cell verbatim, reserved words included.

- `age_years`:

  float; in **years**; `NA` when not a characteristics\[age\] row, or
  the cell was refused (age_refusal says why). The age, a range's
  midpoint, or a bound.

- `age_min_years`:

  float; in **years**; `NA` when as age_years. Youngest the cell allows;
  0 for UpperBound.

- `age_max_years`:

  float; in **years**; `NA` when not an age row; refused; OR precision
  is LowerBound, meaning UNBOUNDED (mzLib's +infinity, which JSON cannot
  carry). Oldest the cell allows.

- `age_precision`:

  string; `NA` when not an age row, or refused. Exact \| Range \|
  LowerBound \| UpperBound (SdrfAgePrecision).

- `age_follows_specification`:

  bool; `NA` when not an age row, or refused. true for nYnMnD / nW
  grammar, false for unambiguous words.

- `age_refusal`:

  string; `NA` when not an age row, or the cell was read. empty \|
  reserved_word \| no_unit \| unreadable. Classified by the bridge;
  mzLib's TryParse answers only true/false.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`sdrf_samples_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md).

- `mzlib_usage_error` (usage):

  as for sdrf validate

- `mzlib_bridge_error` (correctness):

  mzLib cannot read a document

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- A source name is unique only within a document; across a bulk result
  key on (source_index, source_name).

- Age units: a month is 1/12 year, a week 7/365.25, a day 1/365.25, an
  hour 1/(365.25 x 24).

- The age column is matched ignoring case, so characteristics\[Age\] is
  still read; sdrf validate reports the casing.

- A conflicting characteristics\[age\] has no parsed age: the column is
  withheld.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. For many files call
[`sdrf_samples_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.samples`; many files:
  `pymzlib.sdrf.samples_many`

- Rust (mzLibRust): `mzlib::sdrf::samples` with `BulkOptions`; many
  files: `mzlib::sdrf::samples_many`

- R (mzLibR): `sdrf_samples`; many files:
  [`sdrf_samples_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md)

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The refusal classification (empty, reserved_word, no_unit, unreadable)
  is bridge-side. It uses mzLib's own SdrfValidator.ReservedWords, but a
  TryParse overload returning the reason would move it upstream where it
  belongs; propose to mzLib.

- Sample order is first appearance in the document, recovered through
  BySourceName's own dictionary comparer because the dictionary's
  enumeration order is not contractual; an ordered return from mzLib
  would remove that step.

- since.pymzlib is null until released; Rust and R are not yet ported.

- Bulk spellings are the cross-binding decision of 2026-09-23: a
  separate \<verb\>\_many in every binding (R too, not vectorisation).
  Rust and R names here are intended, not yet ported.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics source name, characteristics and factor value columns,
  and the age grammar (Dai et al., Nat Commun 12, 5854, 2021)

## See also

[`sdrf_samples_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md),
[`sdrf_parse_ages`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md),
[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)

## Examples

``` r

samples <- sdrf_samples("sdrf_cohort.sdrf.tsv")
samples
#> <mzlibr_sdrf_samples> sdrf_cohort.sdrf.tsv
#>   6 samples over 12 rows; 1 column(s) withheld as conflicting
#>   4 caveat(s) - read x$caveats once
ages <- samples$records[samples$records$column_name == "characteristics[age]", ]
ages[, c("source_name", "value", "age_years", "age_precision", "age_refusal")]
#>    source_name         value age_years age_precision   age_refusal
#> 6           S1           58Y      58.0         Exact          <NA>
#> 15          S2       40Y-85Y      62.5         Range          <NA>
#> 24          S3         >=90Y      90.0    LowerBound          <NA>
#> 33          S4            63        NA          <NA>       no_unit
#> 42          S5 Not available        NA          <NA> reserved_word
#> 50          S6       40Y-40Y      40.0         Exact          <NA>
samples$records[samples$records$status == "conflicting", c("source_name", "column_name")]
#>    source_name              column_name
#> 54          S6 characteristics[disease]
```
