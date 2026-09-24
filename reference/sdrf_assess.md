# Decide whether an SDRF file describes its samples, or is a valid skeleton

Decide whether each SDRF file's sample half describes an experimental
design - Informative, Partial or Skeleton - with the per-column counts
behind the verdict.

A file generated from a list of data files - reserved words in every
sample column, one replicate number, no factor - passes
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md),
because reserved words are the specification's correct way to say
nothing. For grouping results by biology it is the same as having no
SDRF at all. mzLib's `SdrfSampleInformativeness.Assess` is the gate for
that.

## Usage

``` r
sdrf_assess(path, timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_sdrf_assessment`: `path`; `verdict` - `"Informative"`,
`"Partial"` or `"Skeleton"`; the three checks `factor_value_varies`,
`sample_is_described` and `biological_replicate_varies`; `row_count`
data rows; `caveats`; `column_names`; and `records`, the evidence - one
row per column a check read: `role`, `column_name`, `rows` (data rows),
`filled` (rows with a real answer), `absent` (rows empty or a reserved
word), `distinct_values` (values, compared ignoring case and surrounding
space) and `fill_rate`, the fraction (0-1) `filled / rows`.

## The verdict is how many checks pass

All three is `Informative`, none is `Skeleton`, anything between is
`Partial`. `factor_value_varies`: some `factor value[...]` column holds
two or more real answers. `sample_is_described`: some
`characteristics[...]` column other than organism and biological
replicate holds a real answer - organism is left out because a search
can fill it without a human. `biological_replicate_varies`:
`characteristics[biological replicate]` holds two or more real answers;
`FALSE` also when the column is missing. `Partial` is often legitimate -
read the `caveats`.

## Wraps

Wire verb `sdrf assess`. Generated from the bridge's verb spec
`sdrf.assess.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfSampleInformativeness.Assess`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfSampleInformativeness.cs)
  in `mzLib/Readers/Sdrf/SdrfSampleInformativeness.cs` at mzLib
  `23c2490e`

- [`SdrfColumnCoverage`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfCoverage.cs)
  in `mzLib/Readers/Sdrf/SdrfCoverage.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; default absent. One .sdrf.tsv file. Exactly one of path and
  paths-stdin.

- wire `--paths-stdin`:

  flag; default absent. Many paths on stdin, as for sdrf validate. *Not
  an argument here: sdrf_assess_many() is the many-documents form.*

- wire `--threads`:

  int; in **documents**; default `1`; range `>= 1, or -1`. Documents
  assessed at once; output identical at any value. *Not an argument
  here: only sdrf_assess_many() takes threads.*

- wire `--on-error`:

  string; default `fail`; range `fail | skip`. As for sdrf validate.
  *Not an argument here: only sdrf_assess_many() takes on_error.*

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The absolute path assessed.

- `verdict`:

  string; never `NA`. Informative (all three checks pass), Skeleton
  (none), Partial (between). mzLib's SdrfSampleVerdict member name.

- `factor_value_varies`:

  bool; never `NA`. Some factor value\[...\] column holds \>= 2 distinct
  real answers.

- `sample_is_described`:

  bool; never `NA`. Some characteristics\[...\] column other than
  organism and biological replicate holds a real answer.

- `biological_replicate_varies`:

  bool; never `NA`. characteristics\[biological replicate\] holds \>= 2
  distinct real answers; false when the column is missing.

- `row_count`:

  int; in **rows**; never `NA`. Data rows in the document.

- `column_names`:

  string\[\]; never `NA`. The evidence table's columns.

- `records` (wire `columns`):

  table; never `NA`. One row per column a check read: factor values,
  then sample characteristics, then biological replicate, each in
  mzLib's (ordinal) order.

- `caveats`:

  string\[\]; never `NA`. Four static caveats.

## Columns of `records`

- `role`:

  string; never `NA`. factor_value \| sample_characteristic \|
  biological_replicate: which check read this column.

- `column_name`:

  string; never `NA`. The column as the header spells it.

- `rows`:

  int; in **rows**; never `NA`. Rows, counted once per row per column
  NAME (a repeated column does not inflate it).

- `filled`:

  int; in **rows**; never `NA`. Rows with a real answer in at least one
  of the column's cells.

- `absent`:

  int; in **rows**; never `NA`. Rows that are empty or a reserved word
  (any case): rows - filled.

- `distinct_values`:

  int; in **values**; never `NA`. Different real answers, compared
  case-insensitively after trimming.

- `fill_rate`:

  float; in **fraction (0-1)**; never `NA`. filled / rows; 0 when rows
  is 0. May cross as an integer 0 or 1.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md).

- `mzlib_usage_error` (usage):

  as for sdrf validate (path / paths-stdin / threads / on-error rules,
  missing or repeated files)

- `mzlib_bridge_error` (correctness):

  mzLib cannot read a document

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Partial is often legitimate (a single-condition study has nothing to
  vary); the booleans say which check failed and the caller decides.

- On the 1,236-file bigbio curated corpus at 4f823dcd mzLib measured 300
  Informative, 742 Partial, 194 Skeleton.

- Organism is excluded from sample_is_described because a search can
  fill it without a human describing the samples.

- Looks only at the sample half; says nothing about structure (sdrf
  validate) or cross-file agreement (sdrf lint).

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. For many files call
[`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.assess`; many files:
  `pymzlib.sdrf.assess_many`

- Rust (mzLibRust): `mzlib::sdrf::assess` with `BulkOptions`; many
  files: `mzlib::sdrf::assess_many`

- R (mzLibR): `sdrf_assess`; many files:
  [`sdrf_assess_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- fill_rate is a C# double that serialises 1.0 as 1; bindings must read
  it as a float.

- since.pymzlib is null until released; Rust and R are not yet ported.

- Bulk spellings are the cross-binding decision of 2026-09-23: a
  separate \<verb\>\_many in every binding (R too, not vectorisation).
  Rust and R names here are intended, not yet ported.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics sample columns (characteristics, factor value) (Dai et
  al., Nat Commun 12, 5854, 2021)

## See also

[`sdrf_assess_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md),
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md),
[`sdrf_samples`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)

## Examples

``` r

assessment <- sdrf_assess("sdrf_cohort.sdrf.tsv")
assessment
#> <mzlibr_sdrf_assessment> sdrf_cohort.sdrf.tsv: Informative
#>   factor value varies: yes; sample described: yes; biological replicate varies: yes
#>   4 caveat(s) - read x$caveats once
assessment$records[, c("role", "column_name", "distinct_values", "fill_rate")]
#>                    role                           column_name distinct_values
#> 1          factor_value                 factor value[disease]               2
#> 2 sample_characteristic                  characteristics[age]               5
#> 3 sample_characteristic            characteristics[cell type]               0
#> 4 sample_characteristic              characteristics[disease]               2
#> 5 sample_characteristic        characteristics[organism part]               1
#> 6 sample_characteristic                  characteristics[sex]               2
#> 7  biological_replicate characteristics[biological replicate]               6
#>   fill_rate
#> 1 1.0000000
#> 2 0.8333333
#> 3 0.0000000
#> 4 1.0000000
#> 5 1.0000000
#> 6 1.0000000
#> 7 1.0000000
```
