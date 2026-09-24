# The samples of many SDRF files, in one bridge call

One row per sample, column and position: what each source name's rows
agree on, the columns they disagree about named and withheld, and
characteristics\[age\] read into years.

The samples of many SDRF files, in one bridge call

## Usage

``` r
sdrf_samples_many(paths, threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- paths:

  A character vector of `.sdrf.tsv` paths, in the order to report them.
  Each may be given once.

- threads:

  Documents read at once. `1` by default; `-1` uses every core. The
  result is the same at any value.

- on_error:

  `"fail"`, the default, or `"skip"`, as for
  [`sdrf_validate_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md).

- timeout:

  Seconds to allow for the whole batch, or `NULL`, the default, to wait
  indefinitely.

## Value

An `mzlibr_sdrf_samples_batch`: `file_count` documents given,
`read_count` documents read, `failed_count` documents not read,
`sample_count` samples over every document read, `record_count` rows in
`records`, and `caveats`.

`records` has the columns of
[`sdrf_samples`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md) -
`sample_row_count` in rows; `age_years`, `age_min_years` and
`age_max_years` in years - preceded by `source_index` (1-based, indexing
`files`) and `source_path`. \*\*Key a sample on `source_index` and
`source_name` together\*\*: a source name is unique only within a
document. `files` has one row per input: `path`, `sample_count`,
`row_count`, `conflict_count`, a `problems` list column, and
`error_kind` / `error_message`, `NA` for a document that was read.

## Wraps

Wire verb `sdrf samples` with `--paths-stdin`. Generated from the
bridge's verb spec `sdrf.samples.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfSampleBlock.BySourceName`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfSampleBlock.cs)
  in `mzLib/Readers/Sdrf/SdrfSampleBlock.cs` at mzLib `23c2490e`

- [`SdrfAge.TryParse`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfAge.cs)
  in `mzLib/Readers/Sdrf/SdrfAge.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `paths` (wire `--paths-stdin`):

  flag; default absent. Many paths on stdin, as for sdrf validate.

- `threads`:

  int; in **documents**; default `1`; range `>= 1, or -1`. Documents
  read at once; output identical at any value.

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. As for sdrf validate.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **documents**; never `NA`. Paths given.

- `read_count`:

  int; in **documents**; never `NA`. Documents read.

- `failed_count`:

  int; in **documents**; never `NA`. Documents not read; non-zero only
  under `--on-error` skip.

- `record_count`:

  int; in **rows**; never `NA`. Rows in columns over every document read
  (BULK.md section 2).

- `sample_count`:

  int; in **samples**; never `NA`. Samples over every document read.

- `files`:

  object\[\]; never `NA`. One per input, in input order: {path,
  sample_count, row_count, conflict_count, problems, error}; all but
  path and error null when unread.

- `column_names`:

  string\[\]; never `NA`. source_index, source_path, then the
  single-document columns.

- `records` (wire `columns`):

  table; never `NA`. Rows in input order, then the order within each
  document, whatever `--threads` is.

- `caveats`:

  string\[\]; never `NA`. As for one document.

## Columns of `records`

The first two columns say which input each row came from -
`source_index` is 1-based in R, so it indexes `files` directly - and
rows are grouped by input, in input order, whatever `threads` is.

- `source_index`:

  int; never `NA`. 0-based input position.

- `source_path`:

  string; never `NA`. Absolute path.

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
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.samples`; many files:
  `pymzlib.sdrf.samples_many`

- Rust (mzLibRust): `mzlib::sdrf::samples` with `BulkOptions`; many
  files: `mzlib::sdrf::samples_many`

- R (mzLibR): `sdrf_samples`; many files: `sdrf_samples_many`

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

[`sdrf_samples`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)

## Examples

``` r

batch <- sdrf_samples_many(c("sdrf_cohort.sdrf.tsv", "sdrf_cohort_partner.sdrf.tsv"))
batch
#> <mzlibr_sdrf_samples_batch> sdrf samples
#>   2 inputs: 2 read, 0 failed
#>   62 rows in `records`, 62 records
batch$files[, c("path", "sample_count", "conflict_count")]
#>                           path sample_count conflict_count
#> 1         sdrf_cohort.sdrf.tsv            6              1
#> 2 sdrf_cohort_partner.sdrf.tsv            2              0
```
