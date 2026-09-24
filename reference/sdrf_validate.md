# Check one SDRF file against the specification's structural rules

Check SDRF-Proteomics files against the specification's structural
rules, one row per finding.

Calls mzLib's `SdrfValidator.Validate`. Structure only: required and
recommended columns, column order, casing, malformed names, ragged rows,
integer replicate and fraction columns, reserved-word casing, and the
one hard row rule - `source name` + `assay name` + the `comment` column
for the label must be unique. Controlled-vocabulary accessions are
**not** resolved.

## Usage

``` r
sdrf_validate(path, timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_sdrf_validation`: `path`, `is_valid`, and the counts -
`error_count` and `warning_count` findings of each severity,
`message_count` findings in all, `row_count` data rows in the document
(the header is not a row) - plus `caveats`, `column_names`, and
`records`, a data.frame with one row per finding: `severity` (`"Error"`
or `"Warning"`), `rule` (a stable name such as `"RequiredColumn"`),
`message`, `row_index` (the data row, 1-based), `line_number` (the file
line, 1-based, the header being line 1) and `column_name`. `row_index`
and `line_number` are `NA` for a finding about the whole document;
`column_name` is `NA` when the finding is not about one column.

## Warnings never make a document invalid

`is_valid` is `TRUE` whenever there is no `Error`. mzLib calibrated
every severity against the 1,236-file curated corpus, and a rule that
fired on most curated files was judged wrong rather than the files. And
a valid document may still say nothing about its samples: that is what
[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)
asks.

## Wraps

Wire verb `sdrf validate`. Generated from the bridge's verb spec
`sdrf.validate.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfValidator.Validate`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfValidator.cs)
  in `mzLib/Readers/Sdrf/SdrfValidator.cs` at mzLib `23c2490e`

- [`SdrfValidationResult`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfValidation.cs)
  in `mzLib/Readers/Sdrf/SdrfValidation.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; default absent. One .sdrf.tsv file. Exactly one of path and
  paths-stdin is required.

- wire `--paths-stdin`:

  flag; default absent. Read many paths from stdin, one per line
  (BULK.md section 1). Blank lines are ignored; a repeated path is a
  usage error. Switches the result to the bulk shape below. *Not an
  argument here: sdrf_validate_many() is the many-documents form.*

- wire `--threads`:

  int; in **documents**; default `1`; range `>= 1, or -1`. Documents
  validated at once; -1 is every core. The output is byte-identical at
  any value (tested at 1, 4 and -1). *Not an argument here: only
  sdrf_validate_many() takes threads.*

- wire `--on-error`:

  string; default `fail`; range `fail | skip`. fail: the first failure
  IN INPUT ORDER aborts the call. skip: recorded in files\[i\].error and
  the batch carries on. skip with path is a usage error. *Not an
  argument here: only sdrf_validate_many() takes on_error.*

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The absolute path validated.

- `is_valid`:

  bool; never `NA`. No Error findings. Warnings never make a document
  invalid.

- `error_count`:

  int; in **findings**; never `NA`. Findings of severity Error.

- `warning_count`:

  int; in **findings**; never `NA`. Findings of severity Warning.

- `message_count`:

  int; in **findings**; never `NA`. All findings; the length of every
  column.

- `row_count`:

  int; in **rows**; never `NA`. Data rows in the document.

- `column_names`:

  string\[\]; never `NA`. The table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per finding, in mzLib's
  order.

- `caveats`:

  string\[\]; never `NA`. Four static caveats.

## Columns of `records`

- `severity`:

  string; never `NA`. Error (cannot be reliably consumed) or Warning
  (deviates, still joinable). mzLib's SdrfValidationSeverity member
  name.

- `rule`:

  string; never `NA`. Stable rule id: EmptyHeader, RequiredColumn,
  RecommendedColumn, EmptyColumnName, ColumnNameCase,
  MalformedColumnName, ColumnOrdering, NoRows, RowWidth, EmptyCell,
  ReservedWordCase, NonIntegerValue, RowKeyUniqueness (at 23c2490e).

- `message`:

  string; never `NA`. Human-readable, including the offending value.

- `row_index`:

  int; `NA` when the finding is about the whole document (a missing
  column, an empty header, no rows), not one row. 0-based index into the
  data rows.

- `line_number`:

  int; `NA` when same as row_index: a document-level finding. 1-based
  line in the file: row_index + 2, because the header is line 1.

- `column_name`:

  string; `NA` when the finding is not about one column (RowWidth,
  RowKeyUniqueness, NoRows, EmptyHeader). The column involved; for a
  ragged row's extra cells, "(column N)".

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md).

- `mzlib_usage_error` (usage):

  neither or both of path and paths-stdin; stdin empty under
  paths-stdin; a repeated path

- `mzlib_usage_error` (usage):

  a file not found (under fail, checked before any document is read)

- `mzlib_usage_error` (usage):

  threads 0 or below -1 or not an integer; on-error not fail or skip, or
  given without a value; on-error skip with path

- `mzlib_bridge_error` (correctness):

  mzLib cannot read a document (MzLibException for an empty file); under
  fail, the lowest input position's failure is raised

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- STRUCTURE ONLY: column presence, order, casing, cell shape and row-key
  uniqueness against v1.1.0. Controlled-vocabulary accessions are not
  resolved.

- Severities were calibrated by mzLib against the 1,236-file bigbio
  curated corpus; a rule that fired on most curated files was demoted,
  not the files condemned.

- A document of reserved words validates cleanly. Use sdrf assess for
  whether it says anything, sdrf lint for whether several agree.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. For many files call
[`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.validate`; many files:
  `pymzlib.sdrf.validate_many`

- Rust (mzLibRust): `mzlib::sdrf::validate` with `BulkOptions`; many
  files: `mzlib::sdrf::validate_many`

- R (mzLibR): `sdrf_validate`; many files:
  [`sdrf_validate_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md)

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- since.pymzlib is null until the pyMzLib release carrying
  feat/sdrf-1.0.592 is tagged; Rust and R are not yet ported.

- Bulk spellings are the cross-binding decision of 2026-09-23: a
  separate \<verb\>\_many in every binding (R too, not vectorisation).
  Rust and R names here are intended, not yet ported.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics; mzLib validates against specification v1.1.0 (Dai et
  al., Nat Commun 12, 5854, 2021)

## See also

[`sdrf_validate_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md),
[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md),
[`sdrf_lint`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md)

## Examples

``` r

result <- sdrf_validate("sdrf_skeleton.sdrf.tsv")
result
#> <mzlibr_sdrf_validation> sdrf_skeleton.sdrf.tsv: NOT valid
#>   7 error(s), 2 warning(s) over 2 rows
#>   ! Error RequiredColumn x7
#>   ! Warning RecommendedColumn x2
#>   4 caveat(s) - read x$caveats once
result$records[result$records$severity == "Error", c("rule", "line_number", "column_name")]
#>             rule line_number                                 column_name
#> 1 RequiredColumn          NA                             technology type
#> 2 RequiredColumn          NA                         comment[instrument]
#> 3 RequiredColumn          NA                              comment[label]
#> 4 RequiredColumn          NA             comment[cleavage agent details]
#> 5 RequiredColumn          NA                comment[technical replicate]
#> 6 RequiredColumn          NA                comment[fraction identifier]
#> 7 RequiredColumn          NA comment[proteomics data acquisition method]
```
