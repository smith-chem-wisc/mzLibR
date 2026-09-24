# Validate many SDRF files in one bridge call

Check SDRF-Proteomics files against the specification's structural
rules, one row per finding.

One bridge process for the whole corpus instead of one per file, with
the parallelism inside the bridge where it belongs. The result does not
depend on `threads`: rows are always in input order.

## Usage

``` r
sdrf_validate_many(paths, threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- paths:

  A character vector of `.sdrf.tsv` paths, in the order to report them.
  Each may be given once.

- threads:

  Documents validated at once. `1`, the default, holds one document in
  memory at a time; `-1` uses every core.

- on_error:

  `"fail"`, the default, raises on the first unreadable document in
  input order. `"skip"` records the failure in that document's row of
  `files` and carries on.

- timeout:

  Seconds to allow for the whole batch, or `NULL`, the default, to wait
  indefinitely.

## Value

An `mzlibr_sdrf_validation_batch`: `file_count` documents given,
`read_count` documents validated, `failed_count` documents that could
not be read, `valid_count` documents read with no `Error`,
`message_count` findings over every document read, `record_count` rows
in `records`, and `caveats`.

`records` has the columns of
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)'s,
preceded by `source_index` - the document's 1-based position in `paths`,
so it indexes `files` - and `source_path`. `files` has one row per
input: `path`, `is_valid`, `error_count`, `warning_count`,
`message_count`, `row_count`, and `error_kind` / `error_message`, which
are `NA` for a document that was read. A document that was not read has
`NA` for every count.

## Wraps

Wire verb `sdrf validate` with `--paths-stdin`. Generated from the
bridge's verb spec `sdrf.validate.yaml` (bridge commit `5db922d4cfe1`)
by `scripts/build-man.R`; the spec owns these facts, and all three
bindings render the same ones.

- [`SdrfValidator.Validate`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfValidator.cs)
  in `mzLib/Readers/Sdrf/SdrfValidator.cs` at mzLib `23c2490e`

- [`SdrfValidationResult`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfValidation.cs)
  in `mzLib/Readers/Sdrf/SdrfValidation.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `paths` (wire `--paths-stdin`):

  flag; default absent. Read many paths from stdin, one per line
  (BULK.md section 1). Blank lines are ignored; a repeated path is a
  usage error. Switches the result to the bulk shape below.

- `threads`:

  int; in **documents**; default `1`; range `>= 1, or -1`. Documents
  validated at once; -1 is every core. The output is byte-identical at
  any value (tested at 1, 4 and -1).

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. fail: the first failure
  IN INPUT ORDER aborts the call. skip: recorded in files\[i\].error and
  the batch carries on. skip with path is a usage error.

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

- `valid_count`:

  int; in **documents**; never `NA`. Documents read with no Error.

- `message_count`:

  int; in **findings**; never `NA`. Findings over every document read;
  equals record_count.

- `files`:

  object\[\]; never `NA`. One per input, in input order: {path,
  is_valid, error_count, warning_count, message_count, row_count,
  error}. All but path and error are null when the document was not
  read; error is null when it was, else {kind: usage \| correctness,
  message}.

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

  int; never `NA`. 0-based position of the document in the stdin list.

- `source_path`:

  string; never `NA`. Absolute path of the document.

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
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.validate`; many files:
  `pymzlib.sdrf.validate_many`

- Rust (mzLibRust): `mzlib::sdrf::validate` with `BulkOptions`; many
  files: `mzlib::sdrf::validate_many`

- R (mzLibR): `sdrf_validate`; many files: `sdrf_validate_many`

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

[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)

## Examples

``` r

batch <- sdrf_validate_many(
  c("sdrf_skeleton.sdrf.tsv", "sdrf_cohort.sdrf.tsv", "missing.sdrf.tsv"),
  on_error = "skip"
)
batch
#> <mzlibr_sdrf_validation_batch> sdrf validate
#>   3 inputs: 2 read, 1 failed
#>   11 rows in `records`, 11 records
#>   ! missing.sdrf.tsv: SDRF file not found: 'missing.sdrf.tsv'.
batch$files[, c("path", "is_valid", "error_count", "error_kind")]
#>                     path is_valid error_count error_kind
#> 1 sdrf_skeleton.sdrf.tsv    FALSE           7       <NA>
#> 2   sdrf_cohort.sdrf.tsv     TRUE           0       <NA>
#> 3       missing.sdrf.tsv       NA          NA      usage
```
