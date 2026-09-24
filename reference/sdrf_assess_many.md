# Assess many SDRF files in one bridge call

Decide whether each SDRF file's sample half describes an experimental
design - Informative, Partial or Skeleton - with the per-column counts
behind the verdict.

The gate for a corpus: which of these documents describe an experimental
design?

## Usage

``` r
sdrf_assess_many(paths, threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- paths:

  A character vector of `.sdrf.tsv` paths, in the order to report them.
  Each may be given once.

- threads:

  Documents assessed at once. `1` by default; `-1` uses every core. The
  result is the same at any value.

- on_error:

  `"fail"`, the default, or `"skip"`, as for
  [`sdrf_validate_many`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md).

- timeout:

  Seconds to allow for the whole batch, or `NULL`, the default, to wait
  indefinitely.

## Value

An `mzlibr_sdrf_assessment_batch`: `file_count` documents given,
`read_count` documents assessed, `failed_count` documents not read,
`record_count` rows in `records`, `verdict_counts` (a list counting
documents per verdict: `informative`, `partial`, `skeleton`), and
`caveats`.

`records` has the evidence columns of
[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md) -
`rows`, `filled` and `absent` in rows, `distinct_values` in values,
`fill_rate` a fraction (0-1) - preceded by `source_index` (1-based,
indexing `files`) and `source_path`. `files` has one row per input:
`path`, `verdict`, the three checks, `row_count`, and `error_kind` /
`error_message`, `NA` for a document that was read; an unread document
has `NA` for every other fact.

## Wraps

Wire verb `sdrf assess` with `--paths-stdin`. Generated from the
bridge's verb spec `sdrf.assess.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfSampleInformativeness.Assess`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfSampleInformativeness.cs)
  in `mzLib/Readers/Sdrf/SdrfSampleInformativeness.cs` at mzLib
  `23c2490e`

- [`SdrfColumnCoverage`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfCoverage.cs)
  in `mzLib/Readers/Sdrf/SdrfCoverage.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `paths` (wire `--paths-stdin`):

  flag; default absent. Many paths on stdin, as for sdrf validate.

- `threads`:

  int; in **documents**; default `1`; range `>= 1, or -1`. Documents
  assessed at once; output identical at any value.

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

- `verdict_counts`:

  object; in **documents**; never `NA`. {informative, partial, skeleton}
  over the documents read.

- `files`:

  object\[\]; never `NA`. One per input, in input order: {path, verdict,
  factor_value_varies, sample_is_described, biological_replicate_varies,
  row_count, error}; all but path and error null when unread.

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
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.assess`; many files:
  `pymzlib.sdrf.assess_many`

- Rust (mzLibRust): `mzlib::sdrf::assess` with `BulkOptions`; many
  files: `mzlib::sdrf::assess_many`

- R (mzLibR): `sdrf_assess`; many files: `sdrf_assess_many`

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

[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)

## Examples

``` r

batch <- sdrf_assess_many(
  c("sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv")
)
batch$files[, c("path", "verdict")]
#>                     path     verdict
#> 1   sdrf_cohort.sdrf.tsv Informative
#> 2 sdrf_skeleton.sdrf.tsv    Skeleton
#> 3     PXD000070.sdrf.tsv     Partial
unlist(batch$verdict_counts)
#> informative     partial    skeleton 
#>           1           1           1 
# The documents fit to group results by biology:
batch$files$path[batch$files$verdict != "Skeleton"]
#> [1] "sdrf_cohort.sdrf.tsv" "PXD000070.sdrf.tsv"  
```
