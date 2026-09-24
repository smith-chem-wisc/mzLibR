# Merge several SDRF documents into one analysis table

Merge several SDRF documents into one analysis table, stamping each row
with the document it came from.

Columns are the union of every document's, ordered by SDRF's own block
structure, and a name that repeats is carried at the highest
multiplicity any single document used, so nothing is dropped. A cell a
document did not have is filled with the reserved word
`"not available"`, and a `comment[source document]` column records which
document each row came from.

## Usage

``` r
sdrf_pool(documents, out = NULL, limit = NULL, offset = 0, timeout = 60)
```

## Arguments

- documents:

  A character vector of paths. **Name it** to choose each document's
  provenance label:
  `c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv")`.
  Name every element or none.

- out:

  Write the merged document here as SDRF. The **whole** document is
  written regardless of `limit` and `offset`.

- limit:

  Maximum rows to return. `NULL`, the default, returns all of them.

- offset:

  Rows to skip.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_pooled_sdrf`, which is also an `mzlibr_sdrf`: everything
[`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
returns, plus `document_count`, `paths`, `labels`, and `written` (a list
of `path` and `row_count`, or `NULL` when `out` was not given).
`row_count` counts the rows of the whole pooled document;
`returned_count` the rows returned, starting `offset` rows in.

## Give your documents names

With an unnamed vector, mzLib falls back to
`containing-folder/file-stem`, which depends on where the files happen
to sit - so the same documents pooled from a different directory produce
a different table. The returned `caveats` say so when that fallback was
used.

## The result is an analysis table, not something to deposit

Source name, assay name and label together are unique within one
document, but two experiments may both have a `"Sample 1"`, so a pooled
table will usually violate SDRF's uniqueness rule. Use
[`sdrf_source_documents`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_source_documents.md)
as part of any key.

## Wraps

Wire verb `sdrf pool`. Generated from the bridge's verb spec
`sdrf.pool.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfCollection.Merge`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfCollection.cs)
  in `mzLib/Readers/Sdrf/SdrfCollection.cs` at mzLib `23c2490e`

- [`SdrfDocument.WriteResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfDocument.cs)
  in `mzLib/Readers/Sdrf/SdrfDocument.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `documents` (wire `--stdin`):

  string\[\]; required; range `>= 1 line`. One document per line: path,
  optionally a TAB and a provenance label. Label every line or none; a
  partial set is a usage error. Blank lines are ignored.

- `out`:

  path; default absent. Write the WHOLE merged document here as SDRF,
  whatever offset and limit are. The write is atomic since mzLib 1.0.592
  (#1324): a temp file moved into place.

- `offset`:

  int; in **rows**; default `0`; range `>= 0`. Skip this many merged
  rows in the payload.

- `limit`:

  int; in **rows**; default absent; range `>= 0`. Return at most this
  many merged rows. Absent means every row.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `document_count`:

  int; in **documents**; never `NA`. Documents pooled.

- `paths`:

  string\[\]; never `NA`. The paths, in stdin order.

- `labels`:

  string\[\]; never `NA`. The label stamped for each document: the one
  supplied, or mzLib's containing-folder/file-stem default.

- `columns` (wire `column_names`):

  string\[\]; never `NA`. The union header in SDRF block order; a
  repeated name at the highest multiplicity any one document used;
  comment\[source document\] added.

- `row_count`:

  int; in **rows**; never `NA`. Merged rows before the window.

- `returned_count`:

  int; in **rows**; never `NA`. Rows in rows.

- `offset`:

  int; in **rows**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever rows were left out.

- `rows`:

  string\[\]\[\]; never `NA`. Row-major raw cells, as for sdrf read.

- `written`:

  object; `NA` when out was not given. {path, row_count}; row_count is
  the whole merged document.

- `caveats`:

  string\[\]; never `NA`. Four static caveats, plus a fifth when no
  labels were supplied.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  stdin empty; a blank path; a file not found; some lines labelled and
  some not

- `mzlib_usage_error` (usage):

  offset or limit negative or not an integer

- `mzlib_bridge_error` (correctness):

  mzLib cannot read one of the documents

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- An ANALYSIS table, not a deposit: source name + assay name +
  comment\[label\] is unique per document, not across documents. Key on
  comment\[source document\] too.

- Without labels, provenance depends on where the files sit, so the
  table is not reproducible on another machine.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.pool`

- Rust (mzLibRust): `mzlib::sdrf::pool_with` with `PoolOptions`

- R (mzLibR): `sdrf_pool`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- since.mzlibrust / since.mzlibr are null: ported and merged 2026-09-17
  (mzLibRust#23, mzLibR#19), not yet in a tagged release.

- The Rust example's PoolInput::labelled constructor is written from the
  pool_with signature, not checked against mzLibRust's source.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics, the format pooled (Dai et al., Nat Commun 12, 5854,
  2021)

## See also

[`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md),
[`sdrf_source_documents`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_source_documents.md)

## Examples

``` r

pooled <- sdrf_pool(c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv"),
  limit = 4)
pooled
#> <mzlibr_pooled_sdrf> 2 documents: malaria, colon
#>   38 columns (some names repeat - use sdrf_all())
#>   24 rows in the document, 4 returned
#>   ! truncated - rows were left behind
#>   ! Every cell a source document did not have is filled with the reserved word "not available". Those words are therefore partly an artefact of pooling rather than something an experiment said, and a fill-rate computed over this table is not the fill rate of the originals.
#>   ! Provenance is in the 'comment[source document]' column. It is stamped only where a source document had none, so pooling an already-pooled table keeps the inner labels rather than overwriting them.
#>   ! This is an ANALYSIS table, not something to deposit. source name + assay name + comment[label] is unique within one document, but two experiments may both have a "Sample 1", so the pooled table will usually violate SDRF's uniqueness rule. Use the source-document column as part of any key.
#>   ! Columns are the UNION over all 2 documents, ordered by SDRF's own block structure, and a name that repeats is carried at the highest multiplicity any single document used, so no values are dropped.
sdrf_source_documents(pooled)
#> [1] "malaria" "malaria" "malaria" "malaria"
```
