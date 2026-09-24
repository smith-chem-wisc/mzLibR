# Read one SDRF-Proteomics experimental-design file

Read one SDRF-Proteomics file as its header and rows of raw cells,
keeping repeated column names and ragged rows exactly as written.

Every other reader here answers what a search **found**. SDRF answers
what was **searched** - which sample, which organism part, which
replicate, which instrument settings - which is the half you need to
compare results across experiments.

## Usage

``` r
sdrf_read(path, limit = NULL, offset = 0, timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- limit:

  Maximum rows to return. `NULL`, the default, returns all of them.

- offset:

  Rows to skip.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**Use this, not
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md),
for SDRF.**
[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
joins each row's cells into one semicolon-separated string, and SDRF's
own `NT=...;AC=...` grammar puts semicolons inside cells, so that string
cannot be split back apart.

## Value

An `mzlibr_sdrf`: `path`, `columns` (a character vector that may
repeat), `rows` (a list of character vectors, one per row), `row_count`
(the whole document), `returned_count`, `offset`, `truncated` and
`caveats`.

## The shape is row-major, by necessity

SDRF column names are data, and they **repeat**: 649 files in the
curated corpus carry `comment[modification parameters]` more than once,
up to eight times in one file. A data.frame would force those names
unique. So `columns` keeps the duplicates, `rows` holds the cells, and
position links them. Use
[`sdrf_value`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
for the first cell under a name and
[`sdrf_all`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md)
for every one.

**Rows are ragged.** A row may carry fewer cells than `columns` has
entries - PXD059974 in mzLib's own fixtures has a 46-column header with
17 of its 22 rows carrying 42 cells - and mzLib preserves that so the
file round-trips.
[`sdrf_value`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
gives `NA` for a position a row does not reach.

**Cells are raw strings, never interpreted.** The key=value grammar
arrives exactly as written, because it cannot be told apart from a
`comment[file uri]` cell whose pre-signed URL contains `Signature=` and
`Expires=`.

**A reserved word is a real value.** `"not available"` and
`"not applicable"` mean the experiment stated an absence. `NA` means the
document has no such column. Do not collapse them.

## Reading is not validating

This reads whatever the file holds and makes no claim that it is
correct. mzLib's judgements about a document are their own functions:
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)
checks the specification's structural rules,
[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)
asks whether the sample columns describe a design at all,
[`sdrf_samples`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)
merges each sample's rows with its age parsed, and
[`sdrf_lint`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md)
finds concepts several documents wrote inconsistently.

## Wraps

Wire verb `sdrf read`. Generated from the bridge's verb spec
`sdrf.read.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfDocument.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfDocument.cs)
  in `mzLib/Readers/Sdrf/SdrfDocument.cs` at mzLib `23c2490e`

- [`SdrfHeader`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfHeader.cs)
  in `mzLib/Readers/Sdrf/SdrfHeader.cs` at mzLib `23c2490e`

- [`SdrfRow.Cells`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfRow.cs)
  in `mzLib/Readers/Sdrf/SdrfRow.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A .sdrf.tsv file.

- `offset`:

  int; in **rows**; default `0`; range `>= 0`. Skip this many data rows.
  The whole file is parsed on every call.

- `limit`:

  int; in **rows**; default absent; range `>= 0`. Return at most this
  many data rows. Absent means every row.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The path as given.

- `columns` (wire `column_names`):

  string\[\]; never `NA`. The header, verbatim and in order. A LIST that
  may repeat names (comment\[modification parameters\] up to eight
  times) and may contain empty names; a name identifies a position, not
  a column.

- `row_count`:

  int; in **rows**; never `NA`. Data rows in the whole file, before the
  window.

- `returned_count`:

  int; in **rows**; never `NA`. Rows in rows.

- `offset`:

  int; in **rows**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever rows were left out by offset or limit,
  including an offset past the end.

- `rows`:

  string\[\]\[\]; never `NA`. One list of raw cells per row, row-major.
  RAGGED: a row may be shorter than column_names; a missing position is
  absent, not empty.

- `caveats`:

  string\[\]; never `NA`. Three static caveats, plus a count of short
  rows over the WHOLE document when any row is short.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  path missing or blank; offset or limit negative or not an integer

- `mzlib_bridge_error` (correctness):

  mzLib cannot read the file (MzLibException for an empty or missing
  file)

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- readers read-records also opens .sdrf.tsv, but joins each row's cells
  with ';' into one string that cannot be split back, because the SDRF
  grammar is itself ';'-delimited. Use this verb for SDRF.

- offset and limit window after a full parse.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.read`

- Rust (mzLibRust): `mzlib::sdrf::read_with` with `ReadOptions`

- R (mzLibR): `sdrf_read`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- since.mzlibrust / since.mzlibr are null: the ports (mzLibRust#23,
  mzLibR#19) merged 2026-09-17 but neither repo has cut a release tag
  yet (checked with gh api .../tags, 2026-09-23).

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics, the format read (Dai et al., Nat Commun 12, 5854,
  2021)

## See also

[`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md),
[`sdrf_value`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md),
[`sdrf_all`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md),
[`sdrf_records`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md),
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
doc
#> <mzlibr_sdrf> PXD000070.sdrf.tsv
#>   31 columns (some names repeat - use sdrf_all())
#>   6 rows in the document, 6 returned
#>   ! Cells are RAW STRINGS, never interpreted. The SDRF key=value grammar ("NT=Oxidation;AC=UNIMOD:35") is left exactly as written, because it cannot be told apart from a cell that merely contains '=' and ';' - comment[file uri] routinely carries pre-signed URLs that do.
#>   ! column_names is a LIST, not a set: names may repeat, and their order is part of the document. comment[modification parameters] legitimately appears many times in one file, so a name identifies a POSITION rather than a column.
#>   ! An SDRF reserved word - "not available", "not applicable" - is a real value meaning the experiment stated an absence. It is NOT the same as a column the document does not have, and the two must not be collapsed.
sdrf_value(doc, "characteristics[organism]")
#> [1] "plasmodium falciparum" "plasmodium falciparum" "plasmodium falciparum"
#> [4] "plasmodium falciparum" "plasmodium falciparum" "plasmodium falciparum"
```
