# Restrict an SDRF to the files one search read, and report what could not be placed

Restrict an SDRF to the files one search read, joining each searched
file to the row of its acquisition, and report what the join could not
place.

Answers **which rows of this SDRF did my search read?** A search reads
`X-calib.mzML` while the SDRF names the acquired `X.raw`. Calls mzLib's
`SdrfSearchScope.Restrict`, which joins them **by stem**, ignoring case,
after removing the derivative suffixes MetaMorpheus writes (`-calib`,
`-averaged`, in any order and number). Each kept row is carried cell for
cell, with the searched file's name in a new
`comment[searched data file]` column directly after
`comment[data file]`; nothing is renumbered.

## Usage

``` r
sdrf_restrict(path, searched_files, out = NULL, limit = NULL, offset = 0, timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- searched_files:

  The files the search read, as paths or bare names, a character vector.
  **Name an element by the acquired file it came from** when you know
  it - `c("X-calib.mzML", "X.raw" = "sample_A.mzML")` - and that
  overrides the suffix rule for it; leave the others unnamed. Sent on
  stdin, one line each.

- out:

  Also write the **whole** restricted document here as SDRF, whatever
  `limit` and `offset` are. Must end in `.tsv`; name it `.sdrf.tsv`,
  because mzLib adds that to a name that does not end in it.
  `written$path` names the file written.

- limit:

  Maximum rows to return. `NULL`, the default, returns all of them.

- offset:

  Rows to skip.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**Restrict first, then read a design.**
[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
and
[`sdrf_isobaric_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md)
match `searched_files` by exact name, and key on
`comment[searched data file]` when a document has it, so a restricted
document joins the search's own file names.

## Value

An `mzlibr_restricted_sdrf`, which is also an `mzlibr_sdrf`, so
[`sdrf_value`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md),
[`sdrf_all`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md)
and
[`sdrf_records`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
work on it. `path` is the SDRF restricted, as given; `columns` the
restricted header, which may repeat names; `rows` one character vector
of raw cells per returned row, in the input's row order - a kept row
shorter than the header is padded with empty cells, so unlike
[`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
the rows are never ragged. `row_count` counts the rows kept, in the
**whole** restricted document, whatever `limit` and `offset`;
`returned_count` the rows returned, starting `offset` rows in;
`truncated` whether rows were left behind; `input_row_count` the rows of
the SDRF before the restriction.

`searched_files` lists the searched files as sent, in order, and
`acquired_name_count` counts the files you named an acquired file for.
What the restriction could not place: `searched_without_row`, the
searched files (file names) the SDRF has no row for - reported, never
invented; `dropped_data_files`, the acquired files the SDRF describes
that this search did not read (a row naming no file is dropped too, and
not listed); `ambiguous`, acquisitions two or more searched files
claimed, as `"stem: searched, searched"` - the first in name order kept
the rows, and the others have none and are not in
`searched_without_row`; and `shared_stem`, searched files whose rows
name two or more **different** acquired files (one run listed as `X.raw`
and `X.mzML`), as `"searched: acquired, acquired"`, every such row kept.
`written` is a list (`path`, the file mzLib wrote, and `row_count`, the
rows of the whole document), or `NULL` without `out`. `caveats` says how
the join works and what it cannot tell you.

## Paths, on any platform

The searched files may be full paths; the join uses the file name. A
path written with the other platform's separator (`D:\search\x.mzML` on
Linux) is reduced here too, but the design readers split only the
running platform's separator, so give
[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
and
[`sdrf_isobaric_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md)
paths as this platform writes them
([`file.path()`](https://rdrr.io/r/base/file.path.html) does).

## Wraps

Wire verb `sdrf restrict`. Generated from the bridge's verb spec
`sdrf.restrict.yaml` (bridge commit `45833ca804d6`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfSearchScope.Restrict`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfSearchScope.cs)
  in `mzLib/Readers/Sdrf/SdrfSearchScope.cs` at mzLib `d10017ce`

- [`SdrfSearchScoping`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfSearchScope.cs)
  in `mzLib/Readers/Sdrf/SdrfSearchScope.cs` at mzLib `d10017ce`

- [`SdrfDocument.WriteResults`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfDocument.cs)
  in `mzLib/Readers/Sdrf/SdrfDocument.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `path`:

  path; required. One .sdrf.tsv file, deposited, improved or drafted. A
  file that does not exist is a usage error.

- `searched_files` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line of 1 or 2 tab-separated fields`. The files the search read,
  one per line, as paths or bare names: searched, or searched TAB
  acquired when the caller knows the acquired file a searched one came
  from (MetaMorpheus reads it from an mzML's own sourceFile). An
  acquired name overrides mzLib's suffix rule for that file and may be
  given for some lines and not others; with none, mzLib is given no map
  (its default, acquiredNameOf = null). Blank lines are skipped and a
  UTF-8 BOM is dropped; fields are trimmed. No line, a blank searched
  field, more than two fields, or one searched file given two different
  acquired names (compared ignoring case, as mzLib collapses searched
  files) is a usage error.

- `out`:

  path; default absent; range `.tsv only, any case (PYB-4)`. Also write
  the WHOLE restricted document here, whatever offset and limit are,
  through mzLib's SdrfDocument.WriteResults: atomically (a temp file
  moved into place), CRLF line endings, a final line break. Must end in
  .tsv, case-insensitive, checked before the SDRF is read (bridge thread
  010, PYB-4, as sdrf pool `--out`). A missing folder is created
  (PYB-5). mzLib adds .sdrf.tsv to a name that does not already end in
  it, so searched.tsv is written as searched.tsv.sdrf.tsv; written.path
  says which.

- `offset`:

  int; in **rows**; default `0`; range `>= 0`. Skip this many restricted
  rows in the payload. The file `--out` writes is never windowed.

- `limit`:

  int; in **rows**; default absent; range `>= 0`. Return at most this
  many restricted rows. Absent means every row.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The SDRF path as given.

- `searched_files`:

  string\[\]; never `NA`. The searched field of every stdin line,
  trimmed, in stdin order, as passed to mzLib (duplicates included;
  mzLib collapses those that differ only in case).

- `acquired_name_count`:

  int; in **files**; never `NA`. Searched files given an acquired name;
  0 when none was, and mzLib was then given no map.

- `input_row_count`:

  int; in **rows**; never `NA`. Data rows in the SDRF before the
  restriction.

- `columns` (wire `column_names`):

  string\[\]; never `NA`. The restricted header: the input's, with
  comment\[searched data file\] directly after comment\[data file\] and
  any earlier comment\[searched data file\] removed. A LIST that may
  repeat names, as for sdrf read.

- `row_count`:

  int; in **rows**; never `NA`. Rows kept (SdrfSearchScoping.Document),
  before the window.

- `returned_count`:

  int; in **rows**; never `NA`. Rows in rows.

- `offset`:

  int; in **rows**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever rows were left out by offset or limit,
  including an offset past the end.

- `rows`:

  string\[\]\[\]; never `NA`. One list of raw cells per kept row,
  row-major, as sdrf read shapes a document, in input row order. Never
  ragged: mzLib pads a short row with empty cells.

- `searched_without_row`:

  string\[\]; never `NA`. SdrfSearchScoping.SearchedWithoutRow: searched
  files (file names) the SDRF has no row for. Reported, never invented.

- `dropped_data_files`:

  string\[\]; never `NA`. SdrfSearchScoping.DroppedDataFiles: acquired
  files the SDRF describes that this search did not read, in row order.
  A row naming no file is dropped and not listed.

- `ambiguous`:

  string\[\]; never `NA`. SdrfSearchScoping.Ambiguous: acquisitions two
  or more searched files claimed, as stem: searched, searched. The first
  in ordinal name order kept the rows.

- `shared_stem`:

  string\[\]; never `NA`. SdrfSearchScoping.SharedStem: searched files
  whose rows name two or more DIFFERENT acquired files, as searched:
  acquired, acquired. Every such row is kept. The channel rows of one
  multiplexed file name one file, so are not listed.

- `written`:

  object; `NA` when out was not given. {path, row_count}: path is the
  file mzLib wrote (the out path, plus .sdrf.tsv when it did not end in
  it); row_count is the whole restricted document.

- `caveats`:

  string\[\]; never `NA`. Six static caveats.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  path missing or the file does not exist

- `mzlib_usage_error` (usage):

  stdin has no searched file; a line with a blank searched field or more
  than two fields; one searched file given two different acquired names

- `mzlib_usage_error` (usage):

  `--out` without a value, or not ending in .tsv (checked before the
  SDRF is read)

- `mzlib_usage_error` (usage):

  offset or limit negative or not an integer

- `mzlib_bridge_error` (correctness):

  the SDRF has no comment\[data file\] column (mzLib's
  ArgumentException), or mzLib cannot read it (MzLibException), or the
  `--out` folder cannot be created or written

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- The join is by STEM, case-insensitively, after mzLib removes
  MetaMorpheus's derivative suffixes (-calib, -averaged, in any order
  and number), so X-calib.mzML finds the row naming X.raw. A searched
  TAB acquired line overrides that rule for its file. A row whose stem
  joins nothing is tried once more by an earlier search's
  comment\[searched data file\].

- Nothing else changes: no replicate or fraction is renumbered.
  comment\[searched data file\] goes directly after comment\[data
  file\]; an earlier search's is replaced, because it is this search's
  field.

- A kept row shorter than the header is padded with empty cells, so the
  document is never ragged, unlike what sdrf read returns for the same
  file.

- ambiguous: the other claimants of an acquisition get no row and are
  not in searched_without_row. shared_stem: the document can describe
  one run twice, and only the depositor can say which row is right.

- The design verbs match searched files by exact name. Restrict first:
  sdrf design and sdrf isobaric-design key on comment\[searched data
  file\] when a document has it.

- `--out` writes the whole document with CRLF line endings, whatever the
  input used, and mzLib adds .sdrf.tsv to a name that does not end in
  it. Name the file .sdrf.tsv.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.restrict`

- Rust (mzLibRust): `mzlib::sdrf::restrict_with` with `RestrictOptions`

- R (mzLibR): `sdrf_restrict`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- PUSH-TO-MAINLAND: SdrfDocument.WriteResults returns nothing and adds
  .sdrf.tsv to a name not ending in it, so the bridge reports the path
  by mzLib's own test (ResultFile.CanRead) instead of being told. Ask
  mzLib to return the path written, or to refuse rather than rename.
  sdrf pool `--out` has the same rename and still reports the path as
  given.

- No `--paths-stdin` form: stdin already carries the searched files, and
  one search restricts one SDRF.

- SdrfSampleCopy (#1451, same PR) is not projected: it only feeds
  SdrfBuilder.Build, which is G-sdrf-verbs.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics, the format restricted (Dai et al., Nat Commun 12,
  5854, 2021)

## See also

[`sdrf_searched_data_files`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_searched_data_files.md),
[`sdrf_isobaric_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md),
[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)

## Examples

``` r

# Two calibrated runs, the raw file of one of them as well, a renamed file whose acquisition you
# know, and one the SDRF does not describe:
pool1 <- "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
pool5 <- "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
searched <- c(paste0("D:\\search\\", pool1, "1-calib.mzML"), paste0(pool1, "2-calib.mzML"),
  paste0(pool1, "2.raw"), "sample_A.mzML", "not_in_sdrf.mzML")
names(searched) <- c("", "", "", paste0(pool5, "1.raw"), "")
scoped <- sdrf_restrict("PXD008841.subset.sdrf.tsv", searched, limit = 2)
scoped
#> <mzlibr_restricted_sdrf> PXD008841.subset.sdrf.tsv: 30 of 40 rows kept for 5 searched file(s)
#>   32 columns (some names repeat - use sdrf_all())
#>   30 rows in the document, 2 returned
#>   ! truncated - rows were left behind
#>   ! The join is by STEM, case-insensitively, after mzLib removes MetaMorpheus's derivative suffixes (-calib, -averaged, in any order and number), so X-calib.mzML finds the row naming X.raw. A searched<TAB>acquired line overrides that rule for its file. A row whose stem joins nothing is tried once more by an earlier search's comment[searched data file].
#>   ! Kept rows are carried cell for cell: no replicate or fraction is renumbered. mzLib adds comment[searched data file] (the searched file's name, no folder) directly after comment[data file], and replaces one an earlier search wrote.
#>   ! A kept row shorter than the header is padded with empty cells, so this document is never ragged, unlike what sdrf read returns for the same file.
#>   ! ambiguous: two or more searched files claimed one acquisition; the first in name order kept it, and the others get no row and are not in searched_without_row. shared_stem: one searched file joined rows naming two or more DIFFERENT acquired files (X.raw and X.mzML); every such row is kept, so the document can describe one run twice, and only the depositor can say which is right. A row naming no file is dropped without being listed in dropped_data_files.
#>   ! The design verbs match searched files by EXACT name. Restrict first: sdrf design and sdrf isobaric-design key on comment[searched data file] when a document has it, so the restricted document joins the search's own file names.
#>   ! --out writes the WHOLE restricted document, whatever limit and offset are, through mzLib's SdrfDocument.WriteResults: atomically, with CRLF line endings whatever the input used. mzLib adds .sdrf.tsv to a name that does not already end in it, so written.path is the file actually written.
#>   searched without a row: not_in_sdrf.mzML
#>   acquired but not searched: HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr02.raw
#>   ambiguous: HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02: HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02-calib.mzML, HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02.raw
c(input = scoped$input_row_count, kept = scoped$row_count, returned = scoped$returned_count)
#>    input     kept returned 
#>       40       30        2 
sdrf_searched_data_files(scoped)
#> [1] "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
#> [2] "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
scoped$searched_without_row
#> [1] "not_in_sdrf.mzML"
scoped$dropped_data_files
#> [1] "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr02.raw"
writeLines(strwrap(scoped$ambiguous, 80))
#> HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02:
#> HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02-calib.mzML,
#> HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02.raw
```
