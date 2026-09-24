# The complete file list of a PRIDE Archive project, from its FTP directory tree

List every file a PRIDE project holds by walking its FTP directory tree,
subdirectories included, with approximate sizes.

The authoritative counterpart to
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md).
Where that returns PRIDE's REST manifest - knowingly incomplete for some
projects, omitting for PXD000001 the two largest of 13 files - this
walks the project's FTP directory (subdirectories included) and returns
everything the project holds (mzLib \#1121). Reach for it whenever
completeness or a true project size matters.

## Usage

``` r
pride_list_ftp_files(accession, timeout = 300)
```

## Arguments

- accession:

  A project accession, e.g. `"PXD000001"`.

- timeout:

  Seconds to allow for the whole walk, which spans one request per
  directory.

## Details

This is a **listing** surface only.
[`pride_download`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
and
[`pride_download_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
operate on the REST manifest, so a file that appears \*only\* here - the
whole point of this function - is fetched directly from its `url` with
an ordinary HTTPS client (e.g. `download.file(f$url, f$file_name)`).

## Value

A data.frame with one row per file and columns `relative_path`,
`file_name`, `url`, `approximate_size_bytes`, `approximate_size_mb`, and
`extension`.

Sizes are PRIDE's rounded directory-index values - good for a
project-size estimate but not exact. For the precise transfer size of
one file, issue an HTTP HEAD against its `url`. See
[`pride_approximate_total_size_bytes`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_approximate_total_size_bytes.md).

An unknown accession **raises** `mzlib_project_not_found` (as
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
does): mzLib resolves the project before walking, so a typo fails loudly
rather than returning nothing.

## Wraps

Wire verb `pride ftp-files`. Generated from the bridge's verb spec
`pride.ftp-files.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PrideArchiveClient.GetProjectFilesFromFtpAsync`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideArchiveClient.GetProjectAsync`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideFtpFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideFtpFile.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideFtpFile.cs` at mzLib
  `23c2490e`

- [`PrideArchiveExtensions.TotalApproximateSizeBytes`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs` at
  mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `accession`:

  string; required. The project accession, e.g. PXD000001. The bridge
  checks only that it is not blank. mzLib resolves the project first
  (for its publication date, which locates the FTP directory), so an
  unknown accession fails instead of listing nothing.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `accession`:

  string; never `NA`. The accession, echoed as given.

- `file_count`:

  int; in **files**; never `NA`. Entries in files. 0 only if the
  resolved project's directory listed nothing (the bindings raise their
  not-found error on it).

- `approximate_total_size_bytes`:

  int; in **bytes**; never `NA`. Sum of approximate_size_bytes (mzLib
  TotalApproximateSizeBytes): an estimate over the COMPLETE file list.
  PXD000001 comes to about 1.44 GB.

- `files`:

  object\[\]; never `NA`. One entry per file in walk order (each
  directory's autoindex order, depth first); fields under
  result.columns.

## Columns

- `relative_path`:

  string; never `NA`. Path from the project's FTP root, '/'-separated
  and URL-decoded, e.g. generated/summary.mztab for a file in a
  subdirectory.

- `file_name`:

  string; never `NA`. The last segment of relative_path.

- `url`:

  string; never `NA`. HTTPS URL on ftp.pride.ebi.ac.uk that the file can
  be fetched from with any HTTP client.

- `approximate_size_bytes`:

  int; in **bytes**; never `NA`. PRIDE's autoindex size (e.g. '429M'),
  expanded with 1024-based multipliers and so rounded to about three
  significant figures. 0 when the index showed '-' or an unparseable
  size: 0 means unknown, not empty. For the exact size, HTTP HEAD the
  url and read Content-Length.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  accession missing or blank

- `mzlib_service_unavailable` (service_unavailable):

  the bridge's own classification (Program.ClassifyError; mzLib#1350 is
  not merged): a timeout or cancellation, a socket failure, a request
  that never got a response (refused connection, DNS, TLS), an HTTP 408,
  429 or 5xx on the project lookup or on any directory listing, or a
  body cut off in transit

- `mzlib_bridge_error` (correctness):

  no project has that accession, or it has no publication date to locate
  its directory (MzLibException; the bindings re-map it to their
  not-found error); a directory listing answering any other status;
  nesting deeper than 64 levels (MzLibException: the listing may be
  cyclic)

## Caveats

- This is the complete file list; pride files (PRIDE's REST manifest)
  omits files, for PXD000001 five of 13 including the two largest. The
  price is metadata: no category, checksum or CV locations here.

- Sizes are approximate (rounded autoindex values), so
  approximate_total_size_bytes is an estimate, but an estimate over
  every file, unlike pride files' total_size_bytes.

- A listing surface only: pride download selects from the REST manifest,
  so a file that appears only here cannot be fetched by pride download.
  Fetch it from its url with an ordinary HTTPS client.

- One HTTP request per directory: the walk's cost grows with the tree,
  not the file count.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.pride.list_ftp_files`

- Rust (mzLibRust): `mzlib::pride::list_ftp_files_with`

- R (mzLibR): `pride_list_ftp_files`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- No recorded fixture: the verb needs a pride_PXD000001_ftp_files.json
  recording (live, 13 files) before it can have a checked example, and
  before the bindings' doctests can drop +SKIP.

- An unknown accession crosses as correctness (MzLibException), not
  usage, although it is the caller's mistake. Every binding re-maps it
  to its not-found error; the wire itself does not say 'not found'.

- Directory order is whatever PRIDE's Apache autoindex serves; not
  verified to be stable across requests, so files order is not promised.

- mzlib::pride::list_ftp_files_with takes the timeout as its second
  argument rather than an options struct; no Rust options type exists
  for this verb.

## References

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive and its FTP layout (Perez-Riverol et al., The PRIDE
  database at 20 years: 2025 update, NAR)

## See also

[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
for the rich REST metadata;
[`pride_approximate_total_size_bytes`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_approximate_total_size_bytes.md).

## Examples

``` r

listing <- pride_list_ftp_files("PXD000001")
listing[, c("relative_path", "approximate_size_mb")]
#>             relative_path approximate_size_mb
#> 1              README.txt            0.001638
#> 2                run1.raw          220.200960
#> 3   hidden_from_rest.mzML          449.839104
#> 4 generated/summary.mztab            0.864256
pride_approximate_total_size_bytes(listing)
#> [1] 670905958
```
