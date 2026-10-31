# The complete file list of a PRIDE Archive project, from its FTP directory tree

List every file a PRIDE project holds by walking its FTP directory tree,
subdirectories included, with approximate sizes.

The authoritative counterpart to
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md).
Where that returns PRIDE's REST manifest - knowingly incomplete for some
projects, omitting for PXD000001 its two largest files - this walks the
project's FTP directory (subdirectories included) and returns everything
the project holds (mzLib \#1121). Reach for it whenever completeness or
a true project size matters.

## Usage

``` r
pride_list_ftp_files(accession, timeout = NULL)
```

## Arguments

- accession:

  A project accession, e.g. `"PXD000001"`.

- timeout:

  Seconds to allow for the whole walk, which spans one request per
  directory. `NULL`, the default, adds no limit of mzLibR's own, because
  mzLib already bounds every request: it retries a transient failure
  three times, after 5, 20 and 60 s, and each attempt can take up to
  100 s. An outage that hangs is therefore reported as
  `mzlib_service_unavailable` after about 8 minutes per request. A
  number cuts the whole call off after that many seconds with
  `mzlib_timeout`, which is not `mzlib_service_unavailable`.

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
`pride.ftp-files.yaml` (bridge commit `84c3b984b6b9`) by
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

  int; in **files**; never `NA`. Entries in files. Since mzLib 1.0.594
  (#1412) an empty root listing for a project that exists is retried and
  then thrown (see errors), so 0 means only that the root listed
  subdirectories that held no files. The bindings raise their not-found
  error on 0.

- `approximate_total_size_bytes`:

  int; in **bytes**; never `NA`. Sum of approximate_size_bytes (mzLib
  TotalApproximateSizeBytes): an estimate over the COMPLETE file list.
  PXD000001 comes to about 1.44 GB.

- `files`:

  object\[\]; never `NA`. One entry per file in walk order (each
  directory's autoindex order, depth first); fields under
  result.columns. mzLib does not sort, so the order is PRIDE's and is
  not promised stable across requests.

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

  the bridge's own classification (Program.ClassifyError; since mzLib
  1.0.593 (#1350) PRIDE raises every transport failure as an
  HttpRequestException carrying the cause or the status, and messages
  name the file and host, never the URL): a timeout or cancellation, a
  socket failure, a request that never got a response (refused
  connection, DNS, TLS), an HTTP 408, 429 or 5xx on the project lookup
  or on any directory listing, or a body cut off in transit

- `mzlib_bridge_error` (correctness):

  no project has that accession, or it has no publication date to locate
  its directory (MzLibException; the bindings re-map it to their
  not-found error); a directory listing answering any other status;
  nesting deeper than 64 levels (MzLibException: the listing may be
  cyclic); the project's root listing still empty after mzLib's retries
  (#1412: an HttpRequestException with no status and no cause, which
  mzLib treats as transient but the bridge cannot tell from a refusal,
  so it crosses as correctness until PRIDE-B3); a 403 from
  ftp.pride.ebi.ac.uk after mzLib's retries (mzLib retries it as EBI's
  rate limit, but the bridge classifies by status, so it crosses as
  correctness; PRIDE-B3 asks mzLib to mark such failures)

## Caveats

- Since mzLib 1.0.594 (#1412), PrideArchiveClient retries a transient
  failure, that request alone, after 5, 20 and 60 s (a 429 or 503
  Retry-After of up to 60 s replaces the wait); transient means no HTTP
  status at all, 408, 429, any 5xx, or a 403 from ftp.pride.ebi.ac.uk,
  EBI's rate limit. Each attempt can take up to the client's 100 s
  timeout, so a service that hangs takes up to about 485 s per request
  before the failure is reported. A binding's own timeout shorter than
  that reports a hanging outage as a timeout, not as
  service_unavailable.

- This is the complete file list; pride files (PRIDE's REST manifest)
  omits files, for PXD000001 five of 14 on 2026-10-21, including the two
  largest. The price is metadata: no category, checksum or CV locations
  here.

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

- An unknown accession crosses as correctness (MzLibException), not
  usage, although it is the caller's mistake; the wire has no 'not
  found'. Each binding re-maps it, but not alike: pyMzLib matches
  mzLib's two not-found messages, while mzLibRust and mzLibR re-map
  every MzLibException, so the 64-level cyclic-listing failure also
  reads as not found there.

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
#>                                                           relative_path
#> 1                                                           F063721.dat
#> 2                                                 F063721.dat-mztab.txt
#> 3                                    PRIDE_Exp_Complete_Ac_22134.xml.gz
#> 4                                      PRIDE_Exp_mzData_Ac_22134.xml.gz
#> 5                                PXD000001_community_annotated.sdrf.tsv
#> 6                                                   PXD000001_mztab.txt
#> 7                                                            README.txt
#> 8   TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML
#> 9  TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzXML
#> 10          TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
#> 11            TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw
#> 12                                             erwinia_carotovora.fasta
#> 13                   generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz
#> 14                 generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz
#>    approximate_size_mb
#> 1            20.971520
#> 2             0.305152
#> 3            10.485760
#> 4             9.751757
#> 5             0.004096
#> 6             0.864256
#> 7             0.001638
#> 8           449.839104
#> 9           472.907776
#> 10          243.269632
#> 11          220.200960
#> 12            1.677722
#> 13            5.976883
#> 14            0.103424
pride_approximate_total_size_bytes(listing)
#> [1] 1436359680
```
