# List the files in a PRIDE Archive project

List a PRIDE Archive project's files as PRIDE's REST API publishes them,
with every page already fetched: name, size, checksum, category and
download locations.

Paging is handled for you: however many pages the project spans, you get
one data.frame.

## Usage

``` r
pride_list_files(accession, page_size = NULL, timeout = NULL)
```

## Arguments

- accession:

  A project accession, e.g. `"PXD000001"`. Case and surrounding
  whitespace are normalised. The check is grammatical only: a
  well-formed accession that does not exist costs one live request
  before it fails, because PXD accessions are not fixed-width forever
  and rejecting a valid future one would be worse.

- page_size:

  How many files to request per underlying API call. Affects only how
  the manifest is fetched, never what you get back. `NULL` means mzLib's
  default.

- timeout:

  Seconds to allow for the whole fetch. `NULL`, the default, adds no
  limit of mzLibR's own, because mzLib already bounds every request: it
  retries a transient failure three times, after 5, 20 and 60 s, and
  each attempt can take up to 100 s. An outage that hangs is therefore
  reported as `mzlib_service_unavailable` after about 8 minutes per
  request. A number cuts the whole call off after that many seconds with
  `mzlib_timeout`, which is not `mzlib_service_unavailable`.

## What this manifest is, and is not

**This is what PRIDE's REST API publishes, which is not always
everything in the project.** On 2026-10-21, for PXD000001, the API
listed **9** of the **14** files in the FTP tree, and the five it omits
include the two largest, a 450 MB `.mzML` and the matching 472 MB
`.mzXML`, which are exactly the modern open-format conversions most
people want. (The example below replays a recording from 2026-10-03,
when it listed 8.) The omission is PRIDE's, not mzLib's.

So a manifest that looks short may be short. If completeness matters -
mirroring a project, budgeting a download, proving you analysed
everything - cross-check the FTP directory at
`https://ftp.pride.ebi.ac.uk/pride/data/archive/<year>/<month>/<accession>/`.

## Value

A data.frame with one row per file and the columns `file_name`,
`file_size_bytes`, `size_mb`, `extension`, `category`,
`category_accession`, `checksum`, `https_url`, `downloadable`,
`submission_date`, `publication_date`, `updated_date`,
`project_accession`, and a `locations` list column (see
[`pride_locations`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_locations.md)).

`file_size_bytes` is **the size PRIDE reports, which is not always what
you will transfer**: for compressed files it is frequently the
\*decompressed\* size. In PXD000001 the reported size of
`PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz` is 16,448,103 bytes and the
actual download is 5,984,662 - **2.75x** smaller. The `.mztab.gz`
behaves the same way; the `.xml.gz` does not. PRIDE's own metadata is
inconsistent here, so neither mzLibR nor mzLib can correct it. Treat the
sum as an upper bound on transfer, and see
[`pride_total_size_bytes`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md).

## Wraps

Wire verb `pride files`. Generated from the bridge's verb spec
`pride.files.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PrideArchiveClient.GetProjectFilesAsync`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideArchiveFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveFile.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveFile.cs` at mzLib
  `23c2490e`

- [`PrideArchiveExtensions.TotalSizeBytes`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs` at
  mzLib `23c2490e`

- [`PrideArchiveExtensions.TryGetHttpsDownloadUrl`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs` at
  mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `accession`:

  string; required. The project accession, e.g. PXD000001. Sent to PRIDE
  as given: the bridge checks only that it is not blank (the bindings
  upper-case it and check its shape first, because PRIDE is
  case-sensitive).

- `page_size` (wire `--page-size`):

  int; in **files**; default absent; range `>= 1`. Files requested per
  underlying API call. Changes how many requests the fetch takes, never
  what comes back: mzLib pages until the manifest is complete. Absent,
  the bridge sends none and GetProjectFilesAsync's default applies (100
  at d10017ce). Leave it absent; above PRIDE's server-side cap the
  manifest is complete only while PRIDE reports total_records.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `project_accession` (wire `accession`):

  string; never `NA`. The accession, echoed as given.

- `file_count`:

  int; in **files**; never `NA`. Entries in files. 0 for an unknown or
  private accession: PRIDE answers both with an empty manifest, not an
  error.

- `total_size_bytes`:

  int; in **bytes**; never `NA`. Sum of file_size_bytes over files
  (mzLib TotalSizeBytes). A sum of PRIDE-reported sizes over an
  incomplete manifest; see caveats.

- `files`:

  object\[\]; never `NA`. One entry per file in repository order; fields
  under result.columns.

## Columns

- `file_name`:

  string; never `NA`. The file's name as PRIDE lists it, e.g.
  TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw.

- `file_size_bytes`:

  int; in **bytes**; never `NA`. Size as PRIDE reports it, which is not
  the download length. For compressed files PRIDE often reports the
  DECOMPRESSED size (PXD000001's .mgf.gz: 16448103 reported vs 5984662
  transferred; its .xml.gz: 10,677,205 reported vs 10,668,000). pride
  checksums gives the exact length.

- `checksum`:

  string; never `NA`. The repository's checksum; the empty string when
  PRIDE provides none (never null). It is a SHA-1 where PRIDE gives one
  (40 hex digits; measured 2026-10-21 by hashing PXD010820's
  Maxquant_Output.zip), not the MD5 that pride checksums lists, and it
  is empty for every PXD000001 file.

- `category`:

  string; never `NA`. PRIDE file category value, e.g. RAW, PEAK, SEARCH,
  RESULT, OTHER.

- `category_accession`:

  string; never `NA`. The category's PRIDE CV accession, e.g.
  PRIDE:0000410.

- `https_url`:

  string; `NA` when the file has no location reachable over HTTPS (for
  example Aspera-only); pride download cannot fetch it. Direct HTTPS
  URL: the FTP location scheme-upgraded to HTTPS on PRIDE's FTP host.

- `locations`:

  object\[\]; never `NA`. Every published location as {accession, name,
  value} CV terms, e.g. {PRIDE:0000469, FTP Protocol, ftp://...}. May be
  empty.

- `submission_date`:

  string; never `NA`. ISO-8601 timestamp with offset (DateTimeOffset),
  e.g. 2012-03-13T00:00:00+00:00. Never null: a date PRIDE omitted
  arrives as 0001-01-01T00:00:00+00:00.

- `publication_date`:

  string; never `NA`. As submission_date.

- `updated_date`:

  string; never `NA`. As submission_date.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  accession missing or blank; page-size not an integer, given without a
  value, or \<= 0 (checked by the bridge before mzLib)

- `mzlib_service_unavailable` (service_unavailable):

  the bridge's own classification (Program.ClassifyError; since mzLib
  1.0.593 (#1350) PRIDE raises every transport failure as an
  HttpRequestException carrying the cause or the status, and messages
  name the file and host, never the URL): a timeout or cancellation, a
  socket failure, an HttpRequestException with no status but an inner
  cause (refused connection, DNS, TLS), an HTTP 408, 429 or 5xx, or a
  response body cut off in transit

- `mzlib_bridge_error` (correctness):

  any other HTTP status (e.g. 404), a paging contract violation
  (MzLibException: an identical page while total_records says more
  remain), the MaxPages backstop (an HttpRequestException mzLib composes
  with no status and no inner exception, which ClassifyHttpFailure
  deliberately does not excuse as an outage)

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

- This is PRIDE's REST manifest, which is knowingly incomplete: for
  PXD000001 on 2026-10-21 it listed 9 of the 14 files in the FTP tree,
  omitting five including the two largest (on 2026-10-03, when the
  fixtures were recorded, it listed 8). Use pride ftp-files when
  completeness or the true project size matters.

- An unknown accession and a private project both return file_count 0,
  not an error. The bindings raise their not-found error on an empty
  manifest; the wire does not.

- total_size_bytes is not the number of bytes a download transfers:
  PRIDE reports some compressed files at their decompressed size, and
  the sum covers only the files the manifest lists. The two errors run
  in opposite directions and do not cancel.

- Dates are timestamps with an offset here, unlike pride search, whose
  dates are bare calendar dates.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.pride.list_files`

- Rust (mzLibRust): `mzlib::pride::list_files_with` with `ListOptions`

- R (mzLibR): `pride_list_files`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- An absent PRIDE date crosses as 0001-01-01T00:00:00+00:00, not null
  (PrideArchiveFile's dates are non-nullable DateTimeOffset, and
  ToWireFile does not map default to null as pride search's ToWireDate
  does). Not observed on a real manifest; unverified whether PRIDE ever
  omits one.

## References

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive, the repository this verb lists (Perez-Riverol et al.,
  The PRIDE database at 20 years: 2025 update, NAR)

## See also

[`pride_download_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md),
which is usually what you want next.

## Examples

``` r

files <- pride_list_files("PXD000001")
files[, c("file_name", "category", "size_mb", "downloadable")]
#>                                                     file_name category
#> 1                  PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz    OTHER
#> 2                    PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz     PEAK
#> 3 TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML     PEAK
#> 4                                    erwinia_carotovora.fasta    OTHER
#> 5   TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw      RAW
#> 6                                       F063721.dat-mztab.txt    OTHER
#> 7                          PRIDE_Exp_Complete_Ac_22134.xml.gz   RESULT
#> 8                                                 F063721.dat   SEARCH
#>      size_mb downloadable
#> 1   0.497985         TRUE
#> 2  16.448103         TRUE
#> 3 243.031280         TRUE
#> 4   1.657668         TRUE
#> 5 220.475548         TRUE
#> 6   0.304798         TRUE
#> 7  10.677205         TRUE
#> 8  21.185462         TRUE
```
