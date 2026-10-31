# The MD5 and exact size PRIDE recorded for each file deposited in a project

List the MD5 and exact size PRIDE recorded for each file deposited in a
project: the record a downloaded file is checked against.

Answers **is the file on my disk the file PRIDE holds?** A file is
PRIDE's when its length is `size_bytes` and its MD5 is `md5`.
[`pride_download`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
and
[`pride_download_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
check neither yet, so for files that matter, check them yourself:
[`file.size()`](https://rdrr.io/r/base/file.info.html) and
[`tools::md5sum()`](https://rdrr.io/r/tools/md5sum.html) are all it
takes, as the example shows.

## Usage

``` r
pride_file_checksums(accession, timeout = NULL)
```

## Arguments

- accession:

  A project accession, e.g. `"PXD000001"`. Upper-cased and trimmed
  first, because PRIDE's checksum route is case-sensitive.

- timeout:

  Seconds to allow for the whole call. `NULL`, the default, adds no
  limit of mzLibR's own, because mzLib already bounds every request: it
  retries a transient failure three times, after 5, 20 and 60 s, and
  each attempt can take up to 100 s. An outage that hangs is therefore
  reported as `mzlib_service_unavailable` after about 8 minutes per
  request, and an empty list costs a second request, the project lookup.
  A number cuts the whole call off after that many seconds with
  `mzlib_timeout`, which is not `mzlib_service_unavailable`.

## Details

**This is a third list, not a view of either of the others.** It covers
the files the submitter deposited at the project's root. PRIDE's own
`generated/` files have no row: a file without a row cannot be checked,
and that does not make it invalid. It names files
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
omits (5 of PXD000001's 12 rows on 2026-10-21). And it can name files
that are not deposits at all: PXD015239's lists an NFS temp file,
`.nfs80360b91005fae5300002478`. Look up the files you have by name; do
not read the list as the project's contents.

The `checksum` column of
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
is not this MD5. Where PRIDE fills it in, it is a 40-character SHA-1,
and for every PXD000001 file it is empty.

## Value

A data.frame with one row per file, in file-name order, and the columns
`file_name`, `md5`, `size_bytes` and `project_accession`. **Its row
names are the file names**, so `checksums["README.txt", ]` looks a file
up by the name it was downloaded under; they are case-sensitive, as
PRIDE's names are.

`file_name` is the bare file name, spelled as
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
spells it, never a path. `md5` is 32 hexadecimal characters as PRIDE
wrote them: lower case in every list seen, but not rewritten, so compare
without regard to case
([`tools::md5sum()`](https://rdrr.io/r/tools/md5sum.html) is lower
case). `size_bytes` is the file's exact size in bytes, the length a
download of it comes to, as a double (raw files pass 2 GiB). Compare a
downloaded file's [`file.size()`](https://rdrr.io/r/base/file.info.html)
with it, not with `file_size_bytes` from
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
which can differ: for PXD000001's `PRIDE_Exp_Complete_Ac_22134.xml.gz`
that one is 10,677,205 bytes and this one is 10,668,000, the length the
download has. `project_accession` is the accession the list was fetched
for.

**Zero rows when the project exists but PRIDE holds no checksums for
it**: that is an answer, not an error, though nothing can then be
checked. An unknown accession raises `mzlib_project_not_found`, and so
does a project not yet published, since PRIDE holds no checksums before
publication. A list mzLib cannot read (a new header, a malformed row, a
name listed twice) stays an `mzlib_bridge_error` with `error_type`
`"MzLibException"`: the project exists and its list is broken, so "check
for a typo" would send you looking for the wrong mistake.

## Checking a download

Check both the length and the MD5. The length costs nothing; the MD5
reads the whole file. The length alone does not prove the bytes: on 7
and 8 October 2026, EBI served about 23 of about 70 large transfers at
the right size with its own error text written over part of them (mzLib
\#1449). The PRIDE article shows the check on two damaged copies of a
file.

## Wraps

Wire verb `pride checksums`. Generated from the bridge's verb spec
`pride.checksums.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PrideArchiveClient.GetFileChecksumsAsync`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `d10017ce`

- [`PrideFileChecksum`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/PrideFileChecksum.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideFileChecksum.cs` at mzLib
  `d10017ce`

- [`PrideArchiveClient.GetProjectAsync`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `d10017ce`

## Parameters: units, ranges and defaults

- `accession`:

  string; required. The project accession, e.g. PXD000001. Sent to PRIDE
  as given: the bridge checks only that it is not blank. PRIDE's
  checksum and project routes are case-sensitive (pxd000001 is an
  unknown accession there), so the bindings upper-case it and check its
  shape first.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `project_accession` (wire `accession`):

  string; never `NA`. The accession, echoed as given.

- `file_count`:

  int; in **files**; never `NA`. Entries in files. 0 means the project
  exists and PRIDE holds no checksum for any of its files, so nothing
  can be checked. It is an answer, not a not-found: PRIDE answers an
  unknown accession with an empty list too, so mzLib then asks for the
  project and fails when there is none (see errors).

- `files`:

  object\[\]; never `NA`. One entry per file, sorted by file_name in
  ordinal order. mzLib returns a map keyed by file name, which has no
  order, so the bridge sorts it to send the same list the same way on
  every run. Fields under result.columns.

## Columns

- `file_name`:

  string; never `NA`. The bare file name, as pride files spells
  file_name. Never a path. Unique in the list and case-sensitive:
  README.txt and readme.txt would be two files.

- `md5`:

  string; never `NA`. The file's MD5 as 32 hexadecimal characters, as
  PRIDE wrote it. Lower case in every list seen, but not rewritten, so
  compare it without regard to case, as mzLib's own checks do. Not pride
  files' checksum, which is a different digest (see caveats).

- `size_bytes`:

  int; in **bytes**; never `NA`. The file's exact size in bytes: the
  length a download of it comes to. A 64-bit count; raw files pass 2
  GiB. Compare a downloaded file's length with this, not with pride
  files' file_size_bytes or pride ftp-files' approximate_size_bytes (see
  caveats).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  accession missing, blank, or given without a value (checked by the
  bridge before mzLib)

- `mzlib_service_unavailable` (service_unavailable):

  the bridge's own classification (Program.ClassifyError; since mzLib
  1.0.593 (#1350) PRIDE raises every transport failure as an
  HttpRequestException carrying the cause or the status): a timeout or
  cancellation, a socket failure, an HttpRequestException with no status
  but an inner cause (refused connection, DNS, TLS), an HTTP 408, 429 or
  5xx on the checksum list or on the project lookup, or a body cut off
  in transit

- `mzlib_bridge_error` (correctness):

  no project has that accession, or it is not yet public (MzLibException
  'PRIDE Archive has no project with accession ...': PRIDE answers both
  with an empty list and 404 on the project; the bindings re-map this
  one message to their not-found error); a list mzLib cannot read
  (MzLibException, never retried: a header other than File-Name,
  File-MD5Checksum and File-Size, an empty line before the last row, a
  row without exactly three tab-separated fields or without a name, an
  MD5 that is not 32 hexadecimal characters, a size that is not a whole
  number of bytes, or a name listed twice); any other HTTP status on
  either request (e.g. 404 on the list)

## Caveats

- Since mzLib 1.0.594 (#1412), PrideArchiveClient retries a transient
  failure, that request alone, after 5, 20 and 60 s (a 429 or 503
  Retry-After of up to 60 s replaces the wait); transient means no HTTP
  status at all, 408, 429, any 5xx, or a 403 from ftp.pride.ebi.ac.uk,
  EBI's rate limit. Each attempt can take up to the client's 100 s
  timeout, so a service that hangs takes up to about 485 s per request
  before the failure is reported. An empty list costs a second request,
  the project lookup, with its own retries. A binding's own timeout
  shorter than that reports a hanging outage as a timeout, not as
  service_unavailable.

- A third list, not a view of pride files or pride ftp-files. It covers
  the files the submitter deposited at the project's root. PRIDE's
  generated/ files have no row, so a file without a row cannot be
  checked; that does not make it invalid. For PXD000001 on 2026-10-21 it
  had 12 rows for the 14 files of the FTP tree (the two missing are in
  generated/), and 5 of the 12 named files the REST manifest omits.

- It can list files that are not deposits: PXD015239's includes an NFS
  temp file, .nfs80360b91005fae5300002478, at 2285905384 bytes. Look up
  the files you have by name; do not read the list as the project's
  contents.

- PRIDE holds no checksums before publication: a project not yet public
  is the not-found error here, not an empty list.

- md5 is not pride files' checksum. That field, where PRIDE fills it in,
  is a 40-character SHA-1 (PXD010820's Maxquant_Output.zip, checked
  against the downloaded file on 2026-10-21), and for every PXD000001
  file it is empty. Do not compare the two.

- size_bytes is the exact download length, which neither other listing
  gives. pride files' file_size_bytes can differ: for PXD000001's
  PRIDE_Exp_Complete_Ac_22134.xml.gz it is 10677205 where the download
  is 10668000 bytes (this list's value, and an HTTP HEAD's
  Content-Length on 2026-10-21), and for a .gz file it can be more than
  ten times larger (PXD015239's 190712-VL-BR-LFQ-1.mzid.gz, per mzLib's
  PrideFileChecksum). pride ftp-files' sizes are rounded.

- pride download does not use this list. It calls
  DownloadProjectFilesAsync, which takes the DownloadFileAsync overload
  without a checksum, so no size, MD5 or EBI error-page check reaches it
  (pride thread, PRIDE-B5). Compare a downloaded file's size and MD5
  with its row yourself.

- No total size is sent. mzLib has none over this list, and a sum made
  in the bridge would be the bridge's number, not mzLib's; sum
  size_bytes yourself.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.pride.file_checksums`

- Rust (mzLibRust): `mzlib::pride::file_checksums_with`

- R (mzLibR): `pride_file_checksums`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- mzLib checks that the project exists only when the list's body is
  empty. A body that is only the header returns 0 rows without that
  check, so a header-only answer for an unknown accession would read as
  a project with no checksums. PRIDE has not been seen to send one: it
  answers an unknown accession with an empty body (PXD999999999,
  2026-10-21).

- No public project with an empty list has been found, so the 0-row
  answer is pinned by a stubbed test only, not by a recorded fixture.

## References

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive, whose checksum list this verb reads (Perez-Riverol et
  al., The PRIDE database at 20 years: 2025 update, NAR)

## See also

[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
[`pride_download_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)

## Examples

``` r

checksums <- pride_file_checksums("PXD000001")
nrow(checksums)
#> [1] 12
checksums["README.txt", c("md5", "size_bytes")]
#>                                         md5 size_bytes
#> README.txt 492753137013825025e6729c8798664f       1645

# A downloaded file is PRIDE's when both match. Here, a copy cut short:
copy <- tempfile(fileext = ".txt")
writeBin(charToRaw(strrep("x", 1000)), copy)
file.size(copy) == checksums["README.txt", "size_bytes"]
#> [1] FALSE
tolower(unname(tools::md5sum(copy))) == tolower(checksums["README.txt", "md5"])
#> [1] FALSE
unlink(copy)
```
