# Download exactly the files you selected

Download a PRIDE project's files, all of them or a selection, into a
directory and report the paths written.

The counterpart to
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
and usually the one you want: filter the data.frame however you like,
then pass the rows. As with
[`pride_download`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md),
the bytes are not checked against PRIDE's record; compare a file that
matters with
[`pride_file_checksums`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_file_checksums.md).

## Usage

``` r
pride_download_files(files, destination, overwrite = TRUE, timeout = NULL)
```

## Arguments

- files:

  Rows of a
  [`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  data.frame. They must all come from one project. Rows whose
  `downloadable` is `FALSE` are refused up front rather than failing
  halfway through a multi-gigabyte transfer.

- destination:

  Directory to write into.

- overwrite:

  When `FALSE`, files already present are left alone.

- timeout:

  Seconds to allow, or `NULL` to wait as long as it takes.

## Details

“`r files <- pride_list_files("PXD000001") small <- files[files$size_mb < 5 & files$downloadable, ] pride_download_files(small, "downloads") `“

## Value

A character vector of paths. These say **where each file is**, not what
was transferred just now: with `overwrite = FALSE` a file already
present is left alone and its path is still returned. Do not read
[`length()`](https://rdrr.io/r/base/length.html) of this as work done.

## Wraps

Wire verb `pride download`. Generated from the bridge's verb spec
`pride.download.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PrideArchiveClient.DownloadProjectFilesAsync`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideArchiveClient.DownloadFileAsync`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideArchiveExtensions.WhereCategory`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs` at
  mzLib `23c2490e`

- [`PrideArchiveExtensions.WhereExtension`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs` at
  mzLib `23c2490e`

- [`PrideArchiveExtensions.GetHttpsDownloadUrl`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveExtensions.cs` at
  mzLib `23c2490e`

## Parameters: units, ranges and defaults

- wire `--accession`:

  string; required. The project accession, e.g. PXD000001. The bridge
  checks only that it is not blank. *Not an argument here: taken from
  the rows' project_accession.*

- `destination` (wire `--dest`):

  path; required. Directory to write into; created if absent, but only
  when a file is actually fetched (mzLib creates it after the
  no-overwrite skip). With zero matches, or every file skipped, no
  directory is created, yet destination_directory is still reported.
  Blank is a usage error (it never falls back to the working directory).

- wire `--category`:

  string; default absent; range
  `a PRIDE category value, case-insensitive: RAW, PEAK, SEARCH, RESULT, OTHER, ...`.
  Keep only files of this category (mzLib WhereCategory). Absent keeps
  every category. Given but blank is a usage error, never 'everything'.
  ANDed with ext. *Not an argument here: a selection is already
  filtered, by `[`.*

- wire `--ext`:

  string\[\]; default absent; range
  `comma-separated; a leading dot is optional`. Keep only files whose
  name ends with one of these extensions, case-insensitive (mzLib
  WhereExtension), e.g. ".raw,.mzML". Absent keeps every type. Given but
  naming no extension (blank or only commas) is a usage error. A
  compressed file's extension is .gz. *Not an argument here: a selection
  is already filtered, by `[`.*

- `overwrite` (wire `--no-overwrite`):

  flag; default `FALSE`. Leave a file already present at the destination
  untouched and make no request for it: a cheap resume. Absent, an
  existing file is replaced.

- wire `--names-from-stdin`:

  flag; default `FALSE`. Select exact files instead of filtering: stdin
  carries file names, one per line (newline-delimited; blank lines
  ignored; each line trimmed; exact, case-sensitive match against the
  manifest's file_name). Exclusive with category and ext. Stdin rather
  than argv because argv caps near 32 KB and any other separator could
  occur in a file name. A name containing a line break cannot be
  selected (the bindings refuse one; PRIDE has never published one).
  *Not an argument here: always set: the selection travels on stdin.*

- `files` (wire `--stdin`):

  string\[\]; default absent; range
  `>= 1 non-blank line when names-from-stdin is given`. Read only with
  names-from-stdin: the file names to fetch, one per line. Every name
  must be in the project's REST manifest (pride files). A UTF-8 BOM on
  the first line is stripped (ReadStdinLines) from pyMzLib \#67 on; the
  0.2.0 bridge keeps it, so that first name fails as not in the project.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `accession`:

  string; never `NA`. The accession, echoed as given.

- `destination_directory`:

  string; never `NA`. dest, made absolute.

- `downloaded_count`:

  int; in **files**; never `NA`. Entries in paths: files written PLUS,
  under no-overwrite, files already present and skipped.

- `paths`:

  string\[\]; never `NA`. One path per selected file, in manifest order:
  mzLib's Path.Combine(dest as given, file_name), so relative when dest
  was relative (resolve against destination_directory's parent, not this
  list).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  accession or dest missing or blank; category given but blank; ext
  given but naming no extension; names-from-stdin together with category
  or ext; names-from-stdin with no non-blank line on stdin

- `mzlib_usage_error` (usage):

  names-from-stdin naming a file the project's manifest does not contain
  (the message lists up to five). Raised AFTER the matching files have
  been downloaded: they are on disk when the error arrives.

- `mzlib_service_unavailable` (service_unavailable):

  the bridge's own classification (Program.ClassifyError; since mzLib
  1.0.593 (#1350) PRIDE raises every transport failure as an
  HttpRequestException carrying the cause or the status, and messages
  name the file and host, never the URL): a timeout, a socket or TLS
  failure, an HTTP 408, 429 or 5xx on the manifest or a file, a body
  that delivered nothing for mzLib's BodyStallTimeout (2 min), or a body
  cut off in transit

- `mzlib_bridge_error` (correctness):

  a selected file has no HTTPS location (NotSupportedException;
  Aspera-only); a manifest file name that is not a bare leaf name
  (ArgumentException: refusing to write outside dest); any other HTTP
  status; a local disk failure while writing (a bare IOException, e.g. a
  full disk, is deliberately NOT treated as an outage); a 403 from
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

- Selects from PRIDE's REST manifest (pride files), which omits some
  files; a file listed only by pride ftp-files cannot be downloaded by
  this verb.

- Files are fetched one at a time, in manifest order. Each streams to
  '\<name\>.partial' and is moved into place only when complete, so an
  interrupted transfer never leaves a truncated file under its final
  name; files finished before a failure stay on disk. Since mzLib
  1.0.594 (#1412), a transfer that fails after its retries leaves
  '\<name\>.partial' and '\<name\>.partial.validator' in dest on
  purpose: the next call resumes from them. dest can therefore hold
  these pairs after a failed call.

- No size, MD5 or EBI error-page check reaches this verb. mzLib 1.0.595
  verifies a download against PRIDE's checksums and the announced size,
  and scans a small body for EBI's error text (#1427, \#1449: on
  2026-10-07 and 08 about 23 of about 70 large transfers came back at
  the right size holding an error page), but only on the
  DownloadFileAsync overload that takes a PrideFileChecksum.
  DownloadProjectFilesAsync, which this verb calls, uses the overload
  without one, so such a file is reported as downloaded. Verify large
  files yourself, against the size and MD5 pride checksums lists, until
  the bridge can ask for the check.

- A filter that matches nothing, or an unknown accession, succeeds with
  downloaded_count 0 on the wire. The bindings raise a usage error when
  a category or ext filter matched nothing.

- Downloads can be very large (PXD000001's single .raw is about 220 MB);
  list first and filter or select.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.pride.download`

- Rust (mzLibRust): `mzlib::pride::download` with `DownloadOptions`

- R (mzLibR): `pride_download`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- No recorded fixture, so no checked example: a replay needs a stub
  bridge that also fakes the written files. The envelope keys above are
  read from the object PrideDownloadAsync returns, not from a recording.

- The unmatched-selection check runs after the download, so a typo in
  one name still downloads every other selected file before failing.
  Checking the manifest first would need the bridge to fetch it
  separately from DownloadProjectFilesAsync.

- A zero-match filter is a usage error in the bindings but a success on
  the wire; whether the wire should raise it (as it does for an
  unmatched names-from-stdin selection) is undecided.

## References

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive, the repository files are fetched from (Perez-Riverol et
  al., The PRIDE database at 20 years: 2025 update, NAR)

## Examples

``` r

files <- pride_list_files("PXD000001")
small <- files[files$size_mb < 1 & files$downloadable, ]
small[, c("file_name", "size_mb")]
#>                                    file_name  size_mb
#> 1 PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz 0.497985
#> 6                      F063721.dat-mztab.txt 0.304798
# \donttest{
# Downloads from EBI, so it needs a bridge and the network.
if (nzchar(Sys.getenv("MZLIB_BRIDGE"))) pride_download_files(small, tempdir())
# }
```
