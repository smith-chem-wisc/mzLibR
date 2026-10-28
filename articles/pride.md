# Finding and downloading PRIDE Archive data

The `pride_` functions find PRIDE Archive projects, list their files and
download the ones you choose, using mzLib’s paging and URL resolution.

| question | function | mzLib |
|----|----|----|
| which projects are about my subject? | [`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md) | `PrideArchiveClient.SearchProjectsAsync` |
| what files does a project have, with categories and checksums? | [`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md) | `PrideArchiveClient.GetProjectFilesAsync` |
| every file, including ones the manifest omits | [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md) | `PrideArchiveClient.GetProjectFilesFromFtpAsync` |
| how big is it, before I download? | [`pride_approximate_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_approximate_total_size_bytes.md), [`pride_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md) | `PrideArchiveExtensions.TotalApproximateSizeBytes` |
| download the files I picked | [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md) | `PrideArchiveClient.DownloadFileAsync` |
| download by category or extension | [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md) | `PrideArchiveClient.DownloadProjectFilesAsync` |
| is my download the file PRIDE holds? | [`pride_file_checksums()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_file_checksums.md) | `PrideArchiveClient.GetFileChecksumsAsync` |

## Find a project

[`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md)
is the discovery step: it turns a subject into the accessions every
other function takes, with every page fetched and no accession repeated.

``` r

hits <- pride_search("plasmodium falciparum schizont")
hits[, c("accession", "submission_type", "publication_date")]
#>   accession submission_type publication_date
#> 1 PXD070842        COMPLETE       2026-06-01
#> 2 PXD020210        COMPLETE       2021-01-25
#> 3 PXD020189        COMPLETE       2021-01-25
#> 4 PXD008250         PARTIAL       2018-02-13
#> 5 PXD001684         PARTIAL       2015-04-23
#> 6 PXD000070        COMPLETE       2014-04-24
hits$organisms[[1]]
#> [1] "Homo sapiens (human)"                "Plasmodium falciparum (isolate 3d7)"
hits$matched_fields[[1]]
#> [1] "references" "title"
```

A hit is a search-index projection, not the project’s metadata:
vocabulary fields arrive as display strings with no accessions, a count
of 0 means “not reported”, and `project_file_names` is not the manifest.
Follow the `accession`.

## List, then choose with `[`

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

Choose rows with ordinary subsetting, then pass them to
[`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md).
That says things the `category` and `extensions` filters of
[`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
cannot, such as “under 5 MB” or “everything except the raw file”.

``` r

small <- files[files$size_mb < 5 & files$downloadable, ]
small$file_name
#> [1] "PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz"
#> [2] "erwinia_carotovora.fasta"                  
#> [3] "F063721.dat-mztab.txt"
pride_total_size_bytes(small)
#> [1] 2460451
```

``` r

# Not run: downloads from PRIDE.
pride_download_files(small, "downloads")
```

## Two things the manifest will not tell you

**Sizes can be the decompressed size.** For some compressed files PRIDE
reports the size of the content, not of the download: the `.mgf.gz`
above is listed at 16 MB and transfers 6 MB. Treat
[`pride_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md)
as an upper bound.

**The manifest can be incomplete.** PRIDE’s REST API omits some of
PXD000001’s files, including the two largest.
[`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
walks the project’s FTP tree instead, subdirectories included:

``` r

listing <- pride_list_ftp_files("PXD000001")
c(ftp_tree = nrow(listing), rest_manifest = nrow(files))
#>      ftp_tree rest_manifest 
#>            14             8
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
round(pride_approximate_total_size_bytes(listing) / 1e9, 2)   # GB
#> [1] 1.44
```

A file that appears only there is fetched from its `url` with any HTTPS
client, for example
[`download.file()`](https://rdrr.io/r/utils/download.file.html).

## Filters that match nothing are errors

[`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
with a `category` or `extensions` filter that matches nothing raises,
rather than reporting success with nothing downloaded. A compressed
file’s extension is `.gz`, so `extensions = ".mgf"` matches nothing in
PXD000001:

``` r

# Not run: downloads from PRIDE.
pride_download("PXD000001", "downloads", category = "PEAK", extensions = ".gz")
```

An unknown accession raises `mzlib_project_not_found` rather than
returning nothing:

``` r

tryCatch(pride_list_files("PXD999999999"),
         mzlib_project_not_found = function(e) conditionMessage(e))
#> [1] "PRIDE returned no files for 'PXD999999999'. Either the accession does not exist (check for a typo) or the project is private. PRIDE does not distinguish the two, so neither can mzLibR."
```

An EBI outage raises `mzlib_service_unavailable`. See
[`?mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md),
and the next section for what to do about one.

## Retries, timeouts and outages

**mzLib retries for you.** A PRIDE request that fails in a way that can
pass (no answer, an HTTP 408, 429 or 5xx, or a 403 from EBI’s FTP host,
its rate limit) is retried three times, after 5, 20 and 60 seconds, and
that request alone (mzLib \#1412). A paged manifest retries only the
page that failed, and a download resumes where it stopped.

**mzLibR adds no timeout of its own by default**, on any `pride_`
function. Each attempt is already limited to 100 seconds, so an outage
that hangs is reported as `mzlib_service_unavailable` after about 8
minutes per request. Pass `timeout` to cut a call off sooner; it then
raises `mzlib_timeout`, which is not `mzlib_service_unavailable`:

``` r

# Not run: contacts PRIDE.
pride_list_files("PXD000001", timeout = 60)   # give up after a minute
```

By the time `mzlib_service_unavailable` reaches you, mzLib has already
retried, so retrying again at once rarely helps. Wait minutes, not
seconds:

``` r

# Not run: shows how to handle a PRIDE outage, which a recording cannot produce.
files <- tryCatch(
  pride_list_files("PXD000001"),
  mzlib_service_unavailable = function(e) {
    Sys.sleep(600)   # mzLib retried already; give EBI ten minutes
    pride_list_files("PXD000001")
  }
)
```

Two failures look like a request PRIDE refused, `mzlib_bridge_error`
with `error_type` `"HttpRequestException"`, but are EBI’s: an FTP
listing that came back empty, and a 403 from EBI’s FTP host, each still
failing after mzLib’s retries. Treat those as an outage.

## Interrupted and suspect downloads

Each file streams to a temporary `.partial` name and moves into place
only once complete, so a cancelled or crashed download leaves no
truncated file under its real name. A transfer that fails after mzLib’s
retries leaves `<name>.partial` and `<name>.partial.validator` in the
folder on purpose: the next call for that file resumes from them (mzLib
\#1412). Leave them in place, or delete both to start over.

mzLib can check a download against PRIDE’s checksums, and catch an EBI
error page served at the right size (mzLib \#1427, \#1449). That check
is not yet on the path
[`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
and
[`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
use, so a file can be reported as downloaded when its bytes are wrong.
For files that matter, check each one against
[`pride_file_checksums()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_file_checksums.md)
before you search it, as the next section shows.

## Checking a download against PRIDE’s checksums

**Is the file on your disk the file PRIDE holds?** The download
functions cannot tell you yet. PRIDE records an MD5 and the exact size
of every file a submitter deposits, and
[`pride_file_checksums()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_file_checksums.md)
returns that record. Its row names are the file names, so look a file up
by the name it was downloaded under:

``` r

checksums <- pride_file_checksums("PXD000001")
nrow(checksums)
#> [1] 12
checksums["TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw", c("md5", "size_bytes")]
#>                                                                                        md5
#> TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw c567cecbb9866796e690a4958746f7b6
#>                                                           size_bytes
#> TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw  220475548
```

`md5` is 32 hexadecimal characters, as PRIDE wrote them; compare it
without regard to case. `size_bytes` is the exact length a download of
the file comes to, which neither other listing gives you.
[`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
rounds its sizes, and
[`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
can be wrong:

``` r

gz <- "PRIDE_Exp_Complete_Ac_22134.xml.gz"
c(manifest = files$file_size_bytes[files$file_name == gz], checksums = checksums[gz, "size_bytes"])
#>  manifest checksums 
#>  10677205  10668000
```

The download is 10,668,000 bytes, the checksum list’s number. Because
each size is exact, so is their sum: what the deposited files weigh,
without the `generated/` files, which have no row (below).

``` r

round(sum(checksums$size_bytes) / 1e9, 2)   # GB
#> [1] 1.43
```

### The check

A file is PRIDE’s when its length and its MD5 both match. The length
costs nothing to check. The MD5 reads the whole file, a few seconds per
gigabyte. Check both, because the length alone does not prove the bytes:
on 7 and 8 October 2026, EBI served about 23 of about 70 large transfers
at the right size with its own error text written over part of them
(mzLib \#1449). Base R has all it takes:

``` r

check_download <- function(path, expected) {
  # NA when the file is the one PRIDE recorded, otherwise what is wrong with it.
  size <- file.size(path)
  if (size != expected$size_bytes) {
    return(paste(size, "bytes where PRIDE records", expected$size_bytes))
  }
  md5 <- unname(tools::md5sum(path))
  if (md5 != tolower(expected$md5)) {
    return(paste("MD5", md5, "where PRIDE records", expected$md5))
  }
  NA_character_
}
```

Here it is on two damaged copies of the project’s 1,645-byte
`README.txt`. The first was cut short. The second is the right size with
EBI’s error text inside, and only the MD5 catches it:

``` r

readme <- checksums["README.txt", ]
damaged <- tempfile(fileext = ".txt")
writeBin(charToRaw(strrep("x", 1000)), damaged)
check_download(damaged, readme)
#> [1] "1000 bytes where PRIDE records 1645"
error_page <- "<Error><Code>ConnectionClosedException</Code></Error>"
writeBin(charToRaw(paste0(error_page, strrep(".", 1645 - nchar(error_page)))), damaged)
check_download(damaged, readme)
#> [1] "MD5 6a8f1a7d39089d34bf9bc9b19e553dd4 where PRIDE records 492753137013825025e6729c8798664f"
```

After a real download, check every file it wrote:

``` r

# Not run: downloads files from PRIDE.
paths <- pride_download("PXD000001", "downloads", category = "RAW")
checksums <- pride_file_checksums("PXD000001")
for (path in paths) {
  name <- basename(path)
  if (!name %in% rownames(checksums)) {
    message(name, ": PRIDE has no checksum for it, so it cannot be checked")
    next
  }
  problem <- check_download(path, checksums[name, ])
  if (!is.na(problem)) message(name, ": ", problem, ". Download it again.")
}
```

To download a damaged file again, pass its row to
[`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md),
or call
[`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
with `overwrite = TRUE`, the default. `overwrite = FALSE` would skip it,
because it is already there.

### What the list covers

It covers the files the submitter deposited at the project’s root. It is
a third list, not a view of either of the others, and three things
follow from that.

**PRIDE’s own `generated/` files have no row.** A file without a row
cannot be checked; that does not make it invalid. In PXD000001 those are
the two files under `generated/`, which
[`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
names without their folder:

``` r

listing$relative_path[!listing$relative_path %in% checksums$file_name]
#> [1] "generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz"  
#> [2] "generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz"
files$file_name[!files$file_name %in% checksums$file_name]
#> [1] "PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz"
#> [2] "PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz"
```

**It names files the REST manifest omits,** such as the 450 MB `.mzML`
that only the FTP listing showed above:

``` r

"TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML" %in% checksums$file_name
#> [1] TRUE
```

**It can name files that are not deposits at all.** PXD015239’s list
includes a 2.3 GB NFS temp file, `.nfs80360b91005fae5300002478`. Look up
the files you have by name; do not read the list as the project’s
contents.

### Two traps

**The manifest’s `checksum` column is not this MD5.** Where PRIDE fills
it in, it is a 40-character SHA-1, a different digest of the same bytes.
For every PXD000001 file it is empty:

``` r

unique(files$checksum)
#> [1] ""
```

**Zero rows is an answer.** It means the project exists but PRIDE holds
no checksums for it, so nothing can be checked. An accession PRIDE does
not know raises `mzlib_project_not_found`, and so does a project not yet
published: PRIDE holds no checksums before publication. A list mzLib
cannot read (a new header, a malformed row, a name listed twice) is an
`mzlib_bridge_error` with `error_type` `"MzLibException"`: not your bug,
and mzLib does not retry it, so report it.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive, whose checksum list this verb reads (Perez-Riverol et
  al., The PRIDE database at 20 years: 2025 update, NAR).
