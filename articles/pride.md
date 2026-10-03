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

An EBI outage raises `mzlib_service_unavailable`, which is the one to
retry. See
[`?mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive, the repository files are fetched from (Perez-Riverol et
  al., The PRIDE database at 20 years: 2025 update, NAR).
