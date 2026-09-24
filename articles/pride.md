# Finding and downloading PRIDE Archive data

The `pride_` functions list a PRIDE Archive project’s files and download
the ones you choose, using mzLib’s paging and URL resolution.

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

pride_download_files(small, "downloads")   # needs the network
```

## Two things the manifest will not tell you

**Sizes can be the decompressed size.** For some compressed files PRIDE
reports the size of the content, not of the download: the `.mgf.gz`
above is listed at 16 MB and transfers 6 MB. Treat
[`pride_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md)
as an upper bound.

**The manifest can be incomplete.** PRIDE’s REST API publishes 8 files
for PXD000001; the FTP tree holds 13, including the two largest.
[`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
walks the tree instead:

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

pride_download("PXD000001", "downloads", category = "PEAK", extensions = ".gz")
```

An unknown accession raises `mzlib_project_not_found`, and an EBI outage
raises `mzlib_service_unavailable`, which is the one to retry. See
[`?mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
