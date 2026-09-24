# Sum the approximate sizes of some FTP files

The honest project-size estimate, and the counterpart to
[`pride_total_size_bytes`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md)
with the trade-offs reversed. It sums over the **complete** FTP listing
([`pride_list_ftp_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)),
so no files are missing - but each size is PRIDE's rounded
directory-index value, so the total is an **estimate**, not an exact
byte count. For PXD000001 it lands near the true 1.44 GB, where
[`pride_total_size_bytes`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md)
reports 0.51 GB over the incomplete REST manifest. For the exact bytes
of one file, issue an HTTP HEAD against its `url`.

## Usage

``` r
pride_approximate_total_size_bytes(files)
```

## Arguments

- files:

  A
  [`pride_list_ftp_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
  data.frame, or a subset.

## Value

The summed approximate size in bytes, as a double.

## See also

[`pride_total_size_bytes`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md),
the REST-manifest counterpart.

## Examples

``` r

listing <- pride_list_ftp_files("PXD000001")
pride_approximate_total_size_bytes(listing)
#> [1] 670905958
```
