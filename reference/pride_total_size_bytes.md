# Total size of a set of PRIDE files

Total size of a set of PRIDE files

## Usage

``` r
pride_total_size_bytes(files)
```

## Arguments

- files:

  A
  [`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  data.frame, or a subset.

## Value

A single number of bytes.

## Two reasons this is not the size of the project

**It over-reports compressed files.** PRIDE frequently gives the
\*decompressed\* size - the MGF in PXD000001 reports 16,448,103 bytes
and downloads as 5,984,662, a factor of **2.75**.

**It sums an incomplete manifest.** For PXD000001 this returns **0.51
GB**; the project on disk is **1.44 GB**, because PRIDE's API omits five
files including the two largest (see
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)).

The two errors run in opposite directions and do **not** cancel:
compressed sizes are inflated, whole files are missing entirely.

## Examples

``` r

files <- pride_list_files("PXD000001")
pride_total_size_bytes(files)   # an upper bound on the transfer; see Details
#> [1] 514278049
```
