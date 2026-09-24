# Where a column first sits in an SDRF document

Comparison is exact and case-**sensitive**, matching the SDRF
specification and mzLib.

## Usage

``` r
sdrf_index_of(doc, column)
```

## Arguments

- doc:

  An
  [`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  or
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  result.

- column:

  A column name.

## Value

The 1-based position of the first column with this name, or `NA`.

## See also

[`sdrf_indexes_of`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_indexes_of.md)

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
sdrf_index_of(doc, "source name")
#> [1] 1
sdrf_index_of(doc, "no such column")
#> [1] NA
```
