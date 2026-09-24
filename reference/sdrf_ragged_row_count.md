# How many rows are shorter than the header

How many rows are shorter than the header

## Usage

``` r
sdrf_ragged_row_count(doc)
```

## Arguments

- doc:

  An
  [`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  or
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  result.

## Value

The number of returned rows carrying fewer cells than there are columns.

## Examples

``` r

ragged <- sdrf_read("PXD059974.sdrf.tsv")
sdrf_ragged_row_count(ragged)
#> [1] 17
```
