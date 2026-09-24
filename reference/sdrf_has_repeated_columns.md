# Whether an SDRF document repeats a column name

Whether an SDRF document repeats a column name

## Usage

``` r
sdrf_has_repeated_columns(doc)
```

## Arguments

- doc:

  An
  [`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  or
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  result.

## Value

`TRUE` or `FALSE`. When `TRUE`,
[`sdrf_records`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
loses cells.

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
sdrf_has_repeated_columns(doc)
#> [1] TRUE
```
