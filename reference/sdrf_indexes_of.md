# Every position a column name occupies in an SDRF document

Every position a column name occupies in an SDRF document

## Usage

``` r
sdrf_indexes_of(doc, column)
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

The 1-based positions carrying this name, in document order; empty when
the column is absent.

## See also

[`sdrf_index_of`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_index_of.md)

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
sdrf_indexes_of(doc, "comment[modification parameters]")
#> [1] 20 21 22 23 24 25 26 27
```
