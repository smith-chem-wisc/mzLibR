# The first cell under a column, one per row

`NA` means **the document does not have this column, or the row is too
short to reach it**. It does not mean "empty": the SDRF reserved words
`"not available"` and `"not applicable"` are real values an experiment
chose to write, and they come back as themselves.

## Usage

``` r
sdrf_value(doc, column)
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

A character vector with one element per returned row.

## See also

[`sdrf_all`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md)
for a column that repeats.

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
sdrf_value(doc, "characteristics[organism part]")
#> [1] "human erythrocytes" "human erythrocytes" "human erythrocytes"
#> [4] "human erythrocytes" "human erythrocytes" "human erythrocytes"
```
