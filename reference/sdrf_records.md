# An SDRF document as a data.frame, when its shape allows one

The convenient shape, and a **lossy** one when a column name repeats:
each name becomes one column holding its **last** occurrence a row
reaches, which is what pyMzLib's `records` does.
[`sdrf_has_repeated_columns`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_has_repeated_columns.md)
says whether that applies. A position a row is too short to reach is
`NA`.

## Usage

``` r
sdrf_records(doc)
```

## Arguments

- doc:

  An
  [`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  or
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  result.

## Value

A data.frame with one row per returned row and one column per distinct
name, in order of first appearance. Names are kept verbatim
(`check.names = FALSE`).

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
sdrf_has_repeated_columns(doc)   # TRUE: the data.frame keeps the last of each name
#> [1] TRUE
records <- sdrf_records(doc)
records[, c("source name", "characteristics[organism]")]
#>   source name characteristics[organism]
#> 1    Sample 1     plasmodium falciparum
#> 2    Sample 2     plasmodium falciparum
#> 3    Sample 3     plasmodium falciparum
#> 4    Sample 4     plasmodium falciparum
#> 5    Sample 5     plasmodium falciparum
#> 6    Sample 6     plasmodium falciparum
```
