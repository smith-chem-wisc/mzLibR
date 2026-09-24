# Which document each pooled row came from

Which document each pooled row came from

## Usage

``` r
sdrf_source_documents(pooled)
```

## Arguments

- pooled:

  An
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  result.

## Value

The `comment[source document]` label of each returned row.

## Examples

``` r

pooled <- sdrf_pool(c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv"),
  limit = 4)
sdrf_source_documents(pooled)
#> [1] "malaria" "malaria" "malaria" "malaria"
```
