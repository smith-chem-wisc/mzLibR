# Number of quantified protein groups

Number of quantified protein groups

## Usage

``` r
flashlfq_protein_count(results)
```

## Arguments

- results:

  A
  [`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  result.

## Value

A single count.

## Examples

``` r

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
  match_between_runs = TRUE)
flashlfq_protein_count(quant)
#> [1] 2
```
