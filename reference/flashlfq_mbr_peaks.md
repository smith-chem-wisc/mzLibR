# Exactly the peaks transferred by match-between-runs

Exactly the peaks transferred by match-between-runs

## Usage

``` r
flashlfq_mbr_peaks(results)
```

## Arguments

- results:

  A
  [`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  result.

## Value

The rows of `results$peaks` whose `detection_type` is `"MBR"`.

## Examples

``` r

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
  match_between_runs = TRUE)
flashlfq_mbr_peaks(quant)[, c("file_name", "sequence", "intensity")]
#>   file_name sequence intensity
#> 2     run_4 PEPTIDEK       500
```
