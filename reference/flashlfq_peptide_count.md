# Number of quantified peptides

Distinct peptide sequences, not rows: `peptides` is long, with one row
per peptide per run.

## Usage

``` r
flashlfq_peptide_count(results)
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
flashlfq_peptide_count(quant)
#> [1] 2
```
