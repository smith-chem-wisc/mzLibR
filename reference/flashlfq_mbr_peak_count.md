# Total match-between-runs peaks across every run

Counts transferred **peaks**, not distinct peptides: one peptide rescued
in two runs is two peaks here. For "how many peptides did MBR rescue",
use
[`flashlfq_mbr_rescued_peptide_count`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_rescued_peptide_count.md).

## Usage

``` r
flashlfq_mbr_peak_count(results)
```

## Arguments

- results:

  A
  [`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  result.

## Details

Either way, do not count MBR from `peptides` - it under-counts, badly
and unevenly. Zero unless `match_between_runs` was on.

## Value

A single count.

## Examples

``` r

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
  match_between_runs = TRUE)
flashlfq_mbr_peak_count(quant)
#> [1] 1
```
