# Distinct peptides rescued by match-between-runs

**What this actually computes: the number of distinct modified sequences
among the MBR peaks.** That is stated in code terms on purpose, because
the prose definition a reader supplies - "peptides quantified in at
least one run \*only\* by MBR" - is subtly different and gives a
different number. On the K562 pair this is **140** while the strict
reading is **135**; the five that differ have both an MBR peak and a
zero-intensity MSMS peak in the same run.

## Usage

``` r
flashlfq_mbr_rescued_peptide_count(results)
```

## Arguments

- results:

  A
  [`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  result.

## Details

It equals
[`flashlfq_mbr_peak_count`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_peak_count.md)
only when no peptide was rescued in more than one run.

## Value

A single count.

## Examples

``` r

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
  match_between_runs = TRUE)
flashlfq_mbr_rescued_peptide_count(quant)
#> [1] 1
```
