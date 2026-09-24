# Label-free quantification with FlashLFQ

[`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
takes a search’s identifications and the mzML runs they came from, and
returns FlashLFQ’s peptide, protein and peak tables. The recording used
here is deliberately tiny - two runs, two peptides - so every table fits
on the page.

``` r

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
                           match_between_runs = TRUE)
quant
#> <mzlibr_quant> AllPSMs.psmtsv, 4 identifications
#>   2 runs, 2 peptides, 2 protein groups, 3 peaks
#>   MBR: 1 transfers in peaks, 1 distinct peptides
#>   proteins: 2 NA (could not be resolved), 0 zero (not measured)
```

## Four tidy tables

``` r

quant$spectra_files[, c("file_name", "peak_count", "mbr_peak_count")]
#>   file_name peak_count mbr_peak_count
#> 1     run_3          3              1
#> 2     run_4          2              0
quant$peptides
#>                  sequence base_sequence protein_groups file_name intensity
#> 1                PEPTIDEK      PEPTIDEK         P12345     run_3      1000
#> 2                PEPTIDEK      PEPTIDEK         P12345     run_4       500
#> 3 AC[Carbamidomethyl]DEFR        ACDEFR  P67890;Q11111     run_3      2000
#> 4 AC[Carbamidomethyl]DEFR        ACDEFR  P67890;Q11111     run_4         0
#>   detection_type
#> 1           MSMS
#> 2            MBR
#> 3           MSMS
#> 4    NotDetected
quant$proteins
#>   protein_group gene_name     organism file_name intensity
#> 1        P12345     GENE1 Homo sapiens     run_3      1000
#> 2        P12345     GENE1 Homo sapiens     run_4       500
#> 3        P67890     GENE2 Homo sapiens     run_3        NA
#> 4        P67890     GENE2 Homo sapiens     run_4        NA
quant$peaks
#>   file_name                sequence base_sequence intensity detection_type
#> 1     run_3                PEPTIDEK      PEPTIDEK      1000           MSMS
#> 2     run_4                PEPTIDEK      PEPTIDEK       500            MBR
#> 3     run_3 AC[Carbamidomethyl]DEFR        ACDEFR      2000           MSMS
#>   retention_time num_identifications protein_groups
#> 1           30.1                   1         P12345
#> 2           30.3                   1         P12345
#> 3           22.0                   1  P67890;Q11111
```

`peptides` and `proteins` are long: one row per peptide or protein group
per run, with the run in `file_name`. Intensities are in the
instrument’s units; `retention_time` in `peaks` is the apex, in minutes.

## 0 and `NA` mean different things

A **peptide** intensity of `0` means the peptide was not measured in
that run. A **protein** intensity of `NA` means FlashLFQ could not
resolve a number at all. mzLibR keeps the two apart, so
[`mean()`](https://rdrr.io/r/base/mean.html) over a protein column
returns `NA` rather than a confidently wrong number, and `na.rm = TRUE`
is a choice you make visibly.

## Count transfers from the peaks

With `match_between_runs = TRUE`, count transfers from `peaks`, never
from the peptide roll-up, which drops most of them. On mzLib’s own K562
pair the peaks hold 140 true transfers and the peptide table shows 52.

``` r

flashlfq_mbr_peaks(quant)
#>   file_name sequence base_sequence intensity detection_type retention_time
#> 2     run_4 PEPTIDEK      PEPTIDEK       500            MBR           30.3
#>   num_identifications protein_groups
#> 2                   1         P12345
flashlfq_mbr_peak_count(quant)
#> [1] 1
flashlfq_mbr_rescued_peptide_count(quant)
#> [1] 1
```

## Re-quantify proteins under a new design

[`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
reruns only the protein roll-up, from the `QuantifiedPeptides.tsv`
FlashLFQ wrote, without re-reading any spectra. The design says which
runs are replicates of which sample, and the result is keyed by sample:

``` r

polished <- flashlfq_median_polish("QuantifiedPeptides.tsv",
  design = data.frame(file_name = c("run_3", "run_4"), condition = c("control", "treated"),
                      biological_replicate = 0))
polished$samples
#>       label condition biological_replicate
#> 1 control_1   control                    0
#> 2 treated_1   treated                    0
polished$proteins
#>   protein_group gene_name     organism    sample intensity
#> 1            P1     GENE1 Homo sapiens control_1    3005.6
#> 2            P1     GENE1 Homo sapiens treated_1    6011.3
#> 3            P2     GENE2 Homo sapiens control_1    2262.7
#> 4            P2     GENE2 Homo sapiens treated_1    2368.5
#> 5            P3     GENE3 Homo sapiens control_1        NA
#> 6            P3     GENE3 Homo sapiens treated_1        NA
```

## Before a real run

- Use a MetaMorpheus `.psmtsv` or `.osmtsv`, and mzML runs; convert
  `.raw` first.
- Give `spectra` a `condition` and `biological_replicate` so FlashLFQ
  knows which runs are comparable; match-between-runs needs a balanced
  design.
- `max_threads` defaults to 1 here and to every core in pyMzLib and
  mzLibRust;
  [`?flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  says why.
