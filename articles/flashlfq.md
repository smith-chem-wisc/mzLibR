# Label-free quantification with FlashLFQ

[FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) is the
label-free quantification engine MetaMorpheus uses. Give it a search
result and the mzML runs it came from, and it measures how much of each
peptide and protein is in each run. The whole pipeline is mzLib’s: its
readers read the result file, its converter makes FlashLFQ
identifications, and `FlashLfqEngine` quantifies them. MetaMorpheus is
not involved.

| question | function | mzLib |
|----|----|----|
| how much of each peptide and protein is in each run? | [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md) | `FlashLfqEngine.Run` |
| fill in peptides missing from a run | `flashlfq_quantify(match_between_runs = TRUE)` | FlashLFQ’s match-between-runs |
| group runs into conditions and replicates | `flashlfq_quantify(spectra = <data.frame>)` | `SpectraFileInfo` |
| how many peaks did match-between-runs transfer? | [`flashlfq_mbr_peaks()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_peaks.md), [`flashlfq_mbr_peak_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_peak_count.md) | the peak table |
| re-roll proteins from a peptide table, without the mzML | [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md) | `FlashLfqResults.CalculateProteinResultsMedianPolish` |

The data are real: two K562 runs from mzLib’s FlashLFQ test data and
their MetaMorpheus search, quantified with match-between-runs on.

``` r

runs <- c("20100614_Velos1_TaGe_SA_K562_3.mzML", "20100614_Velos1_TaGe_SA_K562_4.mzML")
quant <- flashlfq_quantify("AllPSMs.psmtsv", runs, match_between_runs = TRUE)
quant
#> <mzlibr_quant> E:\GitClones\mzLib\mzLib\Test\FlashLFQ\TestData\AllPSMs.psmtsv, 594 identifications
#>   2 runs, 354 peptides, 943 protein groups, 647 peaks
#>   MBR: 140 transfers in peaks, 140 distinct peptides
#>   ! the peptide roll-up shows only 21 of those 140 transfers - read peaks, not peptides.
#>     See ?flashlfq_quantify.
#>   proteins: 10 NA (could not be resolved), 1734 zero (not measured)
```

Any quantifiable search result works, not only MetaMorpheus’s: an
MSFragger `psm.tsv` or a DIA-NN `report.tsv` too. mzLib converts
MSFragger’s retention times from seconds to minutes as it reads them
(mzLib \#1116), so FlashLFQ searches the right window.

## Four tidy tables

``` r

quant$spectra_files[, c("file_name", "peak_count", "mbr_peak_count")]
#>                        file_name peak_count mbr_peak_count
#> 1 20100614_Velos1_TaGe_SA_K562_3        340             62
#> 2 20100614_Velos1_TaGe_SA_K562_4        307             78
head(quant$peptides, 4)
#>                                             sequence        base_sequence
#> 1 AHQLVMEGYNWC[Common Fixed:Carbamidomethyl on C]HDR      AHQLVMEGYNWCHDR
#> 2 AHQLVMEGYNWC[Common Fixed:Carbamidomethyl on C]HDR      AHQLVMEGYNWCHDR
#> 3                               GGAAVDPDSGLEHSAHVLEK GGAAVDPDSGLEHSAHVLEK
#> 4                               GGAAVDPDSGLEHSAHVLEK GGAAVDPDSGLEHSAHVLEK
#>         protein_groups                      file_name  intensity detection_type
#> 1 H0YC23;P62714;P67775 20100614_Velos1_TaGe_SA_K562_3 1930193.13           MSMS
#> 2 H0YC23;P62714;P67775 20100614_Velos1_TaGe_SA_K562_4   80695.04           MSMS
#> 3    A0A7I2V3E1;P09874 20100614_Velos1_TaGe_SA_K562_3  946862.10           MSMS
#> 4    A0A7I2V3E1;P09874 20100614_Velos1_TaGe_SA_K562_4  448493.27           MSMS
head(quant$proteins[quant$proteins$intensity > 0 & !is.na(quant$proteins$intensity), ], 4)
#>    protein_group    gene_name     organism                      file_name
#> 65        P50990 primary:CCT8 Homo sapiens 20100614_Velos1_TaGe_SA_K562_3
#> 66        P50990 primary:CCT8 Homo sapiens 20100614_Velos1_TaGe_SA_K562_4
#> 91        Q9UQE7 primary:SMC3 Homo sapiens 20100614_Velos1_TaGe_SA_K562_3
#> 92        Q9UQE7 primary:SMC3 Homo sapiens 20100614_Velos1_TaGe_SA_K562_4
#>     intensity
#> 65 13493446.6
#> 66  3022424.2
#> 91   431063.4
#> 92  1342055.8
head(quant$peaks[, c("file_name", "sequence", "intensity", "detection_type", "retention_time")], 4)
#>                        file_name        sequence intensity detection_type
#> 1 20100614_Velos1_TaGe_SA_K562_3 RVHVTQEDFEMAVAK  502271.2            MBR
#> 2 20100614_Velos1_TaGe_SA_K562_3        IIGVDINK  258493.6            MBR
#> 3 20100614_Velos1_TaGe_SA_K562_3         MLEYALK  807482.4            MBR
#> 4 20100614_Velos1_TaGe_SA_K562_3        FIDVGGYK  309926.9            MBR
#>   retention_time
#> 1       79.93353
#> 2       81.45669
#> 3       82.14552
#> 4       81.95885
```

`peptides` and `proteins` are long: one row per peptide or protein group
per run, with the run in `file_name`. Intensities are in the
instrument’s units; `retention_time` in `peaks` is the apex, in minutes.

Most protein groups read 0 in both runs, and that is a setting, not a
failure. With `use_shared_peptides_for_protein_quant = FALSE`, the
default, only a group’s own peptides count, and most groups here have
none.

## 0 and `NA` mean different things

A **peptide** intensity of `0` means the peptide was not measured in
that run. A **protein** intensity of `NA` means FlashLFQ could not
resolve a number at all. mzLibR keeps the two apart, so
[`mean()`](https://rdrr.io/r/base/mean.html) over a protein column
returns `NA` rather than a confidently wrong number, and `na.rm = TRUE`
is a choice you make visibly.

``` r

c(not_resolved = sum(is.na(quant$proteins$intensity)),
  not_measured = sum(quant$proteins$intensity == 0, na.rm = TRUE))
#> not_resolved not_measured 
#>           10         1734
```

## Count transfers from the peaks

With `match_between_runs = TRUE`, a peptide identified in one run but
missing from another is still quantified in the second, by transferring
the identification across the aligned retention-time axis. Count
transfers from `peaks`, never from the peptide roll-up, which mirrors
FlashLFQ’s `QuantifiedPeptides.tsv` and drops most of them:

``` r

flashlfq_mbr_peak_count(quant)
#> [1] 140
table(quant$peptides$file_name[quant$peptides$detection_type == "MBR"])
#> 
#> 20100614_Velos1_TaGe_SA_K562_4 
#>                             21
head(flashlfq_mbr_peaks(quant)[, c("file_name", "sequence", "intensity")], 3)
#>                        file_name        sequence intensity
#> 1 20100614_Velos1_TaGe_SA_K562_3 RVHVTQEDFEMAVAK  502271.2
#> 2 20100614_Velos1_TaGe_SA_K562_3        IIGVDINK  258493.6
#> 3 20100614_Velos1_TaGe_SA_K562_3         MLEYALK  807482.4
```

None of run 3’s transfers appear in the peptide table: they read
`"NotDetected"` with intensity 0. A peptide by run matrix built from
`peptides` silently drops most of what MBR filled in; build it from
`peaks`.

## Re-quantify proteins under a new design

[`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
reruns only the protein roll-up, from the `QuantifiedPeptides.tsv`
FlashLFQ wrote, without re-reading any spectra. With no design, each run
is its own sample, and the answer is FlashLFQ’s own:

``` r

polished <- flashlfq_median_polish("K562_QuantifiedPeptides.tsv")
polished$samples
#>                            label condition biological_replicate
#> 1 20100614_Velos1_TaGe_SA_K562_3                              0
#> 2 20100614_Velos1_TaGe_SA_K562_4                              1
key <- function(groups, runs) paste(groups, runs)
same <- match(key(quant$proteins$protein_group, quant$proteins$file_name),
              key(polished$proteins$protein_group, polished$proteins$sample))
identical(quant$proteins$intensity, polished$proteins$intensity[same])
#> [1] TRUE
```

Every protein group’s intensity in every run is reproduced exactly. Give
`design` a data.frame with `file_name`, `condition` and
`biological_replicate` to say which runs are replicates of which sample,
and only the roll-up is redone.

## Before a real run

- Use a MetaMorpheus `.psmtsv` or `.osmtsv`, an MSFragger `psm.tsv` or a
  DIA-NN report, and mzML runs; convert `.raw` first.
- Give `spectra` a `condition` and `biological_replicate` so FlashLFQ
  knows which runs are comparable; match-between-runs needs a balanced
  design.
  [`sdrf_design_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md)
  builds that from an SDRF.
- `max_threads` defaults to 1 here and to every core in pyMzLib and
  mzLibRust;
  [`?flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  says why.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1021/acs.jproteome.7b00608](https://doi.org/10.1021/acs.jproteome.7b00608):
  FlashLFQ: the peak-finding, match-between-runs and protein
  quantification this verb runs (Millikin et al., J Proteome Res 2018).
