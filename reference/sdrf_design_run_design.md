# A label-free design as `flashlfq_median_polish()` takes its `design`

One row per run, keyed by `file_name` - the `Intensity_<file_name>`
column of a `QuantifiedPeptides.tsv` - with the four design fields,
0-based.

## Usage

``` r
sdrf_design_run_design(design)
```

## Arguments

- design:

  An
  [`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
  result.

## Value

A data.frame with `file_name`, `condition`, `biological_replicate`,
`technical_replicate` and `fraction`, in SDRF row order.

## See also

[`sdrf_design_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md),
[`flashlfq_median_polish`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)

## Examples

``` r

d <- sdrf_design("PXD067622.sdrf.tsv",
  condition_columns = c("factor value[genotype]", "factor value[treatment]"))
head(sdrf_design_run_design(d), 3)
#>                               file_name                       condition
#> 1 20240830_HF_LC3_MAA_RK_12032_CA_DMSO4 SPRTN-TurboID CA_DMSO (vehicle)
#> 2 20240830_HF_LC3_MAA_RK_12032_CA_DMSO5 SPRTN-TurboID CA_DMSO (vehicle)
#> 3 20240830_HF_LC3_MAA_RK_12032_CA_DMSO6 SPRTN-TurboID CA_DMSO (vehicle)
#>   biological_replicate technical_replicate fraction
#> 1                    0                   0        0
#> 2                    1                   0        0
#> 3                    2                   0        0
```
