# A label-free design as `flashlfq_quantify()` takes its `spectra`

One row per run: `path` (the design's `full_path`) and the four design
fields, 0-based - so the SDRF drives FlashLFQ with no hand-written
design. Pass `searched_files` to
[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
when the SDRF names files without the directory you keep them in.

## Usage

``` r
sdrf_design_spectra(design)
```

## Arguments

- design:

  An
  [`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
  result.

## Value

A data.frame with `path`, `condition`, `biological_replicate`,
`technical_replicate` and `fraction`, in SDRF row order.

## See also

[`sdrf_design_run_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_run_design.md),
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)

## Examples

``` r

d <- sdrf_design("PXD067622.sdrf.tsv",
  condition_columns = c("factor value[genotype]", "factor value[treatment]"))
head(sdrf_design_spectra(d), 3)
#>                                        path                       condition
#> 1 20240830_HF_LC3_MAA_RK_12032_CA_DMSO4.raw SPRTN-TurboID CA_DMSO (vehicle)
#> 2 20240830_HF_LC3_MAA_RK_12032_CA_DMSO5.raw SPRTN-TurboID CA_DMSO (vehicle)
#> 3 20240830_HF_LC3_MAA_RK_12032_CA_DMSO6.raw SPRTN-TurboID CA_DMSO (vehicle)
#>   biological_replicate technical_replicate fraction
#> 1                    0                   0        0
#> 2                    1                   0        0
#> 3                    2                   0        0
```
