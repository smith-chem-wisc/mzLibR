# Differential abundance: limma's moderated t, BH and meta-analysis

Most proteomics differential-abundance analyses end in limma: a linear
model per protein, with each protein’s noisy variance shrunk toward what
all the proteins together say (empirical-Bayes moderation, Smyth 2004).
R has limma. These functions are the same method **computed by mzLib**,
so an analysis in R reports the statistics MetaMorpheus and the Python
and Rust bindings report, from the same code. The first section lets you
check the agreement rather than take it on trust.

| question | function | mzLib | method |
|----|----|----|----|
| which features change with a coefficient, and how surely? | [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md) | `LinearModel.Fit`, `EmpiricalBayes.Moderate` (#1341) | limma `lmFit` + `eBayes(legacy = TRUE)` |
| which of my p-values survive a false discovery rate? | [`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md) | `MultipleTesting.BenjaminiHochberg` | Benjamini and Hochberg 1995 |
| what is one feature’s effect, pooled across studies? | [`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md) | `RandomEffectsMeta.Pool` (#1357) | DerSimonian and Laird 1986, as metafor `rma(method = "DL")` |

[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
already adjusts its own p-values, so you need
[`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
only for p-values from somewhere else, and
[`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
only when you have one estimate per study.

## Check it against limma

mzLib’s test suite carries a frozen limma run: 400 features, 12 samples
with missing values, a design of intercept, age in decades and sex, and
limma’s `eBayes` output for the age coefficient. This package ships the
input and limma’s output in `extdata/stats/`:

``` r

data <- system.file("extdata", "stats", package = "mzLibR")
fit <- stats_fit(file.path(data, "limma_reference_responses.tsv"),
                 file.path(data, "limma_reference_design.tsv"), "age_decades")
fit
#> <mzlibr_moderated_fit> 400 features x 12 samples, testing age_decades
#>   fitted 400, too_few_observations 0, rank_deficient 0
#>   prior df 4.532
#>   ! residual df differ between features: these are eBayes(legacy = TRUE) statistics
fit$coefficient_names
#> [1] "Intercept"   "age_decades" "sex"
```

Compare every feature with what limma said:

``` r

limma <- read.delim(file.path(data, "limma_ebayes_notrend.tsv"))
worst <- function(ours, theirs) max(abs(ours - theirs) / abs(theirs))
c(t = worst(fit$records$t, limma$t_age),
  p_value = worst(fit$records$p_value, limma$p_age),
  bh_adjusted = worst(fit$records$bh_adjusted, limma$bh_age),
  posterior_variance = worst(fit$records$posterior_variance, limma$s2_post),
  df_total = worst(fit$records$df_total, limma$df_total))
#>                  t            p_value        bh_adjusted posterior_variance 
#>       1.611198e-11       3.690533e-13       3.692947e-13       5.493775e-14 
#>           df_total 
#>       2.236334e-14
```

Every relative difference is below 1e-8. With limma installed you can
run it yourself on the same input,
`limma::eBayes(limma::lmFit(y, x), legacy = TRUE)`; this package’s tests
do exactly that, and `bh_adjusted` is
`p.adjust(p_value, method = "BH")`.

**Which limma: the *legacy* estimator.** This matches
`eBayes(legacy = TRUE)`. When features have different residual degrees
of freedom, limma 3.61 and later default to a different prior estimator.
Omitting missing values per feature is exactly what makes residual df
differ, so on most label-free data default limma reports a slightly
different prior and different moderated statistics. The fit says when
your input is such a case:

``` r

fit$residual_df_differ
#> [1] TRUE
```

## A real experiment with a known answer

mzLib’s RNA test data (mzLib \#1388) holds a dilution series - real
mass-spectrometry data where the right answer is known in advance.
Across 27 runs the lncRNA **MALAT1** was loaded at 500, 250 or 125 ng,
while the luciferase transcript **FLuc** stayed at 500 ng. Halving
MALAT1 should read as a log2 fold change of -1, quartering it as -2.

The responses are each oligo’s log2 FlashLFQ intensity minus that run’s
median log2 intensity of the FLuc-only oligos (blank where FlashLFQ
detected nothing), one row per oligo named `<transcript>:<sequence>`.
The design has an intercept, so each coefficient is the change from the
500 ng runs:

``` r

design <- read.delim(file.path(data, "malat_dilution_design.tsv"))
table(malat_250ng = design$malat_250ng, malat_125ng = design$malat_125ng)
#>            malat_125ng
#> malat_250ng  0  1
#>           0  7  8
#>           1 12  0
fit <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
                 file.path(data, "malat_dilution_design.tsv"), c("malat_250ng", "malat_125ng"))
fit$status_counts
#>               fitted too_few_observations       rank_deficient 
#>                  212                  214                   66
```

**Most oligos could not be fitted, and the fit says why rather than
dropping them.** An oligo seen in two or three runs has no residual
variance to estimate (`too_few_observations`); one seen only in the 500
ng runs cannot have a 250 ng effect (`rank_deficient`). Both are rows
with `NA` statistics, and neither enters the Benjamini-Hochberg family.

The question the experiment was designed to answer, as the median log2
fold change per transcript:

``` r

fitted <- fit$records[fit$records$status == "fitted", ]
fitted$transcript <- sub(":.*$", "", fitted$feature)
keep <- fitted$transcript %in% c("MALAT1", "FLuc")
aggregate(estimate ~ transcript + coefficient, data = fitted[keep, ], FUN = median)
#>   transcript coefficient   estimate
#> 1       FLuc malat_125ng -0.0458699
#> 2     MALAT1 malat_125ng -1.0144266
#> 3       FLuc malat_250ng -0.2164214
#> 4     MALAT1 malat_250ng -0.8914901
```

**Half the MALAT1 reads as about -0.9 against a known -1; a quarter
reads as about -1 against a known -2.** The second is compressed by
about half. The statistics did not create that: the estimates are
ordinary least squares, and moderation changes only their standard
errors. FLuc sits near 0 largely by construction, since every run was
normalised to it.

Which oligos pass a 5% false discovery rate?

``` r

hits <- fitted[!is.na(fitted$bh_adjusted) & fitted$bh_adjusted < 0.05, ]
table(hits$coefficient, hits$transcript)
#>              
#>               20mer2 FLuc MALAT1
#>   malat_125ng      2    0      4
#>   malat_250ng      0    1      2
```

**Few, and not all of them true.** Two of the quarter-load calls are
20-mer oligos, loaded at a constant amount. A 5% false discovery rate
bounds the expected share of wrong calls; it does not promise none, and
a known-answer experiment is how you see it happen.

## Your own data

[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
reads two tab-separated files: a **responses** table, one row per
feature and one column per sample, and a **design**, one row per sample,
matched to the responses’ columns by name.
`write.table(x, path, sep = "\t", quote = FALSE)` writes either. Three
things decide whether the answer means what you think:

- **Log-transform first.** The model is additive in whatever you pass;
  with log2 intensities, `estimate` is a log2 fold change.
- **A 0 is an observation.** Many quantification tables write 0 for “not
  measured”. Make those cells `NA`, or they are fitted as a real zero;
  `zero_count` counts them and a caveat says so.
- **The design is used as written.** Include the intercept column
  yourself, and write each comparison as its own column: contrasts are
  not implemented.

## P-values from anywhere

``` r

adjusted <- stats_adjust(c(0.0002, 0.004, 0.019, NA, 0.031, 0.2, NaN, 0.74))
adjusted$tested_count
#> [1] 6
adjusted$records
#>   p_value bh_adjusted
#> 1  0.0002      0.0012
#> 2  0.0040      0.0120
#> 3  0.0190      0.0380
#> 4      NA          NA
#> 5  0.0310      0.0465
#> 6  0.2000      0.2400
#> 7      NA          NA
#> 8  0.7400      0.7400
```

`NA` marks a feature that was not tested: it stays `NA` and is not
counted in m. Do not write 1 for it; that enlarges the family.

## Pooling across studies

[`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
takes one row per study - `feature`, `estimate`, `standard_error` - and
pools each feature with the DerSimonian-Laird random-effects model.
These are metafor’s own reference cases, from mzLib’s test data:

``` r

inputs <- read.delim(file.path(data, "metafor_dl_inputs.tsv"))
pooled <- stats_meta(data.frame(feature = inputs$case, estimate = inputs$yi,
                                standard_error = inputs$sei))
pooled$records[, c("feature", "studies", "estimate", "tau2", "i_squared", "direction_agree",
                   "direction_disagree", "leave_one_out_max_delta")]
#>         feature studies   estimate       tau2 i_squared direction_agree
#> 1   homogeneous       5  0.3096950 0.00000000 0.0000000               5
#> 2 heterogeneous       6  0.2512721 0.07516970 0.8227530               5
#> 3   two_studies       2 -0.2000000 0.06000000 0.4800000               1
#> 4  random_eight       8  0.1069573 0.09113603 0.8123322               6
#>   direction_disagree leave_one_out_max_delta
#> 1                  0              0.01413259
#> 2                  1              0.09481756
#> 3                  1              0.30000000
#> 4                  2              0.07461657
metafor <- read.delim(file.path(data, "metafor_dl_results.tsv"))
worst(pooled$records$estimate[match(metafor$case, pooled$records$feature)], metafor$estimate)
#> [1] 3.892522e-16
```

`direction_agree` and `leave_one_out_max_delta` are the two checks a
reader of a pooled result asks for first: how many studies point the
same way, and how far the answer moves without any one of them. The
interval uses a normal reference and does not widen for few studies, so
treat a pooled estimate from two or three studies with care.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the step-up adjustment (Benjamini and Hochberg, J R Stat Soc B 57:289,
  1995).
- [doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027):
  the empirical-Bayes moderated t this verb computes (Smyth, Stat Appl
  Genet Mol Biol 3:3, 2004).
- [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007): limma,
  whose eBayes(legacy = TRUE) mzLib reproduces to 1e-8 relative (Ritchie
  et al., Nucleic Acids Res 43:e47, 2015).
- [doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456%2886%2990046-2):
  the moment estimator of the between-study variance (DerSimonian and
  Laird, Control Clin Trials 7:177, 1986).
- [doi:10.18637/jss.v036.i03](https://doi.org/10.18637/jss.v036.i03):
  metafor, whose rma(method = ‘DL’) mzLib reproduces to 1e-8 relative
  (Viechtbauer, J Stat Softw 36(3), 2010).
