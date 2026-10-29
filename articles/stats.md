# Differential abundance: limma's moderated t, contrasts, and the tests around it

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
| which features differ by a contrast, within what interval? | `stats_fit(contrasts = )` | `EmpiricalBayes.ModerateContrasts` (#1452) | limma `contrasts.fit` + `topTable(confint = )` |
| what do several studies’ p-values say together? | [`stats_combine_p()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_combine_p.md) | `PValueCombination.Fisher`, `.Stouffer` (#1367) | Fisher 1925, Stouffer 1949 |
| which features rise or fall with a per-sample trait? | [`stats_correlate()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_correlate.md) | `SpearmanCorrelation.Correlate`, `.PartialCorrelate` (#1367, \#1446) | R `cor.test(method = "spearman")`, ppcor `pcor.test` |
| does a 0/1 outcome (detected or not) follow a coefficient? | [`stats_fit_logistic()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_logistic.md) | `LogisticRegression.Fit` (#1367) | R `glm(family = binomial)` per feature |
| which features change, when samples come in groups? | [`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md) | `MixedModel.Fit` (#1367) | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) with a random intercept per group |
| the same, when some samples are only counted (k of n trials)? | [`stats_fit_latent_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md) | `LogitNormalMixedModel.Fit` (#1469) | `lme4::glmer()` and `lmer()` in one model |
| how does a count’s chance of success follow the latent value? | [`stats_fit_link()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md) | `LogitNormalMixedModel.FitLink` (#1469) | maximum likelihood over the cells |
| does a search really control its FDR? | [`stats_entrapment_fdp()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_entrapment_fdp.md) | `EntrapmentFdp.Sweep` (#1385) | Wen et al. 2025, as FDRBench |

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
#>   prior df 4.532 (moments_legacy)
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
#>       1.612212e-11       2.414497e-13       2.414120e-13       3.156851e-14 
#>           df_total 
#>       3.727223e-16
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

### Unequal degrees of freedom: the marginal-likelihood prior

This experiment is the case the note on limma’s legacy estimator warns
about: oligos are missing from different runs, so they have different
residual degrees of freedom, and the default prior is limma’s legacy
estimator:

``` r

c(residual_df_differ = fit$residual_df_differ, estimator = fit$prior$estimator,
  prior_df = round(fit$prior$df, 1))
#> residual_df_differ          estimator           prior_df 
#>             "TRUE"   "moments_legacy"             "19.6"
```

mzLib has a second estimator for exactly this case,
`"marginal_likelihood"` (mzLib \#1368). It fits the prior by maximum
likelihood with **each oligo’s own degrees of freedom**, so an oligo
seen in five runs counts as five runs’ worth of evidence about the
variance, not as an average one. Ask for it by name:

``` r

ml <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
                file.path(data, "malat_dilution_design.tsv"), c("malat_250ng", "malat_125ng"),
                estimator = "marginal_likelihood")
c(estimator = ml$prior$estimator, prior_df = round(ml$prior$df, 1))
#>             estimator              prior_df 
#> "marginal_likelihood"                  "11"
ml_fitted <- ml$records[ml$records$status == "fitted", ]
ml_hits <- ml_fitted[!is.na(ml_fitted$bh_adjusted) & ml_fitted$bh_adjusted < 0.05, ]
table(ml_hits$coefficient, sub(":.*$", "", ml_hits$feature))
#>              
#>               20mer2 MALAT1
#>   malat_125ng      2      4
#>   malat_250ng      0      1
```

**The prior df falls**: once each oligo’s own degrees of freedom count,
the data support less shrinkage toward the common variance. The
estimates are the same least-squares numbers; only the standard errors
move, and with them the calls - compare this table with the one above.
At half load the FLuc call is gone, which is right, since FLuc did not
change.

- **It is mzLib’s own estimator, not limma’s.** Default limma (3.61 and
  later) handles unequal degrees of freedom with a different method
  (`fitFDistUnequalDF1`), so expect small differences from default
  limma.
- **mzLib keeps `"moments_legacy"` the default**, so that
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  reproduces `eBayes(legacy = TRUE)` exactly. Leaving `estimator` out
  always means mzLib’s default, and `prior$estimator` says which one you
  got.
- **A prior df of about 0.05 is a floor, not a finding.** mzLib searches
  d0 from e^-3 upward.
- One prior serves the whole call: every coefficient and every contrast
  is moderated with it.

## Contrasts: any comparison, with a confidence interval

A design’s coefficients answer one kind of question: each level against
the baseline. The questions a scientist asks are often other
comparisons: 500 ng against 250 ng when the baseline is 125 ng, or the
mean of two groups against a third. Each is a **contrast**, a weighted
sum of the coefficients, and `stats_fit(contrasts = )` tests it with the
same moderated t, plus a confidence interval: limma’s `contrasts.fit()`,
`eBayes()` and `topTable(confint = 0.95)`.

You rarely write the weights yourself.
[`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
builds the design from your sample sheet and writes the contrasts beside
it. With no reference level, every pair of levels is a contrast:

``` r

samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
design <- diff_design(samples, factors = list(factor_spec("malat_load")),
                      out = file.path(tempdir(), "malat_pairs"))
stratum <- design$strata[[1]]
stratum$column_names
#> [1] "intercept"        "malat_load=250ng" "malat_load=500ng"
stratum$warnings
#> [1] "No reference level was declared for factor 'malat_load', so every pair of its levels is compared (GR-14). Levels: 125ng, 250ng, 500ng."
labels <- setNames(design$records$label, design$records$contrast_id)
labels
#>                                     c1                                     c2 
#> "malat_load=250ng vs malat_load=125ng" "malat_load=500ng vs malat_load=125ng" 
#>                                     c3 
#> "malat_load=500ng vs malat_load=250ng"
```

mzLib took the first level, 125 ng, as the coding baseline, and wrote
`design_1.tsv` and `contrasts_1.tsv`. Fit the design and test every
contrast in the file. No coefficient is named, so the coefficient table
is empty and every answer is in `contrasts`:

``` r

pairs <- stats_fit(file.path(data, "malat_dilution_log2_vs_fluc.tsv"), stratum$design_file,
                   contrasts = stratum$contrasts_file)
pairs$contrasts
#> <mzlibr_moderated_contrasts> 3 contrast(s) x 492 features, 95% intervals
#>   c1: 0 1 0
#>   c2: 0 0 1
#>   c3:  0 -1  1
c(coefficient_rows = pairs$row_count, contrast_rows = pairs$contrasts$row_count)
#> coefficient_rows    contrast_rows 
#>                0             1476
```

The weights are on `pairs$coefficient_names`. `c1` and `c2` are single
coefficients. **`c3` is not**: it is the 500 ng coefficient minus the
250 ng one, a comparison no column of this design holds. The answers are
known in advance, +1, +2 and +1 in log2, so look at the MALAT1 oligos:

``` r

known <- c(c1 = 1, c2 = 2, c3 = 1)
summary_of <- function(id) {
  r <- pairs$contrasts$records
  r <- r[r$contrast == id & r$status == "fitted" & startsWith(r$feature, "MALAT1:"), ]
  c(known = known[[id]], median = round(median(r$estimate), 2), oligos = nrow(r),
    interval_holds_it = sum(r$confidence_low <= known[[id]] & known[[id]] <= r$confidence_high),
    excludes_zero = sum(!(r$confidence_low <= 0 & 0 <= r$confidence_high)))
}
t(sapply(c("c1", "c2", "c3"), summary_of))
#>    known median oligos interval_holds_it excludes_zero
#> c1     1   0.18    116               103             7
#> c2     2   1.01    116               105            23
#> c3     1   0.89    116               109            16
```

- **Each oligo’s 95% interval is wide.** It holds the known change for
  most oligos, and also holds 0 for most of them: one oligo, measured in
  a handful of runs, cannot tell a doubling from no change. The medians
  across oligos are much more precise than any one interval.
- **The contrasts locate the compression seen above.** 500 ng against
  250 ng reads close to its known +1; 250 ng against 125 ng reads far
  less. The shortfall is between the two lower loads.
- **The intervals are per feature and not adjusted.** For a decision
  across features use `bh_adjusted`, which each contrast has over its
  own fitted features.

Each row of `pairs$contrasts$records` carries what a coefficient row
does, plus `confidence_low` and `confidence_high`, in response units per
contrast unit: a log2 fold change here.

**Your own contrasts.** A contrasts file is a header of a contrast id
column and then exactly the design’s coefficient columns, in any order,
then one row per contrast: its id and one weight per coefficient. To
compare the mean of 250 ng and 125 ng against 500 ng on the design
above:

``` text
contrast           intercept  malat_load=250ng  malat_load=500ng
low_mean_vs_500    0          0.5               -1
```

(125 ng is the baseline, so its coefficient is 0 and needs no weight of
its own.) A column the design lacks, a design column the file lacks, a
weight that is not a number, or a row of all zeros is refused.

The traps:

- **A stratum can skip a contrast id.** Ids are numbered once over the
  whole analysis, so `c1` is the same comparison in every stratum, but a
  stratum without samples on both sides lists it in `not_run` instead.
  Read the labels from
  [`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md),
  as above.
- **With missing values the standard error is exact per feature**,
  computed on the samples that feature was observed in. limma’s
  `contrasts.fit()` approximates it there, so on data with missing
  values these numbers differ from limma’s by design; with complete data
  they agree.
- **A contrast of one coefficient is that coefficient’s test**, number
  for number. Asking for both in one call costs nothing: the model is
  fitted once and the prior is shared.
- `confidence = 0.9` (a fraction) changes the level; leave it out for
  mzLib’s default, and `contrasts$confidence_level` reports the level
  used.
- **It needs the pyMzLib 0.6.0 bridge.** An older bridge ignores
  `--contrasts` and `--estimator`, so
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  stops with an error naming the bridge it needs rather than return a
  fit that did not do what you asked.

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
  yourself.
  [`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
  writes one from a sample sheet, with a contrasts file for any
  comparison that is no single column; see
  [`vignette("diff", package = "mzLibR")`](https://smith-chem-wisc.github.io/mzLibR/articles/diff.md).

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

## Combining p-values across studies

Sometimes all you have from each study is a p-value.
[`stats_combine_p()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_combine_p.md)
turns several independent p-values for the same feature into one, by one
of two rules: `"fisher"`, X = -2 sum ln p against chi-square on 2k
degrees of freedom, driven by the smallest p-values; or `"stouffer"`, Z
= sum w z / sqrt(sum w^2) against the standard normal, which weighs
every study (for example by the square root of its sample size). There
is no default: mzLib prefers neither, so you choose. mzLib’s own
reference cases show the difference:

``` r

cases <- read.delim(file.path(data, "combination_inputs.tsv"))
cases[cases$case == "mixed", ]
#>     case     p weight
#> 21 mixed 0.001    1.5
#> 22 mixed 0.999    2.5
#> 23 mixed 0.500    3.5
#> 24 mixed 0.050    4.5
#> 25 mixed 0.950    5.5
fisher <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p), method = "fisher")
stouffer <- stats_combine_p(data.frame(feature = cases$case, p_value = cases$p, weight = cases$weight),
                            method = "stouffer")
data.frame(feature = fisher$records$feature, studies = fisher$records$studies,
           fisher = signif(fisher$records$p_value, 3), weighted_stouffer = signif(stouffer$records$p_value, 3))
#>   feature studies   fisher weighted_stouffer
#> 1     few       3 8.71e-02          0.150000
#> 2    tiny       4 3.14e-11          0.002280
#> 3     one       1 3.70e-02          0.037000
#> 4    many      12 1.13e-01          0.032300
#> 5   mixed       5 1.91e-02          0.713000
#> 6    tail       3 2.14e-18          0.000836
```

In `mixed`, one p-value of 0.001 sits among four that disagree. Fisher
says the five studies together are significant, because one of them is;
weighted Stouffer says they are not, because the heaviest study points
the other way. Both agree with R’s
[`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) and
[`qnorm()`](https://rdrr.io/r/stats/Normal.html):

``` r

r_answers <- read.delim(file.path(data, "combination_results.tsv"))
c(fisher = worst(fisher$records$p_value, r_answers$fisher_p),
  weighted_stouffer = worst(stouffer$records$p_value, r_answers$weighted_p))
#>            fisher weighted_stouffer 
#>      5.083515e-15      3.632781e-15
```

**The trap: two-sided p-values have no direction.** Both rules combine
*one-sided* p-values. A two-sided p-value - which is what
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
reports - says how surprising an effect is, not which way it went, so
two studies that disagree count as two that agree. Here are three
proteins, each tested in three cohorts (illustrative numbers, written
here, not measured):

``` r

cohorts <- data.frame(
  feature = rep(c("up_everywhere", "flips_sign", "no_effect"), each = 3),
  p_value = c(0.003, 0.02, 0.3, 0.001, 0.004, 0.4, 0.6, 0.45, 0.85),
  weight = c(4, 5, 3, 4, 5, 3, 4, 5, 3),
  effect = c(0.8, 0.5, 0.3, 1.2, -0.9, 0.2, -0.1, 0.1, -0.05))
ignored <- stats_combine_p(cohorts[, c("feature", "p_value", "weight")], method = "stouffer")
signed <- stats_combine_p(cohorts, method = "stouffer")
data.frame(feature = signed$records$feature, direction_ignored = signif(ignored$records$p_value, 2),
           signed = signif(signed$records$p_value, 2))
#>         feature direction_ignored  signed
#> 1 up_everywhere           6.2e-04 8.4e-05
#> 2    flips_sign           9.5e-05 4.3e-01
#> 3     no_effect           6.9e-01 4.4e-01
```

**Ignoring direction, `flips_sign` is the strongest result of the
three**: it rose sharply in one cohort and fell sharply in another. Give
each test its `effect` (only its sign is used) and mzLib first turns
each two-sided p into the one-sided p-value for effect \> 0, so opposite
effects cancel. The combined `p_value` is the **upper tail**: signed, it
tests effect \> 0, and a Stouffer `p_value` near 1 means a consistently
negative effect.

Before you trust a combined p-value: the studies must be independent
(two analyses of the same samples are one study); it is not adjusted
across features (send `records$p_value` to
[`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md));
and it is not an effect size. If every study reports an estimate and a
standard error on one scale,
[`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
pools those and says how much the studies disagree; prefer it.

## Correlating features with a trait

When the samples differ by a number rather than a group - age, dose, the
amount loaded - the question is which features rise or fall with it.
[`stats_correlate()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_correlate.md)
answers it with Spearman’s rank correlation: for each feature, mzLib
ranks the feature’s values and the trait over the samples where both are
known, correlates the ranks, and tests rho = 0, two-sided. It reads
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)’s
responses table and a **trait** table: one row per sample, then one
value column, whose header names the trait.

``` r

correlation <- stats_correlate(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
                               file.path(data, "malat_dilution_load_ng.tsv"))
correlation
#> <mzlibr_correlations> 492 features x 27 samples against malat_ng
#>   exact 4, asymptotic 308, not_estimable 180
r <- correlation$records
r$transcript <- sub(":.*$", "", r$feature)
estimable <- r[r$method != "not_estimable" & r$transcript %in% c("MALAT1", "FLuc", "20mer2"), ]
aggregate(rho ~ transcript, data = estimable, FUN = median)
#>   transcript       rho
#> 1     20mer2 0.4879500
#> 2       FLuc 0.0000000
#> 3     MALAT1 0.4330127
```

MALAT1’s oligos rise with its load, and FLuc’s, the spike every run was
normalised to, do not.

- `method` says how `p_value` was computed: `"exact"` (the full
  permutation distribution, for n of 9 or fewer with no ties, as R’s
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html)), `"asymptotic"`
  (t on n - 2 degrees of freedom, as `cor.test(exact = FALSE)`), or
  `"not_estimable"`, where `rho` and `p_value` are `NA`.
- **Small n cannot be significant.** Three samples have six orderings,
  so even a perfect rank correlation on three runs has p = 1/3.
- **An asymptotic p-value of exactly 0 means rho is exactly -1 or 1**:
  the t-test’s limit, not a measured probability.
- **For 10 to 1289 samples without ties, R’s default
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html) uses an
  Edgeworth series** that mzLib does not, so p-values there differ
  slightly from R’s default output, but not from
  `cor.test(exact = FALSE)`. This package’s tests hold every row to
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html) run in R.

**Partial correlation.** If something else also varies by sample -
batch, run order, sex - give it in a `given` table of the same shape.
mzLib then correlates what is left of the feature’s ranks and the
trait’s ranks after each is regressed on the covariates’ ranks, as
`ppcor::pcor.test(method = "spearman")` does:

``` r

plain <- stats_correlate(file.path(data, "limma_reference_responses.tsv"),
                         file.path(data, "limma_reference_age.tsv"))
partial <- stats_correlate(file.path(data, "limma_reference_responses.tsv"),
                           file.path(data, "limma_reference_age.tsv"),
                           given = file.path(data, "limma_reference_sex.tsv"))
rbind(plain = plain$method_counts, given_sex = partial$method_counts)
#>           exact asymptotic not_estimable
#> plain        22        378             0
#> given_sex     0        400             0
head(data.frame(feature = plain$records$feature, rho = round(plain$records$rho, 2),
                rho_given_sex = round(partial$records$rho, 2)), 5)
#>   feature  rho rho_given_sex
#> 1      g1 0.38          0.37
#> 2      g2 0.93          0.97
#> 3      g3 0.65          0.64
#> 4      g4 0.12          0.12
#> 5      g5 0.43          0.37
```

The partial test is always asymptotic, on n - 2 - k degrees of freedom
for k covariates. Partial out only nuisance variables: a covariate that
follows the trait takes the trait’s signal with it.

## Is a feature’s detection related to the condition?

Missing values in label-free data are rarely random: a protein near the
detection limit is missed more often where there is less of it. Ask the
data first: **does each feature’s chance of being detected depend on the
condition?** That is a logistic regression per feature on a table of 1
(measured) and 0 (not measured).
[`stats_fit_logistic()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_logistic.md)
runs one per feature, as R’s `glm(family = binomial)`, and tests each
coefficient with a Wald z. The dilution series has a detection table: 1
where the responses file has a value, 0 where it is blank.

``` r

detection <- stats_fit_logistic(file.path(data, "malat_dilution_detected.tsv"),
                                file.path(data, "malat_dilution_design.tsv"),
                                c("malat_250ng", "malat_125ng"))
detection$status_counts
#>               fitted too_few_observations       rank_deficient 
#>                  191                    0                    0 
#>            separated        not_converged 
#>                  301                    0
quarter <- detection$records[detection$records$coefficient == "malat_125ng", ]
quarter$transcript <- sub(":.*$", "", quarter$feature)
detected <- quarter[quarter$status == "fitted" & quarter$transcript %in% c("MALAT1", "FLuc"), ]
aggregate(estimate ~ transcript, data = detected, FUN = median)
#>   transcript   estimate
#> 1       FLuc  0.4054651
#> 2     MALAT1 -0.1823216
```

`estimate` is a change in the **log-odds** of detection: for a 0/1
column, the log of the odds ratio against the 500 ng runs.

**A separated feature has no estimate, and that is not the same as no
effect.** When the design splits a feature’s 0s from its 1s, no
maximum-likelihood estimate exists. R’s
[`glm()`](https://rdrr.io/r/stats/glm.html) prints a huge coefficient
and a warning there; mzLib reports `"separated"` and `NA`. Most oligos
here are separated because they were seen in a handful of runs, or none:

``` r

separated <- quarter[quarter$status == "separated", ]
head(sort(table(events = separated$events), decreasing = TRUE), 4)
#> events
#>  0  2  1  3 
#> 64 55 54 37
```

Here separation means “too sparse to say”. In a cleaner experiment the
same status can mark the strongest result in the table - a protein
detected in every treated sample and in no control - so read `events`
and `observed` before you discard one. Against R’s
[`glm()`](https://rdrr.io/r/stats/glm.html), on mzLib’s own reference
input:

``` r

glm_fit <- stats_fit_logistic(file.path(data, "glm_reference_responses.tsv"),
                              file.path(data, "glm_reference_design.tsv"), c("Intercept", "x", "grp"))
glm <- read.delim(file.path(data, "glm_reference_fit.tsv"))
x <- glm_fit$records[glm_fit$records$coefficient == "x", ]
ok <- glm$status == "fitted"
c(estimate = worst(x$estimate[ok], glm$b1[ok]), standard_error = worst(x$standard_error[ok], glm$se1[ok]))
#>       estimate standard_error 
#>   4.203008e-12   1.341800e-07
```

The estimates agree to better than 1e-8 and the standard errors to 1e-6
([`summary.glm()`](https://rdrr.io/r/stats/summary.glm.html) takes its
weights from the step before the estimate). Wald tests are weakest near
separation and in small samples, and no likelihood-ratio test is
offered. A 0 is an outcome: blank only a cell you genuinely do not know.

## When samples come in groups

Runs acquired on different days, replicates from one donor, studies
pooled into one analysis: each group can sit at its own baseline.
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
treats every sample as independent, so a day-to-day shift either
inflates its residual variance or, when the days are unbalanced, leaks
into the effect. A **mixed model** gives every group its own baseline, a
random intercept with variance tau^2, and estimates each coefficient
*within* groups:
[`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md),
as `nlme::lme(y ~ ..., random = ~ 1 | group)` per feature. The dilution
series’ runs were acquired on five days; a **groups** table says which:

``` r

mixed <- stats_fit_mixed(file.path(data, "malat_dilution_log2_vs_fluc.tsv"),
                         file.path(data, "malat_dilution_design.tsv"),
                         file.path(data, "malat_dilution_day.tsv"), c("malat_250ng", "malat_125ng"))
mixed
#> <mzlibr_mixed_model_fit> 492 features x 27 samples in 5 groups (day), reml, reporting malat_250ng, malat_125ng
#>   fitted 210, too_few_observations 214, rank_deficient 66, too_few_groups 2, variance_ratio_at_limit 0, not_converged 0
malat_medians <- function(result) {
  r <- result$records[result$records$status == "fitted" & startsWith(result$records$feature, "MALAT1:"), ]
  aggregate(cbind(estimate, standard_error) ~ coefficient, data = r, FUN = median)
}
malat_medians(fit)
#>   coefficient   estimate standard_error
#> 1 malat_125ng -1.0144266      1.0888865
#> 2 malat_250ng -0.8914901      0.9967881
malat_medians(mixed)
#>   coefficient  estimate standard_error
#> 1 malat_125ng -1.093018      0.6977282
#> 2 malat_250ng -1.048479      0.7064118
```

The design is unchanged: do not add a column per day, the random
intercept replaces it. `estimator` is `"reml"` because none was given
and that is mzLib’s default. Within days, the median standard errors are
smaller: the day-to-day shift no longer counts as noise.

- **A coefficient constant within each group is tested on the groups.**
  Its `df` is nlme’s “containment” df, a whole number per row; with too
  few groups there is none, and `p_value` is `NA` while `t` is reported.
- **`group_variance` 0 is the boundary**: no group effect. nlme cannot
  reach 0 and reports a tiny positive value there instead.
- **Under `"reml"` the log-likelihood is restricted**: compare it only
  between models with the same fixed effects, and use `estimator = "ml"`
  to compare designs.

Against nlme, on mzLib’s own reference input (seven groups of unequal
size):

``` r

lmm <- stats_fit_mixed(file.path(data, "lmm_reference_responses.tsv"),
                       file.path(data, "lmm_reference_design.tsv"),
                       file.path(data, "lmm_reference_groups.tsv"), c("Intercept", "age", "wp"))
nlme <- read.delim(file.path(data, "lmm_reference_fit.tsv"))
nlme <- nlme[nlme$method == "REML", ]
age <- lmm$records[lmm$records$coefficient == "age", ]
c(estimate = worst(age$estimate, nlme$b1), standard_error = worst(age$standard_error, nlme$se1))
#>       estimate standard_error 
#>   1.610138e-06   2.195060e-07
all(age$df == nlme$df1)
#> [1] TRUE
```

mzLib matches nlme to 1e-5 relative (nlme’s own optimizer precision).
Satterthwaite and Kenward-Roger df (lmerTest) are not implemented, so
p-values differ from lmerTest’s.

## When some samples are only counted

| You want to know | Call | mzLib computes it with |
|----|----|----|
| Does a value follow a trait, when some runs measured it and others only counted successes (k of n)? | [`stats_fit_latent_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md) | `LogitNormalMixedModel.Fit` (#1469) |
| How does a count’s chance of success follow the value (the link)? | [`stats_fit_link()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md) | `LogitNormalMixedModel.FitLink` (#1469) |

Some quantities are measured in some runs and only counted in others. A
PTM site’s occupancy is the usual case: where the modified and
unmodified peptides were both quantified, you have a number; where they
were not, you may know only that 3 of 8 covering peptidoforms carried
the modification.
[`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)
has to drop the counted runs.
[`stats_fit_latent_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md)
uses both kinds.

It models a **latent** value per run, such as the log-odds of the site’s
occupancy, with
[`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)’s
structure: fixed effects from your design, a random intercept per group
(an animal, a donor), and run-to-run scatter. A measured run reads the
latent value. A counted run has k of its n trials succeed, each with a
chance g set by the latent value through a **link**: logit g =
intercept + slope x latent + depth_slope x depth. With the identity link
(0, 1), a count samples the same probability a measurement reads. The
model is a binomial generalized linear mixed model (Breslow and Clayton
1993), the kind R users fit with `lme4::glmer` (Bates et al. 2015).

**You choose which cells are counted, and the link.** mzLibR builds
neither. mzLib’s PTM-QTL code does not yet do it either, so the counts
table is yours to write.

mzLib’s own test of this model simulates 16 animals of 4 runs each, with
a trait that is 0 or 1 per animal and a true trait effect of 1.0 on the
latent scale. Every other run is measured; the rest are counted, 5
trials each, through the link (0.4, 0.6). The package ships three of
those simulations, one per feature. The counts table has the responses
table’s layout, each cell `successes/trials`:

``` r

counts <- read.delim(file.path(data, "lnmm_simulated_counts.tsv"), colClasses = "character")
counts[1, 1:5]
#>   feature s1  s2 s3  s4
#> 1   seed1    3/5    2/5
latent <- stats_fit_latent_mixed(file.path(data, "lnmm_simulated_design.tsv"),
  file.path(data, "lnmm_simulated_groups.tsv"), c("Intercept", "trait"), link = c(0.4, 0.6),
  responses = file.path(data, "lnmm_simulated_responses.tsv"),
  counts = file.path(data, "lnmm_simulated_counts.tsv"))
latent
#> <mzlibr_latent_mixed_model_fit> 3 features x 64 samples in 16 groups (animal), reml, reporting Intercept, trait
#>   96 measured and 96 counted cells, through logit g = 0.4 + 0.6 x latent + 0 x depth, overdispersion 0
#>   fitted 3, too_few_observations 0, rank_deficient 0, separated 0, not_converged 0, too_few_groups 0
unlist(latent$link)
#>      intercept          slope    depth_slope overdispersion 
#>            0.4            0.6            0.0            0.0
```

`link` echoes what mzLib used. The depth slope and overdispersion were
not given, so they are mzLib’s defaults. Now fit the measured runs alone
with
[`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md),
and compare:

``` r

measured <- stats_fit_mixed(file.path(data, "lnmm_simulated_responses.tsv"),
  file.path(data, "lnmm_simulated_design.tsv"), file.path(data, "lnmm_simulated_groups.tsv"),
  c("Intercept", "trait"))
ours <- latent$records[latent$records$coefficient == "trait", ]
theirs <- measured$records[measured$records$coefficient == "trait", ]
data.frame(feature = ours$feature,
           latent = sprintf("%+.2f +- %.2f", ours$estimate, ours$standard_error),
           runs = paste(ours$measured, "+", ours$counted),
           measured_only = sprintf("%+.2f +- %.2f", theirs$estimate, theirs$standard_error),
           measured_runs = theirs$observed)
#>   feature        latent    runs measured_only measured_runs
#> 1   seed1 +1.61 +- 0.33 32 + 32 +1.63 +- 0.33            32
#> 2   seed2 +0.73 +- 0.26 32 + 32 +0.77 +- 0.27            32
#> 3   seed3 +1.26 +- 0.28 32 + 32 +1.36 +- 0.28            32
rbind(latent = latent$records$df[latent$records$coefficient == "Intercept"],
      measured_only = measured$records$df[measured$records$coefficient == "Intercept"])
#>               [,1] [,2] [,3]
#> latent          48   48   48
#> measured_only   16   16   16
```

Every trait estimate is within 2 standard errors of the true 1.0
(mzLib’s own test asks for 3):

``` r

all(abs(ours$estimate - 1) < 2 * ours$standard_error)
#> [1] TRUE
```

**The counted runs barely shrink the trait’s standard error, and that is
correct.** The trait is constant within each animal, so it is tested
between animals, on 16 animals minus 2 coefficients, 14 df in both fits.
More runs per animal sharpen each animal’s baseline, but the number of
animals limits how sure you can be about the trait. Effects that vary
within animals gain more: the intercept is tested on 48 df with the
counts, against 16 without.

What each row says, beyond
[`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)’s
columns:

- `estimate` and `standard_error` are on the **latent** scale: natural
  log-odds here, so `exp(estimate)` is an odds ratio of occupancy.
- `measured` and `counted` are the runs of each kind the feature
  contributed. `df` counts both.
- `residual_variance_given` is `TRUE` when sigma^2 came from your
  `residual_variance` table rather than the fit. Give it for a feature
  with too few measured runs to estimate its scatter.
- `status` can be `"separated"`: the counts are split by the design
  (every run of one group all 0, of another all full), so the estimate
  runs to infinity. You get the status, never a p-value.

**Check it against lme4 and glmmTMB.** mzLib’s reference input is 16
simulated features over 23 runs in 8 animals, every run counted, through
the link (-0.4, 0.6, 0.3) with a depth term, and R’s `glmer` fit of each
with 25-point adaptive quadrature. Those fits have no run-to-run
scatter, so each feature is given mzLib’s tiny sigma^2 of 2e-7, as
mzLib’s own comparison does:

``` r

counted <- stats_fit_latent_mixed(file.path(data, "lnmm_reference_design.tsv"),
  file.path(data, "lnmm_reference_groups.tsv"), c("Intercept", "trait", "x"),
  link = c(-0.4, 0.6, 0.3), counts = file.path(data, "lnmm_link_counts.tsv"),
  depth = file.path(data, "lnmm_link_depth.tsv"),
  residual_variance = file.path(data, "lnmm_link_residual_variance.tsv"), estimator = "ml")
glmer <- read.delim(file.path(data, "lnmm_counted_fit.tsv"))
glmer <- glmer[glmer$link == "link" & glmer$method == "ML", ]
sapply(c(Intercept = 0, trait = 1, x = 2), function(j) {
  ours <- counted$records[counted$records$coefficient == c("Intercept", "trait", "x")[j + 1], ]
  max(abs(ours$estimate - glmer[[paste0("b", j)]]) / glmer[[paste0("se", j)]]) < 1e-5
})
#> Intercept     trait         x 
#>      TRUE      TRUE      TRUE
```

Every estimate is within 1e-5 of glmer’s standard error from glmer’s
own. With run-to-run scatter estimated from counts alone, the reference
is `glmmTMB` (Brooks et al. 2017), which integrates by a Laplace
approximation where mzLib uses quadrature, so they agree to about 1%:

``` r

scatter <- stats_fit_latent_mixed(file.path(data, "lnmm_reference_design.tsv"),
  file.path(data, "lnmm_reference_groups.tsv"), c("Intercept", "trait", "x"), link = "identity",
  counts = file.path(data, "lnmm_scatter_counts.tsv"))
tmb <- read.delim(file.path(data, "lnmm_scatter_fit.tsv"))
tmb <- tmb[tmb$method == "REML", ]
trait <- scatter$records[scatter$records$coefficient == "trait", ]
c(group_variance = max(abs(trait$group_variance - tmb$tau2) / tmb$tau2) < 0.03,
  standard_error = max(abs(trait$standard_error - tmb$se1) / tmb$se1) < 0.03)
#> group_variance standard_error 
#>           TRUE           TRUE
```

### Fitting the link

The identity link says a count samples the very probability a
measurement reads. Often it does not: whether a modified form is seen at
all depends on how much of it there is, not only on its share.
[`stats_fit_link()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md)
fits the link from **cells**: each with a count and the normal its
latent value is taken to sit at (a mean and an SD, such as its feature’s
mean and SD over the measured runs, mzLib’s own suggestion). By default
it fits on **detection**: whether any trial succeeded, the event that
makes a covered zero.

mzLib’s own test of `FitLink` simulates 4,000 cells through the link
(-0.5, 0.5, 0.3) with a depth term, and the package ships them:

``` r

link <- stats_fit_link(file.path(data, "link_cells.tsv"), with_depth = TRUE)
str(link[c("cell_count", "with_depth", "with_overdispersion", "detection_only")])
#> List of 4
#>  $ cell_count         : num 4000
#>  $ with_depth         : logi TRUE
#>  $ with_overdispersion: logi FALSE
#>  $ detection_only     : logi TRUE
round(unlist(link$link), 2)
#>      intercept          slope    depth_slope overdispersion 
#>          -0.51           0.51           0.33           0.00
```

The link is recovered from detection alone. `with_overdispersion` and
`detection_only` were left out, so they are mzLib’s defaults, which the
result reports.
[`stats_fit_link()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md)’s
result goes straight into `stats_fit_latent_mixed(link = )`.

Before you trust a latent-scale mixed model:

- **The inputs are on the latent scale.** Send logit occupancy, not a
  fraction or an intensity.
- **The link is a modelling choice.** The identity link assumes a count
  samples the latent value itself; a fitted one depends on the cells you
  chose.
- **Standard errors are at the estimated variances.** They do not
  include those variances’ own uncertainty, so they are not `glmmTMB`’s
  [`vcov()`](https://rdrr.io/r/stats/vcov.html).
- **Under `"reml"`, the fixed effects are integrated out by a Laplace
  approximation,** as `glmmTMB`’s `REML = TRUE` does. With nothing
  counted it is exactly
  [`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)’s
  REML, and the two agree to about 1e-5. Under `"ml"`, estimates and
  variances agree with
  [`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)’s,
  but its standard errors are scaled differently.
- **The integrals are adaptive Gauss-Hermite quadrature** (Liu and
  Pierce 1994): 40 nodes over each animal’s baseline, 12 over each
  count’s run scatter, as `glmer`’s `nAGQ`. The source is
  `mzLib/StatisticalModels/LogitNormalMixedModel.cs` (mzLib \#1469).
- **Not adjusted across features.** Send `p_value` to
  [`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md).

## Does a search really control its FDR?

An entrapment experiment adds peptides that cannot be in the sample to
the search database. Every entrapment the search reports is then a known
false discovery, which lets you check the search’s FDR claim instead of
trusting it.
[`stats_entrapment_fdp()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_entrapment_fdp.md)
computes the three valid estimators of Wen et al. 2025 at each q-value
threshold: `lower_bound` (above your claimed FDR, the search fails to
control it; below, it proves nothing), `combined` (an upper bound: at or
below your claimed FDR, evidence of control) and `paired` (a tighter
upper bound for 1:1 pairing). The “sample” estimator N_E / (r N_T) is
invalid in both directions and is deliberately absent.

``` r

fdp <- stats_entrapment_fdp(file.path(data, "entrapment_reference_targets.tsv"),
                            file.path(data, "entrapment_reference_entrapments.tsv"), ratio = 1)
fdp
#> <mzlibr_entrapment_fdp_estimates> 7 targets, 4 entrapments, ratio 1, 100 of 100 thresholds with a discovery
fdp$records[fdp$records$q_value_threshold %in% c(0.003, 0.01), ]
#>    q_value_threshold target_count entrapment_count lower_bound  combined
#> 3              0.003            2                1   0.3333333 0.6666667
#> 10             0.010            4                3   0.4285714 0.8571429
#>       paired
#> 3  0.3333333
#> 10 0.8571429
```

This is mzLib’s hand-built search: seven targets and four entrapments.
At 1% the search reports 7 discoveries, 3 of them entrapments: far from
a 1% FDR. Give one row per unit you estimate for (precursor, peptide or
protein group), decoys removed, each with its best q-value; mzLib does
not deduplicate. `ratio` is required: mzLib has no default, and `paired`
is reported only at a ratio of exactly 1, so for 1:1 pairing pass 1, not
a ratio computed from database sizes after exclusions.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the step-up adjustment (Benjamini and Hochberg, J R Stat Soc B 57:289,
  1995).
- [doi:10.1111/j.1420-9101.2005.00917.x](https://doi.org/10.1111/j.1420-9101.2005.00917.x):
  the weighted Z (Stouffer) method and why it outperforms Fisher’s for a
  shared effect (Whitlock, J Evol Biol 18:1368, 2005).
- [doi:10.2307/1412159](https://doi.org/10.2307/1412159): the rank
  correlation coefficient rho (Spearman, Am J Psychol 15:72, 1904).
- [doi:10.5351/CSAM.2015.22.6.665](https://doi.org/10.5351/CSAM.2015.22.6.665):
  ppcor, whose pcor.test(method = ‘spearman’) the partial form (–given)
  reproduces (Kim, Commun Stat Appl Methods 22:665, 2015).
- [doi:10.2307/2347111](https://doi.org/10.2307/2347111): the Edgeworth
  series (AS 89) R’s cor.test uses by default for 10 \<= n \< 1290,
  which mzLib does NOT use (Best and Roberts, Appl Stat 24:377, 1975).
- [doi:10.1038/s41592-025-02719-x](https://doi.org/10.1038/s41592-025-02719-x):
  the lower-bound (eq. 2), combined (eq. 1) and paired (eq. 4)
  entrapment estimators mzLib implements, and why the sample estimator
  is invalid (Wen et al., Nat Methods 22:1454, 2025; the method FDRBench
  implements).
- [doi:10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284):
  generalized linear mixed models, the class of model this verb fits: a
  latent linear predictor with a random intercept, observed through a
  binomial link (Breslow and Clayton, J Am Stat Assoc 88:9, 1993).
- [doi:10.1093/biomet/81.3.624](https://doi.org/10.1093/biomet/81.3.624):
  adaptive Gauss-Hermite quadrature, centred at the integrand’s mode and
  scaled by its curvature, which mzLib uses to integrate the counted
  samples (Liu and Pierce, Biometrika 81:624, 1994).
- [doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01):
  lme4’s lmer and glmer (nAGQ 25), the references mzLib’s fits are held
  to (Bates, Maechler, Bolker and Walker, J Stat Softw 67:1, 2015).
- [doi:10.32614/RJ-2017-066](https://doi.org/10.32614/RJ-2017-066):
  glmmTMB, whose Laplace-approximated REML mzLib’s REML follows and
  whose fits with run scatter are a reference (Brooks et al., R Journal
  9:378, 2017).
- [doi:10.1093/biomet/58.3.545](https://doi.org/10.1093/biomet/58.3.545):
  restricted maximum likelihood, the default estimator (Patterson and
  Thompson, Biometrika 58:545, 1971).
- [doi:10.1111/j.2517-6161.1948.tb00014.x](https://doi.org/10.1111/j.2517-6161.1948.tb00014.x):
  the beta-binomial distribution, which a link with overdispersion gives
  a count cell (Skellam, J R Stat Soc B 10:257, 1948).
- [doi:10.1093/comjnl/7.4.308](https://doi.org/10.1093/comjnl/7.4.308):
  the simplex search mzLib maximizes the link’s likelihood with (Nelder
  and Mead, Comput J 7:308, 1965).
- [doi:10.2307/2344614](https://doi.org/10.2307/2344614): generalized
  linear models and their fitting by iteratively reweighted least
  squares, the estimator mzLib implements (Nelder and Wedderburn, J R
  Stat Soc A 135:370, 1972).
- [doi:10.1093/biomet/71.1.1](https://doi.org/10.1093/biomet/71.1.1):
  separation: when the outcomes are split by the design the
  maximum-likelihood estimate does not exist, which is the separated
  status (Albert and Anderson, Biometrika 71:1, 1984).
- [doi:10.2307/2529876](https://doi.org/10.2307/2529876): the
  random-effects linear model y = X beta + Z u + e this verb fits with
  one random intercept (Laird and Ware, Biometrics 38:963, 1982).
- [doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027):
  the empirical-Bayes moderated t this verb computes, for coefficients
  and contrasts (Smyth, Stat Appl Genet Mol Biol 3:3, 2004).
- [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007): limma,
  whose eBayes(legacy = TRUE), and contrasts.fit with topTable(confint),
  mzLib reproduces to 1e-8 relative (Ritchie et al., Nucleic Acids Res
  43:e47, 2015).
- [doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456%2886%2990046-2):
  the moment estimator of the between-study variance (DerSimonian and
  Laird, Control Clin Trials 7:177, 1986).
- [doi:10.18637/jss.v036.i03](https://doi.org/10.18637/jss.v036.i03):
  metafor, whose rma(method = ‘DL’) mzLib reproduces to 1e-8 relative
  (Viechtbauer, J Stat Softw 36(3), 2010).
