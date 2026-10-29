# Differential design and observations: the design matrix and every peptide, by mzLib

Before any differential test there are two tables to get right, and both
are easy to get quietly wrong. The **design**: which coefficient means
what, which level is the baseline, whether two things you want to
separate can be separated by these samples at all. And the
**observations**: what each peptide measured in each biological sample,
with fractions summed and technical replicates averaged correctly, and a
missing value kept as missing rather than written as 0.

[`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
turns the sample sheet you already keep into mzLib’s design matrix,
which equals R’s
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html), and says
which comparisons cannot be estimated and why.
[`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)
turns a FlashLFQ or MetaMorpheus peptide table into every peptide in
every sample, as a log2 intensity or `NA` with the reason. Both write
files that
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
reads.

| question | function | mzLib | method |
|----|----|----|----|
| what is my design matrix, and which comparisons can these samples answer? | [`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md) | `AnalysisDesign.ByStratum`, `.DefaultContrasts`, `.Estimability` (#1450) | treatment coding as R’s [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) (Wilkinson and Rogers 1973) |
| what did each peptide measure in each sample, and if nothing, why? | [`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md) | `FlashLfqObservations.FromPeptideTable`, `ObservationTable` (#1455) | fractions summed, technical replicates averaged in log2 |
| how far apart do the samples sit on the peptides they share? | [`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md) | `SharedPeptideNormalization.Measure`, `.Apply` | median of ratios to a geometric-mean reference (as Anders and Huber 2010) |

Then test with
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md):
it reads the design file and the contrasts file
[`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
writes, and the responses file
[`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)
writes. See
[`vignette("stats", package = "mzLibR")`](https://smith-chem-wisc.github.io/mzLibR/articles/stats.md).

## From a sample sheet to a design

The differential-abundance article analyses a MALAT1 dilution series -
27 runs at 500, 250 or 125 ng - with a hand-made design of 0/1 columns.
Here is the same experiment as a scientist writes it down: one row per
run, and its load.

``` r

diff_data <- system.file("extdata", "diff", package = "mzLibR")
stats_data <- system.file("extdata", "stats", package = "mzLibR")
sheet <- read.delim(file.path(diff_data, "malat_dilution_samples.tsv"))
head(sheet, 3)
#>                                                             run malat_load
#> 1 2026_09_18_1MALAT_500ng_1FLuc_500ng_20mer-2_50ng_147mingrad-1      500ng
#> 2 2026_09_18_1MALAT_500ng_1FLuc_500ng_20mer-2_50ng_147mingrad-2      500ng
#> 3 2026_09_24_1MALAT_500ng_1FLuc_500ng_20mer-2_50ng_147mingrad-1      500ng
table(sheet$malat_load)
#> 
#> 125ng 250ng 500ng 
#>     8    12     7
```

Say which column is a factor and which level is the baseline. Everything
else is mzLib’s:

``` r

design <- diff_design(file.path(diff_data, "malat_dilution_samples.tsv"),
                      factors = list(factor_spec("malat_load", "500ng")),
                      out = file.path(tempdir(), "malat_design"))
stratum <- design$strata[["all"]]
stratum$column_names
#> [1] "intercept"        "malat_load=125ng" "malat_load=250ng"
c(rank = stratum$rank, aliased = length(stratum$aliased_terms), warnings = length(stratum$warnings))
#>     rank  aliased warnings 
#>        3        0        0
design$records[, c("contrast_id", "label", "estimability")]
#>   contrast_id                                label estimability
#> 1          c1 malat_load=125ng vs malat_load=500ng    estimable
#> 2          c2 malat_load=250ng vs malat_load=500ng    estimable
design$records$fitted_weights
#> [[1]]
#> [1] 0 1 0
#> 
#> [[2]]
#> [1] 0 0 1
```

- **The columns** are named `factor=level`. The reference comes first
  and has no column; the other levels follow in ordinal (character)
  order, so `125ng` before `250ng`. That is R’s
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) with the
  same level order.
- **The default contrasts** are each level against the reference. Each
  is one design column here, so `fitted_weights` has a single 1.
- **`out`** wrote `design_1.tsv`, the design table
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  reads, and `contrasts_1.tsv`, each estimable contrast’s weights on the
  same columns, which `stats_fit(contrasts = )` reads.

Fit it. Each contrast against the reference is its design column, tested
by name:

``` r

responses <- file.path(stats_data, "malat_dilution_log2_vs_fluc.tsv")
fit <- stats_fit(responses, stratum$design_file, c("malat_load=250ng", "malat_load=125ng"))
fit$status_counts
#>               fitted too_few_observations       rank_deficient 
#>                  212                  214                   66
malat <- fit$records[fit$records$status == "fitted" & startsWith(fit$records$feature, "MALAT1:"), ]
aggregate(estimate ~ coefficient, data = malat, FUN = median)
#>        coefficient   estimate
#> 1 malat_load=125ng -1.0144266
#> 2 malat_load=250ng -0.8914901
```

The same numbers as the hand-made design, because it is the same model:
this package’s live tests hold the two fits equal, feature by feature.

**Or test the contrasts file itself**, and get a confidence interval
with each answer. Name the coefficients too, and one call returns both
tables:

``` r

both <- stats_fit(responses, stratum$design_file, c("malat_load=250ng", "malat_load=125ng"),
                  contrasts = stratum$contrasts_file)
both$contrasts$weights
#> $c1
#> [1] 0 1 0
#> 
#> $c2
#> [1] 0 0 1
c2 <- both$contrasts$records[both$contrasts$records$contrast == "c2", ]
coefficient <- both$records[both$records$coefficient == "malat_load=250ng", ]
identical(c2$t, coefficient$t)
#> [1] TRUE
row <- c2[c2$status == "fitted" & startsWith(c2$feature, "MALAT1:"), ][1, ]
row[, c("feature", "observed", "estimate", "confidence_low", "confidence_high")]
#>                  feature observed   estimate confidence_low confidence_high
#> 493 MALAT1:AAAAAAAAACUCG        4 -0.7533737      -3.681669        2.174922
```

A contrast of one column is that column’s test, number for number, so
`c2` repeats the `malat_load=250ng` rows and adds the interval. This
oligo was measured in 4 runs, and its 95% interval cannot tell halving
from no change. The contrasts file earns its keep for comparisons that
are no single column, such as 250 ng against 125 ng: declare no
reference level and
[`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
writes every pair, as the differential-abundance article shows.

## Covariates, units and batch

A covariate is a number in a unit. Say the unit, and how the slope
should read: here age in years, as a slope per decade, centred at 50
years. Batch is a fixed block. The sheet is mzLib’s own test design, and
base R’s [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
for it ships beside it:

``` r

design <- diff_design(file.path(diff_data, "design_factor_batch_covariate_samples.tsv"),
                      factors = list(factor_spec("condition", "ctrl")),
                      covariates = list(covariate_spec("age", "years", scale = 10, scaled_unit = "decade",
                                                       centre = 50)),
                      batch = "batch")
stratum <- design$strata[[1]]
stratum$column_names
#> [1] "intercept"         "condition=treated" "batch=b2"         
#> [4] "age (per decade)"
design$records[2, c("label", "covariate_unit", "covariate_scale")]
#>              label covariate_unit               covariate_scale
#> 2 age (per decade)         decade years/10, centred at 50 years
r_matrix <- as.matrix(read.delim(file.path(diff_data, "design_factor_batch_covariate_matrix.tsv"),
                                 check.names = FALSE))
colnames(r_matrix)
#> [1] "(Intercept)"      "conditiontreated" "batchb2"          "I((age - 50)/10)"
max(abs(unname(stratum$matrix) - unname(r_matrix)))
#> [1] 0
```

Every value equals R’s
`model.matrix(~ condition + batch + I((age - 50)/10))`; this package’s
tests also run
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) here and
compare.

**Leave a field out and mzLib chooses.**
`covariate_spec("age", "years")` sends no scale, scaled unit or centre,
so the slope is per year and the centre is the stratum’s mean;
`covariate_centres` reports the centre used. Batch enters the model only
when you name its column, and drops out with a warning in a stratum that
has only one batch.

## When the samples cannot answer

Here every sample of condition A ran in batch 1 and every sample of B in
batch 2:

``` r

design <- diff_design(file.path(diff_data, "design_confounded_samples.tsv"),
                      factors = list(factor_spec("condition", "A")), batch = "batch")
stratum <- design$strata[[1]]
stratum$column_names
#> [1] "intercept"   "condition=B" "batch=2"
c(rank = stratum$rank)
#> rank 
#>    2
stratum$aliased_terms
#> [1] "condition" "batch"
design$records[, c("label", "estimability")]
#>                        label             estimability
#> 1 condition=B vs condition=A not_estimable:confounded
is.null(design$records$fitted_weights[[1]])
#> [1] TRUE
```

**No analysis of these samples can separate the condition from the
batch**, and mzLib says so instead of returning a number. The verdicts
are `"estimable"`; `"not_estimable:confounded"`, where the comparison
lies along a direction two or more terms share, as here; and
`"not_estimable:rank_deficient"`, where it needs a column no sample
informs, such as an interaction cell with no samples.

The design file still holds a design that can be fitted:
`fitted_column_names` leaves out the later of two aliased columns,
because a fit refuses a design without full rank. A contrast that is not
estimable has no fitted weights (`NULL`), so it never reaches the
contrasts file.

## Strata: one design per tissue

When samples are different kinds of thing - liver and muscle - a
comparison is made within each kind. Name the stratum column, and each
stratum gets its own design:

``` r

design <- diff_design(file.path(diff_data, "stratified_samples.tsv"),
                      factors = list(factor_spec("condition", "A")), strata = "organism part")
lapply(design$strata, function(s) s$sample_ids)
#> $`organism part=kidney`
#> [1] "s7" "s8"
#> 
#> $`organism part=liver`
#> [1] "s1" "s3" "s5"
#> 
#> $`organism part=muscle`
#> [1] "s2" "s4" "s6"
design$strata[["organism part=kidney"]]$not_run
#>                        label         reason id
#> 1 condition=B vs condition=A not_in_stratum c1
design$records[, c("stratum", "contrast_id", "label")]
#>                stratum contrast_id                      label
#> 1  organism part=liver          c1 condition=B vs condition=A
#> 2 organism part=muscle          c1 condition=B vs condition=A
```

Kidney has only condition A, so there is nothing to compare there, and
`not_run` says so: it skips `c1`, `condition=B vs condition=A`. **A
contrast id is numbered once over the whole analysis:** `c1` is the same
comparison in liver and in muscle, and a stratum that cannot run it
skips the id. With `out`, stratum *i* writes `design_i.tsv`, numbered in
the order of `strata`, because a stratum key can hold characters a file
name cannot.

## From a peptide table to observations

[`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)
reads the peptide table FlashLFQ (`QuantifiedPeptides.tsv`) or
MetaMorpheus (`AllQuantifiedPeptides.tsv`) wrote, and needs the runs:
which biological sample each file belongs to and, if it applies, its
fraction and technical replicate. **Grouping comes only from what you
give here**; mzLib never reads it out of a file name. This MetaMorpheus
table is mzLib’s test data: three peptides in 18 runs, each run its own
sample:

``` r

peptides <- file.path(diff_data, "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv")
groups <- file.path(diff_data, "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
header <- strsplit(readLines(peptides, n = 1), "\t")[[1]]
files <- sub("-calib$", "", sub("^Intensity_", "", header[startsWith(header, "Intensity_")]))
head(files, 2)
#> [1] "QE-002106_GM1_a" "QE-002107_GM1_b"
runs <- data.frame(file = files, sample = files)
observed <- diff_observations(peptides, runs, basis = "mbr_kept", protein_groups = groups)
observed
#> <mzlibr_observations> 3 peptides x 18 samples, mbr_kept
#>   quantified 8, mbr_transferred 1, ambiguous_peak 0, identified_not_quantified 0, not_detected 45
#>   normalization: none on 1 peptides
```

MetaMorpheus quantifies calibrated files, so its columns end in
`-calib`; give the file name without it, and mzLib matches either. Each
row is one peptide in one sample:

``` r

transferred <- observed$records[observed$records$state == "mbr_transferred", ]
transferred[, c("base_sequence", "sample", "log2_intensity", "unique_protein_group")]
#>    base_sequence          sample log2_intensity unique_protein_group
#> 10     CACASHVAK QE-002118_GM6_a       17.13515               Q16643
observed$records$state[is.na(observed$records$log2_intensity)][1]
#> [1] "not_detected"
```

**`basis` decides whether match-between-runs values count**, and mzLib
has no default: it is your scientific choice. Under `"msms_only"` the
same cell has no value, and still says why:

``` r

only <- diff_observations(peptides, runs, basis = "msms_only", protein_groups = groups)
same <- only$records[only$records$sample == "QE-002118_GM6_a" & only$records$base_sequence == "CACASHVAK", ]
same[, c("log2_intensity", "state")]
#>    log2_intensity           state
#> 10             NA mbr_transferred
```

The protein-group table adds each group’s q-value and flags. It is how a
contaminant is recognised, and mzLib warns about groups the peptides
name that it does not describe:

``` r

observed$protein_groups[observed$protein_groups$is_contaminant | observed$protein_groups$is_decoy, c("name", "genes")]
#>           name  genes
#> 4       P02769    ALB
#> 6 DECOY_P62750 RPL23A
observed$warnings
#> [1] "2 protein group(s) named by peptides are not in the protein-group table, so their q-value and decoy/contaminant status are unknown: Q16643, Q9UI26."
```

What each `state` means:

- `quantified`: a value from MS/MS-identified peaks only.
- `mbr_transferred`: the value includes a match-between-runs peak. Under
  `"msms_only"` it is left out, and `log2_intensity` is `NA`.
- `ambiguous_peak`: a peak was found but withheld as ambiguous. No
  value.
- `identified_not_quantified`: identified by MS/MS, but no peak could be
  quantified. No value.
- `not_detected`: no value.

A sample whose runs disagree takes the most informative state among
them. With fractions, mzLib sums a technical replicate’s fractions on
the linear scale, then averages the sample’s replicates in log2: give
`runs` `fraction` and `technical_replicate` columns.

## Normalization, and the assumption it makes

[`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)
always measures how far each sample sits from the others, and subtracts
it only when you ask. mzLib’s rule: over the peptides with a value in
every sample (decoys, contaminants and entrapments excluded), each
sample’s shift is the median of its log2 value minus the peptide’s mean.
Two runs of FlashLFQ’s E. coli spike-in vignette, which their names say
were spiked at 1x and 2x:

``` r

spike_table <- file.path(diff_data, "EcoliSpikeIn_QuantifiedPeptides.tsv")
spike_runs <- data.frame(file = c("09-04-18_EcoliSpikeInSingleShot1x", "09-04-18_EcoliSpikeInSingleShot2x"),
                         sample = c("spike_1x", "spike_2x"))
spike <- diff_observations(spike_table, spike_runs, basis = "msms_only")
spike$normalization[c("setting", "reference_set_size")]
#> $setting
#> [1] "none"
#> 
#> $reference_set_size
#> [1] 86
spike$normalization$warnings
#> [1] "86 reference peptide(s) have a value in every sample and 86 in at least half, under the floor of 100, so no sample is shifted (none)."
```

**Below mzLib’s floor of 100 shared peptides, nothing is shifted**,
because a median of a few peptides can move a whole sample by itself.
You can lower the floor with `min_peptides`. Here is what that does on
this table:

``` r

shifted <- diff_observations(spike_table, spike_runs, basis = "msms_only", min_peptides = 50,
                             normalize = TRUE, responses_out = file.path(tempdir(), "spike_in_responses.tsv"))
shifted$normalization$setting
#> [1] "shared_peptide_median"
round(shifted$normalization$shifts, 2)
#> spike_1x spike_2x 
#>    -0.47     0.47
by_peptide <- split(spike$records$log2_intensity, spike$records$full_sequence)
both_runs <- Filter(function(v) !anyNA(v), by_peptide)
c(peptides = length(both_runs), median_2x_minus_1x = round(median(sapply(both_runs, diff)), 2))
#>           peptides median_2x_minus_1x 
#>              86.00               0.94
```

**The normalization removed the spike-in.** This slice of the vignette
holds almost only E. coli peptides, so nearly every peptide was doubled.
The peptides measured in both runs sit close to one log2 unit higher in
the 2x run, the doubling the names promise, and the two shifts take all
of it away. Median normalization assumes most peptides do not change
between samples. When most of them do - a table of the spiked organism
alone, a loading series, a pull-down - that assumption is false, and
normalizing deletes the result. **Read the shifts before you subtract
them.**

`responses_out` also wrote the table’s own peptide-by-sample shape,
which
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
reads as its responses:

``` r

shifted$responses_written[c("row_count", "sample_count")]
#> $row_count
#> [1] 243
#> 
#> $sample_count
#> [1] 2
```

Its rows are peptides by full sequence, its columns the samples, and an
empty cell means no value. For more than one stratum, narrow
[`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)
to one stratum’s samples with `samples` (a design file from
[`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
will do): the shifts are measured over every sample of the table, so
measure each stratum on its own.

## What these functions will not do

- **No protein roll-up and no row choice.** mzLib does not yet turn the
  observation table into a protein-by-sample table, choose which rows a
  test uses (decoys, contaminants and q-value cuts included), or
  assemble fits into its differential result table. mzLibR does not do
  it either. Every peptide is a row; filter with `protein_groups`,
  `unique_protein_group` and the result’s `protein_groups` table.
- **Only the default contrasts.**
  [`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
  lists each level against its reference (every pair when there is none)
  and each covariate’s slope. mzLib’s custom contrasts have no wire form
  yet; write your own contrasts file for
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  instead.
- **The peptide table does not record FlashLFQ’s design**, so mzLib
  cannot check that your runs group the files as FlashLFQ did. Give the
  runs the design FlashLFQ quantified with.
- **A table written in a comma-decimal locale does not read correctly.**
  FlashLFQ writes in the machine’s culture, and mzLib reads with a
  decimal point.
- **mzLibR does no arithmetic.** Every number on this page came from
  mzLib.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.2307/2346786](https://doi.org/10.2307/2346786): the
  model-formula notation (main effects and a:b interactions) that R’s
  stats::model.matrix implements and mzLib’s design matrix equals
  (Wilkinson and Rogers, Appl Stat 22:392, 1973).
- [doi:10.1021/acs.jproteome.7b00608](https://doi.org/10.1021/acs.jproteome.7b00608):
  FlashLFQ, whose peptide table and detection types this verb reads
  (Millikin et al., J Proteome Res 2018).
- [doi:10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106):
  the median-of-ratios idea the shared-peptide normalization follows:
  each sample against the log of a geometric-mean reference (Anders and
  Huber, Genome Biol 11:R106, 2010).
