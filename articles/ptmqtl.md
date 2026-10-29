# PTM site occupancy and PTM-QTL: occupancy, pairs and trait effects, by mzLib

A modification site is rarely all-or-nothing. In one run 15% of a serine
is phosphorylated, in the next 40%. That fraction, the site’s
**occupancy**, is what changes with age, dose or genotype, and what two
sites can share. Measuring it from label-free data is easy to get wrong:
a modified form that was identified but not quantified is not a zero, a
site seen only in its modified form is not truly 100%, and two sites
read off one peptide are correlated by construction. mzLib’s
`Quantification.PtmQtl` engine keeps those cases apart.

| question | function | mzLib | method |
|----|----|----|----|
| what fraction of each site carries its modification, run by run? | [`ptmqtl_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md) | `SiteOccupancyCalculator.Calculate` (#1369) | intensity occupancy, as MetaMorpheus’s `IntensityOccupancy_` |
| the same, for a sample measured as several fractions? | `ptmqtl_occupancy(samples = )` | `SiteOccupancyCalculator.Calculate` with a run-to-sample map | fractions summed before the ratio |
| which two modifications sit on one molecule? | `ptmqtl_pairs(type = "physical")` | `PtmPairEngine.Physical` (#1430) | median co-occupancy |
| which two sites rise and fall together across runs? | `ptmqtl_pairs(type = "covarying")` | `PtmPairEngine.CoVarying` (#1430, \#1446) | Spearman’s rho, Benjamini-Hochberg |
| does a site’s occupancy change with a per-run trait? | [`ptmqtl_trait_effects()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md) | `SiteTraitEffects.Fit` (#1446) | a mixed model per site, a random intercept per replicate |

Every number below comes from mzLib’s own hand-computed test cases, so
each is one mzLib’s tests assert.

## The input: one row per peptidoform, per run, per protein

Every function reads a tab-separated table with a header row, mzLib’s
own `PeptidoformObservation` as columns: `run`; `full_sequence` in
MetaMorpheus notation (`PEPS[Common Biological:Phosphorylation on S]K`);
`protein_accession`; `start_residue` and `end_residue`, the peptide’s
first and last residue in that protein, counting from 1; `intensity`,
its MS/MS-quantified intensity in that run, **blank** when it was
identified but not quantified; and optionally `previous_residue`, the
protein residue before the peptide. Four rules decide whether the
numbers mean anything:

- **One row per (run, peptidoform, protein), never one per PSM.** Two
  rows for one peptidoform in one run would count its intensity twice;
  mzLib refuses them.
- **MS/MS identifications only.** A match-between-runs transfer is not
  an identification.
- **Blank, not 0, for “identified but not quantified”.** FlashLFQ writes
  0 for a peptide it did not quantify. Blank those cells.
- **A shared peptide is one row per protein.** Each protein then gets
  its full intensity. MetaMorpheus instead splits a shared peptide
  across protein groups, so its `IntensityOccupancy_` columns can
  differ.

**mzLib does not build this table for you yet.** It has no reader that
turns a MetaMorpheus or FlashLFQ result into these rows, because that
step makes choices: which PSMs count, how a shared peptide is placed.
mzLibR will not make them on mzLib’s behalf. You can assemble the table
from two files mzLibR already reads, the PSMs (for runs, proteins and
positions) and the FlashLFQ peptide table (for intensities). One way,
with your own choices marked:

``` r

# Not run: needs your MetaMorpheus search output.
psms <- readers_read_records("AllPSMs.psmtsv")$records
peptides <- readers_read_quantified_peptides("AllQuantifiedPeptides.tsv")$records

# Intensity per (run, full sequence), MS/MS-detected only. FlashLFQ's 0 means "not quantified".
msms <- peptides[peptides$detection_type == "MSMS" & !is.na(peptides$intensity) & peptides$intensity > 0, ]
intensity <- setNames(msms$intensity, paste(msms$sample_label, msms$sequence, sep = "\r"))

# MetaMorpheus writes one value when every protein shares it, else one per protein.
per_protein <- function(cell, count) {
  parts <- strsplit(cell, "|", fixed = TRUE)[[1]]
  if (length(parts) == 1L) rep(parts, count) else parts
}

# Your choice of PSMs: targets, q-value below 1%, one unambiguous peptidoform.
kept <- psms[psms$decoy_contam_target == "T" & psms$q_value < 0.01 & psms$ambiguity_level == "1", ]
rows <- do.call(rbind, lapply(seq_len(nrow(kept)), function(i) {
  p <- kept[i, ]
  accessions <- strsplit(p$protein_accession, "|", fixed = TRUE)[[1]]
  spans <- per_protein(p$start_and_end_residues_in_parent_sequence, length(accessions))  # "[69 to 84]"
  before <- per_protein(p$previous_residue, length(accessions))
  ends <- do.call(rbind, strsplit(gsub("[][]", "", spans), " to ", fixed = TRUE))
  data.frame(run = p$file_name_without_extension, full_sequence = p$full_sequence,
             protein_accession = accessions, start_residue = as.integer(ends[, 1]),
             end_residue = as.integer(ends[, 2]), previous_residue = ifelse(nchar(before) == 1, before, ""))
}))
rows <- unique(rows)
rows$intensity <- unname(intensity[paste(rows$run, rows$full_sequence, sep = "\r")])
write.table(rows[, c("run", "full_sequence", "protein_accession", "start_residue", "end_residue",
                     "intensity", "previous_residue")],
            "observations.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
```

[`unique()`](https://rdrr.io/r/base/unique.html) keeps one row per (run,
peptidoform, protein) however many PSMs the peptidoform has, and
`na = ""` writes a blank intensity for a peptidoform that was identified
but not quantified. Check that the run names of the two files match:
FlashLFQ labels a calibrated file `<run>-calib`.

## Occupancy, run by run

The cases file holds several of mzLib’s hand-computed cases, each on its
own protein, named for what it tests:

``` r

ptmqtl_data <- system.file("extdata", "ptmqtl", package = "mzLibR")
occ <- ptmqtl_occupancy(file.path(ptmqtl_data, "cases.tsv"))
occ
#> <mzlibr_site_occupancy> 9 sites in 6 runs, 14 cells
#>   quantified 9, floor 1, count_only 3, not_detected 1
r <- occ$records
r[r$protein_accession == "RATIO", c("site_key", "position", "fraction", "unmodified_quantified")]
#>                                           site_key position fraction
#> 3 RATIO:S13:Common Biological:Phosphorylation on S       13     0.15
#>   unmodified_quantified
#> 3                  TRUE
```

`RATIO` has 30 units of the phosphopeptide, 70 of the same peptide
unmodified, and 100 of a missed-cleavage peptide that also covers the
serine: 30 / 200 = 0.15. Every quantified form that covers the position
counts in the denominator, whatever its own sequence.

**Each cell has a state.** Only a `"quantified"` cell has a fraction:

``` r

r[r$protein_accession == "STATES", c("run", "state", "fraction")]
#>           run        state fraction
#> 6   countOnly   count_only       NA
#> 7       floor        floor       NA
#> 8 notDetected not_detected       NA
#> 9       quant   quantified      0.1
```

| state | what was seen | read it as |
|----|----|----|
| `quantified` | the modified form was quantified | a measured fraction |
| `floor` | the modified form was identified but not quantified, while forms covering the site were | low, **not zero** |
| `count_only` | the modified form was identified, but nothing covering the site was quantified | present, amount unknown |
| `not_detected` | forms covering the site were identified (quantified or not), and none carried the modification | a bound below detection, **not zero** |

So never average `fraction` with the `NA`s read as 0: they are censored,
not absent. Each cell also counts the peptidoforms behind it,
identifications rather than intensity: of the `covering_peptidoforms`
identified there, `modified_peptidoforms` carried the modification.

``` r

r[r$protein_accession == "STATES", c("run", "state", "covering_peptidoforms", "modified_peptidoforms")]
#>           run        state covering_peptidoforms modified_peptidoforms
#> 6   countOnly   count_only                     1                     1
#> 7       floor        floor                     2                     1
#> 8 notDetected not_detected                     1                     0
#> 9       quant   quantified                     2                     1
```

**A 1 can be a ceiling.** When a site is seen only in its modified form,
every quantified form covering it carries the modification, so the
fraction is 1. The unmodified form may be there below detection.
`unmodified_quantified` is `FALSE` then:

``` r

r[r$protein_accession == "CEILING", c("fraction", "unmodified_quantified")]
#>   fraction unmodified_quantified
#> 1        1                 FALSE
r[r$protein_accession == "SHARED_POSITION", c("modification", "fraction")]
#>                             modification fraction
#> 4          Common Biological:HexNAc on S      0.3
#> 5 Common Biological:Phosphorylation on S      0.2
```

**Two modifications at one position share a denominator**: on
`SHARED_POSITION`, 20 units carry phosphorylation, 30 HexNAc and 50
neither. **Positions count from 1, in the protein**: a protein
N-terminal modification sits at the residue it is on, residue 2 when the
initiator methionine was removed.
[`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
reads MetaMorpheus’s occupancy file, which puts the N-terminus at
position 0, so do not join the two tables on position.

**What mzLib leaves out.** The cases file also holds rows that make no
site: oxidation and carbamidomethylation (mzLib always skips the
`Common Variable` and `Common Fixed` categories), a peptide N-terminal
modification, and N-acetylation on a peptide starting at residue 2 after
a K. C-terminal modifications are never sites, because the table does
not carry the protein’s length. To keep only some modifications, pass
`modifications` (exact names) or `categories` (for example
`"Common Biological"`).

## A sample measured as several fractions

Fractions of one sample are summed before the ratio, as FlashLFQ sums
them. Give the run-to-sample map as a file:

``` r

fractions <- ptmqtl_occupancy(file.path(ptmqtl_data, "fractions.tsv"),
                              samples = file.path(ptmqtl_data, "fractions_samples.tsv"))
c(runs = fractions$run_count, samples = fractions$sample_count)
#>    runs samples 
#>       3       2
fractions$records[, c("protein_accession", "run", "state", "fraction")]
#>   protein_accession run        state fraction
#> 1      FLOOR_ACROSS   A        floor       NA
#> 2            SUMMED   A   quantified      0.2
#> 3            SUMMED   B not_detected       NA
```

`SUMMED` reads 0.3 in fraction 1 and 0.1 in fraction 2. Sample A is
(30 + 10) / (100 + 100) = 0.2, not the mean of 0.3 and 0.1. On
`FLOOR_ACROSS` the phosphopeptide was identified, unquantified, in one
fraction, and the site was covered in the other: a floor for the sample.
**Map only fractions together.** Technical replicates are separate
measurements of one sample; give each its own name, and let
[`ptmqtl_trait_effects()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md)
treat them as replicates of one biological sample.

## Pairs of modifications

A peptidoform carrying two modifications is direct evidence that both
sit on one molecule:

``` r

physical <- ptmqtl_pairs(file.path(ptmqtl_data, "cases.tsv"), type = "physical")
physical$records[, c("site_a", "site_b", "statistic", "n", "statistic_n", "p_value")]
#>                                              site_a
#> 1 TWO_ON_ONE:K14:Common Biological:Acetylation on K
#>                                                  site_b statistic n statistic_n
#> 1 TWO_ON_ONE:S13:Common Biological:Phosphorylation on S      0.25 2           1
#>   p_value
#> 1      NA
```

`n` is the runs where the doubly modified form was identified (2).
`statistic` is the median **co-occupancy** over the runs where it was
quantified (1): 25 of the 100 units covering both positions carried both
marks. There is no test for physical pairs, so `p_value` and `q` are
`NA`.

Two sites that rise and fall together across runs are co-varying. mzLib
correlates their occupancies with Spearman’s rho, on the runs where both
were quantified, and adjusts the p-values with Benjamini-Hochberg within
each family. `type` is required: mzLib has no default between its two
methods.

``` r

batch <- file.path(ptmqtl_data, "batch.tsv")
plain <- ptmqtl_pairs(batch, type = "covarying")
str(plain[c("covarying_run_count", "min_quantified_fraction", "exclude_ceiling")])
#> List of 3
#>  $ covarying_run_count    : num 12
#>  $ min_quantified_fraction: num 0.7
#>  $ exclude_ceiling        : logi NA
pair_of <- function(result) {
  r <- result$records
  r[startsWith(r$site_a, "PAIR_S") & startsWith(r$site_b, "PAIR_T"), c("statistic", "p_value", "q", "fdr_family", "n")]
}
pair_of(plain)
#>    statistic     p_value           q fdr_family  n
#> 15 0.7622378 0.003950449 0.003950449    A:inter 12
```

`0.7` is mzLib’s default, reported because you did not set it: a site
enters pairs only if it is quantified in at least 70% of the runs.
`exclude_ceiling` is `NA` because you did not choose: mzLib decides it
inside `PtmPairEngine.CoVarying`, and leaves ceilings out. It is `FALSE`
with `keep_ceiling = TRUE`.

A co-varying pair that was not tested says why in `test_status`
(`"not_estimable"`, `"too_few_shared_runs"` or
`"both_sites_mostly_bounds"`; `"tested"` otherwise, `NA` for physical
pairs). `run_level_removed` says whether rho is partial given the run
level; `bounds_included` and `run_covariate_removed` report two mzLib
options (#1468) that
[`ptmqtl_pairs()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_pairs.md)
does not send yet, so here they are `FALSE`.

**The trap: a batch correlates everything.** These twelve runs come in
two batches, and every site reads 0.2 higher in the second. That alone
gives `PAIR_S` and `PAIR_T` a high rho, though within a batch they are
nearly unrelated. `remove_run_level = TRUE` correlates the two sites
given each run’s **level**, the median occupancy of the other sites in
that run:

``` r

given <- ptmqtl_pairs(batch, type = "covarying", remove_run_level = TRUE)
pair_of(given)
#>     statistic  p_value        q fdr_family  n
#> 15 0.02857143 0.933543 0.933543    A:inter 12
c(given$run_level_group, given$min_runs_for_run_level, given$min_sites_for_run_level)
#> [1] "none" "10"   "3"
```

Given the batch, the pair is unrelated. A pair is tested this way only
when it shares at least `min_runs_for_run_level` runs that have a level,
and a run has one only with at least `min_sites_for_run_level` other
sites there. **Remove the run level for a batch or loading shift only**:
if the level follows your design - one condition reads higher
everywhere - removing it removes the effect you are studying. Sites of
different chemistry can have different batch shifts:
`run_level_group = "category"` or `"modification"` gives each group its
own level, and a two-column file (`modification`, `group`) as
`run_level_groups` gives any other grouping.

- **`overlapping` pairs have no `q`.** When one peptide covers both
  sites, their occupancies share intensities and are correlated by
  construction, so mzLib reports them but leaves them out of the
  adjustment.
- **`fdr_family`** is the type (`P` physical, `A` co-varying) and
  whether the two sites are on one protein (`intra`) or two (`inter`).
  `q` is adjusted within it. It is not a target-decoy q-value.
- **`keep_ceiling = TRUE`** keeps runs where a site was seen only
  modified. Two unrelated sites that both read 1 in low-load runs then
  correlate through detection alone; that is why mzLib leaves them out
  by default.

## Does occupancy change with a trait?

[`ptmqtl_trait_effects()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md)
fits, for each site, logit(occupancy) against a per-run trait, with a
random intercept for each biological replicate, so repeated injections
of one animal or donor are not counted as independent evidence. The
example is mzLib’s `TraitDesign`: six replicates at three doses, two
injections each, built so that logit(occupancy) = -2 + 0.5 x dose + a
replicate offset. First compute occupancy, and write it to a file:

``` r

written <- ptmqtl_occupancy(file.path(ptmqtl_data, "trait.tsv"),
                            out = file.path(tempdir(), "trait_occupancy.tsv"))
c(cells_written = written$written$row_count, rows_returned = nrow(written$records))
#> cells_written rows_returned 
#>            24             0
written$state_counts
#>   quantified        floor   count_only not_detected 
#>           23            1            0            0
head(read.delim(file.path(ptmqtl_data, "trait_runs.tsv")), 3)
#>    run replicate dose injection
#> 1 R0_0        R0 0.25         1
#> 2 R0_1        R0 0.25         0
#> 3 R1_0        R1 0.25         1
```

The runs table gives each run its replicate, then the trait, then any
covariates, each column named by its header. Then fit (here on the file
the call above writes, as this package ships it):

``` r

effects <- ptmqtl_trait_effects(file.path(ptmqtl_data, "trait_occupancy.tsv"),
                                file.path(ptmqtl_data, "trait_runs.tsv"))
effects
#> <mzlibr_trait_effects> 2 sites against dose, given injection (min 6 runs, 4 replicates)
#>   fitted 2, below_support 0, too_few_observations 0, rank_deficient 0, too_few_groups 0, variance_ratio_at_limit 0, not_converged 0
effects$records[, c("protein_accession", "status", "effect", "standard_error", "df", "q", "runs", "replicates")]
#>   protein_accession status    effect standard_error df           q runs
#> 1     CENSORED_SITE fitted 0.5428795     0.06507235  3 0.003610195   10
#> 2        TRAIT_SITE fitted 0.5008514     0.05519014  4 0.001634666   12
#>   replicates
#> 1          5
#> 2          6
```

`TRAIT_SITE` recovers the 0.5 it was built with. **`effect` is in logit
units per unit of the trait** - natural-log odds, mzLib’s default scale,
which `effects$scale` reports - so `exp(0.5)` = 1.65: each dose unit
multiplies the odds of the modified form by about 1.65. The degrees of
freedom are between replicates (six replicates, less the intercept and
the dose), because the dose never varies within a replicate.
`CENSORED_SITE` is the same site with one replicate’s two runs replaced
by a floor and a ceiling. Both leave the fit - only quantified,
non-ceiling cells are responses - so it rests on five replicates and
three degrees of freedom.

**Support.** A site needs `min_runs` runs and `min_replicates`
replicates with a response, mzLib’s defaults reported in the result. Ask
for six replicates, and `CENSORED_SITE` is below support:

``` r

strict <- ptmqtl_trait_effects(file.path(ptmqtl_data, "trait_occupancy.tsv"),
                               file.path(ptmqtl_data, "trait_runs.tsv"), min_replicates = 6)
strict$records[, c("protein_accession", "status", "effect", "replicates")]
#>   protein_accession        status    effect replicates
#> 1     CENSORED_SITE below_support        NA          5
#> 2        TRAIT_SITE        fitted 0.5008514          6
```

This answers “does the modified share change?”, not “is the modification
detected more often?”: floors, count-only and undetected cells are left
out, never entered as 0. `q` adjusts over this call’s fitted sites
only - not across datasets, and not a target-decoy q-value.

## What these functions will not do

- **No reader builds the input.** mzLib has none yet, so the table above
  is yours to assemble.
- **No test for physical pairs.** mzLib reports their co-occupancy and
  counts only.
- **No detection model.** For “is the modification detected more
  often?”, fit a logistic model per site with
  [`stats_fit_logistic()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_logistic.md)
  on a 0/1 table of your own.
- **mzLibR does no arithmetic.** Every number on this page came from
  mzLib.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.2307/1412159](https://doi.org/10.2307/1412159): the rank
  correlation coefficient rho of co-varying pairs (Spearman, Am J
  Psychol 15:72, 1904).
- [doi:10.5351/CSAM.2015.22.6.665](https://doi.org/10.5351/CSAM.2015.22.6.665):
  ppcor, whose pcor.test(method = ‘spearman’) the partial correlation
  given the run level reproduces (Kim, Commun Stat Appl Methods 22:665,
  2015).
- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the step-up adjustment behind q, within each fdr_family (Benjamini and
  Hochberg, J R Stat Soc B 57:289, 1995).
- [doi:10.2307/2529876](https://doi.org/10.2307/2529876): the
  random-intercept linear mixed model fitted per site (Laird and Ware,
  Biometrics 38:963, 1982).
