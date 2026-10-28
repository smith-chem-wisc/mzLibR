# Find PTM-PTM pairs within one dataset, with mzLib's PtmPairEngine

PTM-PTM pairs within one dataset: two modifications identified together
on one peptidoform (physical), and two sites whose occupancies rise and
fall together across runs (co-varying), with Benjamini-Hochberg
q-values.

Two kinds of evidence, which you choose between with `type`:

## Usage

``` r
ptmqtl_pairs(observations, type, samples = NULL, modifications = NULL,
  categories = NULL, min_quantified_fraction = NULL, keep_ceiling = FALSE,
  remove_run_level = FALSE, run_level_group = NULL, run_level_groups = NULL,
  timeout = 600)
```

## Arguments

- observations:

  [`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md)'s
  observations table, with the same rules.

- type:

  `"physical"`, `"covarying"` or `"both"` (physical rows first).
  Required: mzLib has no default between its two methods.

- samples:

  [`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md)'s
  run-to-sample table; pairs are then across samples.

- modifications:

  As
  [`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md):
  only these modifications are sites.

- categories:

  As
  [`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md):
  only these categories are sites.

- min_quantified_fraction:

  Co-varying only. A site enters pairs only if quantified in at least
  this fraction of the runs (fraction of the runs, above 0 and at most
  1), rounded up. `NULL`: mzLib's default (0.7 at mzLib `d10017ce`),
  which `min_quantified_fraction` in the result reports.

- keep_ceiling:

  Co-varying only. Keep runs where a site was seen only modified.
  `FALSE` sends nothing, so mzLib's default (leave them out) applies.

- remove_run_level:

  Co-varying only. Correlate given each run's level - the median
  occupancy of the other tested sites of the pair's group - to remove a
  batch or loading shift. **Batch only**: where the level follows the
  design, it removes the effect you study.

- run_level_group:

  With `remove_run_level`: `"category"` or `"modification"`, for one
  level per group of sites. `NULL`: one level for every site, mzLib's
  default. Not with `run_level_groups`.

- run_level_groups:

  With `remove_run_level`: a TSV with a header row, then `modification`,
  `group`, for any grouping, such as chemistry. Every modification among
  the sites needs a group.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

\- **physical** (`PtmPairEngine.Physical`): both modifications were
identified on one peptidoform, so on one molecule. `statistic` is the
median co-occupancy over the runs where both were quantified. There is
no test. - **covarying** (`PtmPairEngine.CoVarying`): the two sites'
occupancies rise and fall together (or oppositely) across runs.
`statistic` is Spearman's rho on the runs where both are quantified,
with a Benjamini-Hochberg `q` per family. Floors, undetected cells and,
by default, ceilings are left out: two unrelated sites both read 1 in
low-load runs, which would correlate them through detection alone.

Occupancy is computed inside by the same mzLib call
[`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md)
makes.

## Value

An `mzlibr_ptm_pairs`:

\- `records`, a data.frame with one row per pair of sites, physical
pairs first: `result_type`, `"physical"` (both on one peptidoform) or
`"covarying"` (across runs); `site_a` and `site_b`, the two sites' keys,
`site_a` first by ordinal order; `same_protein`; `overlapping`, one
peptide covers both (such a co-varying pair is correlated by
construction); `statistic` (physical: fraction (0 to 1), the median
co-occupancy; covarying: unitless correlation (-1 to 1), Spearman's
rho), `NA` when physical pairs were never quantified together, or not
estimable; `p_value`, a fraction (0 to 1), two-sided Spearman, `NA` for
physical (no test) or not estimable; `q`, a fraction (0 to 1),
Benjamini-Hochberg within `fdr_family`, `NA` when not tested (physical,
overlapping, not estimable), not a target-decoy q-value; `fdr_family`,
`P` or `A` then `intra` or `inter` protein, e.g. `"A:inter"`; `n`
(runs), physical: runs identified together, covarying: runs where both
were quantified; `statistic_n` (runs), the runs `statistic` was computed
on; `spearman_method`, `"exact"`, `"asymptotic"` or `"not_estimable"`,
`NA` for physical; and `run_level_removed`. - `observations_file`,
`samples_file` and `run_level_groups_file` (`NA` when not given), as
absolute paths; `type`; `observation_count` rows read; `run_count` runs;
`sample_count` samples, `NA` without `samples`; `run_names`; and, `NA`
for `type = "physical"`: `covarying_run_count` (runs), the runs where
any site was identified or covered, the denominator of
`min_quantified_fraction` (fraction of the runs, 0 to 1), the guard
used; `exclude_ceiling` and `remove_run_level`; `run_level_group`,
`"none"`, `"category"`, `"modification"` or `"file"`, `NA` when the run
level was not removed; `min_runs_for_run_level` (runs) and
`min_sites_for_run_level` (sites), mzLib's thresholds, `NA` when not
removed; `modifications` and `categories`, `NULL` for no filter;
`ignored_columns`; `physical_count` pairs (`NA` for `"covarying"`) and
`covarying_count` pairs (`NA` for `"physical"`); `row_count` rows;
`column_names`. - `caveats`: what each pair type means, and the run
level when it was removed, plus any that apply to this input.

## Wraps

Wire verb `ptmqtl pairs`. Generated from the bridge's verb spec
`ptmqtl.pairs.yaml` (bridge commit `60b1d25968a1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PtmPairEngine.Physical`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/PtmPairs.cs)
  in `mzLib/Quantification/PtmQtl/PtmPairs.cs` at mzLib `d10017ce`

- [`PtmPairEngine.CoVarying`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/PtmPairs.cs)
  in `mzLib/Quantification/PtmQtl/PtmPairs.cs` at mzLib `d10017ce`

- [`PtmPairEngine.DefaultMinQuantifiedFraction`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/PtmPairs.cs)
  in `mzLib/Quantification/PtmQtl/PtmPairs.cs` at mzLib `d10017ce`

- [`PtmPairEngine.MinRunsForRunLevel`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/PtmPairs.cs)
  in `mzLib/Quantification/PtmQtl/PtmPairs.cs` at mzLib `d10017ce`

- [`PtmPairEngine.MinSitesForRunLevel`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/PtmPairs.cs)
  in `mzLib/Quantification/PtmQtl/PtmPairs.cs` at mzLib `d10017ce`

- [`PtmPair`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/PtmPairs.cs)
  in `mzLib/Quantification/PtmQtl/PtmPairs.cs` at mzLib `d10017ce`

- [`SiteOccupancyCalculator.Calculate`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `d10017ce`

- [`SiteOccupancyCalculator.CombineObservations`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `observations`:

  path; required. ptmqtl occupancy's observations table, with the same
  rules: a header row naming run, full_sequence, protein_accession,
  start_residue, end_residue (1-based), intensity (instrument units;
  blank, NA or NaN is identified but not quantified) and optionally
  previous_residue; one row per (run, peptidoform, protein), MS/MS
  identifications only.

- `type`:

  string; required; range `physical | covarying | both`. Which pairs:
  physical (PtmPairEngine.Physical, mzLib's type P), covarying
  (PtmPairEngine.CoVarying, type A), or both, physical rows first.
  Required: they are different evidence, and mzLib has no default
  between its two methods.

- `samples`:

  path; default absent. Optional. ptmqtl occupancy's run-to-sample table
  (a header row, then run, sample). With it the observations are first
  combined per sample by mzLib's CombineObservations, so both pair types
  are computed across samples rather than fractions. Map only fractions
  together. Absent, none is sent: pairs are across runs.

- `modifications`:

  string; default absent; range
  `exact names, tab-separated, none blank`. Optional. As ptmqtl
  occupancy: count only these modifications (Category:IdWithMotif,
  separated by tabs) as sites. Not with categories. Absent, no filter.

- `categories`:

  string; default absent; range
  `category names, tab-separated, none blank`. Optional. As ptmqtl
  occupancy: count only modifications of these categories. Not with
  modifications. Absent, no filter.

- `min_quantified_fraction` (wire `--min-quantified-fraction`):

  float; in **fraction of the runs (0 to 1)**; default
  `mzLib's PtmPairEngine.DefaultMinQuantifiedFraction, 0.7 at the pin`;
  range `> 0 and <= 1`. Co-varying only. A site enters co-varying pairs
  only if its occupancy is quantified (a ceiling counts only with
  keep-ceiling) in at least this fraction of the runs, rounded up: the
  runs (or samples) where any site was identified or covered, reported
  as covarying_run_count. Absent, the bridge sends none and mzLib's
  default applies.

- `keep_ceiling` (wire `--keep-ceiling`):

  flag; default
  `absent: mzLib's excludeCeiling default, true at the pin, so ceilings are left out`.
  Co-varying only. Keep a run where a site was seen only modified (a
  fraction of 1 with unmodified_quantified false). Two unrelated sites
  both read 1 in low-load runs, so ceilings can correlate sites through
  detection alone; mzLib leaves them out unless asked.

- `remove_run_level` (wire `--remove-run-level`):

  flag; default
  `absent: mzLib's removeRunLevel default, false at the pin`. Co-varying
  only. Correlate given each run's level: the partial Spearman
  correlation with the run's median occupancy over the other tested
  sites of the pair's group (the pair's own two left out). Removes a
  batch or loading shift that moves many sites in the same runs. Batch
  only: where the run level follows the design, removing it tests the
  correlation within conditions instead.

- `run_level_group` (wire `--run-level-group`):

  string; default
  `absent: mzLib's runLevelGroup default, null, one group of every site`;
  range `category | modification`. With remove-run-level only. Whose run
  level applies to a site: category, one level per modification category
  (mzLib's ModificationSite.Category); or modification, one level per
  modification name. A pair within one group is correlated given that
  group's level, a pair across two groups given both. Not with
  run-level-groups.

- `run_level_groups` (wire `--run-level-groups`):

  path; default absent. With remove-run-level only. A TSV with a header
  row, then modification, group: one row per modification name, mapping
  it to its group, for example every phosphorylation to phospho. Every
  modification among the sites needs a group. The wire form of mzLib's
  runLevelGroup for any grouping, such as chemistry. Not with
  run-level-group.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `observations_file`:

  string; never `NA`. Absolute path of the observations table.

- `samples_file`:

  string; `NA` when no `--samples`: pairs are across runs. Absolute path
  of the run-to-sample table.

- `type`:

  string; never `NA`. physical, covarying or both, as asked.

- `observation_count`:

  int; in **rows**; never `NA`. Observation rows read.

- `run_count`:

  int; in **runs**; never `NA`. Distinct runs in the observations.

- `sample_count`:

  int; in **samples**; `NA` when no `--samples`. Distinct samples in the
  run-to-sample table.

- `run_names`:

  string\[\]; never `NA`. The runs, or with `--samples` the samples, the
  pairs were computed over, in ordinal order.

- `covarying_run_count`:

  int; in **runs**; `NA` when type physical: no co-varying pairs
  computed. The runs (or samples) where any site was identified or
  covered: the denominator of min_quantified_fraction.

- `min_quantified_fraction`:

  float; in **fraction of the runs (0 to 1)**; `NA` when type physical.
  The guard used: as given, or mzLib's default.

- `exclude_ceiling`:

  bool; `NA` when type physical. Whether ceilings were left out: false
  with `--keep-ceiling`, else mzLib's default.

- `remove_run_level`:

  bool; `NA` when type physical. Whether co-varying pairs were
  correlated given the run level.

- `run_level_group`:

  string; `NA` when the run level was not removed (or type physical).
  none (one group of every site, mzLib's default), category,
  modification, or file (`--run-level-groups`).

- `run_level_groups_file`:

  string; `NA` when no `--run-level-groups`. Absolute path of the
  modification-to-group table.

- `min_runs_for_run_level`:

  int; in **runs**; `NA` when the run level was not removed. mzLib's
  PtmPairEngine.MinRunsForRunLevel: with the run level removed, a pair
  sharing fewer runs with a level is not_estimable.

- `min_sites_for_run_level`:

  int; in **sites**; `NA` when the run level was not removed. mzLib's
  PtmPairEngine.MinSitesForRunLevel: a run whose group has fewer other
  sites with a value there has no level.

- `modifications`:

  string\[\]; `NA` when no `--modifications` filter. The modification
  names the sites were limited to.

- `categories`:

  string\[\]; `NA` when no `--categories` filter. The categories the
  sites were limited to.

- `ignored_columns`:

  string\[\]; never `NA`. Columns of the observations file that were not
  read; empty when every column was read.

- `physical_count`:

  int; in **pairs**; `NA` when type covarying. Physical pairs found.

- `covarying_count`:

  int; in **pairs**; `NA` when type physical. Co-varying pairs tested or
  reported.

- `row_count`:

  int; in **rows**; never `NA`. Pairs in the table.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per pair, physical first, each type in
  mzLib's order (by site_a, then site_b).

- `caveats`:

  string\[\]; never `NA`. What each pair type means for the types
  computed, and the run level when removed; plus one with `--samples`,
  one when columns were not read, and one counting not_estimable pairs
  when there are any.

## Columns of `records`

- `result_type`:

  string; never `NA`. physical (mzLib's PairResultType P: both
  modifications identified on one peptidoform, so on one molecule) \|
  covarying (type A: occupancies correlated across runs).

- `site_a`:

  string; never `NA`. The first site's key, as ptmqtl occupancy's
  site_key; ordinally before site_b. A pair is written once.

- `site_b`:

  string; never `NA`. The second site's key.

- `same_protein`:

  bool; never `NA`. Both sites on the same protein accession.

- `overlapping`:

  bool; never `NA`. Some peptide covers both sites: one span on one
  protein, or one peptide (by base sequence) mapped to both proteins.
  Always true for physical. A co-varying overlapping pair shares
  peptides, so it is correlated by construction, and is left out of its
  family's adjustment.

- `statistic`:

  float; in **physical: fraction (0 to 1), co-occupancy; covarying:
  unitless correlation (-1 to 1)**; `NA` when physical: never quantified
  together in any run; covarying: not_estimable. physical: median over
  runs of the co-occupancy, the intensity of peptidoforms carrying both
  modifications over the intensity of peptidoforms covering both
  positions. covarying: Spearman's rho of the two occupancies across
  runs; with the run level removed, the partial rho given it.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when physical (no test in this
  mzLib), or not_estimable. Two-sided Spearman p-value against rho = 0,
  computed as spearman_method says.

- `q`:

  float; in **fraction (0 to 1)**; `NA` when not tested: physical,
  overlapping, or not_estimable. Benjamini-Hochberg adjusted p within
  fdr_family, over that family's non-overlapping pairs. Not a
  target-decoy q-value.

- `fdr_family`:

  string; never `NA`. The family the adjustment runs over, mzLib's
  PtmPair.FdrFamily verbatim: P or A (physical or co-varying), then
  intra (same protein) or inter, for example A:inter.

- `n`:

  int; in **runs**; never `NA`. physical: runs where a doubly modified
  peptidoform was identified. covarying: runs where both occupancies
  were quantified (with the run level removed, and had a level).

- `statistic_n`:

  int; in **runs**; never `NA`. Runs statistic was computed on.
  physical: runs where the doubly modified form and the forms covering
  both positions were quantified, at most n. covarying: equal to n.

- `spearman_method`:

  string; `NA` when physical: no correlation. exact (the full
  permutation distribution, for 9 or fewer runs without ties, plain
  only) \| asymptotic (t on n - 2 degrees of freedom, or n - 2 - k given
  k run levels: one, or two for a pair across groups, less any level the
  other carries entirely) \| not_estimable (too few runs, a constant
  occupancy, or with the run level removed fewer than
  min_runs_for_run_level runs).

- `run_level_removed`:

  bool; never `NA`. Whether statistic and p_value are the partial
  correlation given the run level. Always false for physical.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  observations missing, blank or not found; type missing or not
  physical, covarying or both; samples or run-level-groups given without
  a value or not found; modifications and categories both given, or
  either blank

- `mzlib_usage_error` (usage):

  min-quantified-fraction, keep-ceiling, remove-run-level,
  run-level-group or run-level-groups with type physical;
  run-level-group or run-level-groups without remove-run-level; both of
  them; run-level-group other than category or modification;
  min-quantified-fraction not above 0 and at most 1, or not a number;
  keep-ceiling or remove-run-level written with a value

- `mzlib_usage_error` (usage):

  a table malformed as ptmqtl occupancy describes; the run-level-groups
  file without two fields per row, with a blank or repeated
  modification, or with no group for a modification among the sites

- `mzlib_usage_error` (usage):

  mzLib refuses the observations, in its own words, as ptmqtl occupancy
  lists

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- physical: both modifications on one peptidoform, so on one molecule.
  statistic is the median co-occupancy over the statistic_n runs where
  both were quantified. There is no test in this mzLib: p_value and q
  are null.

- covarying: Spearman's rho across runs on the runs where both sites are
  quantified. Floors and not_detected cells are not values and are left
  out, and so, unless keep-ceiling, are ceilings.

- A site enters co-varying pairs only if quantified in at least
  min_quantified_fraction of covarying_run_count, rounded up.

- q is Benjamini-Hochberg within fdr_family. An overlapping pair is
  reported with q null: shared peptides correlate it by construction. q
  is not a target-decoy q-value.

- With remove-run-level, a run's level is the median occupancy there of
  the other tested sites of the pair's group; a run with fewer than
  min_sites_for_run_level such sites has no level, and a pair sharing
  fewer than min_runs_for_run_level runs with a level is not_estimable
  (the partial test's t on n - 3 df is too liberal there). Use it for a
  batch or loading shift only: where the run level follows the design,
  it removes the effect under study.

- The observations follow ptmqtl occupancy's rules and exclusions: MS/MS
  identifications only, one row per (run, peptidoform, protein).

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.ptmqtl.pairs`

- Rust (mzLibRust): `mzlib::ptmqtl::pairs_with` with `PairsOptions`

- R (mzLibR): `ptmqtl_pairs`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## References

- [doi:10.2307/1412159](https://doi.org/10.2307/1412159): the rank
  correlation coefficient rho of co-varying pairs (Spearman, Am J
  Psychol 15:72, 1904)

- [doi:10.5351/CSAM.2015.22.6.665](https://doi.org/10.5351/CSAM.2015.22.6.665):
  ppcor, whose pcor.test(method = 'spearman') the partial correlation
  given the run level reproduces (Kim, Commun Stat Appl Methods 22:665,
  2015)

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the step-up adjustment behind q, within each fdr_family (Benjamini and
  Hochberg, J R Stat Soc B 57:289, 1995)

## See also

[`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md),
[`stats_correlate`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_correlate.md)

## Examples

``` r

# Two modifications on one molecule, mzLib's PhysicalPairsCarryCoOccupancy case:
cases <- system.file("extdata", "ptmqtl", "cases.tsv", package = "mzLibR")
found <- ptmqtl_pairs(cases, type = "physical")
found$records[, c("site_a", "statistic", "n", "statistic_n", "p_value")]
#>                                              site_a statistic n statistic_n
#> 1 TWO_ON_ONE:K14:Common Biological:Acetylation on K      0.25 2           1
#>   p_value
#> 1      NA

# Twelve runs in two batches, every site 0.2 higher in the second. The batch alone correlates
# two sites; given each run's level, they are unrelated:
batch <- system.file("extdata", "ptmqtl", "batch.tsv", package = "mzLibR")
plain <- ptmqtl_pairs(batch, type = "covarying")
given <- ptmqtl_pairs(batch, type = "covarying", remove_run_level = TRUE)
pick <- function(r) r$statistic[startsWith(r$site_a, "PAIR_S") & startsWith(r$site_b, "PAIR_T")]
round(c(plain = pick(plain$records), given = pick(given$records)), 2)
#> plain given 
#>  0.76  0.03 
```
