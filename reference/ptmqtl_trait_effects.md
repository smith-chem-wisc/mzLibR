# Fit, per site, the effect of a per-run trait on logit occupancy, with mzLib's SiteTraitEffects

Per site, the effect of a per-run trait on the log odds of occupancy
(natural log or log2), from a random-intercept mixed model per site, so
that repeated runs of one biological replicate are not counted as
independent.

The model, per site, is logit(occupancy) = b0 + b_trait x trait +
covariates + a random intercept per biological replicate + error, fitted
by REML (mzLib's `MixedModel`). The random intercept is the point:
repeated runs of one animal or donor are not counted as independent, so
a trait that is constant within each replicate is tested at the
between-replicate degrees of freedom.

## Usage

``` r
ptmqtl_trait_effects(occupancy, runs, min_runs = NULL, min_replicates = NULL,
  logit_epsilon = NULL, keep_ceiling = FALSE, timeout = 600)
```

## Arguments

- occupancy:

  The table
  [`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md)
  writes with `out` (columns found by header name). A `fraction` cell is
  used as written, as mzLib does for a stored fraction. The peptidoform
  counts the 0.7.0 bridge writes are read if present and not used by the
  fit, so a table written before them reads too.

- runs:

  A TSV with a header row, then one row per run: the run, its biological
  replicate, the trait (its header names it; any unit), then any
  covariates (each header names one). A run with a blank trait is not
  used.

- min_runs:

  Fewest runs with a response for a site to be fitted (runs). `NULL`:
  mzLib's default (6 at mzLib `d10017ce`), which `min_runs` in the
  result reports.

- min_replicates:

  Fewest distinct replicates among them (replicates). `NULL`: mzLib's
  default (4 at mzLib `d10017ce`).

- logit_epsilon:

  Fractions are clamped to `[epsilon, 1 - epsilon]` before the logit
  (fraction, above 0 and below 0.5). `NULL`: mzLib's default (5e-5 at
  mzLib `d10017ce`).

- keep_ceiling:

  Use runs where the site was seen only modified. `FALSE` sends nothing,
  so mzLib's default (leave them out) applies.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

Only quantified cells are responses. Floors, count-only and undetected
cells are left out, never entered as 0, so this answers "does the
modified share change?", not "is the modification detected more often?".
Ceilings are left out by default.

## Value

An `mzlibr_trait_effects`:

\- `records`, a data.frame with one row per site with a usable response,
by site key: `site_key`, `protein_accession`, `position`, `residue` and
`modification`, as
[`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md)
writes them; `status`, one of `"fitted"`, `"below_support"`,
`"too_few_observations"`, `"rank_deficient"`, `"too_few_groups"`,
`"variance_ratio_at_limit"` or `"not_converged"`; `effect` (log odds of
occupancy per trait unit, on the scale `scale` reports; logit occupancy
per trait unit at mzLib's default, `"natural_log_odds"`), the change in
the log odds of occupancy per unit of the trait, so `exp(effect)` is the
odds ratio, `NA` when not fitted; `standard_error` (log odds of
occupancy per trait unit, as `scale`); `df` (degrees of freedom), by the
containment rule; `p_value`, a fraction (0 to 1), two-sided, for the
trait; `q`, a fraction (0 to 1), Benjamini-Hochberg over this call's
fitted sites; `runs`, runs contributing a response; `replicates`,
distinct replicates among them; `median_fraction` (fraction, 0 to 1),
median occupancy over those runs; `covariate_effects` (log odds of
occupancy per covariate unit, as `scale`), a list column, one per
covariate in `covariate_names` order, each `NA` when not fitted;
`replicate_variance` (squared log odds of occupancy, as `scale`),
between replicates; and `residual_variance` (squared log odds of
occupancy, as `scale`), run to run within a replicate. Statistics are
`NA` where `status` is not `"fitted"`. - `occupancy_file` and
`runs_file`, as absolute paths; `trait_name`, the runs table's third
header; `covariate_names`; `cell_count` occupancy rows (cells) read;
`run_count` runs in the runs table; `trait_run_count`, runs with a trait
value (runs); `occupancy_runs_not_used`; `min_runs` (runs) and
`min_replicates` (replicates), the thresholds used, yours or mzLib's
default; `logit_epsilon` (fraction, 0 to 1), the clamp used;
`exclude_ceiling`; `scale`, the unit of every effect, standard error and
variance, `"natural_log_odds"` (mzLib's default, the logit) or
`"log2_odds"` (mzLib \#1464; `NA` from a bridge before pyMzLib 0.7.0);
`status_counts`, sites (sites) per status; `ignored_columns`;
`row_count` rows; `column_names`. - `caveats`: the model and what it
leaves out, plus any that apply to this input.

## What it is not

It says nothing about whether \*detecting\* the modification changes
with the trait: only quantified cells enter. `q` adjusts over this
call's fitted sites only - not across datasets, and not a target-decoy
q-value.

## Wraps

Wire verb `ptmqtl trait-effects`. Generated from the bridge's verb spec
`ptmqtl.trait-effects.yaml` (bridge commit `4e180efe8b59`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SiteTraitEffects.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteTraitEffects.cs)
  in `mzLib/Quantification/PtmQtl/SiteTraitEffects.cs` at mzLib
  `cc67b7a2`

- [`SiteTraitOptions`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteTraitEffects.cs)
  in `mzLib/Quantification/PtmQtl/SiteTraitEffects.cs` at mzLib
  `cc67b7a2`

- [`SiteTraitEffect`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteTraitEffects.cs)
  in `mzLib/Quantification/PtmQtl/SiteTraitEffects.cs` at mzLib
  `cc67b7a2`

- [`OccupancyScale`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteTraitEffects.cs)
  in `mzLib/Quantification/PtmQtl/SiteTraitEffects.cs` at mzLib
  `cc67b7a2`

- [`SiteRunOccupancy`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `cc67b7a2`

- [`MixedModel.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/MixedModel.cs)
  in `mzLib/StatisticalModels/MixedModel.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `occupancy`:

  path; required. The occupancy table ptmqtl occupancy `--out` writes,
  one row per (site, run). Columns are found by header name:
  protein_accession, position (1-based), residue (one letter),
  modification, run, state (quantified, floor, count_only or
  not_detected), modified_intensity and covering_intensity (instrument
  units, finite), unmodified_quantified (true or false), and optionally
  fraction (fraction, 0 to 1; blank for none), which mzLib then uses
  instead of modified_intensity / covering_intensity, as it does for a
  fraction a producer wrote exactly, and optionally
  covering_peptidoforms and modified_peptidoforms (integers; blank for
  not supplied), read into mzLib's cell; the fit does not use them.
  site_key and category are derived by mzLib and not read; any other
  column is listed in ignored_columns. A (site, run) given twice is
  refused, since mzLib's fit would use the first.

- `runs`:

  path; required. A TSV with a header row, then one row per run: the run
  (as the occupancy table names it), its biological replicate (the
  random-intercept group, for example the animal or donor), the trait
  (its header names the trait; any unit, which effect is then per), then
  any covariates (each header names one; every value finite for a run
  with a trait). Columns are by position. A run with a blank, NA or NaN
  trait is not used, as mzLib reads a run absent from the trait. A run
  of the occupancy table missing here is not used either, and is listed
  in occupancy_runs_not_used.

- `min_runs` (wire `--min-runs`):

  int; in **runs**; default
  `mzLib's SiteTraitOptions.MinRuns, 6 at the pin`; range `>= 0`. The
  fewest runs contributing a response for a site to be fitted; fewer is
  below_support. Absent, the bridge sets nothing and mzLib's default
  applies.

- `min_replicates` (wire `--min-replicates`):

  int; in **replicates**; default
  `mzLib's SiteTraitOptions.MinReplicates, 4 at the pin`; range `>= 0`.
  The fewest distinct biological replicates among those runs; fewer is
  below_support. Absent, mzLib's default applies.

- `logit_epsilon` (wire `--logit-epsilon`):

  float; in **fraction (0 to 1)**; default
  `mzLib's SiteTraitOptions.LogitEpsilon, 5e-5 at the pin (half the 4-decimal resolution of stored fractions)`;
  range `> 0 and < 0.5`. Each fraction is clamped to \[epsilon, 1 -
  epsilon\] before the logit, so 0 and 1 have finite logits. Absent,
  mzLib's default applies.

- `keep_ceiling` (wire `--keep-ceiling`):

  flag; default
  `absent: mzLib's SiteTraitOptions.ExcludeCeiling, true at the pin, so ceilings are left out`.
  Use a run where the site was seen only modified (a fraction of 1 with
  unmodified_quantified false) as a response, clamped below 1. mzLib
  leaves such a ceiling out unless asked.

- wire `--scale`:

  string; default
  `absent: mzLib's SiteTraitOptions.Scale, natural_log_odds at the pin`;
  range `natural_log_odds | log2_odds`. The unit the effects and
  variances are reported in (mzLib's OccupancyScale, \#1464).
  natural_log_odds: ln(f / (1 - f)), the logit, so exp(effect) is the
  odds ratio per unit. log2_odds: log2(f / (1 - f)), so an effect of +1
  doubles the odds per unit; mzLib divides the natural-log fit's effects
  and standard errors by ln 2 and its variances by (ln 2)^2, and leaves
  df, p and q unchanged: a change of unit, not a refit. Absent, the
  bridge sets nothing and mzLib's default applies. *Not an argument
  here: not ported yet: the pyMzLib 0.7.0 bridge's option, a later
  mzLibR port; mzLib's default, natural_log_odds, always applies, and
  the result's scale says so.*

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `occupancy_file`:

  string; never `NA`. Absolute path of the occupancy table.

- `runs_file`:

  string; never `NA`. Absolute path of the runs table.

- `trait_name`:

  string; never `NA`. The runs table's third header: the trait.

- `covariate_names`:

  string\[\]; never `NA`. The runs table's headers after the trait, in
  order; empty without covariates. covariate_effects follows this order.

- `cell_count`:

  int; in **cells**; never `NA`. Occupancy rows read.

- `run_count`:

  int; in **runs**; never `NA`. Runs in the runs table.

- `trait_run_count`:

  int; in **runs**; never `NA`. Runs with a trait value: the runs the
  design is built over.

- `occupancy_runs_not_used`:

  string\[\]; never `NA`. Runs of the occupancy table with no trait
  value in the runs table, in ordinal order; their cells were not used.

- `min_runs`:

  int; in **runs**; never `NA`. The support threshold used: as given, or
  mzLib's default.

- `min_replicates`:

  int; in **replicates**; never `NA`. The replicate threshold used: as
  given, or mzLib's default.

- `logit_epsilon`:

  float; in **fraction (0 to 1)**; never `NA`. The clamp used: as given,
  or mzLib's default.

- `exclude_ceiling`:

  bool; never `NA`. Whether ceiling cells were left out: false with
  `--keep-ceiling`, else mzLib's default.

- `scale`:

  string; never `NA`. The scale used, natural_log_odds or log2_odds: as
  given, or mzLib's default. It is the unit of effect, standard_error,
  covariate_effects and the two variances.

- `status_counts`:

  map\<string,int\>; in **sites**; never `NA`. Sites per status, keyed
  by every status value.

- `ignored_columns`:

  string\[\]; never `NA`. Columns of the occupancy file that were not
  read (other than site_key and category); empty when none.

- `row_count`:

  int; in **rows**; never `NA`. Sites in the table.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per site with at least one usable response
  (a quantified cell, not a ceiling unless `--keep-ceiling`, in a run
  with a trait), in ordinal site_key order.

- `caveats`:

  string\[\]; never `NA`. Seven static caveats, one of them naming the
  scale and its units; plus one listing occupancy runs not used, one
  counting sites not fitted by status, and one naming columns not read,
  when they apply.

## Columns of `records`

- `site_key`:

  string; never `NA`. mzLib's ModificationSite.Key, as ptmqtl occupancy
  writes it.

- `protein_accession`:

  string; never `NA`. The site's protein.

- `position`:

  int; never `NA`. 1-based residue position, as ptmqtl occupancy.

- `residue`:

  string; never `NA`. The one-letter residue.

- `modification`:

  string; never `NA`. Category:IdWithMotif.

- `status`:

  string; never `NA`. fitted \| below_support (fewer runs or replicates
  than min_runs or min_replicates) \| too_few_observations \|
  rank_deficient \| too_few_groups \| variance_ratio_at_limit \|
  not_converged (mzLib's MixedModel could not fit it, as stats fit names
  those outcomes). Only fitted rows carry statistics.

- `effect`:

  float; in **log odds of occupancy per trait unit (natural log or log2,
  as scale says)**; `NA` when not fitted (see status). The change in the
  log odds of occupancy per unit of the trait, other covariates held
  fixed. natural_log_odds: exp(effect) is the odds ratio of the modified
  form per trait unit. log2_odds: 2^effect is, so +1 doubles the odds.

- `standard_error`:

  float; in **log odds of occupancy per trait unit (as scale)**; `NA`
  when not fitted. Standard error of effect, on the same scale.

- `df`:

  float; in **degrees of freedom**; `NA` when not fitted. Denominator
  degrees of freedom of the t test, by nlme's containment rule: for a
  trait constant within each replicate, replicates - 1 - the
  between-replicate coefficients.

- `p_value`:

  float; in **fraction (0 to 1)**; `NA` when not fitted. Two-sided
  p-value of the trait coefficient.

- `q`:

  float; in **fraction (0 to 1)**; `NA` when not fitted, so not in the
  family. Benjamini-Hochberg over the sites with a p-value in this call.
  Not a target-decoy q-value.

- `runs`:

  int; in **runs**; never `NA`. Runs contributing a response.

- `replicates`:

  int; in **replicates**; never `NA`. Distinct biological replicates
  among those runs.

- `median_fraction`:

  float; in **fraction (0 to 1)**; never `NA`. Median occupancy over the
  contributing runs, before the clamp.

- `covariate_effects`:

  float\[\]; in **log odds of occupancy per covariate unit (as scale)**;
  never `NA`. One coefficient per covariate, in covariate_names order;
  each element null when the site was not fitted. Empty without
  covariates.

- `replicate_variance`:

  float; in **squared log odds of occupancy (as scale)**; `NA` when not
  fitted. Between-replicate variance of the log odds of occupancy (the
  random intercept).

- `residual_variance`:

  float; in **squared log odds of occupancy (as scale)**; `NA` when not
  fitted. Within-replicate (run-to-run) variance of the log odds of
  occupancy.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  occupancy or runs missing, blank or not found; min-runs or
  min-replicates not an integer or below 0; logit-epsilon not a number,
  or not above 0 and below 0.5; scale not natural_log_odds or log2_odds;
  an option given without a value; keep-ceiling written with a value

- `mzlib_usage_error` (usage):

  the occupancy file empty, header-only, missing a column, with a
  repeated column name, a row of the wrong width, a position that is not
  an integer, a residue that is not one letter, a state that is not one
  of the four, an intensity that is not finite, a fraction that is not a
  number, a peptidoform count that is not blank or an integer,
  unmodified_quantified not true or false, an empty accession,
  modification or run, or a (site, run) twice

- `mzlib_usage_error` (usage):

  the runs file with fewer than three columns, a repeated header after
  the replicate, a row of the wrong width, an empty or repeated run, or
  a trait or covariate that is not a number

- `mzlib_usage_error` (usage):

  mzLib refuses the input, in its own words: a run with a trait and no
  replicate; a run with a trait whose covariates are not all finite; a
  design (intercept, trait, covariates over the runs with a trait) not
  of full column rank, such as a constant trait

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Per site, logit(occupancy) = b0 + b_trait trait + the covariates + a
  random intercept per biological replicate + error, fitted by REML
  (mzLib's MixedModel). Repeated runs of one replicate (technical
  injections) are not counted as independent: a trait constant within
  each replicate is tested at the between-replicate degrees of freedom.

- Only quantified cells are responses; floor, count_only and
  not_detected cells are left out, never entered as 0. This is the
  intensity part of the question: it says nothing about whether
  detecting the modification changes with the trait.

- Ceiling cells are left out unless keep-ceiling.

- Fractions are clamped to \[logit_epsilon, 1 - logit_epsilon\] before
  the logit. effect is per unit of the trait as given: rescale the trait
  to change the unit.

- effect, standard_error, covariate_effects and the two variances are on
  the scale echoed: natural_log_odds (exp(effect) is the odds ratio per
  unit) or log2_odds (an effect of +1 doubles the odds). log2_odds is
  the natural-log fit divided by ln 2 (variances by (ln 2)^2); df,
  p_value and q do not change.

- A site with fewer than min_runs runs or min_replicates replicates is
  below_support, with null statistics; a site with no usable response
  does not appear.

- q adjusts over this call's fitted sites only: not across datasets, and
  not a target-decoy q-value.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.ptmqtl.trait_effects`

- Rust (mzLibRust): `mzlib::ptmqtl::trait_effects_with` with
  `TraitEffectsOptions`

- R (mzLibR): `ptmqtl_trait_effects`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.2307/2529876](https://doi.org/10.2307/2529876): the
  random-intercept linear mixed model fitted per site (Laird and Ware,
  Biometrics 38:963, 1982)

- [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):
  the step-up adjustment behind q (Benjamini and Hochberg, J R Stat Soc
  B 57:289, 1995)

## See also

[`ptmqtl_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md),
[`stats_fit_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)
for the same model on any responses

## Examples

``` r

# mzLib's TraitDesign: six replicates at three doses, two injections each, with
# logit(occupancy) = -2 + 0.5 x dose + a replicate offset. trait_runs.tsv gives each run's
# replicate, dose and an injection covariate; trait_occupancy.tsv is ptmqtl_occupancy(out = )'s.
data <- system.file("extdata", "ptmqtl", package = "mzLibR")
fit <- ptmqtl_trait_effects(file.path(data, "trait_occupancy.tsv"),
  file.path(data, "trait_runs.tsv"))
fit
#> <mzlibr_trait_effects> 2 sites against dose, given injection (min 6 runs, 4 replicates)
#>   fitted 2, below_support 0, too_few_observations 0, rank_deficient 0, too_few_groups 0, variance_ratio_at_limit 0, not_converged 0
fit$records[, c("protein_accession", "status", "effect", "df", "runs", "replicates")]
#>   protein_accession status    effect df runs replicates
#> 1     CENSORED_SITE fitted 0.5428795  3   10          5
#> 2        TRAIT_SITE fitted 0.5008514  4   12          6
```
