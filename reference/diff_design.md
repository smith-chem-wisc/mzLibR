# Build each stratum's design matrix and default contrasts from a sample sheet, with mzLib

Turn a sample sheet and a model (factors with reference levels,
covariates with units, batch, interactions, strata) into mzLib's design
per stratum: the design matrix, its rank and aliased terms, the default
contrasts with their estimability, and the design and contrast tables
stats fit reads.

mzLib's `AnalysisDesign.ByStratum` splits the samples by stratum, and in
each stratum builds one additive model: an intercept, each factor
treatment-coded against its reference level, batch, each covariate as
`(value - centre) / scale`, and any interactions. The matrix equals R's
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) for the
same formula. It then lists `DefaultContrasts` - each level against its
reference, every pair when there is no reference, each covariate's
slope - and judges whether each can be estimated at all from these
samples.

## Usage

``` r
diff_design(samples, factors = NULL, covariates = NULL, interactions = NULL,
  strata = NULL, batch = NULL, out = NULL, timeout = 60)
```

## Arguments

- samples:

  The sample sheet: a TSV with a header row and one row per biological
  sample. First column the sample id; the other columns are available by
  header. A column the model does not name is ignored
  (`unused_columns`).

- factors:

  Condition factors, in column order: a list of
  [`factor_spec`](https://smith-chem-wisc.github.io/mzLibR/reference/factor_spec.md),
  or a character vector of column names (no reference: every pair of
  levels is compared), or a mix in a list.

- covariates:

  Numeric covariates, in column order: a list of
  [`covariate_spec`](https://smith-chem-wisc.github.io/mzLibR/reference/covariate_spec.md),
  each with its stored unit.

- interactions:

  Two-factor interactions, as a list of pairs of factor names, e.g.
  `list(c("condition", "sex"))`.

- strata:

  Columns whose values split the samples into strata, each designed on
  its own (a character vector).

- batch:

  The column holding each sample's batch, fitted as a fixed block.
  `NULL`: no batch.

- out:

  A directory to write `design_<i>.tsv` (for
  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)'s
  `design`) and `contrasts_<i>.tsv` (for its `contrasts`) into, per
  stratum. `NULL`: nothing is written.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_diff_design`:

\- `records`, a data.frame with one row per (stratum, default contrast),
strata in order: `stratum`, the stratum key; `contrast_id`, `c1`, `c2`,
... numbered **once over the whole analysis** (mzLib \#1476), so an id
is the same comparison in every stratum; `label`,
`factor=level vs factor=reference`, or a slope's column; `numerator` and
`denominator`, `factor=level` (`NA` for a slope); `covariate`,
`covariate_unit` (the unit a slope is per, e.g. `decade`) and
`covariate_scale` (e.g. `years/10, centred at 50 years`), `NA` when not
a slope; `estimability`, `"estimable"`, `"not_estimable:confounded"` or
`"not_estimable:rank_deficient"`; and two list columns, `weights`
(contrast weight per design column, unitless), one per entry of the
stratum's `column_names`, and `fitted_weights` (contrast weight per
fitted column, unitless), one per entry of its `fitted_column_names` -
the contrasts file's row - and `NULL` when the contrast is not
estimable. - `strata`, a list of `mzlibr_design_stratum`, named by
stratum key (`"all"` for one stratum, else the stratum columns sorted by
name as `name=value` joined by `;`), each with: `sample_count` samples;
`sample_ids`, the matrix rows, ordinal by sample id; `column_names`,
mzLib's columns (`intercept`, `factor=level` per non-reference level,
`batch=level`, `name (per unit)` for a covariate, `a=level:b=level` for
an interaction); `fitted_column_names`, the columns a fit takes (the
design file's coefficient columns); `rank` (columns); `aliased_terms`,
terms that share a direction with another; `warnings`; `not_run`, a
data.frame of the contrasts this stratum cannot run (mzLib's
`NotRunContrast`): `label`, `reason` (`"not_in_stratum"`) and `id`, the
default contrast's id, which this stratum skips (`NA` when no default
contrast stands for the entry, such as a factor with one level over the
whole analysis); `covariate_centres` (the covariate's stored unit), the
centre subtracted from each covariate; `matrix` and `fitted_matrix`, R
matrices with the sample ids as row names, in coding (0/1 indicators;
covariate columns in scaled units), equal to R's
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html); and
`design_file` and `contrasts_file`, the files written for this stratum,
`NA` without `out` (and `contrasts_file` `NA` when nothing here is
estimable). - `samples_file`, as an absolute path; `sample_count`
samples read, across every stratum; `factors`, `covariates` and
`interactions` as given (a covariate field `NA` was not given, so
mzLib's default applied); `stratum_columns`; `batch_column` (`NA` when
batch is not in the model); `unused_columns`; `stratum_count` strata;
`estimability_counts`, default contrasts (contrasts) per verdict across
strata; `written`, `directory` and `file_count` (files), or `NULL`
without `out`; `row_count` rows; `column_names`. - `caveats`: how the
matrix is coded and what the contrast ids mean, plus counts of contrasts
that cannot be estimated or run.

## When the samples cannot answer

A comparison the samples cannot separate - a condition run in only one
batch - is `"not_estimable:confounded"`, and one that needs a column no
sample informs (an empty interaction cell) is
`"not_estimable:rank_deficient"`. mzLib says so instead of returning a
number. The design file still holds a design that can be fitted:
`fitted_column_names` leaves out the later of two aliased columns.
Contrast ids are numbered once over every sample of the analysis, so
`c1` is the same comparison in every stratum; a stratum without samples
on both sides skips the id and lists it in `not_run`.

## Wraps

Wire verb `diff design`. Generated from the bridge's verb spec
`diff.design.yaml` (bridge commit `45833ca804d6`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`AnalysisDesign.ByStratum`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `cc67b7a2`

- [`NotRunContrast`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/Differential/DesignInputs.cs)
  in `mzLib/Quantification/Differential/DesignInputs.cs` at mzLib
  `cc67b7a2`

- [`AnalysisDesign.Create`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `d10017ce`

- [`AnalysisDesign.DefaultContrasts`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `d10017ce`

- [`AnalysisDesign.Estimability`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `d10017ce`

- [`AnalysisDesign.FittedMatrix`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `d10017ce`

- [`AnalysisDesign.FittedWeights`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `d10017ce`

- [`AnalysisDesign.CovariateCentre`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/AnalysisDesign.cs)
  in `mzLib/Quantification/Differential/AnalysisDesign.cs` at mzLib
  `d10017ce`

- [`AnalysisSample`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/DesignInputs.cs)
  in `mzLib/Quantification/Differential/DesignInputs.cs` at mzLib
  `d10017ce`

- [`FactorSpec`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/DesignInputs.cs)
  in `mzLib/Quantification/Differential/DesignInputs.cs` at mzLib
  `d10017ce`

- [`CovariateSpec`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/DesignInputs.cs)
  in `mzLib/Quantification/Differential/DesignInputs.cs` at mzLib
  `d10017ce`

- [`DesignSpec`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/DesignInputs.cs)
  in `mzLib/Quantification/Differential/DesignInputs.cs` at mzLib
  `d10017ce`

- [`Contrast`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/DesignInputs.cs)
  in `mzLib/Quantification/Differential/DesignInputs.cs` at mzLib
  `d10017ce`

## Parameters: units, ranges and defaults

- `samples`:

  path; required. The sample sheet: a TSV with a header row and one row
  per biological sample. The first column is the sample id (unique, any
  header). Every other column is available to the model by its header; a
  column the model does not name is ignored and listed in
  unused_columns. Cells are trimmed. A factor, stratum or batch cell is
  used as text, exactly as written. A covariate cell is a number in the
  covariate's stored unit; a blank, NA or NaN covariate cell is no
  value, and mzLib refuses that sample by name.

- wire `--stdin`:

  string\[\]; required; range
  `one model term per line; fields tab-separated`. The model, one term
  per line, its fields separated by tabs; an empty optional field is not
  given, so mzLib's default applies. factor NAME, optionally then
  REFERENCE: a condition factor treatment-coded against REFERENCE
  (FactorSpec); with no REFERENCE every pair of its levels is compared
  and mzLib warns. covariate NAME UNIT, optionally then SCALE,
  SCALED_UNIT and CENTRE: a numeric column entered as (value - CENTRE) /
  SCALE (CovariateSpec). UNIT is the stored unit (years); SCALE is
  stored units per scaled unit (10), default 1; SCALED_UNIT names the
  unit a slope is per (decade), default SCALE stored units; CENTRE is in
  the stored unit, default the stratum's mean. interaction FIRST SECOND:
  the two-factor interaction of two named factors. stratum NAME: a
  column whose value splits the samples into strata, each designed on
  its own. batch NAME: the column holding each sample's batch, fitted as
  a fixed block whenever it is named (mzLib's DesignSpec.FitBatch
  default; the bridge never sets it). Factors and covariates enter in
  the order given. A column has one role, and is named once. *Not an
  argument here: the model is the factors, covariates, interactions,
  strata and batch arguments, sent one stdin line per term.*

- `out`:

  path; default absent. A directory (created if absent) to write the
  pivot files into: for stratum i in the order of strata, design_i.tsv
  and, when any contrast is estimable there, contrasts_i.tsv. Files are
  numbered, not named by stratum key, because a key can hold any
  character a sheet cell can; each strata entry names its files. Files
  of the same name are overwritten. Absent, nothing is written.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `samples_file`:

  string; never `NA`. Absolute path of the sample sheet.

- `sample_count`:

  int; in **samples**; never `NA`. Samples read from the sheet, across
  every stratum.

- `factors`:

  object\[\]; never `NA`. The factors given, in order: name and
  reference_level (null when none was given, so every pair of levels is
  compared).

- `covariates`:

  object\[\]; never `NA`. The covariates given, in order: name, unit,
  and scale, scaled_unit and centre as given (each null when not given,
  so mzLib's default applied).

- `interactions`:

  object\[\]; never `NA`. The interactions given, in order, each as
  first and second.

- `stratum_columns`:

  string\[\]; never `NA`. The stratum columns given; empty for one
  stratum.

- `batch_column`:

  string; `NA` when no batch line: batch is not in the model. The batch
  column given.

- `unused_columns`:

  string\[\]; never `NA`. Sheet columns after the sample id that the
  model does not name, in sheet order.

- `stratum_count`:

  int; in **strata**; never `NA`. Designs built: one per distinct
  stratum key, or 1.

- `strata`:

  object\[\]; never `NA`. One object per stratum, in mzLib's order
  (ordinal by stratum key); fields under result.tables.strata.

- `estimability_counts`:

  map\<string,int\>; in **contrasts**; never `NA`. Default contrasts per
  verdict across strata: estimable, not_estimable:confounded,
  not_estimable:rank_deficient.

- `written`:

  object; `NA` when no `--out`: nothing was written. directory
  (absolute) and file_count, the design and contrast files written.

- `row_count`:

  int; in **rows**; never `NA`. Default contrasts across strata: the
  table's row count.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (stratum, default contrast): strata in
  order, each stratum's contrasts in mzLib's order.

- `caveats`:

  string\[\]; never `NA`. Five static caveats, plus one when any
  contrast is not estimable (with counts) and one when any stratum
  cannot run a contrast the model asks for.

## Columns of `records`

- `stratum`:

  string; never `NA`. The stratum key (AnalysisDesign.StratumKey): the
  stratum columns sorted by name, each name=value, joined by ';', with
  %, ; and = in a value escaped as %25, %3B and %3D; all when there is
  no stratum column.

- `contrast_id`:

  string; never `NA`. c1, c2, ... in mzLib's default order, numbered
  once over every sample of the analysis (AnalysisDesign.ByStratum,
  since mzLib \#1476), so an id is the same comparison in every stratum.
  The name a contrasts file gives the contrast.

- `label`:

  string; never `NA`. Readable: factor=level vs factor=reference for a
  level, the covariate's column name for a slope (age (per decade)).

- `numerator`:

  string; `NA` when a slope. The numerator level as factor=level.

- `denominator`:

  string; `NA` when a slope. The reference (or earlier) level as
  factor=level.

- `covariate`:

  string; `NA` when not a slope. For a slope, the covariate.

- `covariate_unit`:

  string; `NA` when not a slope. For a slope, the unit the effect is
  per: SCALED_UNIT, or SCALE stored units (10 years), or the stored unit
  at scale 1.

- `covariate_scale`:

  string; `NA` when not a slope. For a slope, mzLib's statement of how
  the stored covariate became that unit, e.g. years/10, centred at 50
  years.

- `estimability`:

  string; never `NA`. mzLib's verdict: estimable \|
  not_estimable:confounded (the direction it needs is shared by two or
  more terms, e.g. a condition only ever in one batch) \|
  not_estimable:rank_deficient (it needs a column no sample informs,
  e.g. an empty interaction cell). Judged on the columns divided by
  their lengths, so the same in any unit.

- `weights`:

  float\[\]; in **contrast weight per design column (unitless)**; never
  `NA`. One weight per entry of the stratum's column_names: the contrast
  is the weighted sum of the coefficients.

- `fitted_weights`:

  float\[\]; in **contrast weight per fitted column (unitless)**; `NA`
  when not estimable: no fit gives it an estimate.
  AnalysisDesign.FittedWeights: the weights on the stratum's
  fitted_column_names, the row the contrasts file writes. A fit of the
  fitted matrix gives with them the estimate and standard error every
  least-squares fit of the full matrix gives.

## Fields of `strata`

a list of `mzlibr_design_stratum`, named by stratum key.

- `stratum`:

  string; never `NA`. The stratum key, as the table's stratum column.

- `sample_count`:

  int; in **samples**; never `NA`. Samples in the stratum.

- `sample_ids`:

  string\[\]; never `NA`. The rows of matrix and of the design file,
  ordinal by sample id.

- `column_names`:

  string\[\]; never `NA`. mzLib's columns: intercept; factor=level for
  each non-reference level; batch=level; name (per unit) for a
  covariate; a=level:b=level for an interaction.

- `fitted_column_names`:

  string\[\]; never `NA`. AnalysisDesign.FittedColumns by name: every
  column at full rank, else each column kept unless it is redundant with
  those kept before it. The design file's coefficient columns.

- `rank`:

  int; in **columns**; never `NA`. The numerical rank, by LinearModel's
  own full-rank test: the number of fitted columns.

- `aliased_terms`:

  string\[\]; never `NA`. Terms (factor names, batch, covariate names,
  intercept) that share a direction with another term and so cannot be
  separated from it; empty at full rank.

- `warnings`:

  string\[\]; never `NA`. What mzLib left out of the model or changed
  from the specification here, in plain words: a factor with one level,
  a reference level with no sample, one batch, no reference declared.

- `not_run`:

  object\[\]; never `NA`. Contrasts the model asks for that this stratum
  cannot run (mzLib's NotRunContrast), each with id (the default
  contrast's id, the same in every stratum; null when no default
  contrast stands for the entry, such as a factor with one level over
  the whole analysis), label and reason (not_in_stratum).

- `covariate_centres`:

  map\<string,float\>; in **the covariate's stored unit**; never `NA`.
  The centre subtracted from each covariate here: the one given, else
  the stratum's mean.

- `matrix`:

  float\[\]\[\]; in **coding (0/1 indicators; covariate columns in
  scaled units)**; never `NA`. AnalysisDesign.Matrix, one row per
  sample_ids entry and one value per column_names entry. Equals R's
  stats::model.matrix for the same formula and level order.

- `fitted_matrix`:

  float\[\]\[\]; in **coding (0/1 indicators; covariate columns in
  scaled units)**; never `NA`. AnalysisDesign.FittedMatrix: matrix with
  only the fitted columns. It has full column rank, so LinearModel.Fit
  takes it; it is the design file's body.

- `design_file`:

  string; `NA` when no `--out`. Absolute path of design_i.tsv: header
  sample then fitted_column_names, then one row per sample, the fitted
  matrix in round-trip decimal. Exactly what stats fit `--design` reads.

- `contrasts_file`:

  string; `NA` when no `--out`, or no contrast is estimable in this
  stratum. Absolute path of contrasts_i.tsv: header contrast then
  fitted_column_names (the design file's coefficient columns, same
  order), then one row per estimable contrast: contrast_id and
  fitted_weights. Each row is the ContrastWeights(name, weights) mzLib's
  EmpiricalBayes.ModerateContrasts takes for a fit of the design file:
  exactly what stats fit `--contrasts` reads, beside `--design`
  design_i.tsv (pyMzLib 0.6.0).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--samples` missing, blank or not found; `--out` given without a
  directory; the sheet empty, header-only, with an empty or repeated
  column name, an empty or repeated sample id (repeats across strata
  included), or a row whose field count differs from its header; a
  covariate cell that is not a number (other than blank, NA, NaN)

- `mzlib_usage_error` (usage):

  a stdin line that is not factor, covariate, interaction, stratum or
  batch, or has the wrong number of fields; an empty name or unit; a
  SCALE or CENTRE that is not a finite number; a column named twice or
  in two roles; a second batch line; an interaction given twice; a named
  column the sheet does not have

- `mzlib_usage_error` (usage):

  mzLib refuses the input (its ArgumentException, reclassified): a
  sample with no value for a factor or covariate, a non-positive scale,
  batch on some samples only, an interaction naming a factor that is not
  in the model

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Each stratum's matrix is mzLib's AnalysisDesign: an intercept; each
  factor treatment-coded against its reference level, its other levels
  in ordinal order; batch against its first level; each covariate as
  (value - centre) / scale; then any interactions. It equals R's
  stats::model.matrix for the same formula and level order.

- Only mzLib's default contrasts are listed
  (AnalysisDesign.DefaultContrasts): each level against its factor's
  reference, every pair of levels when a factor has no reference, and
  each covariate's slope. With an interaction in the model, a factor's
  contrast is its effect at the reference level of the other factor.

- Contrast ids (c1, c2, ...) are numbered once over every sample of the
  analysis (AnalysisDesign.ByStratum), so an id names the same
  comparison in every stratum. A stratum without samples on both sides
  of a contrast skips its id and lists it in not_run.

- A design file holds only the fitted columns (fitted_column_names):
  mzLib leaves out a column that is redundant with earlier ones, because
  LinearModel.Fit refuses a design without full column rank. A contrast
  that is not estimable has no fitted weights and is not in the
  contrasts file.

- Strata are never columns: each stratum gets its own design, centres
  and fit; only the default contrasts' ids are shared. The individual a
  sample came from is never a fixed column either; it is for a mixed
  model's random intercept.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.diff.design`

- Rust (mzLibRust): `mzlib::diff::design_with` with `DesignOptions`

- R (mzLibR): `diff_design`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Only DefaultContrasts is exposed. mzLib's AnalysisDesign.Between (any
  two levels) and Custom (weights by column name) have no wire form yet;
  with no reference declared, DefaultContrasts already gives every pair.

- AnalysisSample.FromExperimentalDesign is not reachable: mzLib has no
  reader from MetaMorpheus's ExperimentalDesign.tsv to an
  IExperimentalDesign (MetaMorpheus's EngineLayer.ExperimentalDesign has
  one), so this verb starts from a sample sheet.

## References

- [doi:10.2307/2346786](https://doi.org/10.2307/2346786): the
  model-formula notation (main effects and a:b interactions) that R's
  stats::model.matrix implements and mzLib's design matrix equals
  (Wilkinson and Rogers, Appl Stat 22:392, 1973)

## See also

[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md),
which reads the files `out` writes;
[`diff_observations`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)

## Examples

``` r

# The MALAT1 dilution series stats_fit() analyses, described by each run's load rather than by a
# hand-made 0/1 design:
samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
design <- diff_design(samples, factors = list(factor_spec("malat_load", "500ng")),
  out = file.path(tempdir(), "malat_design"))
design
#> <mzlibr_diff_design> 27 samples in 1 stratum(s), 2 default contrast(s)
#>   estimable 2, not_estimable:confounded 0, not_estimable:rank_deficient 0
#>   all: 27 samples, columns intercept, malat_load=125ng, malat_load=250ng
stratum <- design$strata[["all"]]
stratum$fitted_column_names
#> [1] "intercept"        "malat_load=125ng" "malat_load=250ng"
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
#> 
```
