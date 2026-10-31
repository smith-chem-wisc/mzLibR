# Fit the detection link: how a cell's latent value sets the chance that a trial succeeds

Fit the detection link of a latent-scale mixed model (how a cell's
latent value, such as the log-odds of a site's occupancy, sets the
chance that a trial succeeds) by maximum likelihood over cells with a
latent mean and SD and a count.

Answers **how does a count's chance of success follow the latent
value?** mzLib's `LogitNormalMixedModel.FitLink` (mzLib \#1469) finds
the `DetectionLink` - logit g = intercept + slope x latent + depth_slope
x depth - that makes your cells' counts most likely, each cell's latent
value taken as a normal with the mean and SD you give, by maximum
likelihood. By default it fits on **detection**: whether any trial
succeeded, which is what makes a covered zero. Pass the result to
[`stats_fit_latent_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md)
as `link`.

## Usage

``` r
stats_fit_link(cells, with_depth = NULL, with_overdispersion = NULL,
  detection_only = NULL, timeout = 600)
```

## Arguments

- cells:

  A TSV whose header names `latent_mean`, `latent_sd`, `successes` and
  `trials`, and optionally `depth`, in any order and no other column:
  one row per covered cell. The latent mean and SD are in latent units,
  successes of trials is the cell's count, and depth is in your depth
  units. Without a depth column every depth is mzLib's default (0 at the
  pin).

- with_depth:

  Fit the depth slope (mzLib's `withDepth`). `NULL` (the default) sends
  nothing: mzLib's declared default (`FALSE` at the pin), which
  `with_depth` in the result reports.

- with_overdispersion:

  Fit rho, the beta-binomial correlation of the trials in a cell
  (mzLib's `withOverdispersion`). `NULL` (the default): mzLib's declared
  default (`FALSE` at the pin).

- detection_only:

  `TRUE` fits on whether any trial succeeded, `FALSE` on the full counts
  (mzLib's `detectionOnly`). `NULL` (the default): mzLib's declared
  default (`TRUE` at the pin).

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

A slope below 1 says the chance of a success rises more slowly than the
latent: on PTM data, for example, whether a modified form is seen
depends on how much of it there is, not only on its share. **The cells,
and their latent means and SDs, are your choice**: mzLibR does not build
them. mzLib's own example is a covered cell's feature mean and SD over
its measured cells.

## Value

An `mzlibr_link_fit`: `cells_file`, the cells table's absolute path;
`cell_count`, rows read (cells); `depth_given`, whether the cells file
has a `depth` column; `with_depth`, `with_overdispersion` and
`detection_only`, as passed: what you gave, else mzLib's declared
default; `link`, the fitted `mzlibr_detection_link`, ready for
[`stats_fit_latent_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md)'s
`link`, or `NULL` when mzLib's search did not converge - `intercept`
(natural log-odds, the per-trial log-odds at latent 0 and depth 0),
`slope` (natural log-odds per latent unit), `depth_slope` (natural
log-odds per depth unit; exactly 0 when `with_depth` is `FALSE`) and
`overdispersion` (fraction, the intra-cell correlation, 0 to below 1;
exactly 0 when `with_overdispersion` is `FALSE`); and `caveats`, what
the fit is and is not.

## What it is not

A point estimate only: mzLib returns no standard errors. The cells and
their latent normals are modelling choices, and the link fitted from
them is only as good as they are. With `with_depth = TRUE` and no depth
column, every depth is 0 and the data say nothing about `depth_slope`.

## Wraps

Wire verb `stats fit-link`. Generated from the bridge's verb spec
`stats.fit-link.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`LogitNormalMixedModel.FitLink`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`LinkCell`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

- [`DetectionLink`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/StatisticalModels/LogitNormalMixedModel.cs)
  in `mzLib/StatisticalModels/LogitNormalMixedModel.cs` at mzLib
  `cc67b7a2`

## Parameters: units, ranges and defaults

- `cells`:

  path; required; range `>= 1 row`. A TSV with a header row naming the
  columns latent_mean, latent_sd, successes and trials, and optionally
  depth, in any order and no other (so a misspelt optional column cannot
  fall back to the default silently), then one row per covered cell
  (mzLib's LinkCell). latent_mean and latent_sd are the normal the
  cell's latent value is taken to sit at, in latent units (e.g. the mean
  and SD of its feature's measured cells, mzLib's example); successes of
  trials is its count; depth is in the caller's depth unit. Without a
  depth column every cell is constructed without one, so its depth is
  LinkCell's own default (0 at cc67b7a2). mzLib refuses a cell with a
  latent mean, SD or depth that is not finite, an SD below 0, no trials,
  or more successes than trials. Which cells, and their latent normals,
  are the caller's modelling choices.

- `with_depth` (wire `--with-depth`):

  bool; default absent; range `true | false`. mzLib's withDepth: fit the
  depth slope (else it is 0). Absent, the bridge passes FitLink's
  declared default (false at cc67b7a2), read from mzLib, and echoes it.

- `with_overdispersion` (wire `--with-overdispersion`):

  bool; default absent; range `true | false`. mzLib's
  withOverdispersion: fit the beta-binomial intra-cell correlation rho
  (else it is 0, binomial). Absent, FitLink's declared default (false at
  cc67b7a2), echoed.

- `detection_only` (wire `--detection-only`):

  bool; default absent; range `true | false`. mzLib's detectionOnly.
  true: fit on whether any trial succeeded, P(k = 0) against P(k \>= 1),
  the event that makes a covered zero. false: fit on the full count, P(k
  of n). Absent, FitLink's declared default (true at cc67b7a2), echoed.
  mzLib's documentation: on PTM data the full counts follow occupancy
  more weakly than detection does, and a link fitted on them predicts
  covered zeros worse.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `cells_file`:

  string; never `NA`. Absolute path of the cells table.

- `cell_count`:

  int; in **cells**; never `NA`. Rows read.

- `depth_given`:

  bool; never `NA`. Whether the cells file has a depth column.

- `with_depth`:

  bool; never `NA`. withDepth as passed: the option given, else mzLib's
  declared default.

- `with_overdispersion`:

  bool; never `NA`. withOverdispersion as passed: the option given, else
  mzLib's declared default.

- `detection_only`:

  bool; never `NA`. detectionOnly as passed: the option given, else
  mzLib's declared default.

- `link`:

  object; `NA` when mzLib's search did not converge within its iteration
  limit (FitLink returned null). The fitted DetectionLink, ready for
  stats fit-latent-mixed `--link`; fields under result.tables.link.

- `caveats`:

  string\[\]; never `NA`. Three static caveats (one says which
  likelihood was fitted), plus one each when the depth slope or
  overdispersion was not fitted, when the depth slope was fitted without
  a depth column, and when the search did not converge.

## Fields of `link`

an `mzlibr_detection_link`, a named list of these four numbers.

- `link$intercept` (wire `intercept`):

  float; in **natural log-odds**; never `NA`. a: the per-trial log-odds
  at latent 0 and depth 0.

- `link$slope` (wire `slope`):

  float; in **natural log-odds per latent unit**; never `NA`. b: how
  fast the per-trial log-odds follow the latent. Below 1, the chance of
  a success rises more slowly than the latent.

- `link$depth_slope` (wire `depth_slope`):

  float; in **natural log-odds per depth unit**; never `NA`. c: the
  depth term's coefficient; exactly 0 when with_depth is false.

- `link$overdispersion` (wire `overdispersion`):

  float; in **fraction (intra-cell correlation, 0 to \< 1)**; never
  `NA`. rho, the correlation of trials within a cell; exactly 0 when
  with_overdispersion is false.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--cells` missing, blank or not found; the table empty or header-only;
  a header with an empty, repeated, missing or unknown column (each
  named); a row whose field count differs from its header; a
  latent_mean, latent_sd or depth that is not a number (blank, NA or NaN
  pass to mzLib); successes or trials not a whole number;
  `--with-depth`, `--with-overdispersion` or `--detection-only` other
  than true or false, or without a value

- `mzlib_usage_error` (usage):

  mzLib refuses the cells (its ArgumentException, reclassified): a
  latent mean, SD or depth that is not finite, an SD below 0, no trials,
  or successes outside 0 to trials

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Maximum likelihood over the cells, each cell's probability integrated
  over its latent normal by 12-node adaptive Gauss-Hermite quadrature
  and maximized by Nelder-Mead. mzLib returns the point estimate only:
  no standard errors.

- By default the fit is on detection (any success or none), not on how
  many trials succeeded: on PTM data, mzLib's documentation says, the
  full counts follow occupancy more weakly and predict covered zeros
  worse.

- The link is in natural log-odds and is what stats fit-latent-mixed
  `--link` takes. Which cells enter, and their latent means and SDs, are
  the caller's modelling choices; the bridge builds neither.

- With `--with-depth` true and no depth column, every cell's depth is 0
  and the data say nothing about depth_slope.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.fit_link`

- Rust (mzLibRust): `mzlib::stats::fit_link_with` with `FitLinkOptions`

- R (mzLibR): `stats_fit_link`

## Since

Wire protocol 1; pyMzLib 0.7.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.1093/comjnl/7.4.308](https://doi.org/10.1093/comjnl/7.4.308):
  the simplex search mzLib maximizes the link's likelihood with (Nelder
  and Mead, Comput J 7:308, 1965)

- [doi:10.1093/biomet/81.3.624](https://doi.org/10.1093/biomet/81.3.624):
  adaptive Gauss-Hermite quadrature, which integrates each cell's
  probability over its latent normal (Liu and Pierce, Biometrika 81:624,
  1994)

- [doi:10.1111/j.2517-6161.1948.tb00014.x](https://doi.org/10.1111/j.2517-6161.1948.tb00014.x):
  the beta-binomial distribution, the count model with overdispersion
  (Skellam, J R Stat Soc B 10:257, 1948)

## See also

[`stats_fit_latent_mixed`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md)

## Examples

``` r

# mzLib's own test cells: 4,000 cells simulated through the link (-0.5, 0.5, 0.3), which this
# package ships. The link is recovered from detection alone:
cells <- system.file("extdata", "stats", "link_cells.tsv", package = "mzLibR")
fit <- stats_fit_link(cells, with_depth = TRUE)
fit
#> <mzlibr_link_fit> 4000 cells, fitted on detection (any success or none)
#>   logit g = -0.5108 + 0.5065 x latent + 0.3313 x depth, overdispersion 0
round(unlist(fit$link), 2)
#>      intercept          slope    depth_slope overdispersion 
#>          -0.51           0.51           0.33           0.00 
```
