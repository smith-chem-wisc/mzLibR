# Fit each run's map from minutes onto a library's iRT scale, and convert times either way

Fit each run's map from retention time in minutes onto a spectral
library's iRT scale from anchors you choose, and convert the times you
send with them, in either direction.

Answers **where should this library precursor elute in my run?** and
**what is this run's minute on the library's scale?** A spectral library
stores retention time on a run-independent scale, iRT (Escher et al.
2012); a run reports minutes. mzLib's `IrtCalibration.Fit` maps one onto
the other from **anchors**: precursors you identified confidently in the
run, each with its observed apex time in minutes and its library iRT.
The fit is robust: it starts from binned medians and downweights anchors
far from the curve (Tukey bisquare), so some wrong identifications do
not bend it. It is strictly increasing, so it inverts. By default it is
a LOWESS curve (Cleveland 1979) that follows a gradient that bends;
`kind = "linear"` fits one straight line.

## Usage

``` r
library_calibrate_rt(anchors, queries = NULL, kind = NULL, minimum_anchors = NULL,
  bandwidth = NULL, knots = NULL, robustness_iterations = NULL, on_error = "fail",
  timeout = 60)
```

## Arguments

- anchors:

  One row per anchor: a data.frame with columns `group` (any name: the
  run), `rt_minutes` (the observed retention time, in **minutes**) and
  `irt` (the library's indexed retention time, in **iRT**), or a list of
  `list(group, rt_minutes, irt)`. They travel on stdin, one line each.
  mzLib refuses an `NA` or `NaN` time.

- queries:

  Times to convert, in the same shape: each row gives its `group` and
  exactly one of `rt_minutes` (minutes, converted to iRT) or `irt` (iRT,
  converted to minutes), the other `NA`. A data.frame may leave out the
  column it never gives. Each group's queries use that group's fit.
  `NULL` (the default) converts nothing: you get each group's fit and
  its anchors' residuals.

- kind:

  `"linear"` or `"lowess"`, mzLib's `IrtCalibrationKind`. `NULL` (the
  default) sends nothing, so mzLib's default, LOWESS, applies.

- minimum_anchors:

  A group with fewer anchors than this is refused, not fitted (anchors).
  `NULL` (the default) leaves mzLib's default, 20 anchors at the pin.
  mzLib requires at least 3.

- bandwidth:

  LOWESS only: the **fraction of anchors** in each local fit, in (0,
  1\]. Smaller follows the gradient, and its noise, more closely. `NULL`
  (the default) leaves mzLib's default, 0.3.

- knots:

  LOWESS only: how many evenly spaced times between the first and last
  anchor the curve is fitted at (points); between them it is a straight
  line. `NULL` (the default) leaves mzLib's default, 100 points. At
  least 3.

- robustness_iterations:

  Rounds of Tukey bisquare reweighting after the first fit (iterations).
  `NULL` (the default) leaves mzLib's default, 3 iterations. 0 or more.

- on_error:

  `"fail"` (the default) stops at the first group mzLib refuses, with an
  `mzlib_usage_error` naming it; `"skip"` records the refusal in that
  group's row of `groups` and fits the rest. An option mzLib refuses
  fails the call either way.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

Each **group** (a run, typically) is fitted on its own anchors, and the
**queries** you send for it are converted with that fit: give a time in
minutes to get its iRT, or an iRT to get the minutes at which that
precursor should elute in this run.

**You choose the anchors.** mzLib applies no q-value cutoff, and mzLibR
adds none. **Send the queries in the same call as their anchors**: mzLib
keeps the fitted curve's points private, so a calibration cannot be
returned and used later. Refitting the same anchors gives the same map,
because mzLib sorts them first.

## Value

An `mzlibr_rt_calibration`.

`records` is a data.frame with one row per anchor, then one per query,
in the order sent: `group`, as sent; `role`, `"anchor"`, `"to_irt"` (you
gave minutes) or `"to_minutes"` (you gave iRT); `rt_minutes` (min), as
given, or the converted time for a `"to_minutes"` query, `NA` when the
group was not fitted and you gave none; `irt` (iRT), as given, or the
converted value for a `"to_irt"` query, `NA` likewise; `fitted_irt`
(iRT), anchors only, the curve's iRT at the anchor's time, so
`irt - fitted_irt` is its residual, `NA` on a query or when the group
was not fitted; and `is_extrapolated`, whether `rt_minutes` lies outside
the group's anchors, where the map is a straight-line continuation
rather than a fit, `NA` when the group was not fitted.

`groups` is a data.frame, one row per group in order of first
appearance: `group`; `kind`, the shape mzLib fitted (`NA`: not fitted);
`anchor_count`, the anchors it was fitted from (anchors; the anchors
sent, if refused); `query_count`, the times converted (queries);
`residual_sd`, how far the anchors that fit sit from the curve, in **iRT
units**, mzLib's robust residual SD with wrong anchors set aside (divide
by the curve's slope, iRT per minute, to read it in minutes; `NA`: not
fitted); `first_anchor_minutes` and `last_anchor_minutes`, the span the
map is fitted over, in minutes (`NA`: not fitted); and `error_kind`,
`error_type` and `error_message`, why mzLib refused the group under
`on_error = "skip"` (`"usage"`, `"ArgumentException"` and mzLib's
message), `NA` when it was fitted. `groups[groups$group == name, ]` is
one group; `groups[!is.na(groups$error_kind), ]` the groups mzLib
refused.

`kind` is the kind applied, yours or mzLib's default; `minimum_anchors`
(anchors), `bandwidth` (fraction of anchors), `knots` (points) and
`robustness_iterations` (iterations) the options applied, yours or
mzLib's; `on_error` as requested. `group_count` counts the groups sent
(groups), `fitted_count` those mzLib fitted (groups), `failed_count`
those it refused (groups; non-zero only under `"skip"`), `anchor_count`
the anchors sent (anchors), `query_count` the queries sent (queries),
and `row_count` the rows of `records` (rows):
`anchor_count + query_count`. `column_names` lists the columns.
`caveats` says what a calibration can and cannot tell you, plus any that
apply to this call (queries outside the anchors, groups not fitted).

## The traps

**iRT is the scale of one library.** A calibration maps one run onto the
scale of the library the anchors' iRT values came from, so iRT from two
libraries is not comparable. \*\*`residual_sd` is robust\*\*: it
understates the spread when many anchors are wrong, and a LOWESS fit to
few anchors can pass through nearly all of them, so a `residual_sd` near
0 there means an overfit, not a perfect run. With few anchors, prefer
`kind = "linear"`. **Outside the first and last anchor the map is a
straight line**, not a fit: read `is_extrapolated` before trusting a
converted time there.

## Wraps

Wire verb `library calibrate-rt`. Generated from the bridge's verb spec
`library.calibrate-rt.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`IrtCalibration.Fit`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs)
  in `mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs`
  at mzLib `cc67b7a2`

- [`IrtCalibrationOptions`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs)
  in `mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs`
  at mzLib `cc67b7a2`

- [`IrtCalibrationKind`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs)
  in `mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs`
  at mzLib `cc67b7a2`

- [`IrtCalibrationModel.ToIrt`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs)
  in `mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs`
  at mzLib `cc67b7a2`

- [`IrtCalibrationModel.ToRtMinutes`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs)
  in `mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs`
  at mzLib `cc67b7a2`

- [`IrtCalibrationModel.IsExtrapolated`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs)
  in `mzLib/Chromatography/RetentionTimeCalibration/IrtCalibration.cs`
  at mzLib `cc67b7a2`

- [`Irt, RtMinutes`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/MzLibUtil/RetentionTimeUnits.cs)
  in `mzLib/MzLibUtil/RetentionTimeUnits.cs` at mzLib `cc67b7a2`

## Parameters: units, ranges and defaults

- `anchors` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 anchor line; optionally a '--' line and query lines`. Two lists
  of tab-separated lines, group, rt_minutes, irt, split by one line
  holding only '–' (BULK.md section 2). No header; blank lines are
  skipped; a UTF-8 BOM is dropped. ANCHORS come first: one (observed
  retention time in minutes, library iRT) pair per line; mzLib asks for
  the observed apex time of confidently identified precursors
  (IrtCalibration.cs:29-31), and which ones is the caller's choice. A
  blank, NA or NaN anchor value reaches mzLib as NaN, which it refuses.
  QUERIES follow the '–', if any: each gives exactly one of rt_minutes
  (convert to iRT, IrtCalibrationModel.ToIrt) or irt (convert to
  minutes, ToRtMinutes) and leaves the other blank. group is any
  non-empty name, typically the run; each group is fitted by its own
  IrtCalibration.Fit call on its own anchors, and its queries are
  converted by that fit. Groups are reported in order of first
  appearance, over both lists.

- `kind`:

  string; default `lowess`; range `linear | lowess`.
  IrtCalibrationOptions.Kind. lowess (mzLib's default) is a robust
  locally weighted regression for a gradient that bends; linear is one
  straight line, fitted robustly. Absent sends nothing, so mzLib's
  default applies.

- `minimum_anchors` (wire `--minimum-anchors`):

  int; in **anchors**; default `20`; range `>= 3 (mzLib's own check)`.
  IrtCalibrationOptions.MinimumAnchors: a group with fewer anchors is
  refused rather than fitted. Absent means mzLib's default, 20 at the
  pin; the bridge sends no number of its own.

- `bandwidth`:

  float; in **fraction of anchors**; default `0.3`; range
  `(0, 1] (mzLib's own check)`. IrtCalibrationOptions.Bandwidth, LOWESS
  only: the share of a group's anchors in each local fit (at least 3).
  Smaller follows the gradient more closely and noise too. Absent means
  mzLib's default, 0.3 at the pin. mzLib checks it for linear too, and
  does not use it there.

- `knots`:

  int; in **points**; default `100`; range `>= 3 (mzLib's own check)`.
  IrtCalibrationOptions.Knots, LOWESS only: how many evenly spaced times
  between the first and last anchor the curve is fitted at; between them
  it is a straight line. Absent means mzLib's default, 100 at the pin.

- `robustness_iterations` (wire `--robustness-iterations`):

  int; in **iterations**; default `3`; range `>= 0 (mzLib's own check)`.
  IrtCalibrationOptions.RobustnessIterations: rounds of Tukey bisquare
  reweighting after the first fit, each downweighting anchors with large
  residuals. 0 keeps the first fit, which is already weighted from
  binned medians. Absent means mzLib's default, 3 at the pin.

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. fail: the first group
  mzLib refuses (too few anchors, a non-finite one, anchors at one time
  only, or iRT that does not increase with time) stops the call with a
  usage error naming the group. skip: that group's error is set, its
  rows carry the times given and null for the rest, and the other groups
  are fitted. An option mzLib refuses fails the call either way.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `kind`:

  string; never `NA`. linear or lowess: the kind applied, as given or
  mzLib's default.

- `minimum_anchors`:

  int; in **anchors**; never `NA`. The MinimumAnchors applied: as given,
  or mzLib's default.

- `bandwidth`:

  float; in **fraction of anchors**; never `NA`. The Bandwidth applied:
  as given, or mzLib's default.

- `knots`:

  int; in **points**; never `NA`. The Knots applied: as given, or
  mzLib's default.

- `robustness_iterations`:

  int; in **iterations**; never `NA`. The RobustnessIterations applied:
  as given, or mzLib's default.

- `on_error`:

  string; never `NA`. fail or skip, as requested.

- `group_count`:

  int; in **groups**; never `NA`. Distinct groups over both lists.

- `fitted_count`:

  int; in **groups**; never `NA`. Groups mzLib fitted.

- `failed_count`:

  int; in **groups**; never `NA`. Groups mzLib refused; non-zero only
  under on-error skip.

- `anchor_count`:

  int; in **anchors**; never `NA`. Anchor lines, over every group.

- `query_count`:

  int; in **queries**; never `NA`. Query lines, over every group.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns: anchor_count +
  query_count.

- `groups`:

  object\[\]; never `NA`. One per group in order of first appearance:
  {group; kind (linear or lowess, IrtCalibrationModel.Kind; null if not
  fitted); anchor_count (anchors, IrtCalibrationModel.AnchorCount, or
  the lines sent if not fitted); query_count (queries); residual_sd (iRT
  units: mzLib's robust residual SD of the anchors that fit,
  IrtCalibrationModel.ResidualSd; null if not fitted);
  first_anchor_minutes and last_anchor_minutes (min: the span the map is
  fitted over, FirstAnchor and LastAnchor; null if not fitted); error
  ({kind: usage, type: the mzLib exception name, message: mzLib's
  message} under on-error skip, else null)}.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per stdin line, in input order: the
  anchors, then the queries.

- `caveats`:

  string\[\]; never `NA`. Five static caveats, plus one when any query
  falls outside its group's anchors (with the count), and one when any
  group was not fitted.

## Columns of `records`

- `group`:

  string; never `NA`. As given on stdin.

- `role`:

  string; never `NA`. anchor; to_irt (a query that gave rt_minutes); or
  to_minutes (a query that gave irt).

- `rt_minutes`:

  float; in **min**; `NA` when a to_minutes query in a group that was
  not fitted, or an anchor sent blank or NA in one (mzLib refused it).
  Retention time in this run: as given for an anchor or a to_irt query;
  IrtCalibrationModel.ToRtMinutes(irt) for a to_minutes query.

- `irt`:

  float; in **iRT**; `NA` when a to_irt query in a group that was not
  fitted, or an anchor sent blank or NA in one (mzLib refused it).
  Indexed retention time on the anchors' library scale: as given for an
  anchor or a to_minutes query; IrtCalibrationModel.ToIrt(rt_minutes)
  for a to_irt query.

- `fitted_irt`:

  float; in **iRT**; `NA` when a query, whose irt is itself the
  converted value, or a group that was not fitted. Anchors only: the
  fitted map's iRT at the anchor's time,
  IrtCalibrationModel.ToIrt(rt_minutes). irt minus fitted_irt is the
  anchor's residual.

- `is_extrapolated`:

  bool; `NA` when the group was not fitted.
  IrtCalibrationModel.IsExtrapolated(rt_minutes): the time lies before
  the group's first anchor or after its last, where the map is a
  straight-line continuation rather than a fit.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  kind given without a value or not linear or lowess; minimum-anchors,
  knots or robustness-iterations not an integer; bandwidth not a number;
  any of them given without a value; on-error not fail or skip

- `mzlib_usage_error` (usage):

  stdin has no anchor line; a line with other than 3 tab-separated
  fields; an empty group; a value that is not a number (other than
  blank, NA, NaN); a second '–' line; a query with both rt_minutes and
  irt, or with neither; a query value that is not finite

- `mzlib_usage_error` (usage):

  mzLib refuses an option (IrtCalibrationOptions out of range:
  MinimumAnchors below 3, Bandwidth outside (0, 1\], Knots below 3,
  RobustnessIterations below 0; mzLib's ArgumentOutOfRangeException,
  reclassified), whatever on-error says

- `mzlib_usage_error` (usage):

  with on-error fail, mzLib refuses a group's anchors: fewer than
  minimum_anchors (a group with queries and no anchors included), a
  non-finite value, all at one time, or a fit whose iRT does not
  increase with time (mzLib's ArgumentException, reclassified); the
  message starts with the group

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- iRT is the scale of the library the anchors' iRT values came from. A
  calibration maps one run onto that library's scale, so iRT values from
  two libraries are not comparable, and the same run calibrated against
  two libraries gives two maps.

- The anchors are the caller's choice: mzLib asks for the observed apex
  retention time, in minutes, of confidently identified precursors
  paired with their library iRT (IrtCalibration.cs:29-31). The bridge
  applies no q-value cutoff and removes no duplicate. The fit is robust
  to some wrong anchors (binned medians, then Tukey bisquare
  reweighting), not to a wrong choice of them.

- residual_sd is mzLib's robust residual SD, in iRT units, of the
  anchors that fit: residuals beyond 3 robust SDs are set aside before
  it is computed (IrtCalibration.InlierSd). It understates the spread
  when many anchors are wrong, and a LOWESS fit to few anchors can pass
  through nearly all of them, so a residual_sd near 0 there means an
  overfit, not a perfect run. irt minus fitted_irt on the anchor rows
  shows every residual.

- Outside each group's first and last anchor the map continues as a
  straight line, not a fit (IrtCalibrationModel.IsExtrapolated). Read
  is_extrapolated before trusting a converted time there.

- The fitted curve itself does not cross the wire: mzLib keeps its knots
  private and its constructor internal (IrtCalibrationModel), so a
  calibration cannot be returned and reused. Send the queries in the
  same call as their anchors; refitting the same anchors gives the same
  map, because mzLib sorts them first.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.library.calibrate_rt`

- Rust (mzLibRust): `mzlib::library::calibrate_rt` with
  `RtCalibrationOptions`

- R (mzLibR): `library_calibrate_rt`

## Since

Wire protocol 1; pyMzLib 0.7.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The fitted curve cannot be returned: IrtCalibrationModel's knots and
  end slopes are private fields and its constructor is internal
  (IrtCalibration.cs:273-278). A public read-only Knots (and the two end
  slopes) would let the envelope carry each group's curve, at most Knots
  points, for plotting and for reuse without refitting
  (PUSH-TO-MAINLAND).

- IrtCalibration.Line, the two-point line a search uses before any
  anchors exist, is not projected: it is useful only inside a search
  engine. Revisit on demand.

- Calibrated minutes cannot be written back into a .msl:
  MslLibrary.WithCalibratedRetentionTimes takes a slope and intercept,
  and mzLib does not route it through IrtCalibrationModel (mzLib \#1387
  says that is a separate PR).

## References

- [doi:10.1002/pmic.201100463](https://doi.org/10.1002/pmic.201100463):
  the iRT concept: a normalized, run-independent retention time scale
  that a run is mapped onto by a calibration (Escher et al., Proteomics
  12:1111, 2012)

- [doi:10.1080/01621459.1979.10481038](https://doi.org/10.1080/01621459.1979.10481038):
  robust LOWESS with tricube local lines and bisquare reweighting at 6 x
  the median absolute residual, the curve IrtCalibrationKind.Lowess fits
  (Cleveland, J Am Stat Assoc 74:829, 1979)

## See also

[`library_read`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md)
for a library's iRT values,
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
for a search's observed retention times

## Examples

``` r

# Three 24-minute DIA runs of PBMCs, from a DIA-NN report in mzLib's test data. One precursor,
# AAAAAAAAGAFAGR at charge 2, is held out of the anchors, and each run predicts where it should
# elute from its library iRT:
pbmc <- read.delim(system.file("extdata", "calibration", "diann_pbmc_rt.tsv", package = "mzLibR"))
held <- pbmc$precursor_id == "AAAAAAAAGAFAGR2"
anchors <- data.frame(group = pbmc$run[!held], rt_minutes = pbmc$retention_time[!held],
                      irt = pbmc$indexed_retention_time[!held])
queries <- data.frame(group = pbmc$run[held], irt = 11.886096)
calibration <- library_calibrate_rt(anchors, queries, kind = "linear", minimum_anchors = 10)
calibration
#> <mzlibr_rt_calibration> linear: 3 of 3 group(s) fitted from 34 anchors, 3 queries converted
#>   20240316_PBMCs_24min_3Th-1_centroid: linear, 13 anchors over 1.23 to 18.74 min, residual SD 4.07 iRT
#>   20240317_PBMCs_24min_3Th-1_centroid: linear, 11 anchors over 1.261 to 18.73 min, residual SD 3.65 iRT
#>   20240315_PBMCs_24min_3Th-1_centroid: linear, 10 anchors over 1.231 to 18.59 min, residual SD 0.67 iRT
predicted <- calibration$records[calibration$records$role == "to_minutes", ]
# Each prediction is within a quarter of a minute of where the precursor was seen:
data.frame(run = substr(predicted$group, 1, 8), predicted = round(predicted$rt_minutes, 2),
           seen = round(pbmc$retention_time[held], 2))
#>        run predicted seen
#> 1 20240316      9.22 9.21
#> 2 20240317      9.45 9.23
#> 3 20240315      9.17 8.94
```
