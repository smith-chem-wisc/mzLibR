# Estimate a search's false discovery proportion from its entrapment hits

Estimate the false discovery proportion among a search's discoveries
from its entrapment hits, as a lower bound, a combined upper bound and a
paired upper bound, at each of several q-value thresholds.

An entrapment experiment adds peptides that cannot be in the sample
(from another species, or made from the target peptides) to the search
database. Every entrapment the search reports is then a known false
discovery, which lets you check the search's FDR claim instead of
trusting it. mzLib's `EntrapmentFdp.Sweep` computes the three valid
estimators of Wen et al. 2025 (the method FDRBench implements) at each
q-value threshold:

## Usage

``` r
stats_entrapment_fdp(targets, entrapments, ratio, thresholds = NULL, timeout = 600)
```

## Arguments

- targets:

  A TSV whose header names exactly `score` and `q_value`, in any order:
  one row per original-target candidate (precursor, peptide or protein
  group, whichever you are estimating for), decoys removed, each with
  its best q-value. Higher scores are better.

- entrapments:

  A TSV whose header names `score`, `q_value`, `partner_score` and
  `partner_q_value`, and optionally `foreign`: one row per entrapment
  candidate. The partner columns copy the score and q-value of the one
  target this entrapment was made from, when the search identified it;
  leave both blank otherwise. `foreign` is 1 for an entrapment with no
  partner by construction (another species); blank or absent is mzLib's
  default, not foreign. Other columns are refused, so select these
  first.

- ratio:

  r, the entrapment-to-target database size ratio (unitless). Required:
  mzLib has no default. `paired` is reported only when it is exactly 1,
  so for 1:1 pairing pass 1, not a ratio computed from database sizes
  after exclusions.

- thresholds:

  The q-value thresholds to estimate at, as fractions (0.01 for 1%).
  Leave it `NULL` for mzLib's default sweep (0.001 to 0.100 in steps of
  0.001 at mzLib `d10017ce`), echoed in `thresholds`.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

\- `lower_bound` = N_E / (N_T + N_E). Above the FDR you claimed, the
search **fails** to control it. Below, it proves nothing. - `combined` =
N_E (1 + 1/r) / (N_T + N_E). An upper bound: at or below the claimed
FDR, it is **evidence of control**. - `paired`: a tighter upper bound
when each target has exactly one entrapment made from it (`ratio = 1`),
using whether each discovered entrapment's own target was discovered
too.

The "sample" estimator N_E / (r N_T) is invalid in both directions and
is deliberately absent.

## Value

An `mzlibr_entrapment_fdp_estimates`:

\- `records`, a data.frame with one row per q-value threshold with a
discovery: `q_value_threshold` (a q-value, fraction), where a candidate
is discovered at or below it; `target_count`, N_T, targets discovered
(discoveries); `entrapment_count`, N_E, entrapments discovered
(discoveries); and three estimates of the false discovery proportion:
`lower_bound` (fraction of discoveries), a **lower bound** (above your
claimed FDR, the search fails to control it); `combined` (fraction of
discoveries), an **upper bound** (at or below your claimed FDR, evidence
of control; can exceed 1); and `paired` (fraction of discoveries), a
tighter upper bound for 1:1 pairing (can exceed 1), `NA` where it does
not apply: a ratio not exactly 1, or every entrapment foreign. -
`targets_file` and `entrapments_file`, as absolute paths; `ratio`, r as
you gave it (entrapment-to-target database size ratio, unitless);
`thresholds`, the q-value thresholds estimated at, in order (yours or
mzLib's default); `thresholds_given`; `target_candidate_count` and
`entrapment_candidate_count`, the rows of the two files (candidates);
`foreign_count` and `partnered_count` entrapments (candidates);
`threshold_count` thresholds; `empty_threshold_count`, thresholds with
no discoveries (thresholds), which have no row; `row_count` rows;
`column_names`. - `caveats`: what each estimate can and cannot show,
plus any that apply to this input.

## The input traps

Give one row per unit estimated, decoys removed, each with its best
q-value: mzLib does not deduplicate, and it does not check a partner's
copied score against the targets file. The k-fold paired estimator for k
entrapments per target is not implemented.

## Wraps

Wire verb `stats entrapment-fdp`. Generated from the bridge's verb spec
`stats.entrapment-fdp.yaml` (bridge commit `45833ca804d6`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`EntrapmentFdp.Sweep`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EntrapmentFdp.cs)
  in `mzLib/StatisticalModels/EntrapmentFdp.cs` at mzLib `d10017ce`

- [`EntrapmentFdp.DefaultQValueThresholds`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EntrapmentFdp.cs)
  in `mzLib/StatisticalModels/EntrapmentFdp.cs` at mzLib `d10017ce`

- [`ScoredIdentification`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EntrapmentFdp.cs)
  in `mzLib/StatisticalModels/EntrapmentFdp.cs` at mzLib `d10017ce`

- [`ScoredEntrapment`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EntrapmentFdp.cs)
  in `mzLib/StatisticalModels/EntrapmentFdp.cs` at mzLib `d10017ce`

- [`EntrapmentFdpPoint`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/StatisticalModels/EntrapmentFdp.cs)
  in `mzLib/StatisticalModels/EntrapmentFdp.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `targets`:

  path; required. A TSV with a header row naming EXACTLY the columns
  score and q_value, in any order, then one row per original-target
  candidate (mzLib's ScoredIdentification): the unit you estimate the
  FDP of (precursor, peptide or protein group), decoys removed, each
  with its best q-value. score is the search's score, higher is better,
  on the same scale as the entrapments'; q_value is the q-value the
  search reported, a fraction. Every cell must be a finite number. Any
  other column is refused: select the two before writing the file.

- `entrapments`:

  path; required. A TSV with a header row naming the columns score,
  q_value, partner_score and partner_q_value, and optionally foreign, in
  any order, then one row per entrapment candidate (mzLib's
  ScoredEntrapment). partner_score and partner_q_value are a COPY of the
  score and q-value of the one original target this entrapment was made
  from, when the search identified it; leave both blank (or NA) when it
  did not, when it cannot be resolved to one target, or when there is
  none. mzLib does not check them against the targets file. foreign is 1
  or true for an entrapment with no target partner by construction (a
  foreign species), 0 or false otherwise; a blank cell, or no foreign
  column, is mzLib's default (ScoredEntrapment.IsForeign false). Any
  other column is refused, so a misspelt foreign column cannot fall back
  to the default.

- `ratio`:

  float; in **entrapment-to-target database size ratio (unitless)**;
  required; range
  `> 0 and finite (mzLib's rule; anything else is mzLib's refusal, a usage error)`.
  r, the effective ratio of entrapment to target database size, used by
  the combined estimate. REQUIRED: mzLib has no default. paired is
  reported only when r is exactly 1: for 1:1 pairing give 1 itself, not
  a ratio computed from database sizes after exclusions (0.9998 drops
  it).

- `thresholds`:

  float\[\]; in **q-value (fraction, 0 to 1)**; default absent; range
  `comma-separated finite numbers`. The q-value thresholds to estimate
  at, comma-separated, in the order the rows should come in. A candidate
  is discovered at a threshold when its q-value is at or below it.
  Absent, the bridge sends mzLib's EntrapmentFdp.DefaultQValueThresholds
  (0.001, 0.002, ..., 0.100, each i/1000, at d10017ce) and echoes them.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `targets_file`:

  string; never `NA`. Absolute path of the targets table.

- `entrapments_file`:

  string; never `NA`. Absolute path of the entrapments table.

- `ratio`:

  float; in **entrapment-to-target database size ratio (unitless)**;
  never `NA`. r, as given.

- `thresholds`:

  float\[\]; in **q-value (fraction, 0 to 1)**; never `NA`. The
  thresholds estimated at, in order: the `--thresholds` given, or
  mzLib's default.

- `thresholds_given`:

  bool; never `NA`. True when `--thresholds` was given; false when
  thresholds is mzLib's DefaultQValueThresholds.

- `target_candidate_count`:

  int; in **candidates**; never `NA`. Rows read from the targets file.

- `entrapment_candidate_count`:

  int; in **candidates**; never `NA`. Rows read from the entrapments
  file.

- `foreign_count`:

  int; in **candidates**; never `NA`. Entrapments marked foreign.

- `partnered_count`:

  int; in **candidates**; never `NA`. Entrapments with a partner.

- `threshold_count`:

  int; in **thresholds**; never `NA`. The length of thresholds.

- `empty_threshold_count`:

  int; in **thresholds**; never `NA`. Thresholds at which nothing was
  discovered, where every estimate is 0/0, so mzLib leaves them out of
  the table.

- `row_count`:

  int; in **rows**; never `NA`. threshold_count - empty_threshold_count.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per threshold with a discovery, in
  thresholds order.

- `caveats`:

  string\[\]; never `NA`. Five static caveats, plus one each when a
  threshold had no discoveries, when paired is null at every threshold
  (with the ratio and the foreign count), and when no entrapment has a
  partner.

## Columns of `records`

- `q_value_threshold`:

  float; in **q-value (fraction, 0 to 1)**; never `NA`. The threshold.

- `target_count`:

  int; in **discoveries**; never `NA`. N_T: targets with a q-value at or
  below the threshold.

- `entrapment_count`:

  int; in **discoveries**; never `NA`. N_E: entrapments with a q-value
  at or below the threshold.

- `lower_bound`:

  float; in **fraction of discoveries (an FDP estimate)**; never `NA`.
  Eq. 2, N_E / (N_T + N_E). A lower bound on the FDP: above the claimed
  FDR, it shows the search fails to control it; below, it shows nothing.

- `combined`:

  float; in **fraction of discoveries (an FDP estimate; can exceed 1)**;
  never `NA`. Eq. 1, N_E (1 + 1/r) / (N_T + N_E). An upper bound: at or
  below the claimed FDR, it is evidence of control. Not capped at 1.

- `paired`:

  float; in **fraction of discoveries (an FDP estimate; can exceed 1)**;
  `NA` when does not apply: the ratio is not exactly 1, or every
  entrapment is foreign. Eq. 4, (N_E + N\_{E\>=s\>T} + 2 N\_{E\>T\>=s})
  / (N_T + N_E): a tighter upper bound for 1:1 pairing. N\_{E\>=s\>T}
  counts discovered entrapments whose partner was not discovered (no
  partner counts so); N\_{E\>T\>=s} those whose partner was discovered
  but ranked behind them, by lower q-value and then higher score, an
  exact tie counting the entrapment ahead. Not capped at 1.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  targets or entrapments missing, blank or not found; a table empty; a
  header with an empty, repeated, missing or unknown column (each
  named); a row whose field count differs from its header; a score or
  q_value that is not a finite number; a partner with a score but no
  q-value or the reverse, or not finite; a foreign cell other than 0, 1,
  true, false or blank

- `mzlib_usage_error` (usage):

  `--ratio` absent, without a value or not a number; `--thresholds`
  without a value, empty, or with an entry that is not a finite number

- `mzlib_usage_error` (usage):

  mzLib refuses the ratio (its ArgumentOutOfRangeException,
  reclassified): not a positive finite number

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Wen et al. 2025, as FDRBench implements it. lower_bound can only show
  that a search fails to control its FDR; combined is an upper bound,
  evidence that it does; paired is a tighter upper bound for 1:1
  pairing. None is capped at 1.

- A candidate is discovered when its q-value is at or below the
  threshold, and a partner by its own q-value. A discovered pair is
  ranked lower q-value first, then higher score, and an exact tie counts
  against the search.

- Give one row per unit estimated (precursor, peptide or protein group),
  decoys removed, each with its best q-value: mzLib does not
  deduplicate. A partner's score and q-value are copies of its target's
  row, and mzLib does not check them.

- paired is reported only for a ratio of exactly 1 with some entrapment
  not foreign. The sample estimator N_E / (r N_T), invalid in both
  directions, is deliberately absent (mzLib), and the k-fold paired
  estimator for k entrapments per target is not implemented.

- An entrapment generator's exclusions change the inputs in opposite
  directions (mzLib's EntrapmentFdp.ClassifyExclusion): an entrapment
  that is really a target peptide comes out of the entrapments file; a
  target that cannot be traced to one target stays in the targets file
  but is nobody's partner.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.stats.entrapment_fdp`

- Rust (mzLibRust): `mzlib::stats::entrapment_fdp_with` with
  `EntrapmentFdpOptions`

- R (mzLibR): `stats_entrapment_fdp`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## References

- [doi:10.1038/s41592-025-02719-x](https://doi.org/10.1038/s41592-025-02719-x):
  the lower-bound (eq. 2), combined (eq. 1) and paired (eq. 4)
  entrapment estimators mzLib implements, and why the sample estimator
  is invalid (Wen et al., Nat Methods 22:1454, 2025; the method FDRBench
  implements)

## See also

[`stats_adjust`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)

## Examples

``` r

# mzLib's own hand-built search (seven targets, four entrapments, one for each case the paired
# estimate distinguishes), checked against the entrapment project's reference script:
data <- system.file("extdata", "stats", package = "mzLibR")
fdp <- stats_entrapment_fdp(file.path(data, "entrapment_reference_targets.tsv"),
  file.path(data, "entrapment_reference_entrapments.tsv"), ratio = 1)
fdp
#> <mzlibr_entrapment_fdp_estimates> 7 targets, 4 entrapments, ratio 1, 100 of 100 thresholds with a discovery
fdp$records[fdp$records$q_value_threshold == 0.01, ]
#>    q_value_threshold target_count entrapment_count lower_bound  combined
#> 10              0.01            4                3   0.4285714 0.8571429
#>       paired
#> 10 0.8571429
```
