# Compute every site's occupancy in every run, with mzLib's SiteOccupancyCalculator

Intensity-based occupancy of every modification site in every run (or
sample), from a table of MS/MS-identified peptidoforms with their
per-run intensities, with the censoring state of each cell.

For a site (protein, position, modification) in a run, occupancy is the
summed intensity of the quantified peptidoforms that carry the
modification there, over the summed intensity of every quantified
peptidoform that covers the position - the estimator MetaMorpheus writes
as `IntensityOccupancy_`. Each cell keeps its state, so a censored cell
is never mistaken for a 0: a `"floor"` (the modified form was identified
but not quantified) is low but not zero, and `"not_detected"`
(peptidoforms covering the site were identified in the run, quantified
or not, and none carried the modification; mzLib \#1468) is a bound
below what the run could detect. Neither has a `fraction`.

## Usage

``` r
ptmqtl_occupancy(observations, samples = NULL, modifications = NULL, categories = NULL,
  out = NULL, timeout = 600)
```

## Arguments

- observations:

  A TSV with a header row naming its columns, one row per (run,
  peptidoform, protein) identified by MS/MS: `run`, `full_sequence`
  (MetaMorpheus notation), `protein_accession`, `start_residue` and
  `end_residue` (1-based, in the protein), `intensity` (instrument
  units; blank or `NA` for identified but not quantified), and
  optionally `previous_residue` (one letter). A shared peptide is one
  row per protein. Never one row per PSM, and never a match-between-runs
  transfer. mzLib has no reader that builds this table; the ptmqtl
  article shows how to assemble it.

- samples:

  A TSV with a header row, then `run`, `sample`, to sum each sample's
  fractions before the ratio; the `run` column of the result then holds
  the sample. Map only fractions together, never technical replicates.
  Leave it `NULL` for per-run occupancy.

- modifications:

  Count only these modifications as sites: exact names,
  `Category:IdWithMotif`, as a character vector. Not with `categories`.
  `NULL`: no filter.

- categories:

  Count only modifications of these categories (the text before the
  first `:`), as a character vector, e.g. `"Common Biological"`. Not
  with `modifications`.

- out:

  Write the table to this `.tsv` instead of returning it - the file
  [`ptmqtl_trait_effects`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md)
  reads. Then `records` has no rows and `written` says where.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

mzLib always leaves out the categories `Common Variable` and
`Common Fixed`, peptide-terminal modifications and every C-terminal
modification. An N-terminal modification is a site only on a peptide
that starts at residue 1, or at 2 after the initiator Met is removed.

## Value

An `mzlibr_site_occupancy`:

\- `records`, a data.frame with one row per (site, run) where the site
was identified or covered, by protein, then run, then position and
modification: `site_key`, `accession:residue position:modification`, the
key
[`ptmqtl_pairs`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_pairs.md)
and
[`ptmqtl_trait_effects`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md)
use; `protein_accession`; `position`, the 1-based residue (a protein
N-terminal modification at the residue it sits on); `residue`;
`modification`, `Category:IdWithMotif`; `category`, the text before the
first `:`; `run`, the run or with `samples` the sample; `state`, one of
`"quantified"`, `"floor"`, `"count_only"` or `"not_detected"`;
`modified_intensity` (intensity, instrument units), the summed intensity
of the quantified forms carrying the modification here, 0 unless
quantified; `covering_intensity` (intensity, instrument units), of every
quantified form covering the position, 0 for `"count_only"`; `fraction`
(fraction, 0 to 1), modified over covering intensity, `NA` unless
`state` is `"quantified"` - censored, never 0; `unmodified_quantified`,
where `FALSE` with `fraction` 1 is a **ceiling**: the site was seen only
modified; and the peptidoform counts (mzLib \#1468), which count
identifications, not intensity: `covering_peptidoforms` (peptidoforms),
the distinct peptidoforms covering the position identified in the run,
quantified or not - the n of a cell's k of n - and
`modified_peptidoforms` (peptidoforms), those of them carrying this
modification here - the k, 0 for `"not_detected"`. Both are `NA` only
from a table that did not supply them; mzLib's calculation always
does. - `observations_file` and `samples_file` (`NA` without `samples`),
as absolute paths; `observation_count` rows read; `run_count` runs;
`sample_count` samples, `NA` without `samples`; `run_names`;
`site_count` sites; `state_counts`, cells per state; `modifications` and
`categories`, the filter given, `NULL` for none; `ignored_columns`,
columns that were not read (a misspelled `previous_residue` shows up
here); `row_count` rows, returned or written; `column_names`; `written`
(`path`, `format`, `row_count`), or `NULL`. - `caveats`: what the
numbers are and are not, plus any that apply to this input.

## Positions

`position` is 1-based, unlike
[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md),
whose N-terminus is 0. A peptide shared by several proteins counts in
full for each; MetaMorpheus gives it a proportional share per protein
group, so its `IntensityOccupancy_` can differ.

## Wraps

Wire verb `ptmqtl occupancy`. Generated from the bridge's verb spec
`ptmqtl.occupancy.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SiteOccupancyCalculator.Calculate`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `cc67b7a2`

- [`PeptidoformObservation`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `d10017ce`

- [`SiteRunOccupancy`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a2/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `cc67b7a2`

- [`ModificationSite`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `d10017ce`

- [`OccupancyState`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/PtmQtl/SiteOccupancy.cs)
  in `mzLib/Quantification/PtmQtl/SiteOccupancy.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `observations`:

  path; required. A TSV with a header row naming its columns, in any
  order, one row per (run, peptidoform, protein) identified by MS/MS:
  run; full_sequence, the peptidoform in MetaMorpheus notation (a
  modification as \[Category:IdWithMotif\] after its residue, a protein
  N-terminal one before the first residue); protein_accession;
  start_residue and end_residue, the peptide's first and last residue in
  the protein, 1-based; intensity, the peptidoform's MS/MS-quantified
  intensity in that run in instrument units (blank, NA or NaN means
  identified but not quantified); and optionally previous_residue, the
  protein residue before the peptide (one letter, blank when unknown). A
  peptide shared by several proteins is one row per protein. One row per
  PSM is refused, since it would count an intensity several times. A
  match-between-runs transfer must not be a row. A 0 intensity is read
  as quantified at zero, and adds nothing to either sum. Any other
  column is not read and is listed in ignored_columns. mzLib has no
  reader that builds this table from MetaMorpheus or FlashLFQ output;
  the pyMzLib ptmqtl guide shows how.

- `samples`:

  path; default absent. Optional. A TSV with a header row, then run,
  sample: two columns, one row per run. With it, each peptidoform's
  intensity is summed over a sample's runs before the ratio (mzLib's
  runToSample overload), and the run column of the result holds the
  sample. Every run of the observations must be listed; a run listed
  here but absent there is ignored, as mzLib ignores it. Map only
  fractions of one sample together, never technical replicates. Absent,
  the bridge sends none and occupancy is per run.

- `modifications`:

  string; default absent; range
  `exact names, tab-separated, none blank`. Optional. Count only these
  modifications as sites: exact names as full_sequence writes them,
  Category:IdWithMotif, separated by tabs (a name holds spaces and
  commas, never a tab). mzLib's includeModification filter, applied
  after its own exclusions. Not with categories. Absent, the bridge
  sends no filter: every modification mzLib does not exclude is a site.

- `categories`:

  string; default absent; range
  `category names, tab-separated, none blank`. Optional. Count only
  modifications of these categories, the text before the first ':' of
  the name (mzLib's ModificationSite.Category), for example Common
  Biological. Separated by tabs. Not with modifications. Absent, no
  filter.

- `out`:

  path; default absent; range `.tsv only, any case (PYB-1)`. Write the
  table to this path as plain tab-separated text instead of returning
  columns: one header row (column_names), null as an empty cell, numbers
  in invariant culture at round-trip precision, booleans as true or
  false, nothing quoted. This is the file ptmqtl trait-effects reads as
  occupancy. Parent directories are created. Must differ from
  observations. Must end in .tsv, any case; any other extension, or
  none, is a usage error before any input is read.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `observations_file`:

  string; never `NA`. Absolute path of the observations table.

- `samples_file`:

  string; `NA` when no `--samples`: occupancy is per run. Absolute path
  of the run-to-sample table.

- `observation_count`:

  int; in **rows**; never `NA`. Observation rows read.

- `run_count`:

  int; in **runs**; never `NA`. Distinct runs in the observations.

- `sample_count`:

  int; in **samples**; `NA` when no `--samples`. Distinct samples in the
  run-to-sample table.

- `run_names`:

  string\[\]; never `NA`. The distinct values of the run column, in
  ordinal order: the runs, or with `--samples` the samples, that have at
  least one cell.

- `site_count`:

  int; in **sites**; never `NA`. Distinct sites (protein, position,
  modification) in the table.

- `state_counts`:

  map\<string,int\>; in **cells**; never `NA`. Cells per state:
  quantified, floor, count_only, not_detected.

- `modifications`:

  string\[\]; `NA` when no `--modifications` filter. The modification
  names the sites were limited to.

- `categories`:

  string\[\]; `NA` when no `--categories` filter. The categories the
  sites were limited to.

- `ignored_columns`:

  string\[\]; never `NA`. Columns of the observations file that are not
  fields of mzLib's observation, so were not read; empty when every
  column was read.

- `row_count`:

  int; in **rows**; never `NA`. Cells in the table, whether returned or
  written.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order; also the header
  `--out` writes.

- `records` (wire `columns`):

  table; `NA` when the table was written to `--out` instead. One row per
  (site, run) where the site was identified or covered, in mzLib's
  order: by protein, then run, then position and modification.

- `written`:

  object; `NA` when `--out` was not given. {path, format: tsv,
  row_count} of the file `--out` wrote.

- `caveats`:

  string\[\]; never `NA`. Seven static caveats; one more with
  `--samples`, one when columns were not read, and one counting the
  censored (floor or not_detected) cells when there are any.

## Columns of `records`

- `site_key`:

  string; never `NA`. mzLib's ModificationSite.Key, accession:residue
  position:modification, for example P1:S13:Common
  Biological:Phosphorylation on S. The key ptmqtl pairs and ptmqtl
  trait-effects name sites by.

- `protein_accession`:

  string; never `NA`. The protein the site is on, as given.

- `position`:

  int; never `NA`. 1-based residue position in the protein. A protein
  N-terminal modification is at the residue it sits on (1, or 2 after
  the initiator Met is removed). readers read-occupancy counts
  differently (0 is the N-terminus): do not join the two on position.

- `residue`:

  string; never `NA`. The one-letter residue at position.

- `modification`:

  string; never `NA`. The modification as full_sequence names it,
  Category:IdWithMotif.

- `category`:

  string; never `NA`. The text before the first ':' of modification
  (mzLib's ModificationSite.Category); empty when the name has none.

- `run`:

  string; never `NA`. The run, or with `--samples` the sample.

- `state`:

  string; never `NA`. quantified (the modified form was quantified: a
  fraction in (0, 1\]) \| floor (the modified form was identified but
  not quantified while forms covering the site were: low but NOT zero)
  \| count_only (the modified form was identified, but nothing covering
  the site was quantified) \| not_detected (a covered-zero, mzLib's
  DEF-OCC-COVERED-ZERO: peptidoforms covering the site were identified
  in this run and none carried the modification; a bound below what the
  run could detect, not a measured 0; since mzLib \#1468 the covering
  forms need not be quantified).

- `modified_intensity`:

  float; in **intensity (instrument units)**; never `NA`. Summed
  intensity of the quantified peptidoforms carrying this modification at
  this position; 0 unless state is quantified.

- `covering_intensity`:

  float; in **intensity (instrument units)**; never `NA`. Summed
  intensity of every quantified peptidoform covering the position,
  modified or not; 0 for count_only.

- `fraction`:

  float; in **fraction (0 to 1)**; `NA` when state is not quantified:
  the occupancy is censored (floor), a bound (not_detected) or unknown
  (count_only), never 0. modified_intensity / covering_intensity.

- `unmodified_quantified`:

  bool; never `NA`. Whether some quantified form covering the position
  does not carry this modification there (the unmodified form, or
  another modification at the same position). false with fraction 1 is a
  ceiling: the site was seen only modified.

- `covering_peptidoforms`:

  int; in **peptidoforms**; `NA` when not supplied (Calculate always
  supplies it; a stored table may not). Distinct peptidoforms covering
  the position identified in this run (or sample), quantified or not:
  mzLib's SiteRunOccupancy.CoveringPeptidoforms (#1468). It counts
  identifications, not intensity: the n of a cell's k of n.

- `modified_peptidoforms`:

  int; in **peptidoforms**; `NA` when not supplied, as
  covering_peptidoforms. Of covering_peptidoforms, those carrying this
  modification at the position (SiteRunOccupancy.ModifiedPeptidoforms):
  the k. 0 for not_detected. A count_only cell with k = n was seen only
  modified; with k \< n both forms were identified.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  observations missing, blank or not found; samples given without a
  value or not found; out given without a value, not ending in .tsv, or
  the observations path; modifications and categories both given, either
  without a value or with a blank name

- `mzlib_usage_error` (usage):

  the observations file empty, header-only, missing a required column,
  with a repeated column name, a row whose field count differs from its
  header, a start_residue or end_residue that is not an integer, an
  intensity that is not a number (other than blank, NA, NaN), or a
  previous_residue longer than one letter; the samples file without two
  fields per row, with a blank or repeated run

- `mzlib_usage_error` (usage):

  mzLib refuses the observations, in its own words: an empty run,
  sequence or accession; start_residue below 1 or a span that does not
  fit the sequence's length; an infinite or negative intensity; one
  (run, peptidoform, protein, start) given twice; a run with no sample,
  or an empty sample

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Intensity occupancy, as MetaMorpheus writes IntensityOccupancy\_
  (QuantProject DEF-OCC-INT): for a site in a run, the summed intensity
  of the quantified peptidoforms carrying the modification at that
  position over the summed intensity of every quantified peptidoform
  covering it. Supply MS/MS identifications only: a match-between-runs
  transfer must not be a row.

- fraction is null unless state is quantified. A floor (identified, not
  quantified, while the site was covered) is low but not zero, and
  not_detected is a bound below detection: neither may be averaged in as
  0.

- A fraction of 1 with unmodified_quantified false is a ceiling: the
  unmodified form may be present below detection.

- position is 1-based, with a protein N-terminal modification at the
  residue it sits on; readers read-occupancy puts the N-terminus at 0.

- mzLib always excludes the categories Common Variable and Common Fixed,
  peptide-terminal modifications, and every C-terminal modification (the
  protein's length is not supplied). An N-terminal modification is a
  site only on a peptide starting at residue 1, or at 2 when
  previous_residue is M or blank. Whether it is protein- or
  peptide-terminal is read from mzLib's modification registry; a name
  the registry does not know is treated as protein N-terminal.

- A peptide shared by several proteins counts in full for each;
  MetaMorpheus gives it a proportional share per protein group, so its
  IntensityOccupancy\_ can differ.

- With samples, fractions are summed before the ratio. Map only
  fractions of one sample together; technical replicates keep their own
  run.

- covering_peptidoforms and modified_peptidoforms count identifications,
  not intensity: the k of n of a cell with no fraction. They tell a
  count_only cell seen only modified from one where both forms were
  identified, which is what ptmqtl pairs `--include-bounds` ranks as a
  high bound or leaves out.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.ptmqtl.occupancy`

- Rust (mzLibRust): `mzlib::ptmqtl::occupancy_with` with
  `OccupancyOptions`

- R (mzLibR): `ptmqtl_occupancy`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## See also

[`ptmqtl_pairs`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_pairs.md),
[`ptmqtl_trait_effects`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md),
[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)

## Examples

``` r

# mzLib's own hand-computed cases, one protein each:
cases <- system.file("extdata", "ptmqtl", "cases.tsv", package = "mzLibR")
occ <- ptmqtl_occupancy(cases)
occ
#> <mzlibr_site_occupancy> 9 sites in 6 runs, 14 cells
#>   quantified 9, floor 1, count_only 3, not_detected 1
occ$records[occ$records$protein_accession == "STATES", c("run", "state", "fraction")]
#>           run        state fraction
#> 6   countOnly   count_only       NA
#> 7       floor        floor       NA
#> 8 notDetected not_detected       NA
#> 9       quant   quantified      0.1
# RATIO is 30 / (30 + 70 + 100): a missed-cleavage form covering the site counts below the line.
occ$records[occ$records$protein_accession == "RATIO", c("modification", "fraction")]
#>                             modification fraction
#> 3 Common Biological:Phosphorylation on S     0.15
```
