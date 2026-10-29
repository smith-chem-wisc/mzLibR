# Every peptide in every biological sample, as log2 intensity or NA with the reason, by mzLib

Build mzLib's observation table from a stored FlashLFQ or MetaMorpheus
peptide table: every peptide in every biological sample as log2
intensity or no value with the reason, with the shared-peptide
normalization shifts, and optionally the peptide-by-sample table stats
fit reads.

mzLib's `FlashLfqObservations.FromPeptideTable` reads a stored FlashLFQ
or MetaMorpheus peptide table, maps each run's detection type to a
state, sums each technical replicate's fractions on the linear scale and
averages a sample's replicates in log2. Then
`SharedPeptideNormalization.Measure` measures how far each sample sits
from the others on the peptides they share; `normalize = TRUE` subtracts
it.

## Usage

``` r
diff_observations(peptides, runs, basis, protein_groups = NULL, samples = NULL,
  normalize = FALSE, min_peptides = NULL, out = NULL, responses_out = NULL,
  timeout = 600)
```

## Arguments

- peptides:

  A FlashLFQ `QuantifiedPeptides.tsv` or a MetaMorpheus
  `AllQuantifiedPeptides.tsv`.

- runs:

  Every run of the table, as a data.frame with `file` and `sample`
  columns, and optionally `fraction` and `technical_replicate` (whole
  numbers; `NA` or absent sends none, so mzLib's 0). `file` is the run's
  name as the `Intensity_` columns write it, without extension; a column
  with `-calib` after it matches too. The result's own `runs` is
  accepted as well.

- basis:

  `"msms_only"` (no match-between-runs value enters) or `"mbr_kept"`.
  Required: mzLib has no default.

- protein_groups:

  A MetaMorpheus `AllQuantifiedProteinGroups.tsv`, for each group's
  q-value and decoy, contaminant and entrapment flags. `NULL`: none
  described.

- samples:

  A TSV whose first column (after a header row) lists the samples to
  keep, e.g. one stratum's design file from
  [`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md).
  `NULL`: every sample.

- normalize:

  Subtract each sample's shift. The shifts are measured and reported
  either way.

- min_peptides:

  The floor on the normalization's reference set (peptides). `NULL`: not
  sent, so mzLib's own floor, 100 at mzLib `d10017ce`.

- out:

  A `.tsv` path to write the long table to instead of returning its
  rows.

- responses_out:

  A `.tsv` path to also write the peptide-by-sample log2 table that
  [`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  reads as `responses`.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_observations`:

\- `records`, a data.frame with one row per (peptide, sample), peptides
ordinal by full sequence: `full_sequence`, the engine's modified
sequence, the key; `base_sequence`; `protein_groups`, a list column of
every group the peptide table names (`UNDEFINED` = in no surviving
group); `unique_protein_group`, the one group it is unique to, `NA` when
shared, `UNDEFINED` or none; `sample`; `log2_intensity` (log2 of
FlashLFQ intensity units), fractions summed, replicates averaged in
log2, less the shift with `normalize = TRUE`, and `NA` for **no value**
(see `state`), never 0; and `state`, one of `"quantified"`,
`"mbr_transferred"`, `"ambiguous_peak"`, `"identified_not_quantified"`
or `"not_detected"`. Empty with `out`. - `peptides_file`,
`protein_groups_file` and `samples_file` as absolute paths, the last two
`NA` when not given; `basis` and `normalize` as given; `min_peptides`
(peptides) as given, `NA` when mzLib's floor applied; `run_count` runs;
`runs`, a data.frame of the runs as mzLib holds them (`file_name`,
`sample_id`, `fraction`, `technical_replicate`), ordinal by file name;
`sample_count` samples; `samples`; `samples_without_values`;
`peptide_count` peptides; `protein_group_count` protein groups;
`protein_groups`, a data.frame of those groups (`name`, `q_value` a
fraction 0 to 1 or `NA` if not finite, `is_decoy`, `is_contaminant`,
`is_entrapment`, `genes` or `NA`); `state_counts`, rows (cells) per
state; `warnings`; `written` (`path`, `format`, `row_count`) or `NULL`;
`responses_written` (`path`, `format`, `row_count` peptides,
`sample_count` samples) or `NULL`; `row_count` rows, `peptide_count` x
`sample_count`; `column_names`. - `normalization`, mzLib's
`NormalizationResult`, measured whether or not it was applied:
`setting`, `"shared_peptide_median"` (the peptides with a value in every
sample), `"shared_peptide_median_half"` (in at least half; used when too
few are in all) or `"none"` (too few either way: no sample shifted);
`reference_set_size` (peptides) in the set used;
`minimum_reference_set_size` (peptides), the floor used; `shifts`, each
sample's shift in log2 units, `NA` for a sample with no shift; and
`warnings`. - `caveats`: how values are combined, what the basis does,
what mzLib cannot check about the table, and what the normalization
assumes.

## The traps

Median normalization assumes most peptides do not change between
samples. When most of them do - a table of the spiked organism alone, a
loading series, a pull-down - normalizing deletes the result: read
`normalization$shifts` before you subtract them. The peptide table does
not record the design FlashLFQ quantified with, so give the runs that
design. mzLib does not yet roll peptides up to proteins or choose a
test's rows, so neither does this.

## Wraps

Wire verb `diff observations`. Generated from the bridge's verb spec
`diff.observations.yaml` (bridge commit `aea14604d5f7`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`FlashLfqObservations.FromPeptideTable`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/FlashLFQ/Differential/FlashLfqObservations.cs)
  in `mzLib/FlashLFQ/Differential/FlashLfqObservations.cs` at mzLib
  `d10017ce`

- [`FlashLfqObservations.ReadProteinGroups`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/FlashLFQ/Differential/FlashLfqObservations.cs)
  in `mzLib/FlashLFQ/Differential/FlashLfqObservations.cs` at mzLib
  `d10017ce`

- [`ObservationTable`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/ObservationTable.cs)
  in `mzLib/Quantification/Differential/ObservationTable.cs` at mzLib
  `d10017ce`

- [`ObservationTable.ForSamples`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/ObservationTable.cs)
  in `mzLib/Quantification/Differential/ObservationTable.cs` at mzLib
  `d10017ce`

- [`ObservationRun`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/ObservationTable.cs)
  in `mzLib/Quantification/Differential/ObservationTable.cs` at mzLib
  `d10017ce`

- [`QuantBasis`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/ObservationTable.cs)
  in `mzLib/Quantification/Differential/ObservationTable.cs` at mzLib
  `d10017ce`

- [`ObservationState`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/ObservationTable.cs)
  in `mzLib/Quantification/Differential/ObservationTable.cs` at mzLib
  `d10017ce`

- [`SharedPeptideNormalization.Measure`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/SharedPeptideNormalization.cs)
  in `mzLib/Quantification/Differential/SharedPeptideNormalization.cs`
  at mzLib `d10017ce`

- [`SharedPeptideNormalization.Apply`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/SharedPeptideNormalization.cs)
  in `mzLib/Quantification/Differential/SharedPeptideNormalization.cs`
  at mzLib `d10017ce`

- [`SharedPeptideNormalization.DefaultMinimumPeptides`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Quantification/Differential/SharedPeptideNormalization.cs)
  in `mzLib/Quantification/Differential/SharedPeptideNormalization.cs`
  at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `peptides`:

  path; required. A FlashLFQ QuantifiedPeptides.tsv, or a MetaMorpheus
  AllQuantifiedPeptides.tsv (mzLib picks the reader by that name
  ending). Its Intensity\_ and Detection Type\_ columns are the runs. An
  IsoTracker table is refused.

- `basis`:

  string; required; range `msms_only | mbr_kept`. Which
  match-between-runs values enter (mzLib's QuantBasis, by its machine
  name). msms_only: no transferred value enters a sum. mbr_kept: the
  transferred values FlashLFQ kept enter. Required: mzLib has no
  default, because it is a scientific choice.

- `runs` (wire `--stdin`):

  string\[\]; required; range
  `>= 1 line; file, sample, then optionally fraction and technical replicate, tab-separated`.
  The experimental design's runs, one per line: the file as the table's
  Intensity\_ columns name it, without its extension (a column with
  -calib after that name matches too); the biological sample it belongs
  to; then optionally its fraction and its technical replicate
  (integers; left out, ObservationRun's default 0). Every run the table
  has must be given, and only those. Grouping comes only from these
  lines: a file name is never parsed for it.

- `protein_groups` (wire `--protein-groups`):

  path; default absent. A MetaMorpheus AllQuantifiedProteinGroups.tsv,
  read by FlashLfqObservations.ReadProteinGroups for each group's
  q-value and its decoy, contaminant and entrapment flags. Absent, no
  group is described: mzLib then warns, and no peptide can be kept out
  of the normalization reference as a decoy or contaminant.

- `samples`:

  path; default absent. A TSV with a header row whose first column lists
  the samples to keep (ObservationTable.ForSamples), such as a sample
  sheet or a design file from diff design. Narrow to one stratum this
  way: the normalization is measured over every sample of the table it
  is given. Absent, every sample is kept.

- `normalize`:

  flag; default `FALSE`. Subtract each sample's shift from its log2
  values (SharedPeptideNormalization.Apply). The shifts are measured and
  reported either way.

- `min_peptides` (wire `--min-peptides`):

  int; in **peptides**; default absent; range `>= 1`. The floor on
  either normalization reference set (Measure's minimumPeptides).
  Absent, the bridge sends none and mzLib's
  SharedPeptideNormalization.DefaultMinimumPeptides applies: 100 at this
  pin. Below the floor on the every-sample set mzLib tries the
  half-samples set, and below it there too no sample is shifted (setting
  none).

- `out`:

  path; default absent; range `a .tsv path`. Write the long table here
  (header, then one line per row; protein_groups ';'-joined; an empty
  cell for null) instead of returning its rows. Absent, the rows are
  returned.

- `responses_out` (wire `--responses-out`):

  path; default absent; range `a .tsv path`. Also write the table's own
  peptide-by-sample shape: header full_sequence then the samples, one
  row per peptide, each cell ObservationTable.Log2Intensity in
  round-trip decimal (after `--normalize`, if given) or empty for no
  value. Exactly what stats fit `--responses` reads. Peptides only:
  mzLib does not yet roll peptides up to proteins or choose a test's
  rows.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `peptides_file`:

  string; never `NA`. Absolute path of the peptide table.

- `protein_groups_file`:

  string; `NA` when no `--protein-groups`. Absolute path of the
  protein-group table.

- `samples_file`:

  string; `NA` when no `--samples`: every sample was kept. Absolute path
  of the sample list.

- `basis`:

  string; never `NA`. msms_only or mbr_kept, as given.

- `normalize`:

  bool; never `NA`. Whether the shifts were subtracted (`--normalize`).

- `min_peptides`:

  int; in **peptides**; `NA` when not given: mzLib's
  DefaultMinimumPeptides applied
  (normalization.minimum_reference_set_size says which).
  `--min-peptides` as given.

- `run_count`:

  int; in **runs**; never `NA`. Runs in the table (after `--samples`).

- `runs`:

  object\[\]; never `NA`. ObservationTable.Runs, ordinal by file name:
  file_name, sample_id, fraction and technical_replicate, as mzLib holds
  them.

- `sample_count`:

  int; in **samples**; never `NA`. Biological samples: the distinct
  samples of the runs (after `--samples`).

- `samples`:

  string\[\]; never `NA`. The samples, ordinal by id: every sample, with
  or without values.

- `samples_without_values`:

  string\[\]; never `NA`. Samples in which no peptide has a value under
  this basis (a lost run).

- `peptide_count`:

  int; in **peptides**; never `NA`. Peptides (distinct full sequences);
  the table has peptide_count x sample_count rows.

- `protein_group_count`:

  int; in **protein groups**; never `NA`. Groups described by
  `--protein-groups`; 0 without it.

- `protein_groups`:

  object\[\]; never `NA`. Each described group: name (MetaMorpheus's,
  member accessions joined by \|), q_value (fraction 0 to 1; null if not
  finite), is_decoy, is_contaminant, is_entrapment, genes (null when the
  table has none).

- `state_counts`:

  map\<string,int\>; in **cells**; never `NA`. Rows per state, keyed by
  every state in mzLib's order: quantified, mbr_transferred,
  ambiguous_peak, identified_not_quantified, not_detected.

- `normalization`:

  object; never `NA`. mzLib's NormalizationResult, measured whether or
  not `--normalize` applied it; fields under
  result.tables.normalization.

- `warnings`:

  string\[\]; never `NA`. ObservationTable.Warnings: what the reader
  should know that did not stop the build, e.g. protein groups the
  peptides name that the protein-group table does not describe.

- `written`:

  object; `NA` when no `--out`: the rows are in columns. path
  (absolute), format (tsv) and row_count of the long table written.

- `responses_written`:

  object; `NA` when no `--responses-out`. path (absolute), format (tsv),
  row_count (peptides) and sample_count of the responses table written.

- `row_count`:

  int; in **rows**; never `NA`. peptide_count x sample_count, returned
  or written.

- `column_names`:

  string\[\]; never `NA`. The table's columns, in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (peptide, sample): peptides ordinal by
  full sequence, each peptide's samples in order. Empty with `--out`.

- `caveats`:

  string\[\]; never `NA`. Six static caveats, plus one without
  `--protein-groups` and one when the normalization setting is none.

## Columns of `records`

- `full_sequence`:

  string; never `NA`. The engine's full (modified) sequence: the
  peptide's key, unique in the table.

- `base_sequence`:

  string; never `NA`. The unmodified sequence.

- `protein_groups`:

  string\[\]; never `NA`. Every group the peptide table's Protein Groups
  cell names, distinct and ordinal-sorted; UNDEFINED is MetaMorpheus's
  placeholder for a peptide in no surviving group; empty when the cell
  is blank.

- `unique_protein_group`:

  string; `NA` when shared (two or more groups), only UNDEFINED, or no
  group. The one group the peptide is unique to
  (ObservedPeptide.UniqueProteinGroup). Sharing counts against every
  group written, whatever its q-value.

- `sample`:

  string; never `NA`. The biological sample.

- `log2_intensity`:

  float; in **log2 of FlashLFQ intensity units**; `NA` when no value:
  state says why. log2 of the sample's intensity: each technical
  replicate's fractions summed on the linear scale, then the replicates'
  log2 values averaged; minus the sample's shift with `--normalize`.
  Never 0 for missing.

- `state`:

  string; never `NA`. Where the value came from, or why there is none
  (mzLib's ObservationState; a sample without a value takes the most
  informative state among its runs): quantified (MS/MS-identified peaks
  only) \| mbr_transferred (a match-between-runs peak is in the value;
  under msms_only the value is left out and log2_intensity is null) \|
  ambiguous_peak (a peak withheld as ambiguous; no value) \|
  identified_not_quantified (identified by MS/MS, no peak; no value) \|
  not_detected (no value).

## Fields of `normalization`

- `normalization$setting` (wire `setting`):

  string; never `NA`. Which reference set the shifts were measured on:
  shared_peptide_median (peptides with a value in every sample with
  values) \| shared_peptide_median_half (in at least half of them and at
  least 2; used when the every-sample set is under the floor) \| none
  (both under the floor: no sample shifted).

- `normalization$reference_set_size` (wire `reference_set_size`):

  int; in **peptides**; never `NA`. Peptides in the reference set used,
  or in the last set tried when setting is none. Peptides of decoy,
  contaminant or entrapment groups never enter it.

- `normalization$minimum_reference_set_size` (wire
  `minimum_reference_set_size`):

  int; in **peptides**; never `NA`. The floor used: `--min-peptides`,
  else mzLib's DefaultMinimumPeptides.

- `normalization$shifts` (wire `shifts`):

  map\<string,float\>; in **log2 units**; never `NA`. Each sample's
  shift, subtracted from its values with `--normalize`: the median, over
  the reference peptides it has, of its value minus the peptide's mean
  over the samples that have it. On the half set mzLib repeats that
  until the shifts settle and centres them to average 0 within each
  group of linked samples. A value is null when undefined (no reference
  peptide has a value in that sample, or setting none); that sample is
  not shifted.

- `normalization$warnings` (wire `warnings`):

  string\[\]; never `NA`. Which samples have no shift and why, and
  whether the half set did not settle or fell into unlinked groups.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--peptides` or `--basis` missing; `--basis` not msms_only or
  mbr_kept; `--protein-groups`, `--samples`, `--out` or
  `--responses-out` given without a value; a file not found, or not the
  table type named; `--out` or `--responses-out` not ending in .tsv;
  `--min-peptides` below 1 or not an integer

- `mzlib_usage_error` (usage):

  no runs on stdin; a run line with other than 2 to 4 fields, an empty
  file or sample, or a fraction or technical replicate that is not an
  integer; the `--samples` file empty, header-only or with an empty
  first cell

- `mzlib_usage_error` (usage):

  mzLib refuses the input (its ArgumentException, reclassified): two
  runs with one file name or one place (sample, replicate, fraction); a
  table column no run matches, or a run no column matches; IsoTracker
  output; a blank, negative or non-finite intensity or an unreadable
  detection type; a peptide listed twice; a group name holding ';'; a
  `--samples` id the table does not have

- `mzlib_bridge_error` (correctness):

  never expected

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- A missing value is null, never 0, and state says why. Fractions of one
  technical replicate are summed on the linear scale, then a sample's
  technical replicates are averaged on the log2 scale (mzLib's
  ObservationTable).

- Under msms_only a match-between-runs value enters no sum, but its cell
  still says mbr_transferred; under mbr_kept it enters, and the cell
  says mbr_transferred.

- The peptide table does not record the design FlashLFQ quantified with,
  so mzLib cannot check that these runs group the files as FlashLFQ did.
  FlashLFQ withholds every fraction of a sample whose strongest fraction
  is ambiguous; give the runs the design it quantified with, or a value
  it withheld can become a partial sum.

- FlashLFQ writes its peptide table in the culture of the machine it ran
  on and mzLib reads it as invariant (a decimal point), so a table
  written where the decimal separator is a comma does not read
  correctly.

- Every peptide is a row: decoy, contaminant and entrapment peptides,
  and peptides of groups at any q-value, included. mzLib does not yet
  choose the rows a differential test uses, or roll peptides up to
  proteins, so `--responses-out` is peptide by sample only. Filter rows
  with protein_groups and unique_protein_group.

- The normalization shifts are always measured
  (SharedPeptideNormalization.Measure) over every sample of the table,
  after `--samples`; they are subtracted only with `--normalize`.
  Measure one stratum at a time. A sample whose shift is null is not
  shifted.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.diff.observations`

- Rust (mzLibRust): `mzlib::diff::observations_with` with
  `ObservationsOptions`

- R (mzLibR): `diff_observations`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## References

- [doi:10.1021/acs.jproteome.7b00608](https://doi.org/10.1021/acs.jproteome.7b00608):
  FlashLFQ, whose peptide table and detection types this verb reads
  (Millikin et al., J Proteome Res 2018)

- [doi:10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106):
  the median-of-ratios idea the shared-peptide normalization follows:
  each sample against the log of a geometric-mean reference (Anders and
  Huber, Genome Biol 11:R106, 2010)

## See also

[`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md),
[`stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md),
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)

## Examples

``` r

# Two runs of FlashLFQ's E. coli spike-in vignette, which their names say were spiked at 1x and
# at 2x:
peptides <- system.file("extdata", "diff", "EcoliSpikeIn_QuantifiedPeptides.tsv",
  package = "mzLibR")
runs <- data.frame(
  file = c("09-04-18_EcoliSpikeInSingleShot1x", "09-04-18_EcoliSpikeInSingleShot2x"),
  sample = c("spike_1x", "spike_2x"))
table <- diff_observations(peptides, runs, basis = "msms_only")
table
#> <mzlibr_observations> 243 peptides x 2 samples, msms_only
#>   quantified 319, mbr_transferred 43, ambiguous_peak 9, identified_not_quantified 6, not_detected 109
#>   normalization: none on 86 peptides
table$state_counts
#>                quantified           mbr_transferred            ambiguous_peak 
#>                       319                        43                         9 
#> identified_not_quantified              not_detected 
#>                         6                       109 
table$normalization[c("setting", "reference_set_size")]
#> $setting
#> [1] "none"
#> 
#> $reference_set_size
#> [1] 86
#> 
```
