# Quantify a search's peptides across mzML runs with FlashLFQ

Quantify a search's peptides and proteins across mzML runs with
FlashLFQ, returning per-run peptide and protein intensities and every
chromatographic peak.

Quantify a search's peptides across mzML runs with FlashLFQ

## Usage

``` r
flashlfq_quantify(psms, spectra, normalize = FALSE, ppm_tolerance = 10,
  isotope_ppm_tolerance = 5, integrate = FALSE, match_between_runs = FALSE,
  mbr_ppm_tolerance = 10, mbr_q_value_threshold = 0.05,
  use_shared_peptides_for_protein_quant = FALSE, bayesian_protein_quant = FALSE,
  use_pep_q_value = FALSE, max_threads = 1, output_directory = NULL, timeout = NULL)
```

## Arguments

- psms:

  Path to a PSM result file. \*\*Use a MetaMorpheus `.psmtsv` or
  `.osmtsv`.\*\* Every run named in the file must have a matching mzML
  in `spectra`; FlashLFQ matches identifications to runs by base file
  name.

- spectra:

  The mzML runs. Either a character vector of paths, or a data.frame
  with a `path` column and optional `condition`, `biological_replicate`,
  `technical_replicate` and `fraction` columns. Omitted design fields
  get the defaults MetaMorpheus uses with no experimental-design file.

  **mzML only.** `.raw` and `.d` are rejected up front; convert them
  first.

- normalize:

  Whether to normalise between runs.

- ppm_tolerance:

  Mass tolerance for peak-finding, in ppm.

- isotope_ppm_tolerance:

  Mass tolerance for isotope-envelope matching, in ppm.

- integrate:

  Whether to integrate peak areas rather than take apex intensity.

- match_between_runs:

  Whether to transfer identifications between runs.

  MBR needs a **complete, balanced design** to work properly, and
  `mbr_q_value_threshold` is its FDR control - without it roughly 80% of
  transfers are false. Set `condition` and `biological_replicate` on
  `spectra` so FlashLFQ knows which runs are comparable.

  \*\*Whatever you do, count transfers from `peaks`, never from
  `peptides`.\*\* See the return value.

- mbr_ppm_tolerance:

  Mass tolerance for match-between-runs, in ppm.

- mbr_q_value_threshold:

  FDR threshold for accepting a transfer, as a q-value.

- use_shared_peptides_for_protein_quant:

  Whether peptides shared between protein groups contribute to protein
  quantification.

  `FALSE` is the default and it is the main reason protein intensities
  come back as **0**: on mzLib's K562 pair, **847** of 943 protein
  groups are 0 in both runs, mostly because their only evidence is
  shared peptides.

- bayesian_protein_quant:

  Whether to use the Bayesian protein quantification model.

- use_pep_q_value:

  Whether to use PEP q-values rather than q-values for filtering.

- max_threads:

  Worker threads, or `-1` for one per core. **Defaults to 1 here, where
  the bridge, pyMzLib and mzLibRust default to -1.**

  One thread is reproducible by construction, and it is what this
  package has always shipped. The multithreaded nondeterminism that made
  it the default here (the peptide roll-up dropping MBR intensities,
  smith-chem-wisc/mzLib#1111) was fixed by mzLib#1155, which the pinned
  bridge includes, so `-1` is expected to give the same answer faster;
  that has not yet been re-measured on mzLib's K562 pair, which is why
  the default has not moved.

- output_directory:

  Where FlashLFQ should write its TSV output, or `NULL` to write none.

- timeout:

  Seconds to allow, or `NULL` to wait as long as it takes.

## Value

An `mzlibr_quant`: `psm_file`, `identification_count`, `parameters`,
`output_directory`, and four tidy data.frames - `spectra_files`,
`peptides`, `proteins` and `peaks`.

`spectra_files` has one row per run, with `peak_count` and
`mbr_peak_count` in peaks. `peptides` and `proteins` are long - one row
per peptide or protein group per run - with `intensity` in the
instrument's intensity units. `peaks` has one row per chromatographic
peak: `intensity` in instrument units, `retention_time` of the apex in
minutes, and `num_identifications`, the identifications that explain the
peak (more than one means ambiguous).

## Read peaks, not the peptide roll-up

`peptides` is FlashLFQ's roll-up and **it drops most match-between-runs
transfers.** On mzLib's own K562 pair there are **140** true transfers
in `peaks` and the peptide table shows **52** - a 63% under-count.
Worse, it is not evenly spread: per run the peaks give run_3 **62** and
run_4 **78**, while the roll-up gives run_3 **0** and run_4 52. Read
only the roll-up and MBR appears not to have worked at all in half the
experiment.

"Peptides quantified in both runs" is **257** from `peaks` and **169**
from `peptides`.

Nor can you reproduce the roll-up by pivoting `peaks` yourself: where a
run has several peaks for one peptide the roll-up reports one rather
than their sum.

## `proteins$file_name` is a sample, not a file

FlashLFQ measures peptides in **files** but resolves proteins across
**samples**, grouping runs by condition and biological replicate before
the median-polish roll-up. `flashlfq_quantify()` gives every run its own
sample, so `proteins$file_name` and `peptides$file_name` carry the same
run base names and the two frames join cleanly on it.

The distinction only bites if runs are ever grouped into replicates:
`proteins$file_name` would then hold a sample label
(`"condition_replicate"`) covering several runs, while
`peptides$file_name` stayed per run. The column name is kept for
symmetry between the two frames; read it as "the thing this intensity
was measured over".

## What 0 and NA mean, and which is rare

`peptides$intensity` is **0** when the peptide was not measured in that
run. It is never NA.

`proteins$intensity` is **NA** when FlashLFQ could not resolve a number
at all - its median-polish produced NaN. Arithmetic propagates it, so
[`mean()`](https://rdrr.io/r/base/mean.html) on a protein column returns
NA rather than a confidently wrong number, and `na.rm = TRUE` is a
choice you make visibly. mzLibR never applies it on your behalf.

**NA is the rare outcome and 0 is the common one.** On the K562 pair,
**2** proteins are NA and **847** are 0 in both runs. "No usable number"
is 849; "could not be resolved" is 2.

## Wraps

Wire verb `quant flashlfq`. Generated from the bridge's verb spec
`quant.flashlfq.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`FileReader.ReadQuantifiableResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`MzLibExtensions.MakeIdentifications`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/ResultsReading/MzLibExtensions.cs)
  in `mzLib/FlashLFQ/ResultsReading/MzLibExtensions.cs` at mzLib
  `23c2490e`

- [`FlashLfqParameters`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/FlashLfqParameters.cs)
  in `mzLib/FlashLFQ/FlashLfqParameters.cs` at mzLib `23c2490e`

- [`FlashLfqEngine.Run`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/FlashLfqEngine.cs)
  in `mzLib/FlashLFQ/FlashLfqEngine.cs` at mzLib `23c2490e`

- [`FlashLfqResults.WriteResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/FlashLFQResults.cs)
  in `mzLib/FlashLFQ/FlashLFQResults.cs` at mzLib `23c2490e`

- [`SpectraFileInfo`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/MassSpectrometry/SpectraFileInfo.cs)
  in `mzLib/MassSpectrometry/SpectraFileInfo.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `psms`:

  path; required. A quantifiable PSM result file mzLib can read: a
  MetaMorpheus .psmtsv/.osmtsv or an MSFragger psm.tsv (retention time
  converted to minutes by the reader since mzLib#1116). Every run it
  names must have an mzML on stdin, matched by base file name.

- `spectra` (wire `--stdin`):

  string\[\]; required; range `>= 1 line`. The mzML runs, one per line,
  tab-separated: path\[TAB condition\[TAB biorep\[TAB techrep\[TAB
  fraction\]\]\]\]. Blank lines are ignored. Omitted fields default as
  MetaMorpheus does with no design: condition "", biorep = the line's
  0-based index, techrep 0, fraction 0. Replicate and fraction numbers
  are non-negative integers. Base names must be unique.

- `normalize`:

  flag; default `FALSE`. Normalize intensities across runs (FlashLFQ
  Normalize).

- `ppm_tolerance` (wire `--ppm`):

  float; in **ppm**; default `10`; range `> 0`. Mass tolerance for MS1
  peak-finding (PpmTolerance), in ppm.

- `isotope_ppm_tolerance` (wire `--isotope-ppm`):

  float; in **ppm**; default `5`; range `> 0`. Mass tolerance for
  isotope-envelope matching (IsotopePpmTolerance), in ppm.

- `integrate`:

  flag; default `FALSE`. Integrate peak area instead of taking the apex
  intensity (Integrate). FlashLFQ recommends off.

- `match_between_runs` (wire `--mbr`):

  flag; default `FALSE`. Match-between-runs: quantify a peptide in a run
  where it was not identified by transferring the identification
  (MatchBetweenRuns). Needs a complete, balanced design.

- `mbr_ppm_tolerance` (wire `--mbr-ppm`):

  float; in **ppm**; default `10`; range `> 0`. Mass tolerance for MBR
  transfers (MbrPpmTolerance), in ppm.

- `mbr_q_value_threshold` (wire `--mbr-q`):

  float; in **q-value**; default `0.05`; range `0 to 1`. q-value below
  which an MBR transfer is accepted (MbrQValueThreshold).

- `use_shared_peptides_for_protein_quant` (wire `--shared-peptides`):

  flag; default `FALSE`. Let peptides shared between protein groups
  contribute to protein quant (UseSharedPeptidesForProteinQuant).

- `bayesian_protein_quant` (wire `--bayesian`):

  flag; default `FALSE`. Also run FlashLFQ's Bayesian protein
  fold-change engine (BayesianProteinQuant). Its results are not on the
  wire; they reach only out's BayesianProteinQuant.tsv.

- `use_pep_q_value` (wire `--use-pep-q`):

  flag; default `FALSE`. Filter identifications on PEP q-value instead
  of q-value when building them (MakeIdentifications usePepQValue).

- `max_threads` (wire `--threads`):

  int; in **threads**; default `-1`; range
  `any integer; -1 means one fewer than the logical cores`. FlashLFQ
  MaxThreads. mzLib resolves it in the engine constructor: -1, or any
  value \>= the core count, becomes (cores - 1); anything that leaves
  \<= 0 becomes 1. parameters.max_threads echoes the RESOLVED value, so
  it depends on the machine. See caveats for what the thread count does
  to answers at this pin.

- `output_directory` (wire `--out`):

  path; default absent. Directory (created if absent) where FlashLFQ
  also writes QuantifiedPeaks.tsv, QuantifiedPeptides.tsv,
  QuantifiedProteins.tsv and, with bayesian, BayesianProteinQuant.tsv.
  Absent or blank writes nothing.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `psm_file`:

  string; never `NA`. Absolute path of the PSM file.

- `identification_count`:

  int; in **identifications**; never `NA`. Identifications FlashLFQ was
  given (after the q-value filter mzLib's converter applies).

- `parameters`:

  object; never `NA`. The FlashLFQ parameters used, by mzLib name:
  normalize, ppm_tolerance (ppm), isotope_ppm_tolerance (ppm),
  integrate, match_between_runs, mbr_ppm_tolerance (ppm),
  mbr_q_value_threshold (q-value),
  use_shared_peptides_for_protein_quant, bayesian_protein_quant,
  max_threads (the resolved thread count).

- `spectra_files`:

  object\[\]; never `NA`. One per stdin line, in stdin order; fields
  under result.tables.spectra_files.

- `peptide_count`:

  int; in **peptides**; never `NA`. Entries in peptides: distinct
  modified sequences quantified.

- `protein_count`:

  int; in **protein groups**; never `NA`. Entries in proteins.

- `peptides`:

  object\[\]; never `NA`. One per modified sequence, mirroring
  QuantifiedPeptides.tsv; fields under result.tables.peptides.

- `proteins`:

  object\[\]; never `NA`. One per protein group, mirroring
  QuantifiedProteins.tsv; fields under result.tables.proteins.

- `peaks`:

  object\[\]; never `NA`. Every chromatographic peak, run by run in
  stdin order, mirroring QuantifiedPeaks.tsv; fields under
  result.tables.peaks. The complete surface for match-between-runs.

- `output_directory`:

  string; `NA` when out was not given (or was blank): no files were
  written. Absolute path of out.

## Fields of `spectra_files`

- `file_name`:

  string; never `NA`. Run base name (no directory or extension): the key
  of every intensities map.

- `full_path`:

  string; never `NA`. The mzML path as given on stdin.

- `condition`:

  string; never `NA`. Sample-group label; "" when none was given.

- `biological_replicate`:

  int; never `NA`. 0-based.

- `technical_replicate`:

  int; never `NA`. 0-based.

- `fraction`:

  int; never `NA`. 0-based.

- `peak_count`:

  int; in **peaks**; never `NA`. Chromatographic peaks quantified in
  this run.

- `mbr_peak_count`:

  int; in **peaks**; never `NA`. Of those, peaks with detection type MBR
  (transferred into this run). 0 unless mbr.

## Fields of `peptides`

- `sequence`:

  string; never `NA`. Full modified sequence as FlashLFQ renders it: the
  quantified identity.

- `base_sequence`:

  string; never `NA`. Bare residues.

- `protein_groups`:

  string; never `NA`. Protein group names, ';'-joined, distinct.

- `intensity` (wire `intensities`):

  map\<string,float\>; in **intensity (instrument units)**; never `NA`.
  Run base name -\> intensity; 0 means not quantified in that run. A
  value is null only if FlashLFQ produced a non-finite number (the
  bindings read it as 0). Where a peptide has several peaks in a run
  this is ONE of them, not their sum.

- `detection_type` (wire `detection_types`):

  map\<string,string\>; never `NA`. Run base name -\> FlashLFQ
  DetectionType name: MSMS, MBR, MSMSIdentifiedButNotQuantified,
  MSMSAmbiguousPeakfinding, NotDetected, ...

## Fields of `proteins`

- `protein_group`:

  string; never `NA`. Protein group name.

- `gene_name`:

  string; never `NA`. As the PSM file gave it; may be "".

- `organism`:

  string; never `NA`. As the PSM file gave it; may be "".

- `intensity` (wire `intensities`):

  map\<string,float\>; in **intensity (instrument units)**; never `NA`.
  Run base name -\> protein intensity. A value is null when median
  polish could not resolve one (NaN: a degenerate peptide matrix, e.g.
  too few peptides per run); 0 means not measured. Keyed per run even
  when a design groups runs into samples.

## Fields of `peaks`

- `file_name`:

  string; never `NA`. Run base name.

- `base_sequence`:

  string; never `NA`. Of the peak's first identification; "" if it has
  none.

- `sequence`:

  string; never `NA`. Modified sequence of the first identification; ""
  if none.

- `intensity`:

  float; in **intensity (instrument units)**; `NA` when FlashLFQ
  produced a non-finite intensity for the peak. Apex intensity, or the
  integrated area under integrate.

- `detection_type`:

  string; never `NA`. DetectionType name; MBR marks a transferred peak.

- `retention_time`:

  float; in **min**; `NA` when the peak has no apex (or its time is
  non-finite). Apex retention time.

- `num_identifications`:

  int; in **identifications**; never `NA`. Identifications that explain
  the peak; \> 1 means ambiguous.

- `protein_groups`:

  string; never `NA`. Protein groups of all its identifications,
  ';'-joined, distinct.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  psms missing or blank; the PSM file not found; nothing on stdin; a
  stdin line with no path; a run that is not .mzML; an mzML not found;
  two runs sharing a base name; a design field that is not a
  non-negative integer

- `mzlib_usage_error` (usage):

  the PSM file names runs with no mzML on stdin (checked before any
  spectra are read; the message lists up to five); ppm, isotope-ppm,
  mbr-ppm or mbr-q not a number; threads not an integer

- `mzlib_bridge_error` (correctness):

  the PSM file is of a type mzLib cannot read as quantifiable results
  (MzLibException), or is malformed; an mzML fails to parse; the
  FlashLFQ engine fails

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- THREADS AT THIS PIN: mzLib#1111 (peptide roll-up nondeterministically
  dropping MBR intensities at MaxThreads -1, flipping protein groups
  between null and a number) is fixed by mzLib#1155, which is in the
  pin: PEP training rows are now built in a fixed order, and its commit
  reports identical peptides and proteins across runs. The bindings'
  'set threads to 1 to reproduce' warnings predate the fix. The bridge
  project has not re-measured the K562 case at -1 on this pin; until it
  does, 1 remains the conservative choice for published numbers.

- mzML only: .raw and .d are rejected up front; convert them first.

- For match-between-runs read peaks, not peptides: the peptide roll-up
  (as in QuantifiedPeptides.tsv) reports far fewer MBR transfers than
  happened (K562 pair: 140 MBR peaks, 52 MBR peptide entries).

- Peptide intensity 0 means not quantified; protein intensity null means
  unquantifiable. Neither is a measured zero.

- Only FlashLFQ settings listed here are reachable; the rest (isotopes
  required, MBR RT window, donor criterion, IsoTracker, Bayesian
  sampling) run at mzLib's defaults.

- stdout carries only the envelope: the engine's console output is
  redirected to stderr for the run.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.flashlfq.quantify`

- Rust (mzLibRust): `mzlib::flashlfq::quantify_with` with
  `QuantifyOptions`

- R (mzLibR): `flashlfq_quantify`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Re-measure the K562 MBR pair at threads -1 on the 23c2490e bridge
  (G-flashlfq-stale-caveats) before any binding deletes its \#1111
  warning; the caveat above states what mzLib#1155 claims, not a bridge
  measurement.

- Default divergence: the wire and pyMzLib and mzLibRust default threads
  to -1, mzLibR defaults max_threads to 1 and warns on anything else.
  The spec states the wire default; whether the bindings should agree is
  a cross-binding decision.

- flashlfq_small.json is hand-shaped (paths '/abs/...', 4
  identifications), not recorded from a live bridge; a live recording
  would make the example executable against real files.

- ppm tolerances and mbr-q are not range-checked by the bridge (a
  negative or zero tolerance reaches FlashLFQ); the ranges above are
  FlashLFQ's meaningful ones, not enforced ones. An option given with no
  value is silently defaulted by the Arguments parser.

- A PSM file of an unsupported or non-quantifiable type crosses as
  correctness (mzLib's MzLibException from
  FileReader.ReadQuantifiableResultFile), though it is arguably the
  caller's mistake.

- Most wire params project under FlashLFQ names in every binding (psms,
  stdin -\> spectra, ppm -\> ppm_tolerance, isotope-ppm -\>
  isotope_ppm_tolerance, mbr -\> match_between_runs, mbr-ppm -\>
  mbr_ppm_tolerance, mbr-q -\> mbr_q_value_threshold, shared-peptides
  -\> use_shared_peptides_for_protein_quant, bayesian -\>
  bayesian_protein_quant, use-pep-q -\> use_pep_q_value, threads -\>
  max_threads, out -\> output_directory); pyMzLib's spec lint will need
  PYTHON_DEVIATIONS entries for them.

## References

- [doi:10.1021/acs.jproteome.7b00608](https://doi.org/10.1021/acs.jproteome.7b00608):
  FlashLFQ: the peak-finding, match-between-runs and protein
  quantification this verb runs (Millikin et al., J Proteome Res 2018)

## See also

[`flashlfq_mbr_peaks`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_peaks.md),
[`flashlfq_mbr_rescued_peptide_count`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_rescued_peptide_count.md)

## Examples

``` r

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
  match_between_runs = TRUE)
quant
#> <mzlibr_quant> AllPSMs.psmtsv, 4 identifications
#>   2 runs, 2 peptides, 2 protein groups, 3 peaks
#>   MBR: 1 transfers in peaks, 1 distinct peptides
#>   proteins: 2 NA (could not be resolved), 0 zero (not measured)
quant$spectra_files[, c("file_name", "peak_count", "mbr_peak_count")]
#>   file_name peak_count mbr_peak_count
#> 1     run_3          3              1
#> 2     run_4          2              0
# Count transfers from the peaks, never from the peptide roll-up.
flashlfq_mbr_peaks(quant)
#>   file_name sequence base_sequence intensity detection_type retention_time
#> 2     run_4 PEPTIDEK      PEPTIDEK       500            MBR           30.3
#>   num_identifications protein_groups
#> 2                   1         P12345
```
