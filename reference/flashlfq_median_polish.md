# Roll a FlashLFQ peptide table up to protein intensities, under a new design

Roll a FlashLFQ QuantifiedPeptides.tsv up to protein intensities per
sample with FlashLFQ's median polish, under a new experimental design,
without re-reading any spectra.

The second half of
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
on its own. Given a `QuantifiedPeptides.tsv` FlashLFQ already wrote, it
rebuilds FlashLFQ's peptide and protein graph and runs the same
median-polish protein quantification, without re-reading any spectra.
Reach for it to re-quantify proteins under a different experimental
design, or with shared peptides toggled, without paying for peak-finding
again.

## Usage

``` r
flashlfq_median_polish(peptides, design = NULL, use_shared_peptides = FALSE,
  output_directory = NULL, timeout = NULL)
```

## Arguments

- peptides:

  Path to a FlashLFQ `QuantifiedPeptides.tsv` - the file
  [`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  writes into `output_directory`. It must carry `Sequence`,
  `Base Sequence`, `Protein Groups` and one `Intensity_<run>` column per
  run; a blank intensity cell is read as 0 (not measured).

- design:

  The experimental design: a data.frame with one row per run, a
  `file_name` column matching each `Intensity_<file_name>` column, and
  optional `condition`, `biological_replicate`, `technical_replicate`
  and `fraction` columns (0-based whole numbers); or a character vector
  of run names. **Median polish groups runs by condition and biological
  replicate, so this is how you say which runs are replicates of which
  sample** - the whole reason to re-run it. A design must name every run
  in the table and only those. `NULL`, the default, makes each run its
  own biological replicate with a blank condition, which is what
  FlashLFQ assumes when it writes the file with no design.

- use_shared_peptides:

  Whether peptides shared between protein groups contribute. `FALSE` by
  default, when a group with only shared peptides quantifies to 0.

- output_directory:

  Where FlashLFQ should also write `QuantifiedProteins.tsv`, or `NULL`
  to write nothing. Its columns are labelled by the same rule as
  `samples$label`.

- timeout:

  Seconds to allow, or `NULL` to wait as long as it takes.

## Value

An `mzlibr_median_polish`: `peptides_file`, `parameters`,
`output_directory` (`NA` when nothing was written), `peptide_count` -
distinct peptides read from the table - and `protein_count`, in protein
groups; `samples`, a data.frame with one row per sample (`label`,
`condition`, `biological_replicate`, 0-based); and `proteins`, long, one
row per protein group per sample, with `protein_group`, `gene_name`,
`organism`, `sample` (a `samples$label`) and `intensity` in the
instrument's intensity units.

## What 0 and NA mean

A protein `intensity` is **NA** where median polish could not resolve a
number - a degenerate peptide matrix - and **0** where the protein was
not measured in that sample. They are different facts, the same two
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
reports. A sample of several runs (fractions, technical replicates)
reports their sum.

## Sample labels

A sample is labelled by its run name when no design groups runs, and
`"<condition>_<n>"` once one does, where `n` is the biological replicate
plus one: `control_1`. Runs are names, never opened as files: nothing
but the peptide table is read, so peak-level facts - match-between-runs
peaks, retention times - are not available here. Use
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
for those.

## Wraps

Wire verb `quant median-polish`. Generated from the bridge's verb spec
`quant.median-polish.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`FlashLfqResults.CalculateProteinResultsMedianPolish`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/FlashLFQResults.cs)
  in `mzLib/FlashLFQ/FlashLFQResults.cs` at mzLib `23c2490e`

- [`FlashLfqResults.WriteResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/FlashLFQResults.cs)
  in `mzLib/FlashLFQ/FlashLFQResults.cs` at mzLib `23c2490e`

- [`ProteinGroup`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/FlashLFQ/ProteinGroup.cs)
  in `mzLib/FlashLFQ/ProteinGroup.cs` at mzLib `23c2490e`

- [`SpectraFileInfo`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/MassSpectrometry/SpectraFileInfo.cs)
  in `mzLib/MassSpectrometry/SpectraFileInfo.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `peptides`:

  path; required. A FlashLFQ QuantifiedPeptides.tsv. Columns are found
  by header name: Sequence, Base Sequence and Protein Groups are
  required; Gene Names (or Gene Name) and Organism are optional; one
  Intensity\_\<run\> per run, each optionally with a Detection
  Type\_\<run\>. A blank intensity cell is 0 (not measured).

- `design` (wire `--stdin`):

  string\[\]; default absent; range
  `none, or exactly one line per Intensity_ run`. The experimental
  design, one run per line, tab-separated: run_name\[TAB condition\[TAB
  biorep\[TAB techrep\[TAB fraction\]\]\]\], run_name matching an
  Intensity\_\<run_name\> column (a name, never opened as a file).
  Omitted fields: condition "", biorep = the line's 0-based index,
  techrep 0, fraction 0. With no lines at all, every Intensity\_ column
  is its own biological replicate with a blank condition, as FlashLFQ
  assumes with no design file. Median polish groups runs by condition
  and biological replicate, so this is how a caller says which runs are
  replicates of which sample.

- `use_shared_peptides` (wire `--shared-peptides`):

  flag; default `FALSE`. Let peptides shared between protein groups
  contribute (UseSharedPeptidesForProteinQuant). Off: a group with only
  shared peptides quantifies to 0.

- `output_directory` (wire `--out`):

  path; default absent. Directory (created if absent) where FlashLFQ
  also writes QuantifiedProteins.tsv. Absent or blank writes nothing.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `peptides_file`:

  string; never `NA`. Absolute path of the peptide table.

- `parameters`:

  object; never `NA`. {use_shared_peptides_for_protein_quant}.

- `samples`:

  object\[\]; never `NA`. One per (condition, biological replicate),
  ordered by condition (ordinal) then biorep, as the engine orders them;
  fields under result.tables.samples. Their labels are the keys of every
  protein's intensities.

- `peptide_count`:

  int; in **peptides**; never `NA`. Distinct sequences read from the
  table (rows with a Sequence).

- `protein_count`:

  int; in **protein groups**; never `NA`. Entries in proteins.

- `proteins`:

  object\[\]; never `NA`. One per protein group, ordered by
  protein_group (ordinal); fields under result.columns.

- `output_directory`:

  string; `NA` when out was not given (or was blank): nothing was
  written. Absolute path of out.

## Columns of `proteins`

- `protein_group`:

  string; never `NA`. Protein group name; a ';'-joined Protein Groups
  cell is split into one group per name.

- `gene_name`:

  string; never `NA`. From the table's gene column, paired to the group
  by position when the counts match, else the whole ';'-joined cell; ""
  when the table has no gene column.

- `organism`:

  string; never `NA`. As gene_name, from the Organism column.

- `intensity` (wire `intensities`):

  map\<string,float\>; in **intensity (instrument units)**; never `NA`.
  Sample label -\> protein intensity. A value is null when median polish
  could not resolve one (NaN: a degenerate peptide matrix); 0 means not
  measured in that sample. A sample of several runs (fractions,
  technical replicates) reports their sum, as FlashLFQ's own output
  does.

## Fields of `samples`

- `label`:

  string; never `NA`. The run name when no design groups runs (every
  condition blank or Default, and one fraction); otherwise
  '\<condition\>\_\<biorep + 1\>', e.g. control_1.

- `condition`:

  string; never `NA`. The sample's condition; "" with no design.

- `biological_replicate`:

  int; never `NA`. 0-based (the label adds 1).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  peptides missing or blank, or not found; the file empty; a required
  column (Sequence, Base Sequence, Protein Groups) missing; no
  Intensity\_ columns; a sequence on two rows; an intensity cell that is
  not a number; a Detection Type\_ cell that is not a FlashLFQ
  DetectionType name

- `mzlib_usage_error` (usage):

  a design line with no run name, a run named twice, or a
  replicate/fraction that is not a non-negative integer; design runs
  with no Intensity\_ column; Intensity\_ runs the design does not
  mention (a design must name every run and only those)

- `mzlib_bridge_error` (correctness):

  FlashLFQ's reconstruction or median polish fails

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Detection types: where the table has no Detection Type\_ column for a
  run (or the cell is blank), a positive intensity is treated as MSMS
  and 0 as NotDetected. That is this verb's only inference; a present
  but unrecognised value is an error, not a guess.

- Protein intensity null means unquantifiable, 0 means not measured;
  neither is a measured zero.

- Runs are names, not files: nothing but the peptide table is opened,
  and peak-level information (MBR peaks, retention times) is not
  available here. Use quant flashlfq for those.

- With out, QuantifiedProteins.tsv labels its columns by the same rule
  as samples\[\].label at this pin (mzLib#1129 is in 23c2490e). Before
  \#1129 the file's rule was inverted for unfractionated data; the
  values always agreed.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.flashlfq.median_polish`

- Rust (mzLibRust): `mzlib::flashlfq::median_polish_with` with
  `MedianPolishOptions`

- R (mzLibR): `flashlfq_median_polish`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Not yet in mzLibRust or mzLibR. Intended names mirror each binding's
  flashlfq module: Rust mzlib::flashlfq::median_polish (defaults) and
  median_polish_with + MedianPolishOptions {design, use_shared_peptides,
  output_directory, timeout} beside quantify_with + QuantifyOptions (the
  DesignEntry builder in the example is a sketch); R
  flashlfq_median_polish(peptides, design = NULL, use_shared_peptides =
  FALSE, output_directory = NULL, timeout = NULL) beside
  flashlfq_quantify, design as a data frame. The checker cannot verify
  them.

- The bridge's own comment in MedianPolishToWire and pyMzLib's
  median_polish docstring still say the TSV labels disagree until
  mzLib#1129 lands; it is in the pin, so both are stale. The agreement
  stated in the caveat is from reading ProteinGroup.cs at 23c2490e, not
  from a recorded `--out` run.

- median_polish_small.json is hand-shaped ('/abs/...'), not a live
  recording.

- The design's techrep and fraction are parsed but only fraction feeds
  the label rule; how FlashLFQ's median polish uses techrep is mzLib's
  and was not traced.

- The wire params peptides/shared-peptides/out project as
  peptides/use_shared_peptides/output_directory in pyMzLib; the spec
  lint will need PYTHON_DEVIATIONS entries (stdin -\> design,
  shared-peptides -\> use_shared_peptides, out -\> output_directory).

## References

- [doi:10.1021/acs.jproteome.7b00608](https://doi.org/10.1021/acs.jproteome.7b00608):
  FlashLFQ, whose protein quantification this verb reruns (Millikin et
  al., J Proteome Res 2018)

## See also

[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)

## Examples

``` r

polished <- flashlfq_median_polish("QuantifiedPeptides.tsv",
  design = data.frame(file_name = c("run_3", "run_4"), condition = c("control", "treated"),
                      biological_replicate = 0))
polished
#> <mzlibr_median_polish> QuantifiedPeptides.tsv
#>   4 peptides, 3 protein groups, 2 samples: control_1, treated_1
#>   proteins: 2 NA (could not be resolved), 0 zero (not measured)
polished$samples
#>       label condition biological_replicate
#> 1 control_1   control                    0
#> 2 treated_1   treated                    0
polished$proteins
#>   protein_group gene_name     organism    sample intensity
#> 1            P1     GENE1 Homo sapiens control_1    3005.6
#> 2            P1     GENE1 Homo sapiens treated_1    6011.3
#> 3            P2     GENE2 Homo sapiens control_1    2262.7
#> 4            P2     GENE2 Homo sapiens treated_1    2368.5
#> 5            P3     GENE3 Homo sapiens control_1        NA
#> 6            P3     GENE3 Homo sapiens treated_1        NA
```
