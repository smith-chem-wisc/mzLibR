# Read a label-free experimental design out of an SDRF, or every reason it cannot be read

Read a label-free experimental design out of an SDRF, in MetaMorpheus's
terms, or every reason it was refused.

Calls mzLib's `SdrfLabelFreeDesign.Read`. It runs every check
MetaMorpheus's own design validator runs and **refuses rather than
repairs**: a factor written `not available`, two files claiming the same
replicate and fraction, a document with several factor columns and no
declaration - each is a refusal, and all of them are reported at once.
It relabels in only two ways, both recorded in `notes`: study-wide
biological replicate numbers are ranked within each condition
(`22 -> 1`), and rows for files the search does not read are dropped.

## Usage

``` r
sdrf_design(path, condition_columns = NULL, searched_files = NULL, out = NULL,
  timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- condition_columns:

  The `factor value[...]` columns that make up the condition, by exact
  (case-sensitive) name, in the order to join them with `_`. `NULL` uses
  the document's only factor column; a document with several is then
  refused, because joining all of them could split a condition on a
  nuisance factor.

- searched_files:

  The files the search will read, as paths or bare names. Each must be
  named **exactly** (case and extension) by one SDRF row; rows for other
  files are dropped and reported in `notes`, and the paths you give
  become `full_path`. `NULL` takes the SDRF's own file names.

- out:

  Also write MetaMorpheus's `ExperimentalDesign.tsv` (1-based) here.
  Must end in `.tsv`. Written only when the design is valid.
  MetaMorpheus finds the file only when it is named
  `ExperimentalDesign.tsv` and sits beside the spectra.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**A refusal is an answer, not an error.** When `is_valid` is `FALSE`,
`refusals` lists every reason and `records` has no rows, so a design
MetaMorpheus would reject cannot be used by mistake. mzLib refuses
because MetaMorpheus skips quantification with only a warning when its
design file is invalid.

## Value

An `mzlibr_sdrf_design`: `path`; `is_valid`; `file_key_column`, the
column runs were keyed on (`NA` when the SDRF names no file, itself a
refusal); `condition_columns` used, in join order;
`condition_columns_declared` and `searched_files_given`; `file_count`
files (runs) in the design, 0 when refused; `refusals` and `notes` in
mzLib's words; `report`, mzLib's account of all of it; `written`
(`path`, `file_count`), or `NULL`; `caveats`; and `records`, one row per
run in SDRF row order: `full_path`, `file_name` (the run key,
`Intensity_<file_name>`), `condition`, and `biological_replicate`,
`technical_replicate` and `fraction`, **0-based**, as the quant
functions take them. `ExperimentalDesign.tsv` writes the same numbers
plus one.

## From SDRF to FlashLFQ

[`sdrf_design_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md)
gives the design as
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
takes `spectra`, and
[`sdrf_design_run_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_run_design.md)
as
[`flashlfq_median_polish`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
takes `design`, so the SDRF drives the quantification with no
hand-written design. Both refuse a refused design.

Label-free only: an isobaric (TMT, iTRAQ) SDRF needs a channel design
this does not produce. See
[`isobaric_kits`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
for the channels themselves.

## Wraps

Wire verb `sdrf design`. Generated from the bridge's verb spec
`sdrf.design.yaml` (bridge commit `e76157831b15`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfLabelFreeDesign.Read`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Readers/Sdrf/SdrfLabelFreeDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfLabelFreeDesign.cs` at mzLib `0a808fec`

- [`SdrfLabelFreeDesign.WriteExperimentalDesignTsv`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Readers/Sdrf/SdrfLabelFreeDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfLabelFreeDesign.cs` at mzLib `0a808fec`

- [`SdrfLabelFreeDesignOptions`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Readers/Sdrf/SdrfLabelFreeDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfLabelFreeDesign.cs` at mzLib `0a808fec`

## Parameters: units, ranges and defaults

- `path`:

  path; required. One .sdrf.tsv file. A file that does not exist is a
  usage error.

- `condition_columns` (wire `--condition-columns`):

  string; default absent. The factor value\[...\] columns the condition
  is built from, by exact (case-sensitive) name, separated by TAB
  characters, in join order; the condition is their values joined with
  "\_". Absent: mzLib uses the document's only factor column, and
  refuses a document with several. Given but blank, or with a blank
  name, is a usage error.

- `searched_files` (wire `--searched-files-stdin`):

  flag; default absent. stdin holds the files the search will read, one
  per line, as paths or bare names. Each must be named EXACTLY (case and
  extension) by one row; rows for other files are dropped and reported
  in notes; the searched paths replace the SDRF's names in full_path. An
  empty list is mzLib's refusal, not a usage error.

- `out`:

  path; default absent; range `.tsv only, any case (PYB-1)`. Also write
  MetaMorpheus's ExperimentalDesign.tsv (1-based) here, through mzLib's
  writer. Must end in .tsv, case-insensitive (PYB-1), checked before the
  SDRF is read. Written only when is_valid; a refused design writes
  nothing and written is null.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The SDRF path as given.

- `is_valid`:

  bool; never `NA`. No refusals. False is an answer, not an error: the
  table is then empty.

- `file_key_column`:

  string; `NA` when the SDRF has neither comment\[searched data file\]
  nor comment\[data file\], so no row names a file (also a refusal). The
  column files were keyed on: comment\[searched data file\] when
  present, else comment\[data file\] (mapping row MAP-12).

- `condition_columns`:

  string\[\]; never `NA`. The factor columns the condition was built
  from, in join order; empty when none could be chosen.

- `condition_columns_declared`:

  bool; never `NA`. Whether `--condition-columns` was given.

- `searched_files_given`:

  bool; never `NA`. Whether `--searched-files-stdin` was given.

- `file_count`:

  int; in **files**; never `NA`. Rows in the table; 0 when refused.

- `refusals`:

  string\[\]; never `NA`. Every reason the design was refused, all at
  once, in mzLib's words. Empty when valid.

- `notes`:

  string\[\]; never `NA`. Every relabelling on the way to a valid
  design: biological replicates ranked within a condition, with the old
  -\> new mapping (MAP-33), and rows dropped because the search does not
  read their file.

- `report`:

  string; never `NA`. mzLib's human-readable account
  (SdrfLabelFreeDesign.Report): counts, key column, condition columns,
  notes, refusals.

- `column_names`:

  string\[\]; never `NA`. The table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per file, in SDRF row
  order.

- `written`:

  object; `NA` when `--out` was not given, or the design was refused and
  nothing was written. {path, file_count} of the ExperimentalDesign.tsv
  written.

- `caveats`:

  string\[\]; never `NA`. Six static caveats.

## Columns of `records`

- `full_path`:

  string; never `NA`. The file as the SDRF's key column names it, or the
  matching searched path when `--searched-files-stdin` was given.

- `file_name`:

  string; never `NA`. full_path's file name without its extension
  (SpectraFileInfo.FilenameWithoutExtension): the run key the quant
  verbs use.

- `condition`:

  string; never `NA`. The condition columns' values joined with \_.

- `biological_replicate`:

  int; never `NA`. 0-based, ranked within the condition
  (ExperimentalDesign.tsv's Biorep minus 1).

- `technical_replicate`:

  int; never `NA`. 0-based (Techrep minus 1).

- `fraction`:

  int; never `NA`. 0-based (Fraction minus 1).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  path missing or the file does not exist

- `mzlib_usage_error` (usage):

  `--out` without a value, or not ending in .tsv (checked before the
  SDRF is read)

- `mzlib_usage_error` (usage):

  `--condition-columns` without a value, or with a blank name

- `mzlib_bridge_error` (correctness):

  mzLib cannot parse the SDRF at all, or the `--out` directory cannot be
  written

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Label-free only. An isobaric (TMT, iTRAQ) SDRF needs a channel design
  this verb does not produce.

- is_valid = false is an answer: refusals lists every reason at once and
  the table is empty. mzLib refuses rather than repairs because
  MetaMorpheus skips quantification with only a warning when its design
  file is invalid.

- The table is 0-based, as SpectraFileInfo and the quant verbs take it.
  ExperimentalDesign.tsv is 1-based; mzLib converts when it writes.

- notes is the only record of a renumbering; read it.

- A factor holding not available or not applicable is refused: it would
  pair samples nobody said were alike.

- MetaMorpheus finds the design only as a file named
  ExperimentalDesign.tsv beside the spectra; `--out` writes any .tsv
  path.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.design`

- Rust (mzLibRust): `mzlib::sdrf::design_with` with `DesignOptions`

- R (mzLibR): `sdrf_design`

## Since

Wire protocol 1; pyMzLib 0.3.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- No `--paths-stdin` form. The condition columns and searched files are
  per-deposit facts; a batch would need either one shared declaration
  (only sensible with none, i.e. each document's only factor column) or
  a per-document list. Decide before adding sdrf design_many.

- Not yet ported: mzLibRust and mzLibR names above are proposals until
  their ports land (parity order).

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics, the format read (Dai et al., Nat Commun 12, 5854,
  2021)

## See also

[`sdrf_design_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md),
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md),
[`sdrf_samples`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)

## Examples

``` r

both <- c("factor value[genotype]", "factor value[treatment]")
d <- sdrf_design("PXD067622.sdrf.tsv", condition_columns = both)
d
#> <mzlibr_sdrf_design> PXD067622.sdrf.tsv: 24 runs in 8 condition(s)
#>   files keyed on: comment[data file]
#>   condition from: factor value[genotype] + factor value[treatment]
head(d$records[, c("file_name", "condition", "biological_replicate")], 4)
#>                               file_name                               condition
#> 1 20240830_HF_LC3_MAA_RK_12032_CA_DMSO4         SPRTN-TurboID CA_DMSO (vehicle)
#> 2 20240830_HF_LC3_MAA_RK_12032_CA_DMSO5         SPRTN-TurboID CA_DMSO (vehicle)
#> 3 20240830_HF_LC3_MAA_RK_12032_CA_DMSO6         SPRTN-TurboID CA_DMSO (vehicle)
#> 4  20240830_HF_LC3_MAA_RK_12032_CA_FA22 SPRTN-TurboID CA_formaldehyde 1 mM, 1 h
#>   biological_replicate
#> 1                    0
#> 2                    1
#> 3                    2
#> 4                    0

# A refusal lists every reason at once - here one per run, as the treatment is "not available":
refused <- sdrf_design("PXD049018.sdrf.tsv", condition_columns = both)
refused$is_valid
#> [1] FALSE
length(refused$refusals)
#> [1] 20
writeLines(refused$refusals[[1]])
#> Line 2 (MSB67868ABand_01.raw): 'factor value[treatment]' is 'not available'. A condition cannot be built from an unknown factor; fill it in, or leave the column out of the declared condition columns to pool these rows.

# Study-wide replicate numbers are ranked within each condition, and the mapping is kept:
ranked <- sdrf_design("PXD067622_studywide.sdrf.tsv", condition_columns = both)
writeLines(ranked$notes[[2]])
#> Condition 'SPRTN-TurboID CA_formaldehyde 1 mM, 1 h': biological replicates renumbered 22 -> 1, 23 -> 2, 24 -> 3.
```
