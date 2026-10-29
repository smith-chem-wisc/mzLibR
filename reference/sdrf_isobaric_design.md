# Read an isobaric (TMT, TMTpro, iTRAQ, DiLeu) design out of an SDRF, or every reason it cannot be read

Read an isobaric (TMT, TMTpro, iTRAQ, DiLeu) experimental design out of
an SDRF, one row per file and annotated channel, in MetaMorpheus's
terms, or every reason it was refused.

Answers **can this SDRF drive an isobaric quantification, and which
sample is in which channel of which file?** Calls mzLib's
`SdrfIsobaricDesign.Read`. It runs MetaMorpheus's own checks on a TMT
design - one plex, fraction and technical replicate per file; a channel
described once per plex; no sample twice at one fraction in a plex - and
**refuses rather than repairs**, reporting every failure at once.

## Usage

``` r
sdrf_isobaric_design(path, kit = NULL, plex_column = NULL,
  plex_file_name_pattern = NULL, single_plex = NULL, condition_columns = NULL,
  searched_files = NULL, out = NULL, timeout = 60)
```

## Arguments

- path:

  Path to a `.sdrf.tsv` file.

- kit:

  The kit the search used, named as
  [`isobaric_kits`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
  takes it (`"TMT10"`, `"tmt18"`, `"TMT6-plex"`, `"TMT11 on K"`). `NULL`
  sends none, and mzLib then refuses: the kit is never guessed from the
  labels.

- plex_column:

  The column each row's plex is read from, by its exact name.

- plex_file_name_pattern:

  A .NET regular expression matched, ignoring case, against each row's
  file name; the plex is its group named `plex`, else its first group,
  e.g. `"TMT_?pool(?<plex>\\d+)"` (one backslash in the pattern, two in
  R source).

- single_plex:

  The name of the one plex every file belongs to.

- condition_columns:

  The `factor value[...]` columns that make up the condition, exactly as
  for
  [`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md).
  `NULL` uses the document's only factor column.

- searched_files:

  The files the search will read, exactly as for
  [`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md):
  each must be named **exactly** by some row. A search that read
  `X-calib.mzML` for an SDRF naming `X.raw` matches nothing;
  [`sdrf_restrict`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)
  the SDRF first.

- out:

  A **folder**, created when missing. mzLib writes MetaMorpheus's
  `TmtDesign.txt` (1-based) inside it, replacing one already there, but
  only when the design is valid. MetaMorpheus reads it only from the
  folder of the first spectra file.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**It guesses nothing:** not the kit, because one kit's labels can be a
subset of another's, and not the plex, which you declare with exactly
one of `plex_column`, `plex_file_name_pattern` and `single_plex`. Giving
none, or several, is mzLib's refusal.

Channel labels are read by mzLib's `SdrfIsobaricLabels`, the rule
[`isobaric_kits`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)'s
`sdrf_label` writes: for TMT10, `TMT131` (PRIDE's spelling of the last
channel) names `131N` (mzLib \#1467); for TMT11, where `131N` and `131C`
are both channels, it is still refused.

**The numbers are 1-based, as the SDRF writes them.**
[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)'s
are 0-based. That split is mzLib's, and mzLibR converts neither.

## Value

An `mzlibr_sdrf_isobaric_design`. `is_valid` is `TRUE` when there are no
refusals. **A refusal is an answer, not an error**: `refusals` then
lists every reason at once, in mzLib's words, and `records` has no rows.

`kit` is the kit name exactly as you gave it, and `kit_type` mzLib's
name for the kit it resolved to (`"TMT10"`, `"TMT18"`, `"iTRAQ4"`, ...),
as
[`isobaric_kits`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
names it; both `NA` when no kit was given. `kit_channels` lists every
channel of the kit in reporter m/z order (a channel here and absent from
`records` is one the SDRF does not annotate), `NULL` without a kit.
`plex_column`, `plex_file_name_pattern` and `single_plex` echo the plex
source you gave, `NA` for the others; `plex_source` is where mzLib read
the plex from, as its report prints it, `"(not declared)"` when no
single source could be used. `plex_ids` is a named vector, plex name to
its number, 1 to N in case-insensitive name order (the order
MetaMorpheus numbers plexes in), empty when refused; `plex_count` counts
the plexes, 0 when refused. `file_key_column` is the column files were
keyed on - `comment[searched data file]` when the SDRF has it, else
`comment[data file]` - or `NA` when it has neither. `condition_columns`
lists the columns the condition was built from;
`condition_columns_declared` and `searched_files_given` say whether you
passed them. `file_count` counts the files in the design and
`record_count` the rows of `records`, one per file and annotated
channel; both 0 when refused. `notes` says what was dropped or filled on
the way; `report` is mzLib's account of all of it, to print. `written`
is a list (`path` of the `TmtDesign.txt`, `file_count`, `record_count`),
or `NULL` when `out` was not given or the design was refused. `caveats`
says what the design does and does not cover.

`records` is a data.frame, files in SDRF row order, then each file's
channels in reporter m/z order: `full_path` (the searched path you
passed for the file, else the file name the SDRF gives); `plex`, as the
declared source gives it; `plex_id`, the plex's number in `plex_ids`;
`sample_name`, the SDRF `source name` of the sample in the channel;
`channel`, as the kit spells it (`"127N"`); `condition`, the condition
columns' values joined with `_`; and `biological_replicate` (as written,
never ranked), `fraction` (1 when the SDRF has no fraction column) and
`technical_replicate` (1 when it has no technical-replicate column),
**all 1-based**.

## Restrict first, then read the design

A search reads `X-calib.mzML` where the SDRF names the acquired `X.raw`,
so `searched_files` of the search's own names match nothing.
[`sdrf_restrict`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)
joins the two and writes the searched name into
`comment[searched data file]`, and this function keys on that column
when a document has it.

## Wraps

Wire verb `sdrf isobaric-design`. Generated from the bridge's verb spec
`sdrf.isobaric-design.yaml` (bridge commit `45833ca804d6`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfIsobaricDesign.Read`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs` at mzLib `d10017ce`

- [`SdrfIsobaricDesign.WriteTmtDesign`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs` at mzLib `d10017ce`

- [`SdrfIsobaricDesign.PlexIds`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs` at mzLib `d10017ce`

- [`SdrfIsobaricDesignOptions`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs` at mzLib `d10017ce`

- [`SdrfIsobaricFile`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs` at mzLib `d10017ce`

- [`SdrfIsobaricChannel`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs)
  in `mzLib/Readers/Sdrf/SdrfIsobaricDesign.cs` at mzLib `d10017ce`

- [`IsobaricMassTag.TryGetTagType`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/Modifications/IsobaricMassTag.cs)
  in `mzLib/Omics/Modifications/IsobaricMassTag.cs` at mzLib `d10017ce`

- [`IsobaricMassTag.TryGetIsobaricMassTag`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/Modifications/IsobaricMassTag.cs)
  in `mzLib/Omics/Modifications/IsobaricMassTag.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `path`:

  path; required. One .sdrf.tsv file. A file that does not exist is a
  usage error.

- `kit`:

  string; default absent. The kit the search used, resolved exactly as
  isobaric kits `--kit` resolves it: mzLib's
  IsobaricMassTag.TryGetTagType, whole name, case-insensitive,
  optionally followed by " on \<motif\>" (TMT10, tmt10, TMT6-plex, TMT11
  on K), never a substring. It is SdrfIsobaricDesignOptions.Tag, which
  supplies the channel list and their reporter m/z order. A blank or
  unknown name is a usage error naming the known kits. Absent: mzLib is
  given no tag and refuses ("No isobaric tag was given"), because the
  kit is never guessed from the labels.

- `plex_column` (wire `--plex-column`):

  string; default absent. SdrfIsobaricDesignOptions.PlexColumn: the
  column each row's plex is read from, by its exact name. A column not
  in the SDRF, or a cell that is empty or not available, is mzLib's
  refusal. Given with no value is a usage error. mzLib needs exactly one
  of plex-column, plex-file-name-pattern and single-plex; none, or
  several, is mzLib's refusal, not a usage error.

- `plex_file_name_pattern` (wire `--plex-file-name-pattern`):

  string; default absent. SdrfIsobaricDesignOptions.PlexFileNamePattern:
  a .NET regular expression matched, ignoring case, against each row's
  file name (the key column's, comment\[searched data file\] else
  comment\[data file\]). The plex is its group named plex, else its
  first group. A pattern that does not parse, has no group, or does not
  match a file name is mzLib's refusal. Given with no value is a usage
  error.

- `single_plex` (wire `--single-plex`):

  string; default absent. SdrfIsobaricDesignOptions.SinglePlex: the name
  of the one plex every file belongs to. Given with no value is a usage
  error.

- `condition_columns` (wire `--condition-columns`):

  string; default absent. SdrfIsobaricDesignOptions.ConditionColumns,
  exactly as sdrf design takes it: the factor value\[...\] columns the
  condition is built from, by exact name, separated by TAB characters,
  in join order; the condition is their values joined with "\_". Absent:
  mzLib uses the document's only factor column, and refuses a document
  with several. Given but blank, or with a blank name, is a usage error.

- `searched_files` (wire `--searched-files-stdin`):

  flag; default absent. SdrfIsobaricDesignOptions.SearchedFiles, exactly
  as sdrf design takes it: stdin holds the files the search will read,
  one per line, as paths or bare names. Each must be named EXACTLY (case
  and extension) by some row; rows for other files are dropped and
  reported in notes; the searched paths become full_path. A search that
  read X-calib.mzML for an SDRF naming X.raw matches no row: run sdrf
  restrict first, which writes the searched name into comment\[searched
  data file\], the column this reader keys on first.

- `out`:

  path; default absent; range `a folder; an existing file is refused`. A
  FOLDER, created when missing (bridge thread 010, PYB-5). mzLib's
  SdrfIsobaricDesign.WriteTmtDesign writes MetaMorpheus's TmtDesign.txt
  inside it, replacing one already there: 1-based, one row per file and
  annotated channel, header SdrfIsobaricDesign.TmtDesignHeader. The file
  name is MetaMorpheus's
  (GlobalVariables.TmtExperimentalDesignFileName); mzLib has no constant
  for it at the pin, so the bridge holds the one copy. Written only when
  is_valid; a refused design writes nothing and creates no folder, and
  written is null. A path that is an existing file, or no value, is a
  usage error, checked before the SDRF is read.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. The SDRF path as given.

- `is_valid`:

  bool; never `NA`. No refusals. False is an answer, not an error: the
  table is then empty.

- `kit`:

  string; `NA` when `--kit` was not given (mzLib then refuses). The
  `--kit` name exactly as given.

- `kit_type`:

  string; `NA` when `--kit` was not given. The IsobaricMassTagType
  member the kit resolved to (TMT6, TMT10, TMT11, TMT16, TMT18, iTRAQ4,
  iTRAQ8, diLeu4, diLeu12 at d10017ce), as isobaric kits names it.

- `kit_channels`:

  string\[\]; `NA` when `--kit` was not given. Every channel label of
  the kit, in reporter m/z order (IsobaricMassTag.ChannelLabels). A
  channel here and absent from the table is one the SDRF does not
  annotate.

- `plex_column`:

  string; `NA` when `--plex-column` was not given. The plex column as
  given.

- `plex_file_name_pattern`:

  string; `NA` when `--plex-file-name-pattern` was not given. The
  pattern as given.

- `single_plex`:

  string; `NA` when `--single-plex` was not given. The single plex as
  given.

- `plex_source`:

  string; never `NA`. SdrfIsobaricDesign.PlexSource: where the plex was
  read from, as mzLib's report prints it; (not declared) when no single
  source could be used.

- `plex_ids`:

  object; never `NA`. SdrfIsobaricDesign.PlexIds: plex name to its
  IsobaricQuantSampleInfo.PlexId, 1..N in case-insensitive name order
  (MetaMorpheus's ToMzLibDesign rule), listed in id order. Empty when
  refused.

- `plex_count`:

  int; in **plexes**; never `NA`. Entries in plex_ids; 0 when refused.

- `file_key_column`:

  string; `NA` when the SDRF has neither comment\[searched data file\]
  nor comment\[data file\], so no row names a file (also a refusal). The
  column files were keyed on: comment\[searched data file\] when
  present, else comment\[data file\] (MAP-12).

- `condition_columns`:

  string\[\]; never `NA`. The factor columns the condition was built
  from, in join order; empty when none could be chosen.

- `condition_columns_declared`:

  bool; never `NA`. Whether `--condition-columns` was given.

- `searched_files_given`:

  bool; never `NA`. Whether `--searched-files-stdin` was given.

- `file_count`:

  int; in **files**; never `NA`. Files in the design
  (SdrfIsobaricDesign.Files); 0 when refused.

- `record_count`:

  int; in **rows**; never `NA`. Rows in the table: one per file and
  annotated channel; 0 when refused.

- `refusals`:

  string\[\]; never `NA`. Every reason the design was refused, all at
  once, in mzLib's words. Empty when valid.

- `notes`:

  string\[\]; never `NA`. What was dropped or filled on the way to a
  valid design: rows dropped because the search does not read their
  file, and a file with no row of its own for one of its plex's
  channels.

- `report`:

  string; never `NA`. mzLib's human-readable account
  (SdrfIsobaricDesign.Report): counts and kit, key column, plex source,
  condition columns, sample type, notes, refusals.

- `column_names`:

  string\[\]; never `NA`. The table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per file and annotated
  channel: files in SDRF row order, then each file's channels in
  reporter m/z order.

- `written`:

  object; `NA` when `--out` was not given, or the design was refused and
  nothing was written. {path, file_count, record_count} of the
  TmtDesign.txt written; path is the `--out` folder joined with
  TmtDesign.txt.

- `caveats`:

  string\[\]; never `NA`. Eight static caveats.

## Columns of `records`

- `full_path`:

  string; never `NA`. SdrfIsobaricFile.FilePath: the matching searched
  path when `--searched-files-stdin` was given, else the file name the
  key column gives.

- `plex`:

  string; never `NA`. SdrfIsobaricFile.Plex: the plex as the declared
  source gives it (a column cell, a pattern's group, or the single
  plex).

- `plex_id`:

  int; never `NA`. The plex's entry in plex_ids, 1-based: the PlexId
  mzLib's quantification design gives it.

- `sample_name`:

  string; never `NA`. SdrfIsobaricChannel.SampleName: the SDRF source
  name.

- `channel`:

  string; never `NA`. SdrfIsobaricChannel.Label: the channel as the kit
  spells it (127N), whatever form the SDRF's comment\[label\] used
  (TMT127N, NT=TMT127N;AC=...).

- `condition`:

  string; never `NA`. SdrfIsobaricChannel.Condition: the condition
  columns' values joined with \_.

- `biological_replicate`:

  int; never `NA`. SdrfIsobaricChannel.BiologicalReplicate: 1-based, as
  the SDRF writes it, never ranked or shifted.

- `fraction`:

  int; never `NA`. SdrfIsobaricFile.Fraction: 1-based, as written; 1
  when the SDRF has no comment\[fraction identifier\] column (mzLib's
  SdrfDesignRules).

- `technical_replicate`:

  int; never `NA`. SdrfIsobaricFile.TechnicalReplicate: 1-based, as
  written; 1 when the SDRF has no comment\[technical replicate\] column.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  path missing or the file does not exist

- `mzlib_usage_error` (usage):

  `--kit` without a value, or a name mzLib cannot resolve (the same
  message isobaric kits `--kit` gives)

- `mzlib_usage_error` (usage):

  `--plex-column`, `--plex-file-name-pattern` or `--single-plex` without
  a value; `--condition-columns` without a value or with a blank name

- `mzlib_usage_error` (usage):

  `--out` without a value, or naming an existing file (checked before
  the SDRF is read)

- `mzlib_bridge_error` (correctness):

  mzLib cannot parse the SDRF at all, or the `--out` folder cannot be
  created or written

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- 1-BASED, as written. fraction, technical_replicate and
  biological_replicate are the SDRF's own numbers, as mzLib's
  SdrfIsobaricDesign keeps them and TmtDesign.txt holds them. sdrf
  design is 0-based (mzLib's SpectraFileInfo). The split is mzLib's, and
  the bridge converts neither table: never mix the two without adding or
  removing 1.

- is_valid = false is an answer: refusals lists every reason at once and
  the table is empty. A missing kit, no plex source, two plex sources
  and a bad pattern are refusals, not errors.

- Nothing is guessed. The kit is only the one named, because one kit's
  labels can be a subset of another's (TMT11 with 131C unannotated reads
  exactly like TMT10). The plex is only the source declared;
  comment\[sample preparation batch\] is often not a plex (in PXD061609
  it holds three fractionation arms of one plex).

- Rows are the annotated channels only. mzLib's quantification design
  (ToExperimentalDesign) carries every kit channel, an unannotated one
  with biological replicate 0 and no sample; that design is not on the
  wire.

- A plex is one channel-to-sample map: every file of a plex carries all
  of its channels, so fractions share their samples.

- Channel labels are read through mzLib's SdrfIsobaricLabels, one rule
  for readers and writers. For TMT10, TMT131 (PRIDE's spelling of the
  last channel) names 131N (since mzLib \#1467); TMT131 is still refused
  for TMT11, where 131N and 131C are both channels.

- Sample type is not read (mzLib MAP-08): every annotated channel is a
  study sample, a reference or pool channel included.

- `--out` names a folder. MetaMorpheus reads TmtDesign.txt only from the
  folder of the first spectra file, and refuses a folder that also holds
  ExperimentalDesign.tsv. mzLib's writer replaces an existing
  TmtDesign.txt; MetaMorpheus's own `--sdrf`Design never does.

- TMT16 is mzLib's TMTpro 16-plex. MetaMorpheus searches it with its
  TMT18 label; a TMT16 design is read as is by such a run, with 134C and
  135N empty.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.isobaric_design`

- Rust (mzLibRust): `mzlib::sdrf::isobaric_design_with` with
  `IsobaricDesignOptions`

- R (mzLibR): `sdrf_isobaric_design`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- PUSH-TO-MAINLAND: TmtDesign.txt is MetaMorpheus's name
  (GlobalVariables.TmtExperimentalDesignFileName), and mzLib publishes
  no constant for it at d10017ce, so the bridge holds the one copy
  (Sdrf.TmtDesignFileName). Ask mzLib for
  SdrfIsobaricDesign.TmtDesignFileName, or a WriteTmtDesign overload
  that takes the folder, and read the name from there.

- No `--paths-stdin` form, for the reason sdrf design gives: the kit,
  the plex source, the condition columns and the searched files are
  per-deposit facts.

- mzLib's quantification design (ToExperimentalDesign: every kit channel
  with its reporter m/z, an unannotated one with biological replicate 0)
  is not projected. Add it when a bridge verb quantifies reporter ions
  from it.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics, the format read (Dai et al., Nat Commun 12, 5854,
  2021)

## See also

[`sdrf_restrict`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md),
[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md),
[`isobaric_kits`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)

## Examples

``` r

# PXD008841 writes no plex column; the plex survives only in the file names, spelled two ways
# (TMTpool1, TMT_pool5), so one pattern reads both:
d <- sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10",
  plex_file_name_pattern = "TMT_?pool(?<plex>\\d+)")
d
#> <mzlibr_sdrf_isobaric_design> PXD008841.subset.sdrf.tsv: 4 files, 2 plex(es), kit TMT10, 40 channel rows (1-based numbers)
#>   plex from: file-name pattern 'TMT_?pool(?<plex>\d+)' (declared)
#>   files keyed on: comment[data file]
d$plex_ids
#> 1 5 
#> 1 2 
d$records[10, c("channel", "sample_name", "plex_id", "fraction", "biological_replicate")]
#>    channel sample_name plex_id fraction biological_replicate
#> 10    131N        pool       1        1                    1

# Declare no plex source and mzLib refuses rather than guess one:
refused <- sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10")
refused$is_valid
#> [1] FALSE
writeLines(strwrap(refused$refusals[[1]], 90))
#> No plex source was declared. Name the column that holds the plex, give a file-name
#> pattern whose group captures it, or name the single plex. The plex is never inferred
#> (MAP-09): 'comment[sample preparation batch]' is often not a plex.
```
