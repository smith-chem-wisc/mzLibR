# Read PTM site occupancy from a MetaMorpheus protein-group table

Read the PTM site occupancy of a MetaMorpheus
AllQuantifiedProteinGroups.tsv as one row per group, sample group, basis
and modified site.

MetaMorpheus writes each protein group's modification-site occupancy
into two cells per sample group - by PSM count and by intensity - in a
compact text form. This reads them as one row per group, sample group,
basis and modified site.

## Usage

``` r
readers_read_occupancy(path, limit = NULL, offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to a MetaMorpheus `AllQuantifiedProteinGroups.tsv`.

- limit:

  Maximum protein groups to read - groups, not rows; a group with no
  modified sites gives no rows. `NULL`, the default, reads all of them.

- offset:

  Protein groups to skip - groups, not rows.

- out:

  Write the table to this path as tab-separated text and return only a
  summary.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**Read each basis for what it is exact in.** On `basis == "count"` rows
`numerator` and `denominator` are exact PSM counts and `fraction` is
rounded; on `basis == "intensity"` rows `fraction` is exact and the two
intensities are rounded.

## Value

An `mzlibr_occupancy_records`. `record_count` counts the protein groups
in the whole file; `returned_count` the groups read, starting `offset`
groups in; `row_count` the rows in `records`. `truncated_cell_count`
counts occupancy cells MetaMorpheus cut short, in cells.

`records` has one row per group, sample group, basis and site:
`fraction` from 0 to 1, and `numerator` and `denominator` in PSMs on
count rows or in intensity on intensity rows. `position` is the 1-based
residue; `cell_is_truncated` marks a row from a cut cell.
`rows_not_read` is `NA`: rows are not counted for this table.

## entity_index is not an accession index

An entity with no modified sites is dropped from the cell, so
`entity_index` cannot be zipped with the accessions in
`protein_group_name`.

## Wraps

Wire verb `readers read-occupancy`. Generated from the bridge's verb
spec `readers.read-occupancy.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`ModificationOccupancyCell.Parse`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Omics/BioPolymerGroup/ModificationOccupancyCell.cs)
  in `mzLib/Omics/BioPolymerGroup/ModificationOccupancyCell.cs` at mzLib
  `23c2490e`

- [`SampleGroupMeasurement.CountOccupancy`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs)
  in
  `mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A MetaMorpheus AllQuantifiedProteinGroups.tsv.

- `offset`:

  int; in **groups**; default `0`; range `>= 0`. Skip this many protein
  groups (a group with no modified sites gives no rows). A window, not a
  cursor: mzLib parses the whole file on every call, so paging re-reads
  it.

- `limit`:

  int; in **groups**; default absent; range `>= 0`. Return at most this
  many protein groups (a group with no modified sites gives no rows).
  Absent means no limit; there is deliberately no default cap.

- `out`:

  path; default absent. Write the table to this path as tab-separated
  text instead of returning columns. Must differ from path. Parent
  directories are created.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. Absolute path of the input.

- `file_type`:

  string; never `NA`. The mzLib SupportedFileType that was dispatched.

- `reader`:

  string; never `NA`. The mzLib reader class that parsed the file.

- `sample_labels`:

  string\[\]; never `NA`. The file's sample-group labels in header
  order.

- `truncated_cell_count`:

  int; in **cells**; never `NA`. Occupancy cells the writer cut short or
  replaced with 'Output too long for Excel', including those left with
  no complete site (which give no rows).

- `record_count`:

  int; in **groups**; never `NA`. Groups in the whole file, before the
  window.

- `rows_not_read`:

  int; in **rows**; `NA` when not counted for this table. Data rows that
  did not become records.

- `retention_time_unit`:

  string; `NA` when the table has no time column. Always null for this
  table.

- `caveats`:

  string\[\]; never `NA`. What this view cannot be trusted to mean for
  THIS file, each citing the mzLib source.

- `column_names`:

  string\[\]; never `NA`. Column order.

- `absent_fields`:

  string\[\]; never `NA`. Columns the view defines that this file's
  format has no column or source for; null in every row (BULK.md section
  4). Empty when every column has a source.

- `failed_fields`:

  string\[\]; never `NA`. 'field: ExceptionType' for each column whose
  read threw on some returned rows; those cells are null.

- `excluded_fields`:

  object\[\]; never `NA`. {field, type, reason, verb} for each field
  with no column shape; verb names the command that carries it, or null.

- `error`:

  object; `NA` when always, for a single path: a file that cannot be
  read fails the call instead. {kind, type, message}; non-null only in a
  files\[\] entry under on-error skip.

- `returned_count`:

  int; in **groups**; never `NA`. Groups returned in columns (the unit
  offset and limit count in); 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns: more than
  returned_count in a long table, where one record gives several rows; 0
  when written to out.

- `offset`:

  int; in **groups**; never `NA`. The offset applied, in groups.

- `truncated`:

  bool; never `NA`. True whenever records were left out by offset or
  limit.

- `records` (wire `columns`):

  table; `NA` when the table was written to out instead. Column name to
  per-row values, in column_names order.

- `output`:

  object; `NA` when out was not given. {path, format: "tsv", row_count}
  when out was given.

## Columns of `records`

- `protein_group_name`:

  string; never `NA`. The group, as in read-protein-groups.

- `sample_label`:

  string; never `NA`. The sample group whose cell this site came from.

- `basis`:

  string; never `NA`. count (CountOccupancy\_ cell, from PSM counts) or
  intensity (IntensityOccupancy\_ cell).

- `entity_index`:

  int; never `NA`. 0-based index among the cell's entities THAT HAVE
  SITES; not zippable with the group's accessions.

- `position`:

  int; never `NA`. Position as written: 0 is the N-terminus, residues
  count from 1, length + 1 is the C-terminus (OccupancySite.Position).

- `is_n_terminus`:

  bool; never `NA`. position == 0.

- `modification`:

  string; never `NA`. ModificationIdWithMotif, e.g. 'Oxidation on M';
  may contain commas and brackets.

- `fraction`:

  float; in **fraction (0 to 1)**; never `NA`. The fraction as printed:
  two decimals for count cells, four for intensity cells.

- `numerator`:

  float; in **PSMs (basis count) or intensity (basis intensity)**; never
  `NA`. Modified count, or modified intensity to four significant
  digits.

- `denominator`:

  float; in **PSMs (basis count) or intensity (basis intensity)**; never
  `NA`. Total count, or total intensity to four significant digits.

- `cell_is_truncated`:

  bool; never `NA`. The cell this site came from was cut short; its
  complete sites are still rows.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`readers_read_occupancy_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy_many.md).

- `mzlib_usage_error` (usage):

  path missing or blank; file or directory not found; mzLib does not
  recognise the file type (the message points at readers formats)

- `mzlib_usage_error` (usage):

  offset, limit or out given with no value; offset or limit \< 0;
  non-integer value; out resolves to the same file as path

- `mzlib_usage_error` (usage):

  the file is not a MetaMorpheus protein-group table

- `mzlib_bridge_error` (correctness):

  mzLib cannot parse the table (MzLibException). A single malformed cell
  does NOT fail the read: it is named in failed_fields and gives no rows

- `mzlib_usage_error` (usage):

  paths-stdin with path, offset or limit; paths-stdin with a value; no
  paths on stdin; a repeated path; threads 0 or below -1; on-error not
  fail or skip; on-error skip with a single path; out naming an input

- `mzlib_usage_error` (usage):

  with on-error fail, an input that is missing or not this view: the
  usage error of that input, prefixed 'Input \<i\> (\<path\>)'

- `mzlib_bridge_error` (correctness):

  with on-error fail, mzLib failing to parse an input: that input's own
  error type, message prefixed 'Input \<i\> (\<path\>)'

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Trust numerator/denominator on basis count rows and fraction on basis
  intensity rows: each cell rounds the other
  (ModificationOccupancyCell).

- entity_index cannot be zipped with protein_group_name's accessions: an
  entity with no sites is dropped from the cell.

- A malformed cell is refused by mzLib (FormatException), named in
  failed_fields as count_occupancy or intensity_occupancy, and gives no
  rows.

- A truncated cell keeps its complete sites with cell_is_truncated true;
  truncated_cell_count counts every cut cell, including those with no
  rows.

- offset and limit count GROUPS; a group with no modified sites gives no
  rows.

- Many files: `--paths-stdin` reads the list in ONE process, `--threads`
  at a time, into one long table with source_index and source_path
  first. The output does not depend on `--threads` (a tested property);
  the default of 1 is a memory choice, since every reader materialises
  the whole file.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. For many files call
[`readers_read_occupancy_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_occupancy`; many files:
  `pymzlib.readers.read_occupancy_many`

- Rust (mzLibRust): `mzlib::readers::read_occupancy_with` with
  `ReadOptions`; many files: `mzlib::readers::read_occupancy_many`

- R (mzLibR): `readers_read_occupancy`; many files:
  [`readers_read_occupancy_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy_many.md)

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- Peptide-level occupancy: ModificationOccupancyCell also parses peptide
  tables' cells, but QuantifiedPeptideFromTsv at the pin exposes none,
  so this verb reads protein groups only.

- numerator and denominator are floats for both bases; a count-basis
  value crosses as a whole number (1, not 1.0) and bindings must read it
  as a float.

## See also

[`readers_read_occupancy_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy_many.md),
[`readers_read_protein_groups`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md)

## Examples

``` r

sites <- readers_read_occupancy("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
sites
#> <mzlibr_occupancy_records> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv (MetaMorpheusQuantifiedProteinGroups)
#>   18 samples; 95 rows for 6 protein groups
#>   6 records in the file, 6 returned
#>   ! Trust numerator/denominator on basis 'count' rows and fraction on basis 'intensity' rows. MetaMorpheus prints a count cell's fraction to two decimals and an intensity cell's numerator and denominator to four significant digits (ModificationOccupancyCell). Exact intensities are in read-protein-groups.
#>   ! entity_index is the position among the cell's entities THAT HAVE SITES: an entity with none is dropped from the cell, so it cannot be zipped with the group's accessions by position (ModificationOccupancyCell).
#>   ! position counts 0 as the N-terminus and residues from 1; the C-terminus is length + 1 and cannot be recognised without the sequence, so only is_n_terminus is given (OccupancySite).
#>   ! A cell the writer cut short, or replaced with 'Output too long for Excel', keeps its complete sites with cell_is_truncated true; a cut cell with no complete site gives no rows. truncated_cell_count counts both, so a short table is never mistaken for a complete one.
#>   ! offset and limit count GROUPS; each group yields one row per sample group, basis and site, and a group with no modified sites yields none.
counted <- sites$records[sites$records$basis == "count", ]
head(counted[, c("protein_group_name", "position", "modification", "numerator", "denominator")])
#>    protein_group_name position                  modification numerator
#> 1              P68363      329              Deamidation on N         1
#> 3              P68363      329              Deamidation on N         1
#> 5              P68363      329              Deamidation on N         1
#> 7              P68363      329              Deamidation on N         2
#> 9              P68363      329              Deamidation on N         1
#> 11             P05141       52 N6,N6,N6-trimethyllysine on K         1
#>    denominator
#> 1            2
#> 3            2
#> 5            2
#> 7            4
#> 9            4
#> 11           1
```
