# Read the PTM site occupancy of many protein-group tables, in one bridge call

Read the PTM site occupancy of a MetaMorpheus
AllQuantifiedProteinGroups.tsv as one row per group, sample group, basis
and modified site.

The many-files form of
[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md):
one bridge process for the whole list, `threads` files at a time, rows
grouped by file in input order.

## Usage

``` r
readers_read_occupancy_many(paths, out = NULL, threads = 1, on_error = "fail",
  timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to MetaMorpheus
  `AllQuantifiedProteinGroups.tsv` files.

- out:

  Write the long table here as tab-separated text and return only a
  summary.

- threads:

  Files read at once. `1`, the default, or `-1` for one per core. The
  table does not depend on it.

- on_error:

  `"fail"`, the default, or `"skip"` to record an unreadable file in
  `files` and read the rest.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_read_batch`: `records` begins with `source_index` (1-based,
into `paths`) and `source_path`, then the columns of
[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md) -
`fraction` from 0 to 1, `numerator` and `denominator` in PSMs or
intensity by basis. `files` has one row per input, including its
`truncated_cell_count`. `file_count`, `read_count` and `failed_count`
count files; `record_count` counts groups over the files read and
`returned_count` the groups returned; `row_count` counts rows.

## Wraps

Wire verb `readers read-occupancy` with `--paths-stdin`. Generated from
the bridge's verb spec `readers.read-occupancy.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`ModificationOccupancyCell.Parse`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Omics/BioPolymerGroup/ModificationOccupancyCell.cs)
  in `mzLib/Omics/BioPolymerGroup/ModificationOccupancyCell.cs` at mzLib
  `23c2490e`

- [`SampleGroupMeasurement.CountOccupancy`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs)
  in
  `mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `paths` (wire `--paths-stdin`):

  flag; default absent. Read MANY files: newline-delimited paths on
  stdin (blank lines ignored; a repeated path is a usage error).
  Exclusive with path; offset and limit are refused with it.

- `threads`:

  int; in **files**; default `1`; range `>= 1, or -1 for one per core`.
  Files read at once. The output is byte-identical at any value (files
  are emitted in input order); 1 is the default because each file in
  flight is held whole in memory.

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. fail: the first
  unreadable input stops the batch with its usage or correctness error,
  the message starting 'Input \<i\> (\<path\>)'. skip: the failure goes
  in that input's files\[i\].error and the rest are read. skip with a
  single path is a usage error.

- `out`:

  path; default absent. Write the long table here, one file at a time in
  input order, so memory holds at most threads files. A batch that stops
  on an error deletes its partial table. Must differ from every input.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **files**; never `NA`. Paths given.

- `read_count`:

  int; in **files**; never `NA`. Files read.

- `failed_count`:

  int; in **files**; never `NA`. Files not read; non-zero only under
  on-error skip.

- `record_count`:

  int; in **groups**; never `NA`. Groups summed over the files read
  (BULK.md section 2).

- `returned_count`:

  int; in **groups**; never `NA`. Groups returned in columns, summed
  over the files; 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns (more than
  returned_count for a long view); 0 when written to out.

- `on_error`:

  string; never `NA`. fail or skip, as requested. `--threads` is
  deliberately not echoed: the envelope is byte-identical at every
  thread count.

- `column_names`:

  string\[\]; never `NA`. source_index, source_path, then the
  single-file columns.

- `records` (wire `columns`):

  table; `NA` when the table was written to out instead. One long table:
  rows grouped by input in input order, each file's rows in its own
  order, whatever `--threads` is.

- `output`:

  object; `NA` when out was not given. {path, format: "tsv", row_count};
  written one file at a time.

- `files`:

  object\[\]; never `NA`. One per input, in input order: the per-file
  block (path, file_type, reader, sample_labels, truncated_cell_count,
  record_count, rows_not_read, retention_time_unit, caveats,
  column_names, absent_fields, failed_fields, excluded_fields, error).
  An unread input has error set and every fact it could not establish
  null or empty.

## Columns of `records`

The first two columns say which input each row came from -
`source_index` is 1-based in R, so it indexes `files` directly - and
rows are grouped by input, in input order, whatever `threads` is.

- `source_index`:

  int; never `NA`. 0-based position of the row's input in the stdin
  list.

- `source_path`:

  string; never `NA`. Absolute path of the row's input.

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
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_occupancy`; many files:
  `pymzlib.readers.read_occupancy_many`

- Rust (mzLibRust): `mzlib::readers::read_occupancy_with` with
  `ReadOptions`; many files: `mzlib::readers::read_occupancy_many`

- R (mzLibR): `readers_read_occupancy`; many files:
  `readers_read_occupancy_many`

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

[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)

## Examples

``` r

if (FALSE) { # \dontrun{
# Your own files: no recording of a many-file read of this table exists to replay.
batch <- readers_read_occupancy_many(
  c("search_1/AllQuantifiedProteinGroups.tsv", "search_2/AllQuantifiedProteinGroups.tsv"),
  threads = 2, on_error = "skip")
batch$files[, c("path", "record_count", "error_message")]
} # }
```
