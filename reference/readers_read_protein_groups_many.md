# Read many protein-group tables into one long table, in one bridge call

Read a MetaMorpheus AllQuantifiedProteinGroups.tsv as one row per
protein group per sample group, with each sample group's spectral count
and intensity as columns.

The many-files form of
[`readers_read_protein_groups`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md):
the whole list goes to ONE bridge process, which reads `threads` files
at once and returns one table, rows grouped by file in input order
whatever `threads` is. Never loop the one-file function over files
instead.

## Usage

``` r
readers_read_protein_groups_many(paths, out = NULL, threads = 1, on_error = "fail",
  timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to MetaMorpheus
  `AllQuantifiedProteinGroups.tsv` files.

- out:

  Write the long table here as tab-separated text, one file at a time,
  and return only a summary.

- threads:

  Files read at once. `1`, the default, holds one whole file in memory
  at a time; `-1` uses one per core. The table does not depend on it.

- on_error:

  `"fail"`, the default, stops at the first file that cannot be read,
  with its error. `"skip"` records the failure in that file's row of
  `files` and reads the rest.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_read_batch`. `records` is the long table: `source_index`, the
1-based position of the row's file in `paths`, and `source_path`, then
the same columns as
[`readers_read_protein_groups`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md) -
`q_value` a fraction, `spectral_count` in PSMs, `intensity` in the
instrument's intensity units. `files` is a data.frame with one row per
input: its `path`, `file_type`, `reader`, `record_count` in groups,
`sample_labels`, `caveats`, `absent_fields`, `failed_fields`,
`excluded_fields` and, for an input that could not be read,
`error_kind`, `error_type` and `error_message`. `file_count`,
`read_count` and `failed_count` count files; `record_count` counts
groups over the files read and `returned_count` the groups returned;
`row_count` counts rows.

## Wraps

Wire verb `readers read-protein-groups` with `--paths-stdin`. Generated
from the bridge's verb spec `readers.read-protein-groups.yaml` (bridge
commit `5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these
facts, and all three bindings render the same ones.

- [`ProteinGroupFromTsvFile.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/ResultFiles/ProteinGroupFromTsvFile.cs)
  in
  `mzLib/Readers/InternalResults/ResultFiles/ProteinGroupFromTsvFile.cs`
  at mzLib `23c2490e`

- [`ProteinGroupFromTsv.SampleGroups`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs)
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
  block (path, file_type, reader, sample_labels, record_count,
  rows_not_read, retention_time_unit, caveats, column_names,
  absent_fields, failed_fields, excluded_fields, error). An unread input
  has error set and every fact it could not establish null or empty.

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

  string; never `NA`. The members' accessions '\|'-joined in accession
  order (ProteinGroupFromTsv.ProteinGroupName); the join key to
  read-records.

- `gene`:

  string; `NA` when the file has no Gene column (then also in
  absent_fields) or the cell is blank. Gene name(s) as written.

- `organism`:

  string; `NA` when the file has no Organism column (then also in
  absent_fields) or the cell is blank. Organism as written.

- `decoy_contaminant_target`:

  string; never `NA`. T, D, C, or entrapment ET/ED, as written.

- `is_decoy`:

  bool; never `NA`. decoy_contaminant_target contains D (mzLib's
  IsDecoy).

- `is_contaminant`:

  bool; never `NA`. decoy_contaminant_target is C.

- `is_entrapment`:

  bool; never `NA`. decoy_contaminant_target starts with E.

- `q_value`:

  float; in **fraction (0 to 1)**; never `NA`. The group's q-value; the
  table is unfiltered and this is what filtering reads.

- `sample_label`:

  string; `NA` when the file has no per-sample columns: the group is one
  row and this is in absent_fields. The sample group this row measures.

- `spectral_count`:

  int; in **PSMs**; `NA` when the cell is blank, or the file has no
  SpectralCount\_ columns (then in absent_fields); an isobaric
  intensity-only label has none. PSMs from this sample group assigned to
  the group (SpectralCount\_).

- `intensity`:

  float; in **intensity (instrument units)**; `NA` when the cell was
  blank: not quantified in this sample group, which is not zero. The
  group's intensity in this sample group (Intensity\_).

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

  the file is not a MetaMorpheus protein-group table (the message names
  what it is and points at read-records)

- `mzlib_bridge_error` (correctness):

  mzLib cannot parse the table (MzLibException)

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

- The table is UNFILTERED as MetaMorpheus writes it: decoys,
  contaminants and groups above 1% FDR are all rows
  (ProteinGroupFromTsv.cs:15). Filter on q_value and
  decoy_contaminant_target before counting.

- Long format: one row per group per sample group. offset and limit
  count GROUPS, and returned_count counts rows.

- sample_label is the header label verbatim; condition, replicate and
  channel cannot be recovered from it (ProteinGroupFromTsv.cs:103).

- A blank intensity is null, not zero (ProteinGroupFromTsv.cs:114). An
  isobaric file has one counting label and one intensity label per
  channel, so a row can carry only spectral_count or only intensity.

- No member is a leading protein: members are sorted by accession
  (ProteinGroupFromTsv.cs:17).

- Occupancy cells are excluded_fields with verb readers read-occupancy;
  the group's other fields are in read-records, joined on
  protein_group_name.

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

- Python (pyMzLib): `pymzlib.readers.read_protein_groups`; many files:
  `pymzlib.readers.read_protein_groups_many`

- Rust (mzLibRust): `mzlib::readers::read_protein_groups_with` with
  `ReadOptions`; many files: `mzlib::readers::read_protein_groups_many`

- R (mzLibR): `readers_read_protein_groups`; many files:
  `readers_read_protein_groups_many`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- The group-level columns carried beside each sample (identity, q-value,
  decoy label) are a choice made here; coverage, masses and counts stay
  in read-records. Revisit if every binding ends up joining for the same
  field.

- Only the MetaMorpheus 1.1.11 fixture exists (label-free); no isobaric
  or 1.2+ (BioPolymer ...) header table has been read through this verb
  yet.

## See also

[`readers_read_protein_groups`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md)

## Examples

``` r

if (FALSE) { # \dontrun{
# Your own files: no recording of a many-file read of this table exists to replay.
batch <- readers_read_protein_groups_many(
  c("search_1/AllQuantifiedProteinGroups.tsv", "search_2/AllQuantifiedProteinGroups.tsv"),
  threads = 2, on_error = "skip")
batch$files[, c("path", "record_count", "error_message")]
} # }
```
