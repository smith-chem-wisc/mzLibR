# Read many files, each into its format's own fields, in one bridge call

Read any file mzLib recognises into a table of that format's own record
fields, naming every field that could not become a column.

The many-files form of
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md).
Every file must have the same mzLib record type, because one table has
one column set: a list mixing record types is refused, naming the
groups, before any file is parsed. Split such a list by `file_type` from
[`readers_identify_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify_many.md)
first.

## Usage

``` r
readers_read_records_many(paths, out = NULL, threads = 1, on_error = "fail",
  timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to any files mzLib recognises.

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

An `mzlibr_read_batch`. `records` begins with `source_index` (1-based,
into `paths`) and `source_path`, then the format's own columns. `files`
has one row per input: `record_type`, `record_count` in records,
`skipped_count` in items for mzIdentML, `caveats`, `excluded_fields`,
`failed_fields` and the error columns. `file_count`, `read_count` and
`failed_count` count files; `record_count` and `returned_count` count
records; `row_count` counts the rows of `records`.

## Wraps

Wire verb `readers read-records` with `--paths-stdin`. Generated from
the bridge's verb spec `readers.read-records.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`FileReader.ReadResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`IResultFile.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/IResultFile.cs)
  in `mzLib/Readers/BaseClasses/IResultFile.cs` at mzLib `23c2490e`

- [`MzIdentMLResultFile.SkippedMatches`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/ExternalResults/ResultFiles/MzIdentMLResultFile.cs)
  in `mzLib/Readers/ExternalResults/ResultFiles/MzIdentMLResultFile.cs`
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

  int; in **records**; never `NA`. Records summed over the files read
  (BULK.md section 2).

- `returned_count`:

  int; in **records**; never `NA`. Records returned in columns, summed
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
  block (path, file_type, reader, record_type, views, skipped_count,
  skipped, record_count, rows_not_read, retention_time_unit, caveats,
  column_names, absent_fields, failed_fields, excluded_fields, error).
  An unread input has error set and every fact it could not establish
  null or empty.

## Columns of `records`

The columns are **per format**: this file's own mzLib record fields,
named in `column_names`. Every cell follows these rules:

- A cell is null when mzLib's value is null, when a double is NaN or
  infinite (JSON cannot carry it), or when the property threw (then
  named in failed_fields).

- -1 crosses unchanged: in a format's own columns it is often a real
  measurement. The one exception is retention_time and monoisotopic_mass
  on IQuantifiableRecord types (psmtsv, osmtsv, MsFraggerPsm,
  DiaNnReport), where mzLib documents -1 as absent; those are null,
  matching read-results.

- Units are the writing tool's and are not normalised: retention time is
  minutes for psmtsv/osmtsv/MsFraggerPsm/DiaNnReport and undeclared for
  the rest. read-records carries no retention_time_unit; use a typed
  view when units matter.

- Enums cross as their member names, dates as ISO 8601 strings, and a
  list of scalars as ONE ';'-joined string, on the wire and in out
  alike; split it yourself.

- A column named in absent_fields is null in every row, whatever default
  mzLib filled in.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  path missing or blank; file or directory not found

- `mzlib_usage_error` (usage):

  mzLib does not recognise the file type (the message points at readers
  formats)

- `mzlib_usage_error` (usage):

  offset, limit or out given with no value; offset or limit \< 0;
  non-integer value; out resolves to the same file as path

- `mzlib_bridge_error` (correctness):

  the file is recognised but mzLib fails to parse it; the error type is
  mzLib's own (e.g. HeaderValidationException, MzLibException)

- `mzlib_bridge_error` (correctness):

  a vendor format that cannot load on this platform surfaces as its real
  failure (Bruker/timsTOF native DLLs are Windows x64 only)

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

- `mzlib_usage_error` (usage):

  with paths-stdin, inputs of more than one record type: refused before
  any file is parsed, naming the groups, whatever on-error says

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Columns are this format's own fields under mzLib's names. Two formats
  do not share columns, and a same-named column need not mean the same
  thing; use a typed view (read-results, read-features, read-matches,
  read-spectra) for comparable columns.

- No field is silently dropped: anything that cannot be a column is in
  excluded_fields, whose verb names where it IS carried. At mzLib
  1.0.592 that includes the most important values of three new readers,
  which are dictionaries:
  MetaMorpheusQuantifiedProteinGroups.sample_groups
  (read-protein-groups, and read-occupancy for the occupancy cells),
  FlashLFQQuantifiedPeptide.samples (read-quantified-peptides) and
  MzIdentML/MzIdentMLGz.scores (read-matches `--scores`).

- record_count counts what mzLib parsed, not lines. What mzLib drops is
  reported beside it: rows_not_read counts malformed psmtsv/MSFragger
  lines, and skipped_count/skipped list the mzIdentML items it could not
  represent.

- Since mzLib 1.0.592 (#1346), pro_forma on psmtsv/osmtsv from
  MetaMorpheus 1.1.11 or earlier is computed from full_sequence instead
  of null; ambiguous or unconvertible rows stay null. Measured cost:
  12.9 s -\> 15.6 s on a 271,551-row psmtsv.

- Since mzLib 1.0.592 (#1345), mbr_score on FlashLFQQuantifiedPeak is
  null, not 0, when a peak has no score, and current-format peaks tables
  (no MBR Score column) read, with seven more columns.

- offset/limit window after a full parse; memory and time scale with the
  file, not the window.

- pro_forma is never absent at mzLib 1.0.592: when a psmtsv has no
  ProForma column it is computed from Full Sequence, which the reader
  REQUIRES (SpectrumMatchFromTsv.cs:194), so a file without both fails
  to read instead.

- Many files: `--paths-stdin` needs every input to have the same record
  type (the long table has one column set); a mixed list is a usage
  error naming the groups. Otherwise as every bulk verb: one process,
  `--threads` at a time, output identical at any thread count.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_records`; many files:
  `pymzlib.readers.read_records_many`

- Rust (mzLibRust): `mzlib::readers::read_records_with` with
  `ReadOptions`; many files: `mzlib::readers::read_records_many`

- R (mzLibR): `readers_read_records`; many files:
  `readers_read_records_many`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The three mzLib 1.0.592 fixtures (mzid_gz, mm_protein_groups,
  mm_peptides) exist only on pyMzLib#61 until it merges; mzLibRust and
  mzLibR have not copied them.

- absent_fields for read-records is judged only for record types mapped
  with CsvHelper \[Name\] attributes (FlashLFQ peaks, MSFragger,
  MsPathFinderT, TopFD, the \#1347 tables, ...). psmtsv/osmtsv parse
  their own header dictionary, so their optional columns are not judged
  yet; an empty list there means no basis, not nothing absent.

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

## See also

[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)

## Examples

``` r

batch <- readers_read_records_many(
  c("PXD078927_msgf_1_1_0.mzid", "PXD078927_msgf_1_1_0.mzid.gz"),
  threads = 2
)
batch
#> <mzlibr_read_batch> readers read-records
#>   2 inputs: 2 read, 0 failed (on_error = "fail")
#>   24 rows in `records`, 24 records
table(batch$records$source_index)
#> 
#>  1  2 
#> 12 12 
```
