# Read many peptide tables into one long table, in one bridge call

Read a FlashLFQ QuantifiedPeptides.tsv (or MetaMorpheus
AllQuantifiedPeptides.tsv) as one row per peptide per sample, with the
sample's intensity and detection type as columns.

The many-files form of
[`readers_read_quantified_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md):
one bridge process for the whole list, `threads` files at a time, rows
grouped by file in input order.

## Usage

``` r
readers_read_quantified_peptides_many(paths, out = NULL, threads = 1, on_error = "fail",
  timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to FlashLFQ or MetaMorpheus peptide
  tables.

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

An `mzlibr_read_batch`: `records`, the long table, begins with
`source_index` (1-based, into `paths`) and `source_path`, then the
columns of
[`readers_read_quantified_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md) -
`intensity` in the instrument's intensity units, `retention_time` in
minutes. `files` has one row per input with its per-file facts.
`file_count`, `read_count` and `failed_count` count files;
`record_count` counts peptides over the files read and `returned_count`
the peptides returned; `row_count` counts rows.

## Wraps

Wire verb `readers read-quantified-peptides` with `--paths-stdin`.
Generated from the bridge's verb spec
`readers.read-quantified-peptides.yaml` (bridge commit `5db922d4cfe1`)
by `scripts/build-man.R`; the spec owns these facts, and all three
bindings render the same ones.

- [`QuantifiedPeptideFile.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/ResultFiles/QuantifiedPeptideFile.cs)
  in
  `mzLib/Readers/InternalResults/ResultFiles/QuantifiedPeptideFile.cs`
  at mzLib `23c2490e`

- [`QuantifiedPeptideFromTsv.Samples`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/IndividualResultRecords/QuantifiedPeptideFromTsv.cs)
  in
  `mzLib/Readers/InternalResults/IndividualResultRecords/QuantifiedPeptideFromTsv.cs`
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

  int; in **peptides**; never `NA`. Peptides summed over the files read
  (BULK.md section 2).

- `returned_count`:

  int; in **peptides**; never `NA`. Peptides returned in columns, summed
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

- `sequence`:

  string; never `NA`. Full sequence with modifications
  (QuantifiedPeptideFromTsv.Sequence).

- `base_sequence`:

  string; never `NA`. Unmodified sequence.

- `peak_order`:

  int; `NA` when not IsoTracker output: no Peak Order column (then in
  absent_fields). IsoTracker's peak order.

- `protein_groups`:

  string; `NA` when no Protein Groups column, or a blank cell. Protein
  groups the peptide maps to, ';'-joined as written.

- `gene_names`:

  string; `NA` when no Gene Names column, or a blank cell. Gene names as
  written.

- `organism`:

  string; `NA` when no Organism column, or a blank cell. Organism as
  written.

- `sample_label`:

  string; `NA` when the file has no per-sample columns (then in
  absent_fields). The sample this row measures.

- `intensity`:

  float; in **intensity (instrument units)**; `NA` when the cell was
  blank. Intensity\_ as written. FlashLFQ writes 0 for a peptide it did
  not quantify, and 0 crosses as 0: read detection_type.

- `detection_type`:

  string; `NA` when the cell was blank. MSMS, MBR, NotDetected, ...
  verbatim (Detection Type\_).

- `retention_time`:

  float; in **min**; `NA` when not IsoTracker output (then in
  absent_fields) or a blank cell. RetentionTime (min)\_ for this sample.

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

  the file is not a FlashLFQ peptide table

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

- intensity 0 is NOT a measured zero: FlashLFQ writes a literal 0 for a
  peptide it did not quantify in a sample, and mzLib keeps it
  (QuantifiedPeptideFromTsv.cs:14). Filter on detection_type before a
  mean or a log.

- retention_time is IsoTracker output only
  (QuantifiedPeptideFromTsv.cs:61); plain FlashLFQ and MetaMorpheus
  tables name it in absent_fields, and peak_order likewise.

- Long format: offset and limit count PEPTIDES, and returned_count
  counts rows.

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

- Python (pyMzLib): `pymzlib.readers.read_quantified_peptides`; many
  files: `pymzlib.readers.read_quantified_peptides_many`

- Rust (mzLibRust): `mzlib::readers::read_quantified_peptides_with` with
  `ReadOptions`; many files:
  `mzlib::readers::read_quantified_peptides_many`

- R (mzLibR): `readers_read_quantified_peptides`; many files:
  `readers_read_quantified_peptides_many`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- quant median-polish still parses QuantifiedPeptides.tsv with a private
  reader (Quantification.cs); switching it to mzLib's
  QuantifiedPeptideFile is out of this verb's scope.

- No IsoTracker fixture exists, so retention_time and peak_order have
  only been seen absent.

## See also

[`readers_read_quantified_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)

## Examples

``` r

if (FALSE) { # \dontrun{
# Your own files: no recording of a many-file read of this table exists to replay.
batch <- readers_read_quantified_peptides_many(
  c("search_1/AllQuantifiedPeptides.tsv", "search_2/AllQuantifiedPeptides.tsv"),
  threads = 2, on_error = "skip")
batch$files[, c("path", "record_count", "error_message")]
} # }
```
