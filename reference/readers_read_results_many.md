# Read many result files into the uniform record view, in one bridge call

Read a result file through the cross-format quantifiable view: sequence,
retention time, charge, theoretical mass, decoy flag and protein groups,
the same columns for every format that offers it.

The many-files form of
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md).

## Usage

``` r
readers_read_results_many(paths, out = NULL, threads = 1, on_error = "fail",
  timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to files offering the `"quantifiable"`
  view.

- out:

  Write the long table here as tab-separated text, one file at a time,
  and return only a summary.

- threads:

  Files read at once. `1`, the default, holds one whole file in memory
  at a time; `-1` uses one per core. The table does not depend on it.

- on_error:

  `"fail"`, the default, stops at the first file that cannot be read.
  `"skip"` records the failure in that file's row of `files` and reads
  the rest.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_read_batch`. `records` begins with `source_index` (1-based,
into `paths`) and `source_path`, then the columns of
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md):
`retention_time` in each file's own `retention_time_unit` - minutes for
all four formats at mzLib 1.0.592 - `charge_state` a charge,
`monoisotopic_mass` in Da. Convert times with
[`readers_retention_time_in_minutes`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md),
which applies each file's unit to its rows. `files` has one row per
input: `record_count` in records, `rows_not_read` in rows,
`retention_time_unit`, `caveats`, `absent_fields` and the error columns.
`file_count`, `read_count` and `failed_count` count files;
`record_count` and `returned_count` count records over the files read;
`row_count` counts the rows of `records`.

## Wraps

Wire verb `readers read-results` with `--paths-stdin`. Generated from
the bridge's verb spec `readers.read-results.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`FileReader.ReadQuantifiableResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`IQuantifiableRecord`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/IQuantifiableRecord.cs)
  in `mzLib/Readers/BaseClasses/IQuantifiableRecord.cs` at mzLib
  `23c2490e`

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
  block (path, file_type, reader, record_count, rows_not_read,
  retention_time_unit, caveats, column_names, absent_fields,
  failed_fields, excluded_fields, error). An unread input has error set
  and every fact it could not establish null or empty.

## Columns of `records`

The first two columns say which input each row came from -
`source_index` is 1-based in R, so it indexes `files` directly - and
rows are grouped by input, in input order, whatever `threads` is.

- `source_index`:

  int; never `NA`. 0-based position of the row's input in the stdin
  list.

- `source_path`:

  string; never `NA`. Absolute path of the row's input.

- `file_name`:

  string; never `NA`. The spectra file the record came from, as the
  format writes it (not a join key across formats).

- `base_sequence`:

  string; never `NA`. Unmodified sequence.

- `full_sequence`:

  string; never `NA`. Modified sequence; on an ambiguous psmtsv row
  every candidate '\|'-joined.

- `retention_time`:

  float; in **min**; `NA` when mzLib's -1 'absent' sentinel, or a
  non-finite value. Retention time, in retention_time_unit.

- `charge_state`:

  int; in **charge**; never `NA`. Precursor charge.

- `monoisotopic_mass`:

  float; in **Da**; `NA` when mzLib's -1 'absent' sentinel. THEORETICAL
  peptide mass, not the observed precursor.

- `is_decoy`:

  bool; `NA` when the format cannot report decoys: MSFragger and DIA-NN,
  where is_decoy is also in absent_fields. Decoy flag where the format
  carries one.

- `protein_accessions`:

  string; never `NA`. ';'-joined accessions.

- `gene_name`:

  string; never `NA`. ';'-joined gene names.

- `organism`:

  string; never `NA`. ';'-joined organisms.

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

  the file has no quantifiable view (the message names the views it does
  have)

- `mzlib_bridge_error` (correctness):

  mzLib recognises the file and fails to parse it

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

- No q-value, PEP or score: nothing here is FDR-filtered.

- is_decoy is absent for MSFragger and DIA-NN, which write no
  target/decoy label; null there means unknown, not target.

- monoisotopic_mass is the theoretical mass in every format; the
  observed precursor is not in this view.

- psmtsv: monoisotopic_mass keeps the first candidate of an ambiguous
  row while full_sequence keeps them all (SpectrumMatchFromTsv.cs:119,
  :194).

- A malformed psmtsv row is dropped silently by mzLib; rows_not_read
  counts them.

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

- Python (pyMzLib): `pymzlib.readers.read_results`; many files:
  `pymzlib.readers.read_results_many`

- Rust (mzLibRust): `mzlib::readers::read_results_with` with
  `ReadOptions`; many files: `mzlib::readers::read_results_many`

- R (mzLibR): `readers_read_results`; many files:
  `readers_read_results_many`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

## See also

[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)

## Examples

``` r

if (FALSE) { # \dontrun{
# Your own files: no recording of a many-file read of this view exists to replay.
batch <- readers_read_results_many(c("run_1/AllPSMs.psmtsv", "run_2/AllPSMs.psmtsv"),
  threads = 2)
table(batch$records$source_path)
} # }
```
