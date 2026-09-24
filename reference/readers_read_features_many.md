# Read MS1 features from many files, in one bridge call

Read deconvolved MS1 features through the cross-format ms1_features
view: m/z, charge, retention-time range, apex intensity and isotope
count.

The many-files form of
[`readers_read_features`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md).

## Usage

``` r
readers_read_features_many(paths, out = NULL, threads = 1, on_error = "fail",
  timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to `_ms1.feature` or Dinosaur
  `.feature.tsv` files.

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
into `paths`) and `source_path`, then the columns of
[`readers_read_features`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md):
`mz` in m/z, `charge`, `retention_time_start` and `retention_time_end`
in each file's own `retention_time_unit`, and `intensity` in the
instrument's intensity units. `files` has one row per input, with its
`retention_time_unit` - `"unknown"` for `_ms1.feature`. `file_count`,
`read_count` and `failed_count` count files; `record_count` and
`returned_count` count features; `row_count` counts the rows of
`records`.

## Wraps

Wire verb `readers read-features` with `--paths-stdin`. Generated from
the bridge's verb spec `readers.read-features.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`IMs1FeatureFile.GetMs1Features`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/IMs1FeatureFile.cs)
  in `mzLib/Readers/BaseClasses/IMs1FeatureFile.cs` at mzLib `23c2490e`

- [`Ms1Feature.GetSingleChargeFeatures`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/ExternalResults/IndividualResultRecords/Ms1Feature.cs)
  in
  `mzLib/Readers/ExternalResults/IndividualResultRecords/Ms1Feature.cs`
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

  int; in **features**; never `NA`. Features summed over the files read
  (BULK.md section 2).

- `returned_count`:

  int; in **features**; never `NA`. Features returned in columns, summed
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

- `mz`:

  float; in **m/z**; never `NA`. Monoisotopic m/z of this single-charge
  feature.

- `charge`:

  int; in **charge**; never `NA`. Charge state.

- `retention_time_start`:

  float; in **retention_time_unit**; never `NA`. Feature start, in
  retention_time_unit.

- `retention_time_end`:

  float; in **retention_time_unit**; never `NA`. Feature end, in
  retention_time_unit.

- `intensity`:

  float; in **intensity (instrument units)**; `NA` when the file has no
  apex intensity (FLASHDeconv schema, then in absent_fields) or this
  feature has none. Apex intensity (Dinosaur intensityApex; TopFD
  Apex_intensity).

- `number_of_isotopes`:

  int; `NA` when \_ms1.feature: mzLib never sets it (then in
  absent_fields). Isotopes in the feature.

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

  the file has no ms1_features view

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

- \_ms1.feature rows are expanded to one per charge in \[ChargeStateMin,
  ChargeStateMax\] (Ms1Feature.cs:84); Dinosaur is one row per file
  line.

- intensity is the per-charge APEX intensity (Ms1Feature.cs:86); a
  FLASHDeconv file has none and names intensity in absent_fields instead
  of reporting mzLib's zero.

- number_of_isotopes is absent for every \_ms1.feature
  (Ms1Feature.cs:91).

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

- Python (pyMzLib): `pymzlib.readers.read_features`; many files:
  `pymzlib.readers.read_features_many`

- Rust (mzLibRust): `mzlib::readers::read_features_with` with
  `ReadOptions`; many files: `mzlib::readers::read_features_many`

- R (mzLibR): `readers_read_features`; many files:
  `readers_read_features_many`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- retention_time_start/\_end declare their unit as the envelope's
  retention_time_unit, not a fixed unit, because it genuinely varies per
  file.

## See also

[`readers_read_features`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md)

## Examples

``` r

if (FALSE) { # \dontrun{
# Your own files: no recording of a many-file read of this view exists to replay.
batch <- readers_read_features_many(c("a.feature.tsv", "b.feature.tsv"), threads = 2)
batch$files[, c("path", "retention_time_unit")]
} # }
```
