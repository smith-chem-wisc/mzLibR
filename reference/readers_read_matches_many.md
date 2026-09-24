# Read identifications from many files, in one bridge call

Read identifications through the cross-format spectral_match view -
scan, sequences, accession, decoy flag, modifications and the q-value,
rank and threshold a format records - optionally with each engine's
scores as long rows.

The many-files form of
[`readers_read_matches`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md).

## Usage

``` r
readers_read_matches_many(paths, scores = FALSE, out = NULL, threads = 1,
  on_error = "fail", timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to MsPathFinderT, Casanovo or mzIdentML
  files.

- scores:

  Make the table long by score - one row per match and engine score - as
  [`readers_read_matches`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  does.

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
[`readers_read_matches`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md):
`q_value` a fraction from 0 to 1, and with `scores = TRUE` `score_value`
in the engine's own units, named by `score_name`. `files` has one row
per input, with `skipped_count` in items for mzIdentML. `file_count`,
`read_count` and `failed_count` count files; `record_count` and
`returned_count` count matches; `row_count` counts rows, more than
matches with `scores = TRUE`.

## Wraps

Wire verb `readers read-matches` with `--paths-stdin`. Generated from
the bridge's verb spec `readers.read-matches.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`ISpectralMatch`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/ISpectralMatch.cs)
  in `mzLib/Readers/BaseClasses/ISpectralMatch.cs` at mzLib `23c2490e`

- [`MzIdentMLResultFile.SkippedMatches`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/ExternalResults/ResultFiles/MzIdentMLResultFile.cs)
  in `mzLib/Readers/ExternalResults/ResultFiles/MzIdentMLResultFile.cs`
  at mzLib `23c2490e`

- [`MzIdentMLRecord.Scores`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/ExternalResults/IndividualResultRecords/MzIdentMLRecord.cs)
  in
  `mzLib/Readers/ExternalResults/IndividualResultRecords/MzIdentMLRecord.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `scores`:

  flag; default `FALSE`. Make the table long by score: one row per match
  and engine score, adding match_index, score_name and score_value. Only
  mzIdentML records scores (#1306); other formats keep one row per match
  with the two columns in absent_fields. offset and limit still count
  matches.

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

  int; in **matches**; never `NA`. Matches summed over the files read
  (BULK.md section 2).

- `returned_count`:

  int; in **matches**; never `NA`. Matches returned in columns, summed
  over the files; 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns (more than
  returned_count for a long view); 0 when written to out.

- `on_error`:

  string; never `NA`. fail or skip, as requested. `--threads` is
  deliberately not echoed: the envelope is byte-identical at every
  thread count.

- `scores_included`:

  bool; never `NA`. Whether `--scores` was given.

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
  block (path, file_type, reader, skipped_count, skipped, record_count,
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

- `file_name_without_extension`:

  string; never `NA`. The spectra file of the match.

- `one_based_scan_number`:

  int; never `NA`. Scan number; see caveats for mzIdentML index= and
  Casanovo.

- `base_sequence`:

  string; never `NA`. Unmodified sequence.

- `full_sequence`:

  string; never `NA`. Modified sequence.

- `accession`:

  string; `NA` when the property threw (then in failed_fields). Protein
  accession(s).

- `is_decoy`:

  bool; `NA` when every format but MsPathFinderT (then in
  absent_fields). Decoy flag.

- `modifications`:

  string; never `NA`. 'position:name' pairs, ';'-joined, 1 = N-terminus.

- `modification_count`:

  int; never `NA`. Modifications on the match.

- `q_value`:

  float; in **fraction (0 to 1)**; `NA` when the format has no q-value
  (in absent_fields) or this mzIdentML item reports none. PSM-level
  q-value (mzIdentML MS:1002354 or a child; MsPathFinderT QValue).

- `rank`:

  int; `NA` when not mzIdentML (in absent_fields). The item's rank;
  filter on rank == 1.

- `pass_threshold`:

  bool; `NA` when not mzIdentML (in absent_fields). The item's
  passThreshold.

- `match_index`:

  int; never `NA`. 0-based position of the match among the file's
  matches. *Present only with `scores`.*

- `score_name`:

  string; `NA` when the match has no scores, or the format records none
  (in absent_fields). Engine score name as written, e.g.
  MS-GF:SpecEValue. *Present only with `scores`.*

- `score_value`:

  float; in **engine-defined (see score_name)**; `NA` when as
  score_name. The score's value. *Present only with `scores`.*

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

  the file has no spectral_match view

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

- Nothing is FDR-filtered. q_value is filled only by mzIdentML and an
  MsPathFinderT file with a QValue column; an \_IcTarget/\_IcDecoy file
  has none and names it absent (mzLib would report 0).

- mzIdentML: every SpectrumIdentificationItem is a row, lower ranks and
  failed thresholds included (MzIdentMLResultFile.cs:177). Filter on
  rank == 1 and pass_threshold.

- mzIdentML: is_decoy is absent because the isDecoy attribute defaults
  to false when omitted (MzIdentMLResultFile.cs:173).

- mzIdentML: items mzLib cannot represent are skipped, not failed, and
  listed in skipped (MzIdentMLResultFile.cs:123).

- MsPathFinderT: is_decoy comes from an XXX protein-name prefix
  (MsPathFinderTResult.cs:92); accession assumes a UniProt header and
  lands in failed_fields otherwise.

- Casanovo: one_based_scan_number is the mzTab index plus one
  (CasanovoMzTabFile.cs:116).

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

- Python (pyMzLib): `pymzlib.readers.read_matches`; many files:
  `pymzlib.readers.read_matches_many`

- Rust (mzLibRust): `mzlib::readers::read_matches_with` with
  `MatchOptions`; many files: `mzlib::readers::read_matches_many`

- R (mzLibR): `readers_read_matches`; many files:
  `readers_read_matches_many`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- readers_matches_mzid_scores.json (`--scores`) is replayed through the
  pyMzLib conftest rather than listed here: check_verbs.py excludes
  present_when columns from its column check, so a `--scores` fixture
  cannot be a spec example yet.

- q_value, rank and pass_threshold join a view whose other columns are
  Readers.ISpectralMatch members; they are read from the record types
  that carry them. If mzLib adds them to the interface, the per-format
  reads go away.

## See also

[`readers_read_matches`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)

## Examples
