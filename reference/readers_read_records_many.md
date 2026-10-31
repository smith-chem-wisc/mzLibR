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
`84c3b984b6b9`) by `scripts/build-man.R`; the spec owns these facts, and
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

  path; default absent; range `.tsv only, any case (PYB-1)`. Write the
  long table here, one file at a time in input order, so memory holds at
  most threads files. A batch that stops on an error deletes its partial
  table. Must differ from every input. Must end in .tsv, any case
  (PYB-1); any other extension, or none, is a usage error before any
  input is read, and no extension is appended.

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
  DiaNnReport, DiaNnReportParquet), where mzLib documents -1 as absent;
  those are null, matching read-results.

- Units are the writing tool's and are not normalised: retention time is
  minutes for psmtsv/osmtsv/MsFraggerPsm/DiaNnReport/DiaNnReportParquet
  and undeclared for the rest. read-records carries no
  retention_time_unit; use a typed view when units matter.

- Enums cross as their member names, dates as ISO 8601 strings, and a
  list of scalars as ONE ';'-joined string, on the wire and in out
  alike; split it yourself. The one exception is an SDRF row, whose
  header and cells lists are excluded rather than joined (PYB-2): SDRF
  cells contain ';' themselves.

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

- `mzlib_usage_error` (usage):

  out (with path or paths-stdin) whose extension is not .tsv, any case,
  including none (PYB-1); checked before any input is opened

- `mzlib_usage_error` (usage):

  the answer is too large to return as one JSON document (PYB-3 part 1;
  measured: past the ~1.07 G-character .NET string ceiling, not RAM).
  Message: 'this read (N records) is too large to return as one
  document; use `--limit`/`--offset` or `--out`', N the records this
  answer holds (its returned_count, or record_count where the envelope
  has none)

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

- An SDRF reads with ZERO columns and a real record_count: SdrfRow's
  only fields, header and cells, are lists of strings, and a ';'-join
  could not be split back because SDRF cells contain ';'
  (NT=Oxidation;AC=UNIMOD:35). Both are in excluded_fields with reason
  "use sdrf read" (PYB-2, pyMzLib \#76); read the document with sdrf
  read. With out, the file is an empty header line and one empty line
  per record.

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
  REQUIRES (SpectrumMatchFromTsv.cs:197), so a file without both fails
  to read instead.

- Since mzLib 1.0.595 (#1371), psmtsv and osmtsv records carry
  entrapment_fraction: the share of the PSM that is an entrapment
  discovery (DecoyContaminantTargetLabel.EntrapmentFraction), from 0
  to 1. A shared T\|ET row counts half when its sequences differ and 0
  when both proteins carry the same peptide, so count entrapment
  discoveries for an FDP estimate with this column, not is_entrapment.
  It is null where mzLib returns NaN: the Decoy/Contaminant/Target and
  Full Sequence columns name different numbers of candidates, or the
  full sequence is blank or MetaMorpheus's Output too long for Excel
  placeholder beside an entrapment target.

- Since mzLib 1.0.595 (#1371), is_decoy, is_contaminant and
  is_entrapment read every \|-joined parent label: true when any parent
  is D, C or E (DecoyContaminantTargetLabel.IsDecoy, IsContaminant,
  IsEntrapment). On a protein-group table T\|C is now a contaminant and
  T\|ET entrapment; on a psmtsv the E in Output too long for Excel no
  longer reads as entrapment. Labels MetaMorpheus writes read as before.

- Since mzLib 1.0.594 (#1403), a computed pro_forma names
  UniProt-sourced modifications by UNIMOD accession:
  BottomUpExample.psmtsv's first row reads YPIEH\[UNIMOD:34\]GIVTNWDDMEK
  where it read YPIEH\[UniProt:Tele-methylhistidine on H\]GIVTNWDDMEK. A
  pro_forma column the file carries is passed through as written.

- A DIA-NN 2.x report.parquet (DiaNnReportParquet, mzLib \#1384) reads
  into the TSV report's record type, DiaNnPrecursor, so the columns are
  the same; those for fields 2.x no longer writes (File.Name,
  PG.Quantity, MS2.Scan and others, DiaNnParquetReportFile.cs:16) hold
  mzLib's not-reported default: null, except ms2_scan_number, which
  reads 0. absent_fields does not name them (see absent_fields). Decoy
  rows (`--report-decoys`) are skipped by mzLib, so record_count counts
  targets.

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

- absent_fields for read-records is judged only for record types mapped
  with CsvHelper \[Name\] attributes (FlashLFQ peaks, MSFragger,
  MsPathFinderT, TopFD, the \#1347 tables, ...). psmtsv/osmtsv parse
  their own header dictionary, so their optional columns are not judged
  yet; an empty list there means no basis, not nothing absent.

- PYB-3 part 2: refusing a too-large answer BEFORE serialising is not
  built (pyMzLib main b18f013, Program.cs SuccessJson); part 1 catches
  the failure after serialising. A threshold per record will not work:
  dataRepo's re-run on 0.4.0 (dataRepo to pyMzLib 018, 2026-10-04)
  measured 725 and 805 characters per record in two psmtsv files from
  one producer. It proposed source size instead, 0.63 to 0.71 output
  characters per source byte, which puts the limit at about 1.51 to 1.70
  GB of source text. That is two files of one format, and no threshold
  has been agreed.

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
