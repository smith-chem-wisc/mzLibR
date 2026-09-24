# Read any file mzLib recognises, into that format's own fields

Read any file mzLib recognises into a table of that format's own record
fields, naming every field that could not become a column.

The exhaustive verb: if
[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
succeeds on a path, this reads it. All 31 file types, including the 14
that belong to no cross-format view at all - TopPIC, Crux, MSFragger's
peptide and protein tables, the FlashDeconv formats, SDRF - which no
other `readers_` function can touch.

## Usage

``` r
readers_read_records(path, limit = NULL, offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to any file mzLib recognises. A Bruker `.d` directory is also
  accepted.

- limit:

  Maximum records to return. `NULL`, the default, returns all of them.

- offset:

  Records to skip. A window, not a cursor - see
  [`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md).

- out:

  Write a **tab-separated** table here and return only a summary. Read
  it with [`read.delim()`](https://rdrr.io/r/utils/read.table.html), not
  [`read.csv()`](https://rdrr.io/r/utils/read.table.html).

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**For SDRF, use
[`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
instead.** This function joins each SDRF row's cells into one
semicolon-separated string, and SDRF's `NT=...;AC=...` grammar puts
semicolons inside cells, so the joined string cannot be split back
apart.

## Value

An `mzlibr_native_records`. `records` is a data.frame of this format's
own fields, or `NULL` when `out` was given. `record_count` counts the
records in the whole file, `returned_count` the records returned,
starting `offset` records in. `column_names`, `record_type`, `views`,
`excluded_fields` and `failed_fields` describe the table.

## The columns are not uniform, by design

They are this format's own mzLib record fields, under mzLib's names in
`snake_case`: a TopPIC file gives 36 columns, a Crux file 23, an
experiment annotation 5. Read `column_names`, and use
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md),
[`readers_read_features`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md)
or
[`readers_read_matches`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
when you need columns that mean the same thing across formats.

Because the names are mzLib's own, they are **cross-referenceable
against the mzLib source**: a column called `e_value` is
`ToppicPrsm$EValue`, and `record_type` names the class to look in.

## Nothing is silently dropped

`excluded_fields` is a data.frame of the fields that could **not**
become columns, each with the reason - a nested object or a dictionary
has no faithful column shape, and inventing one would mean publishing a
schema mzLib does not have. A column that simply vanished would be
indistinguishable from a field the format does not have.

`failed_fields` names the fields that **raised** while being read, with
the exception type. Several mzLib properties are computed and assume a
UniProt-style FASTA header - Crux's and MsPathFinderT's `accession` are
both `protein_id` split on a pipe - so on other databases they throw.
Those cells arrive `NA` rather than taking the whole read down, but a
failure must not look like missing data.

Note that mzLib's documented `-1` "absent" sentinel is **not** mapped to
`NA` here, unlike in
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md).
In a format's own columns `-1` is frequently a real measurement - a mass
difference, a delta, TopPIC's `feature_score` - and nulling those would
destroy data.

## Wraps

Wire verb `readers read-records`. Generated from the bridge's verb spec
`readers.read-records.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`FileReader.ReadResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`IResultFile.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/IResultFile.cs)
  in `mzLib/Readers/BaseClasses/IResultFile.cs` at mzLib `23c2490e`

- [`MzIdentMLResultFile.SkippedMatches`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/ExternalResults/ResultFiles/MzIdentMLResultFile.cs)
  in `mzLib/Readers/ExternalResults/ResultFiles/MzIdentMLResultFile.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. Any of the file types readers formats lists (36 at
  mzLib 1.0.592), including spectra files and the Bruker .d directories.

- `offset`:

  int; in **records**; default `0`; range `>= 0`. Skip this many
  records. A window, not a cursor: mzLib parses the whole file on every
  call, so paging re-reads it. For a large file use out in one call.

- `limit`:

  int; in **records**; default absent; range `>= 0`. Return at most this
  many records. Absent means no limit; there is deliberately no default
  cap.

- `out`:

  path; default absent. Write the selected window as a tab-separated
  table here instead of returning columns (header = column_names). Must
  differ from path. Parent directories are created.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. Absolute path of the input.

- `file_type`:

  string; never `NA`. The SupportedFileType mzLib dispatched, not a
  guess from the extension.

- `reader`:

  string; never `NA`. The mzLib reader class that parsed it.

- `record_type`:

  string; never `NA`. The mzLib record class the columns are properties
  of (e.g. ToppicPrsm, MzIdentMLRecord, ProteinGroupFromTsv), so a
  column is cross-referenceable against the mzLib source.

- `views`:

  string\[\]; never `NA`. The cross-format views this file also offers;
  same values as readers formats.

- `record_count`:

  int; in **records**; never `NA`. Records mzLib parsed from the whole
  file, before the window. Records, not lines: see caveats.

- `returned_count`:

  int; in **records**; never `NA`. Records in columns; 0 when written to
  out.

- `offset`:

  int; in **records**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever records were left out by offset or
  limit.

- `excluded_fields`:

  table; never `NA`. Every property of record_type that could not become
  a column, as {field, type, reason, verb}: nested objects, lists of
  composites, and dictionaries (reason starts "a dictionary"). verb
  names the command that carries the field (sample_groups: readers
  read-protein-groups; samples: readers read-quantified-peptides;
  mzIdentML scores: readers read-matches) or is null. Empty when every
  property projects.

- `error`:

  object; `NA` when always, for a single path: a file that cannot be
  read fails the call instead. {kind, type, message}; non-null only in a
  files\[\] entry under on-error skip.

- `failed_fields`:

  string\[\]; never `NA`. "field: ExceptionType" for every property that
  threw while being read on a returned record; those cells are null.

- `column_names`:

  string\[\]; never `NA`. The record type's projectable public
  properties in snake_case, base class first, in declaration order.

- `records` (wire `columns`):

  table; `NA` when the table was written to out instead. Column name to
  per-record values, in column_names order.

- `output`:

  object; `NA` when out was not given. {path, format: "tsv", row_count}
  when out was given; row_count is the window written, not the whole
  file.

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

## On the wire but not projected yet

The bridge sends these, and this version of mzLibR does not return them
yet:

- `rows_not_read`:

  arrives with the mzLib 1.0.592 port.

- `retention_time_unit`:

  arrives with the mzLib 1.0.592 port.

- `caveats`:

  arrives with the mzLib 1.0.592 port.

- `absent_fields`:

  arrives with the mzLib 1.0.592 port.

- `skipped_count`:

  arrives with the mzLib 1.0.592 port.

- `skipped`:

  arrives with the mzLib 1.0.592 port.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` belongs
to the verb's many-files form.

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
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_records`; many files:
  `pymzlib.readers.read_records_many`

- Rust (mzLibRust): `mzlib::readers::read_records_with` with
  `ReadOptions`; many files: `mzlib::readers::read_records_many`

- R (mzLibR): `readers_read_records`

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

[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md),
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)

## Examples

``` r

prsms <- readers_read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv", limit = 3)
prsms
#> <mzlibr_native_records> E:\CodeReview\pymzlib-readers-wt\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv (ToppicPrsm)
#>   36 columns from ToppicPrsm - no cross-format view
#>   4 records in the file, 3 returned
#>   ! truncated - records were left behind
#>   1 field(s) could not become columns: alternative_identifications
head(prsms$column_names)
#> [1] "file_name_without_extension" "file_path"                  
#> [3] "prsm_id"                     "spectrum_id"                
#> [5] "dissociation_type"           "one_based_scan_number"      
prsms$records[, c("one_based_scan_number", "base_sequence", "e_value")]
#>   one_based_scan_number base_sequence e_value
#> 1                   259        GYASDS  1e+300
#> 2                   270         THIGY  1e+300
#> 3                   354      SAECTKRF  1e+300
prsms$excluded_fields
#>                         field                      type
#> 1 alternative_identifications List<AlternativeToppicId>
#>                                                    reason
#> 1 a list of composite values has no faithful column shape
```
