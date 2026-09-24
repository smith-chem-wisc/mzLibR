# Read identifications, in the cross-format `spectral_match` view

Read identifications through the cross-format spectral_match view -
scan, sequences, accession, decoy flag, modifications and the q-value,
rank and threshold a format records - optionally with each engine's
scores as long rows.

Six file types offer it: MsPathFinderT's targets, decoys and combined
results, Casanovo's `.mztab`, and mzIdentML `.mzid` and `.mzid.gz` - the
format most search engines can export. These are the identification
formats that share no \*file\*-level interface, so
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
cannot reach them.

## Usage

``` r
readers_read_matches(path, limit = NULL, offset = 0, scores = FALSE, out = NULL,
  timeout = NULL)
```

## Arguments

- path:

  Path to an MsPathFinderT `_IcTarget.tsv` / `_IcDecoy.tsv` /
  `_IcTDA.tsv`, a Casanovo `.mztab`, or an mzIdentML `.mzid` /
  `.mzid.gz`.

- limit:

  Maximum matches to return. `NULL`, the default, returns all of them.

- offset:

  Matches to skip.

- scores:

  Make the table long by score: one row per match and engine score,
  adding `match_index`, `score_name` and `score_value`. `limit` and
  `offset` still count matches. Only mzIdentML records scores; for other
  formats the two score columns are `NA` and named in `absent_fields`.

- out:

  Write a tab-separated table here and return only a summary.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_match_records`. `record_count` counts the matches in the
whole file and `returned_count` the matches returned, starting `offset`
matches in; `row_count` counts the rows in `records` - more rows than
matches with `scores = TRUE`. For mzIdentML, `skipped_count` counts the
identification items mzLib did not turn into rows, and `skipped` lists
each with its reason; both are `NA`/`NULL` for formats that keep no such
list.

`records` is a data.frame with `file_name_without_extension`,
`one_based_scan_number`, `base_sequence`, `full_sequence`, `accession`,
`is_decoy`, `modifications`, `modification_count`, `q_value` - a
fraction from 0 to 1 - `rank` and `pass_threshold`; with `scores = TRUE`
also `match_index`, `score_name` and `score_value`, in whatever units
the engine's score has. A column this format has no source for is `NA`
in every row and is named in `absent_fields`. `rows_not_read`, the rows
that did not become matches, is not counted for this view and is `NA`.

## Nothing here is FDR-filtered

Every match is a row, as the file lists it: for mzIdentML that includes
lower ranks and items that failed the engine's threshold. `q_value` is
filled only by mzIdentML and by an MsPathFinderT file that carries a
`QValue` column; `rank` and `pass_threshold` only by mzIdentML. Filter
before you report - for mzIdentML, on `rank == 1` and `pass_threshold`.

## Three is_decoy traps, all reported in caveats

**MsPathFinderT** infers decoys from the protein \*name\* - mzLib
reports a decoy when the name starts with `XXX`. A database whose decoys
carry a different prefix reads entirely as targets.

**Casanovo** is de novo and writes no target/decoy label at all, so
`is_decoy` is `NA` and in `absent_fields` - never a `FALSE` that would
read as "target".

**mzIdentML**'s `isDecoy` attribute defaults to false when a writer
omits it, so mzLib cannot tell "target" from "not recorded"; `is_decoy`
is `NA` and in `absent_fields` here too.

## Wraps

Wire verb `readers read-matches`. Generated from the bridge's verb spec
`readers.read-matches.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

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

- `path`:

  path; required. An MsPathFinderT \_IcTarget/\_IcDecoy/\_IcTda.tsv, a
  Casanovo .mztab, or an mzIdentML .mzid / .mzid.gz.

- `offset`:

  int; in **matches**; default `0`; range `>= 0`. Skip this many
  matches. A window, not a cursor: mzLib parses the whole file on every
  call, so paging re-reads it.

- `limit`:

  int; in **matches**; default absent; range `>= 0`. Return at most this
  many matches. Absent means no limit; there is deliberately no default
  cap.

- `scores`:

  flag; default `FALSE`. Make the table long by score: one row per match
  and engine score, adding match_index, score_name and score_value. Only
  mzIdentML records scores (#1306); other formats keep one row per match
  with the two columns in absent_fields. offset and limit still count
  matches.

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

- `skipped_count`:

  int; in **items**; `NA` when the format keeps no skip list (every
  format but mzIdentML). mzIdentML items mzLib did not represent
  (crosslinks, unresolvable modifications, substitutions);
  record_count + skipped_count is the items in the file (#1313).

- `skipped`:

  object\[\]; `NA` when as skipped_count.
  {spectrum_identification_item_id, spectrum_id, reason} per skipped
  item.

- `record_count`:

  int; in **matches**; never `NA`. Matches in the whole file, before the
  window.

- `rows_not_read`:

  int; in **rows**; `NA` when never counted for this view. Data rows
  that did not become records.

- `retention_time_unit`:

  string; `NA` when the view has no time column. Always null for this
  view.

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

- `scores_included`:

  bool; never `NA`. Whether `--scores` made the table long by score.

- `returned_count`:

  int; in **matches**; never `NA`. Matches returned in columns (the unit
  offset and limit count in); 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns: more than
  returned_count in a long table, where one record gives several rows; 0
  when written to out.

- `offset`:

  int; in **matches**; never `NA`. The offset applied, in matches.

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
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`readers_read_matches_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches_many.md).

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
any work. For many files call
[`readers_read_matches_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_matches`; many files:
  `pymzlib.readers.read_matches_many`

- Rust (mzLibRust): `mzlib::readers::read_matches_with` with
  `MatchOptions`; many files: `mzlib::readers::read_matches_many`

- R (mzLibR): `readers_read_matches`; many files:
  [`readers_read_matches_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches_many.md)

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

[`readers_read_matches_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches_many.md),
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)

## Examples

``` r

matches <- readers_read_matches("PXD078927_msgf_1_1_0.mzid", limit = 3)
matches
#> <mzlibr_match_records> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\DataFiles\PXD078927_msgf_1_1_0.mzid (MzIdentML)
#>   12 records in the file, 3 returned
#>   ! truncated - records were left behind
#>   ! NOTHING here is FDR-filtered. q_value is the only confidence field this view carries, and only mzIdentML and an MsPathFinderT file with a QValue column fill it (absent_fields says when it is empty). Every one of these formats records scores this view does not expose; read-records has them. Filter before you report.
#>   ! is_decoy is null for this format. mzIdentML's isDecoy attribute is optional and defaults to false, and mzLib reports a decoy only when every peptide evidence says so (MzIdentMLResultFile.cs:173), so false cannot be told apart from 'not stated'. read-records carries mzLib's boolean for a caller who knows the writer sets it.
#>   ! Every SpectrumIdentificationItem is a row, not only the matches the submitter accepted: lower-ranked candidates and items that fail the threshold are here too (MzIdentMLResultFile.cs:177). Filter on rank == 1 and pass_threshold before counting identifications.
#>   ! one_based_scan_number is parsed from the nativeID (MzIdentMLResultFile.cs:159). 'scan=N' gives N, but 'index=N', which peak-list input carries, is a zero-based position in the file and gives N + 1, not an instrument scan number. -1 means the nativeID had neither.
#>   ! Items mzLib cannot represent as one linear match are skipped, not failed: crosslinks, modifications without a resolvable UNIMOD accession, substitutions, and two modifications on one residue (MzIdentMLResultFile.cs:123). They are not rows; skipped_count and skipped name each one and why, so record_count plus skipped_count is the number of items in the file.
#>   ! The engine's own scores (for example MS-GF:SpecEValue) have no common name across search engines (MzIdentMLResultFile.cs:179). Pass scores=true for them as long rows, one per match and score. q_value is null on an item that reports none.
#>   ! accession joins every protein the item's peptide evidence names with '|' (MzIdentMLRecord.cs:45), the same character a UniProt header uses inside one accession, so the cell cannot be split back into proteins reliably.
#>   absent from this file (NA in every row): is_decoy
matches$records[, c("one_based_scan_number", "base_sequence", "q_value", "rank")]
#>   one_based_scan_number  base_sequence   q_value rank
#> 1                 14316   HSNLNDATYQRT 0.0000000    1
#> 2                 14316 KASAGQISVQPTFS 0.6588785    2
#> 3                 14316   TRQYTADNLNSH 0.7383230    3

# One row per match and engine score:
scored <- readers_read_matches("PXD078927_msgf_1_1_0.mzid", limit = 1, scores = TRUE)
scored$records[, c("match_index", "score_name", "score_value")]
#>   match_index        score_name  score_value
#> 1           0    MS-GF:RawScore 1.150000e+02
#> 2           0 MS-GF:DeNovoScore 1.230000e+02
#> 3           0  MS-GF:SpecEValue 3.041116e-14
#> 4           0      MS-GF:EValue 8.393480e-12
#> 5           0      MS-GF:QValue 0.000000e+00
#> 6           0   MS-GF:PepQValue 0.000000e+00
#> 7           0      IsotopeError 1.000000e+00

# Casanovo writes no target/decoy label, so is_decoy is NA rather than a false FALSE.
denovo <- readers_read_matches("Casanovo_5.0.0.mztab")
denovo$absent_fields
#> [1] "is_decoy"       "q_value"        "rank"           "pass_threshold"
```
