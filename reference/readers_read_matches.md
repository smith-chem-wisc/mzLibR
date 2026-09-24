# Read identifications, in the cross-format `spectral_match` view

Read identifications through the cross-format spectral_match view -
scan, sequences, accession, decoy flag, modifications and the q-value,
rank and threshold a format records - optionally with each engine's
scores as long rows.

Four file types offer it: MsPathFinderT's targets, decoys and combined
results, and Casanovo's `.mztab`. These are the identification formats
that share no \*file\*-level interface, so
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
cannot reach them.

## Usage

``` r
readers_read_matches(path, limit = NULL, offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to an MsPathFinderT `_IcTarget.tsv` / `_IcDecoy.tsv` /
  `_IcTDA.tsv`, or a Casanovo `.mztab`.

- limit:

  Maximum matches to return. `NULL`, the default, returns all of them.

- offset:

  Matches to skip.

- out:

  Write a tab-separated table here and return only a summary.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_match_records`. `record_count` counts the matches in the
whole file and `returned_count` the matches returned, starting `offset`
matches in.

`records` is a data.frame with `file_name_without_extension`,
`one_based_scan_number`, `base_sequence`, `full_sequence`, `accession`,
`is_decoy` (`NA` where the format records no target/decoy label),
`modifications` and `modification_count`.

## Nothing here is FDR-filtered, and there is nothing to filter on

mzLib's `ISpectralMatch` carries identity fields only. Every one of
these formats records an E-value or q-value somewhere;
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
will give you those columns. Filter before you report.

## Two is_decoy traps, both reported in caveats

**MsPathFinderT** infers decoys from the protein \*name\* - mzLib
reports a decoy when the name starts with `XXX`. A database whose decoys
carry a different prefix reads entirely as targets.

**Casanovo** is de novo and writes no target/decoy label at all. mzLib
leaves the field at its default `FALSE` and never assigns it, so `FALSE`
would mean \*unknown\*; it arrives as `NA` instead - the rule
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
already applies to MSFragger.

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

- wire `--scores`:

  flag; default `FALSE`. Make the table long by score: one row per match
  and engine score, adding match_index, score_name and score_value. Only
  mzIdentML records scores (#1306); other formats keep one row per match
  with the two columns in absent_fields. offset and limit still count
  matches. *Not an argument here: arrives with the mzLib 1.0.592 port.*

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

- `record_count`:

  int; in **matches**; never `NA`. Matches in the whole file, before the
  window.

- `caveats`:

  string\[\]; never `NA`. What this view cannot be trusted to mean for
  THIS file, each citing the mzLib source.

- `column_names`:

  string\[\]; never `NA`. Column order.

- `error`:

  object; `NA` when always, for a single path: a file that cannot be
  read fails the call instead. {kind, type, message}; non-null only in a
  files\[\] entry under on-error skip.

- `returned_count`:

  int; in **matches**; never `NA`. Matches returned in columns (the unit
  offset and limit count in); 0 when written to out.

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

## On the wire but not projected yet

The bridge sends these, and this version of mzLibR does not return them
yet:

- `reader`:

  arrives with the mzLib 1.0.592 port.

- `skipped_count`:

  arrives with the mzLib 1.0.592 port.

- `skipped`:

  arrives with the mzLib 1.0.592 port.

- `rows_not_read`:

  arrives with the mzLib 1.0.592 port.

- `retention_time_unit`:

  arrives with the mzLib 1.0.592 port.

- `absent_fields`:

  arrives with the mzLib 1.0.592 port.

- `failed_fields`:

  arrives with the mzLib 1.0.592 port.

- `excluded_fields`:

  arrives with the mzLib 1.0.592 port.

- `scores_included`:

  arrives with the mzLib 1.0.592 port.

- `row_count`:

  arrives with the mzLib 1.0.592 port.

- `q_value`:

  arrives with the mzLib 1.0.592 port.

- `rank`:

  arrives with the mzLib 1.0.592 port.

- `pass_threshold`:

  arrives with the mzLib 1.0.592 port.

- `match_index`:

  arrives with the mzLib 1.0.592 port.

- `score_name`:

  arrives with the mzLib 1.0.592 port.

- `score_value`:

  arrives with the mzLib 1.0.592 port.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` belongs
to the verb's many-files form.

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
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_matches`; many files:
  `pymzlib.readers.read_matches_many`

- Rust (mzLibRust): `mzlib::readers::read_matches_with` with
  `MatchOptions`; many files: `mzlib::readers::read_matches_many`

- R (mzLibR): `readers_read_matches`

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

[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)

## Examples

``` r

matches <- readers_read_matches("Casanovo_5.0.0.mztab")
matches
#> <mzlibr_match_records> E:\CodeReview\pymzlib-readers-wt\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\Casanovo_5.0.0.mztab (CasanovoMzTab)
#>   5 records in the file, 5 returned
#>   ! There is no score, E-value or q-value in this view, so NOTHING here is FDR-filtered. Readers.ISpectralMatch carries identity fields only; every one of these formats records confidence in columns this view does not expose. read-records has them. Filter before you report.
#>   ! is_decoy is null for this format. Casanovo is de novo and writes no target/decoy label; mzLib's record leaves the field at its default false and never assigns it (CasanovoMzTabRecord.cs:84), so false would mean 'unknown', not 'target'.
#>   ! one_based_scan_number is the mzTab spectrum INDEX plus one, not necessarily the instrument's scan number (CasanovoMzTabFile.cs:116). When Casanovo was run on an MGF the two are unrelated, so do not join this against a raw file on scan number.
#>   ! full_sequence and modifications are resolved by matching Casanovo's mass shifts against mzLib's modification dictionary (CasanovoMzTabFile.cs:124), not read from named annotations — Casanovo writes none. An empty value therefore means the peptide is unmodified, but a populated one is mzLib's interpretation of a mass, not the search engine's own call.
# Casanovo writes no target/decoy label, so is_decoy is NA rather than a false FALSE.
matches$records[, c("one_based_scan_number", "base_sequence", "is_decoy")]
#>   one_based_scan_number   base_sequence is_decoy
#> 1                     1      AGAHLQGGAK       NA
#> 2                     3       RGTGVENVK       NA
#> 3                     4 LDDPKEEDEEKEEGK       NA
#> 4                     5     HQGVMVGMGQK       NA
#> 5                     7         RQEFEMK       NA
```
