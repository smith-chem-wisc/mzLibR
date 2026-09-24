# Read a result file into the uniform record view

Read a result file through the cross-format quantifiable view: sequence,
retention time, charge, theoretical mass, decoy flag and protein groups,
the same columns for every format that offers it.

Only the three file types offering the `"quantifiable"` view can be read
this way — check
[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
first, or catch the error, which names the views the file does have.

## Usage

``` r
readers_read_results(path, limit = NULL, offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to a MetaMorpheus `.psmtsv` / `.osmtsv` or an MSFragger
  `psm.tsv`.

- limit:

  Maximum records to return. `NULL`, the default, returns all of them.

  **There is no default row limit**, deliberately. A result file can
  carry a million rows, and truncating by default would mean the
  ordinary call returns a table that looks complete and is not.
  `truncated` reports whether anything was left behind.

- offset:

  Records to skip.

  **This is a window, not a cursor.** mzLib materialises the whole file
  on every call — its readers look lazy and are not — so paging re-reads
  and re-parses the file once per page. A loop over pages is quadratic.
  For a large file use `out` instead, in one call.

- out:

  Write the records to this path as a **tab-separated** table and return
  only a summary, instead of carrying them back in the envelope. The
  intended path for large files.

  Tab-separated, not comma-separated, because these fields contain
  commas: MSFragger's mapped proteins are a comma-separated list inside
  a single field. Read it with `read.delim(path)`, not `read.csv(path)`.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely. A large file
  legitimately takes a while.

## Value

An `mzlibr_result_records`. `record_count` counts the records in the
**whole file** regardless of `limit` and `offset`; `returned_count`
counts the records that came back, starting `offset` records in;
`rows_not_read` counts data rows that did not become records.

`records` is a data.frame of the record view, one row per record, or
`NULL` when `out` was given. `retention_time` is in minutes for all four
formats at mzLib 1.0.592, and `retention_time_unit` says so per file;
`charge_state` is a charge; `monoisotopic_mass` is a neutral mass in Da.
`is_decoy` is `NA` where the format records no decoy label.

## Two fields to read before you trust the table

`rows_not_read` — data rows that did not become records. **mzLib drops a
malformed row silently**, so a non-zero value here means the file is
partly unreadable and your table is incomplete. `NA` when the count
could not be established.

`caveats` — what the uniform view cannot be trusted to mean for this
format, each citing the mzLib source it came from. Empty for some
formats and not for others. This is where you learn that, e.g., TopPIC
retention times are seconds while MetaMorpheus's and MSFragger's are
minutes.

Relatedly, `retention_time_unit` is per format and mzLib does **not**
normalise it. Convert with
[`readers_retention_time_in_minutes`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
rather than by hand.

## Wraps

Wire verb `readers read-results`. Generated from the bridge's verb spec
`readers.read-results.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`FileReader.ReadQuantifiableResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`IQuantifiableRecord`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/IQuantifiableRecord.cs)
  in `mzLib/Readers/BaseClasses/IQuantifiableRecord.cs` at mzLib
  `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A MetaMorpheus .psmtsv/.osmtsv, an MSFragger psm.tsv
  or a DIA-NN report.tsv.

- `offset`:

  int; in **records**; default `0`; range `>= 0`. Skip this many
  records. A window, not a cursor: mzLib parses the whole file on every
  call, so paging re-reads it.

- `limit`:

  int; in **records**; default absent; range `>= 0`. Return at most this
  many records. Absent means no limit; there is deliberately no default
  cap.

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

  int; in **records**; never `NA`. Records in the whole file, before the
  window.

- `rows_not_read`:

  int; in **rows**; `NA` when the count is not meaningful for this
  format (only the psmtsv family and MSFragger are one line per record).
  Data rows that did not become records: mzLib drops a malformed line
  silently.

- `retention_time_unit`:

  string; never `NA`. The unit of retention_time for THIS format:
  minutes for all four at mzLib 1.0.592, unknown when mzLib gives no
  basis.

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

  int; in **records**; never `NA`. Records returned in columns (the unit
  offset and limit count in); 0 when written to out.

- `offset`:

  int; in **records**; never `NA`. The offset applied, in records.

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

## On the wire but not projected yet

The bridge sends these, and this version of mzLibR does not return them
yet:

- `reader`:

  arrives with the mzLib 1.0.592 port.

- `absent_fields`:

  arrives with the mzLib 1.0.592 port.

- `failed_fields`:

  arrives with the mzLib 1.0.592 port.

- `excluded_fields`:

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
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_results`; many files:
  `pymzlib.readers.read_results_many`

- Rust (mzLibRust): `mzlib::readers::read_results_with` with
  `ReadOptions`; many files: `mzlib::readers::read_results_many`

- R (mzLibR): `readers_read_results`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

## See also

[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md),
[`readers_retention_time_in_minutes`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)

## Examples

``` r

psms <- readers_read_results("FraggerPsm_FragPipev21.1_psm.tsv", limit = 2)
psms
#> <mzlibr_result_records> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv (MsFraggerPsm)
#>   5 records in the file, 2 returned
#>   ! truncated - records were left behind
#>   retention_time_unit: minutes
#>   ! is_decoy is null for this format: MSFragger's psm.tsv carries no target/decoy column, so mzLib cannot report decoy status (MsFraggerPsm.cs:231) and the field crosses as null. Null means 'unknown', not 'target' - do not filter this format on is_decoy == false.
#>   ! monoisotopic_mass is the THEORETICAL peptide mass (MsFraggerPsm.cs:233, CalculatedPeptideMass), not the observed precursor mass. The psmtsv formats report the theoretical mass here too, so the two are consistent - but neither is what the instrument measured.
#>   ! file_name is the full 'Spectrum File' path including its .pep.xml extension, whereas the psmtsv formats report a bare base name. The field is not a join key across formats.
psms$records[, c("base_sequence", "charge_state", "retention_time")]
#>   base_sequence charge_state retention_time
#> 1       KPVGAAK            2     0.03233000
#> 2     KPAAAAGAK            2     0.04907667
readers_retention_time_in_minutes(psms)
#> [1] 0.03233000 0.04907667
```
