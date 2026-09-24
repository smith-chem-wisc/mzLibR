# Read a FlashLFQ peptide table, one row per peptide per sample

Read a FlashLFQ QuantifiedPeptides.tsv (or MetaMorpheus
AllQuantifiedPeptides.tsv) as one row per peptide per sample, with the
sample's intensity and detection type as columns.

Reads FlashLFQ's `QuantifiedPeptides.tsv`, or the
`AllQuantifiedPeptides.tsv` MetaMorpheus writes with it, long: one row
per peptide per sample, with that sample's intensity and detection type
as columns.

## Usage

``` r
readers_read_quantified_peptides(path, limit = NULL, offset = 0, out = NULL,
  timeout = NULL)
```

## Arguments

- path:

  Path to a FlashLFQ `QuantifiedPeptides.tsv` or MetaMorpheus
  `AllQuantifiedPeptides.tsv`
  ([`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
  reports `FlashLFQQuantifiedPeptide`).

- limit:

  Maximum peptides to return - peptides, not rows; each gives one row
  per sample. `NULL`, the default, returns all of them.

- offset:

  Peptides to skip - peptides, not rows.

- out:

  Write the table to this path as tab-separated text and return only a
  summary.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**An intensity of 0 is not a measurement.** FlashLFQ writes a literal 0
for a peptide it did not quantify in a sample, and mzLib reads it as 0;
`detection_type` - `"MSMS"`, `"MBR"`, `"NotDetected"` and so on - is
what tells the two apart.

## Value

An `mzlibr_quantified_peptide_records`. `record_count` counts the
peptides in the whole file; `returned_count` the peptides returned,
starting `offset` peptides in; `row_count` the rows in `records`,
peptides times samples. `sample_labels` are the file's sample labels in
header order.

`records` has one row per peptide per sample: `intensity` in the
instrument's intensity units (0 is not a measurement - see above);
`detection_type`; and `retention_time` in minutes, written only by
IsoTracker. For any other file `retention_time` and `peak_order` are
`NA` in every row and `absent_fields` names them. `rows_not_read` is
`NA`: every line is a peptide, or the file fails to read.

## Wraps

Wire verb `readers read-quantified-peptides`. Generated from the
bridge's verb spec `readers.read-quantified-peptides.yaml` (bridge
commit `5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these
facts, and all three bindings render the same ones.

- [`QuantifiedPeptideFile.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/ResultFiles/QuantifiedPeptideFile.cs)
  in
  `mzLib/Readers/InternalResults/ResultFiles/QuantifiedPeptideFile.cs`
  at mzLib `23c2490e`

- [`QuantifiedPeptideFromTsv.Samples`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/IndividualResultRecords/QuantifiedPeptideFromTsv.cs)
  in
  `mzLib/Readers/InternalResults/IndividualResultRecords/QuantifiedPeptideFromTsv.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A FlashLFQ QuantifiedPeptides.tsv or MetaMorpheus
  AllQuantifiedPeptides.tsv (readers identify reports
  FlashLFQQuantifiedPeptide).

- `offset`:

  int; in **peptides**; default `0`; range `>= 0`. Skip this many
  peptides (each gives one row per sample). A window, not a cursor:
  mzLib parses the whole file on every call, so paging re-reads it.

- `limit`:

  int; in **peptides**; default absent; range `>= 0`. Return at most
  this many peptides (each gives one row per sample). Absent means no
  limit; there is deliberately no default cap.

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

- `sample_labels`:

  string\[\]; never `NA`. The file's sample labels in header order,
  verbatim.

- `record_count`:

  int; in **peptides**; never `NA`. Peptides in the whole file, before
  the window.

- `rows_not_read`:

  int; in **rows**; `NA` when not counted for this table. Data rows that
  did not become records.

- `retention_time_unit`:

  string; never `NA`. Always "minutes": the per-sample header is
  RetentionTime (min)\_.

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

- `returned_count`:

  int; in **peptides**; never `NA`. Peptides returned in columns (the
  unit offset and limit count in); 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns: more than
  returned_count in a long table, where one record gives several rows; 0
  when written to out.

- `offset`:

  int; in **peptides**; never `NA`. The offset applied, in peptides.

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
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`readers_read_quantified_peptides_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides_many.md).

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
any work. For many files call
[`readers_read_quantified_peptides_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_quantified_peptides`; many
  files: `pymzlib.readers.read_quantified_peptides_many`

- Rust (mzLibRust): `mzlib::readers::read_quantified_peptides_with` with
  `ReadOptions`; many files:
  `mzlib::readers::read_quantified_peptides_many`

- R (mzLibR): `readers_read_quantified_peptides`; many files:
  [`readers_read_quantified_peptides_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides_many.md)

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

[`readers_read_quantified_peptides_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides_many.md),
[`readers_read_protein_groups`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md),
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)

## Examples

``` r

peptides <- readers_read_quantified_peptides("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv",
  limit = 1)
peptides$absent_fields   # this file has no retention times: they are NA, not missing rows
#> [1] "peak_order"     "retention_time"
head(peptides$records[, c("sequence", "sample_label", "intensity", "detection_type")])
#>                                                                                                             sequence
#> 1 [Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AC[Common Fixed:Carbamidomethyl on C]ASHVAK
#> 2 [Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AC[Common Fixed:Carbamidomethyl on C]ASHVAK
#> 3 [Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AC[Common Fixed:Carbamidomethyl on C]ASHVAK
#> 4 [Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AC[Common Fixed:Carbamidomethyl on C]ASHVAK
#> 5 [Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AC[Common Fixed:Carbamidomethyl on C]ASHVAK
#> 6 [Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AC[Common Fixed:Carbamidomethyl on C]ASHVAK
#>            sample_label intensity detection_type
#> 1 QE-002106_GM1_a-calib         0    NotDetected
#> 2 QE-002107_GM1_b-calib         0    NotDetected
#> 3 QE-002108_GM1_c-calib         0    NotDetected
#> 4 QE-002110_GM2_a-calib         0    NotDetected
#> 5 QE-002111_GM2_b-calib         0    NotDetected
#> 6 QE-002112_GM2_c-calib         0    NotDetected
```
