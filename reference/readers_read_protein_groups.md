# Read a MetaMorpheus protein-group table, one row per group per sample group

Read a MetaMorpheus AllQuantifiedProteinGroups.tsv as one row per
protein group per sample group, with each sample group's spectral count
and intensity as columns.

MetaMorpheus writes `AllQuantifiedProteinGroups.tsv` wide: one intensity
column and one spectral-count column per sample group. This reads it
long - one row per protein group per sample group - so the per-sample
values are ordinary columns and nothing arrives as a map.

## Usage

``` r
readers_read_protein_groups(path, limit = NULL, offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to a MetaMorpheus `AllQuantifiedProteinGroups.tsv`
  ([`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
  reports `MetaMorpheusQuantifiedProteinGroups`).

- limit:

  Maximum protein groups to return - groups, not rows; each gives one
  row per sample group. `NULL`, the default, returns all of them.

- offset:

  Protein groups to skip - groups, not rows. A window, not a cursor:
  mzLib parses the whole file on every call.

- out:

  Write the table to this path as tab-separated text and return only a
  summary.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**The table is unfiltered**, as MetaMorpheus writes it: decoys,
contaminants and groups above 1% FDR are all rows. Filter on `q_value`
and `decoy_contaminant_target` before you count anything.

## Value

An `mzlibr_protein_group_records`. `record_count` counts the protein
groups in the whole file; `returned_count` the groups returned, starting
`offset` groups in; `row_count` the rows in `records`, which is groups
times sample groups. `sample_labels` are the file's sample-group labels
in header order, verbatim.

`records` has one row per group per sample group: `q_value` a fraction
from 0 to 1; `spectral_count` in PSMs; and `intensity` in the
instrument's intensity units, `NA` where the cell was blank - not
quantified in that sample group, which is not zero. `gene` and
`organism` are `NA` where the file has no such column, and
`absent_fields` then names them. `rows_not_read` is `NA`: every line is
a group, or the file fails to read.

## Sample labels are verbatim

A label such as `QE-002106_GM1_a-calib` is what MetaMorpheus wrote in
the header. Condition, replicate and channel cannot be recovered from it
reliably, so none is parsed out: map the labels to your design yourself.

## The occupancy columns are elsewhere

MetaMorpheus's per-sample PTM occupancy cells have no column shape and
are listed in `excluded_fields`, whose `verb` column names
`readers read-occupancy`: use
[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md),
which reads them one row per site.

## Wraps

Wire verb `readers read-protein-groups`. Generated from the bridge's
verb spec `readers.read-protein-groups.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`ProteinGroupFromTsvFile.LoadResults`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/ResultFiles/ProteinGroupFromTsvFile.cs)
  in
  `mzLib/Readers/InternalResults/ResultFiles/ProteinGroupFromTsvFile.cs`
  at mzLib `23c2490e`

- [`ProteinGroupFromTsv.SampleGroups`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs)
  in
  `mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsv.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A MetaMorpheus AllQuantifiedProteinGroups.tsv (readers
  identify reports MetaMorpheusQuantifiedProteinGroups).

- `offset`:

  int; in **groups**; default `0`; range `>= 0`. Skip this many protein
  groups (each gives one row per sample group). A window, not a cursor:
  mzLib parses the whole file on every call, so paging re-reads it.

- `limit`:

  int; in **groups**; default absent; range `>= 0`. Return at most this
  many protein groups (each gives one row per sample group). Absent
  means no limit; there is deliberately no default cap.

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

  string\[\]; never `NA`. The file's sample-group labels in header
  order, verbatim (e.g. QE-002106_GM1_a-calib). Condition and replicate
  cannot be recovered from them.

- `record_count`:

  int; in **groups**; never `NA`. Groups in the whole file, before the
  window.

- `rows_not_read`:

  int; in **rows**; `NA` when not counted for this table: one line is
  one group, and every group is read or the file fails. Data rows that
  did not become records.

- `retention_time_unit`:

  string; `NA` when the table has no time column. Always null for this
  table.

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

  int; in **groups**; never `NA`. Groups returned in columns (the unit
  offset and limit count in); 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns: more than
  returned_count in a long table, where one record gives several rows; 0
  when written to out.

- `offset`:

  int; in **groups**; never `NA`. The offset applied, in groups.

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

- `protein_group_name`:

  string; never `NA`. The members' accessions '\|'-joined in accession
  order (ProteinGroupFromTsv.ProteinGroupName); the join key to
  read-records.

- `gene`:

  string; `NA` when the file has no Gene column (then also in
  absent_fields) or the cell is blank. Gene name(s) as written.

- `organism`:

  string; `NA` when the file has no Organism column (then also in
  absent_fields) or the cell is blank. Organism as written.

- `decoy_contaminant_target`:

  string; never `NA`. T, D, C, or entrapment ET/ED, as written.

- `is_decoy`:

  bool; never `NA`. decoy_contaminant_target contains D (mzLib's
  IsDecoy).

- `is_contaminant`:

  bool; never `NA`. decoy_contaminant_target is C.

- `is_entrapment`:

  bool; never `NA`. decoy_contaminant_target starts with E.

- `q_value`:

  float; in **fraction (0 to 1)**; never `NA`. The group's q-value; the
  table is unfiltered and this is what filtering reads.

- `sample_label`:

  string; `NA` when the file has no per-sample columns: the group is one
  row and this is in absent_fields. The sample group this row measures.

- `spectral_count`:

  int; in **PSMs**; `NA` when the cell is blank, or the file has no
  SpectralCount\_ columns (then in absent_fields); an isobaric
  intensity-only label has none. PSMs from this sample group assigned to
  the group (SpectralCount\_).

- `intensity`:

  float; in **intensity (instrument units)**; `NA` when the cell was
  blank: not quantified in this sample group, which is not zero. The
  group's intensity in this sample group (Intensity\_).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`readers_read_protein_groups_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups_many.md).

- `mzlib_usage_error` (usage):

  path missing or blank; file or directory not found; mzLib does not
  recognise the file type (the message points at readers formats)

- `mzlib_usage_error` (usage):

  offset, limit or out given with no value; offset or limit \< 0;
  non-integer value; out resolves to the same file as path

- `mzlib_usage_error` (usage):

  the file is not a MetaMorpheus protein-group table (the message names
  what it is and points at read-records)

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

- The table is UNFILTERED as MetaMorpheus writes it: decoys,
  contaminants and groups above 1% FDR are all rows
  (ProteinGroupFromTsv.cs:15). Filter on q_value and
  decoy_contaminant_target before counting.

- Long format: one row per group per sample group. offset and limit
  count GROUPS, and returned_count counts rows.

- sample_label is the header label verbatim; condition, replicate and
  channel cannot be recovered from it (ProteinGroupFromTsv.cs:103).

- A blank intensity is null, not zero (ProteinGroupFromTsv.cs:114). An
  isobaric file has one counting label and one intensity label per
  channel, so a row can carry only spectral_count or only intensity.

- No member is a leading protein: members are sorted by accession
  (ProteinGroupFromTsv.cs:17).

- Occupancy cells are excluded_fields with verb readers read-occupancy;
  the group's other fields are in read-records, joined on
  protein_group_name.

- Many files: `--paths-stdin` reads the list in ONE process, `--threads`
  at a time, into one long table with source_index and source_path
  first. The output does not depend on `--threads` (a tested property);
  the default of 1 is a memory choice, since every reader materialises
  the whole file.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. For many files call
[`readers_read_protein_groups_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_protein_groups`; many files:
  `pymzlib.readers.read_protein_groups_many`

- Rust (mzLibRust): `mzlib::readers::read_protein_groups_with` with
  `ReadOptions`; many files: `mzlib::readers::read_protein_groups_many`

- R (mzLibR): `readers_read_protein_groups`; many files:
  [`readers_read_protein_groups_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups_many.md)

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- The group-level columns carried beside each sample (identity, q-value,
  decoy label) are a choice made here; coverage, masses and counts stay
  in read-records. Revisit if every binding ends up joining for the same
  field.

- Only the MetaMorpheus 1.1.11 fixture exists (label-free); no isobaric
  or 1.2+ (BioPolymer ...) header table has been read through this verb
  yet.

## See also

[`readers_read_protein_groups_many`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups_many.md),
[`readers_read_occupancy`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md),
[`readers_read_quantified_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)

## Examples

``` r

groups <- readers_read_protein_groups("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
  limit = 1)
groups
#> <mzlibr_protein_group_records> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv (MetaMorpheusQuantifiedProteinGroups)
#>   18 samples; 18 rows for 1 protein groups
#>   6 records in the file, 1 returned
#>   ! truncated - records were left behind
#>   ! The table is UNFILTERED, as MetaMorpheus writes it: decoys, contaminants and groups above 1% FDR are all rows (ProteinGroupFromTsv.cs:15). Filter on q_value and decoy_contaminant_target before counting or comparing groups.
#>   ! sample_label is the column label verbatim, e.g. QE-002106_GM1_a-calib. Condition, replicate and channel cannot be recovered from it (ProteinGroupFromTsv.cs:103); map labels to your design yourself, from an SDRF or your own sample sheet.
#>   ! intensity null means the cell was blank: the group was not quantified in that sample group, which is not zero (ProteinGroupFromTsv.cs:114). An isobaric file has one counting label and one intensity label per channel, so a row can carry only spectral_count or only intensity (ProteinGroupFromTsv.cs:104).
#>   ! No member is a leading protein: protein_group_name lists the members sorted by accession, so the first accession is only the one that sorts first (ProteinGroupFromTsv.cs:17).
#>   ! offset and limit count GROUPS; each group yields one row per sample group.
head(groups$records[, c("protein_group_name", "q_value", "sample_label", "intensity")])
#>   protein_group_name q_value          sample_label intensity
#> 1             P68363       0 QE-002106_GM1_a-calib   1838317
#> 2             P68363       0 QE-002107_GM1_b-calib   4550193
#> 3             P68363       0 QE-002108_GM1_c-calib   5226934
#> 4             P68363       0 QE-002110_GM2_a-calib   5452923
#> 5             P68363       0 QE-002111_GM2_b-calib  14198578
#> 6             P68363       0 QE-002112_GM2_c-calib        NA
```
