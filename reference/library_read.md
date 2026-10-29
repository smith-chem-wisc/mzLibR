# Read an mzLib spectral library (`.msl`): one row per precursor, or per fragment ion

Read an mzLib spectral library (.msl) as one row per precursor, or with
fragments one row per fragment ion, with mzLib's not-available sentinels
as nulls.

Answers **what is in this spectral library?** For each precursor - a
peptide or proteoform at one charge - the m/z, retention time and ion
mobility a search expects to see it at, its labels and provenance, and
on request the fragment ions it expects in the MS2 spectrum.
MetaMorpheus and mzLib write libraries in mzLib's binary `.msl` format;
this reads one with mzLib's own reader, `MslFileTypeHandler.Open` and
`MslLibrary.GetEntry`.

## Usage

``` r
library_read(path, fragments = FALSE, index_only_threshold = NULL, limit = NULL,
  offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to a `.msl` file. Any other extension is an `mzlib_usage_error`:
  text libraries (`.msp`) are not read here.

- fragments:

  One row per **fragment ion** instead of one per precursor, with the
  precursor's columns repeated on each row. `FALSE` by default: a
  library has tens of fragments per precursor, so pair it with `limit`
  or `out` for a large one.

- index_only_threshold:

  A library this many **bytes** or larger opens index-only; a smaller
  one loads whole. `NULL` (the default) sends nothing, so mzLib's own
  default applies: `MslFileTypeHandler.DefaultIndexOnlyThresholdBytes`,
  1073741824 bytes (1 GiB) in the pinned build. `0` opens every library
  index-only. The rows are the same either way.

- limit:

  Maximum **precursors** to return; with `fragments` each gives one row
  per fragment ion. `NULL`, the default, returns all of them.

- offset:

  **Precursors** to skip, in `precursor_index` order.

- out:

  Write a tab-separated table here and return only a summary. The path
  must end in `.tsv` (any case); any other is an `mzlib_usage_error`
  before a read. The file keeps the wire's 0-based `precursor_index`.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait.

## Details

Rows come in the library's stored order (`precursor_index`), not by m/z:
sort on `precursor_mz` for that. mzLib's "not available" values - an ion
mobility of 0.0, a NaN q-value, an NCE of 0, an empty protein
accession - arrive as `NA`.

Large libraries open **index-only**: precursors in memory, fragment
blocks left on disk. mzLib decides by file size (1 GiB by default in the
pinned build); pass `index_only_threshold` to move the line. The rows
are the same either way.

## Value

An `mzlibr_library_records`. `record_count` counts the precursors in the
whole library, before `offset` and `limit`; `returned_count` the
precursors returned (precursors, not rows: the unit `limit` and `offset`
count in), 0 with `out`; `row_count` the rows of `records`, one per
fragment ion with `fragments`; `offset` the precursors skipped;
`truncated` whether precursors were left behind.

The library's facts: `reader` (`"MslLibrary"`); `format_version`, the
file's `.msl` format version; `is_compressed` (a compressed library
always loads whole) and `is_index_only`; `target_count`, precursors that
are not decoys, entrapment targets included; `decoy_count`, decoy
precursors, entrapment decoys included; `proteoform_count`, precursors
whose `molecule_type` is `"Proteoform"`; and `min_precursor_mz` and
`max_precursor_mz`, in m/z, `NA` when the library has no precursors.
`fragments_included` says whether the rows are fragment ions, and
`index_only_threshold` the threshold applied, in bytes: yours, or
mzLib's default.

`file_type` is always `NA` (`.msl` is not one of mzLib's
`SupportedFileType` members); `rows_not_read` is always `NA`: rows are
not counted for a binary library, because mzLib checks the whole file
before returning any precursor; `retention_time_unit` is `NA` at the top
level, because the unit varies by precursor - it is a column.
`absent_fields` names the columns this library has no source for, so
every value is `NA` (`ion_mobility` when no precursor has one,
`protein_accession` and `protein_name` without a protein table,
`gene_name` without gene names); `failed_fields` the columns whose read
threw on some rows; `excluded_fields` is always empty; `caveats` says
what the table cannot be trusted to mean. `output` says where the table
was written, or is `NULL`.

`records` is a data.frame, or `NULL` with `out`. One row per precursor:
`precursor_index` (**1-based in R**, the position in the library and the
row order), `full_sequence` (mzLib bracket notation), `base_sequence`,
`precursor_mz` in m/z, `charge_state` (a charge), `retention_time` in
iRT or minutes **per row** and `retention_time_unit` (`"iRT"` or
`"minutes"`), `ion_mobility` in 1/K0 (V s/cm2; `NA` where mzLib stores
0.0), `is_decoy`, `is_entrapment`, `decoy_contaminant_target` (`"T"`,
`"D"`, `"ET"` or `"ED"`), `protein_accession`, `protein_name` and
`gene_name` (`NA` where mzLib stores an empty string), `is_proteotypic`,
`q_value` (a fraction, 0 to 1; `NA` where mzLib stores NaN, a predicted
library), `elution_group_id`, `molecule_type`, `source_type`
(`"Predicted"`, `"Empirical"` or `"EmpiricalRefined"`),
`dissociation_type` and `nce` (an NCE, a whole number; `NA` where mzLib
stores 0).

With `fragments = TRUE`, one row per fragment ion: those columns, then
`fragment_count` (fragment ions of the precursor), `fragment_mz` in m/z,
`fragment_intensity` (relative intensity within the precursor, 0 to 1),
`product_type`, `secondary_product_type` and `secondary_fragment_number`
(`NA` for a terminal ion), `fragment_number`, `residue_position` (`NA`
for a custom-loss fragment, whose stored position mzLib reads back as
0), `fragment_charge` (a charge), `neutral_loss` in Da (0 is no loss),
`is_internal`, `is_diagnostic`, `exclude_from_quant` and `annotation`. A
precursor with no fragment ions is one row with those cells `NA`.

## Three things to know before you use the numbers

**Retention time has a unit per row.** A library stores either an
indexed retention time (iRT) or run minutes for each precursor
(`MslLibraryEntry.RtIsCalibrated`), and one library can hold both. Read
`retention_time_unit` on each row; there is no conversion between them
without a calibration against a run, which
[`library_calibrate_rt`](https://smith-chem-wisc.github.io/mzLibR/reference/library_calibrate_rt.md)
fits.

\*\*"Not available" is `NA`, not a number.\*\* A text library converted
to `.msl` (from `.msp`) has no NCE, dissociation type, q-value, ion
mobility or protein, so those are `NA` (and `dissociation_type` is
`"Unknown"`).

**The library is unfiltered.** Decoys (`D`, `ED`) and entrapment
precursors (`ET`, `ED`) are rows; filter on `is_decoy` or
`decoy_contaminant_target` before counting targets.

## Wraps

Wire verb `library read`. Generated from the bridge's verb spec
`library.read.yaml` (bridge commit `aea14604d5f7`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`MslFileTypeHandler.Open`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/SpectralLibrary/Mslfiletypehandler.cs)
  in `mzLib/Readers/SpectralLibrary/Mslfiletypehandler.cs` at mzLib
  `d10017ce`

- [`MslFileTypeHandler.IsMslFile`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/SpectralLibrary/Mslfiletypehandler.cs)
  in `mzLib/Readers/SpectralLibrary/Mslfiletypehandler.cs` at mzLib
  `d10017ce`

- [`MslLibrary.Load`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/SpectralLibrary/MslLibrary.cs)
  in `mzLib/Readers/SpectralLibrary/MslLibrary.cs` at mzLib `d10017ce`

- [`MslLibrary.LoadIndexOnly`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/SpectralLibrary/MslLibrary.cs)
  in `mzLib/Readers/SpectralLibrary/MslLibrary.cs` at mzLib `d10017ce`

- [`MslLibrary.GetEntry(int, bool)`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Readers/SpectralLibrary/MslLibrary.cs)
  in `mzLib/Readers/SpectralLibrary/MslLibrary.cs` at mzLib `d10017ce`

- [`MslLibraryEntry`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SpectralLibrary/MslFragmentIon.cs)
  in `mzLib/Omics/SpectralLibrary/MslFragmentIon.cs` at mzLib `d10017ce`

- [`MslFragmentIon`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SpectralLibrary/MslFragmentIon.cs)
  in `mzLib/Omics/SpectralLibrary/MslFragmentIon.cs` at mzLib `d10017ce`

- [`DecoyContaminantTargetLabel.For(bool, bool, bool)`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/BioPolymer/DecoyContaminantTargetLabel.cs)
  in `mzLib/Omics/BioPolymer/DecoyContaminantTargetLabel.cs` at mzLib
  `d10017ce`

## Parameters: units, ranges and defaults

- `path`:

  path; required; range `.msl, any case (MslFileTypeHandler.IsMslFile)`.
  An mzLib binary spectral library. Any other extension is a usage error
  before the file is opened: text libraries (.msp) are not read by this
  verb.

- `fragments`:

  flag; default `FALSE`. One row per fragment ion instead of one per
  precursor: the precursor's columns repeated on each row, then the
  fragment's. A precursor with no fragment ions gives one row whose
  fragment cells are null. Given with a value is a usage error.

- `index_only_threshold` (wire `--index-only-threshold`):

  int; in **bytes**; default `1073741824`; range `>= 0`. Passed to
  MslFileTypeHandler.Open: a library this size or larger opens
  index-only (precursors in memory, fragment blocks read from disk on
  demand), a smaller one is loaded whole. Absent means mzLib's own
  default, MslFileTypeHandler.DefaultIndexOnlyThresholdBytes, 1073741824
  bytes (1 GiB) at the pin; the bridge sends no value of its own. 0
  opens every library index-only. The rows are the same in both modes;
  only memory and speed differ.

- `offset`:

  int; in **precursors**; default `0`; range `>= 0`. Skip this many
  precursors, in precursor_index order. A window, not a cursor: the
  library is opened again on every call.

- `limit`:

  int; in **precursors**; default absent; range `>= 0`. Return at most
  this many precursors; with fragments each gives one row per fragment
  ion. Absent means no limit.

- `out`:

  path; default absent; range `.tsv only, any case (PYB-1)`. Write the
  table to this path as tab-separated text instead of returning columns.
  Must differ from path. Parent directories are created. Must end in
  .tsv, any case (PYB-1); any other extension, or none, is a usage error
  before the library is opened, and no extension is appended.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. Absolute path of the library.

- `file_type`:

  string; `NA` when always: .msl is not an mzLib SupportedFileType
  (mzLib reads it with MslLibrary, outside FileReader), so there is no
  type name to give. The mzLib SupportedFileType, for parity with the
  reader verbs.

- `reader`:

  string; never `NA`. Always MslLibrary, the mzLib class that read the
  file.

- `format_version`:

  int; never `NA`. The file's .msl format version
  (MslFileHeader.FormatVersion): 1 to 3 at the pin, where 3 allows
  compression (MslFormat.CurrentVersion).

- `is_compressed`:

  bool; never `NA`. The fragment section is zstd-compressed
  (MslLibrary.IsCompressed). A compressed library always loads whole.

- `is_index_only`:

  bool; never `NA`. mzLib opened the library index-only
  (MslLibrary.IsIndexOnly): at or above index_only_threshold, and not
  compressed.

- `target_count`:

  int; in **precursors**; never `NA`. Precursors that are not decoys,
  entrapment targets included (MslLibrary.TargetCount).

- `decoy_count`:

  int; in **precursors**; never `NA`. Decoy precursors, entrapment
  decoys included (MslLibrary.DecoyCount).

- `proteoform_count`:

  int; in **precursors**; never `NA`. Precursors whose molecule_type is
  Proteoform (MslLibrary.ProteoformCount).

- `min_precursor_mz`:

  float; in **m/z**; `NA` when the library has no precursors. Smallest
  precursor m/z (MslLibrary.MinPrecursorMz).

- `max_precursor_mz`:

  float; in **m/z**; `NA` when the library has no precursors. Largest
  precursor m/z (MslLibrary.MaxPrecursorMz).

- `record_count`:

  int; in **precursors**; never `NA`. Precursors in the whole library,
  before the window (MslLibrary.PrecursorCount).

- `rows_not_read`:

  int; in **rows**; `NA` when always: the file is binary and mzLib
  validates it whole (magic, version, footer, precursor count, CRC-32)
  before returning any entry, so every precursor is read or the call
  fails. Data rows that did not become records.

- `retention_time_unit`:

  string; `NA` when always: the unit varies by precursor, so it is the
  retention_time_unit column. Always null for this verb.

- `caveats`:

  string\[\]; never `NA`. What this table cannot be trusted to mean,
  each citing the mzLib source.

- `column_names`:

  string\[\]; never `NA`. Column order.

- `absent_fields`:

  string\[\]; never `NA`. Columns this library has no source for, from
  its file flags (MslFormat.FileFlag\*): ion_mobility when no precursor
  has one, protein_accession and protein_name when there is no protein
  table, gene_name when there are no gene names. Null in every row
  (BULK.md section 4).

- `failed_fields`:

  string\[\]; never `NA`. 'field: ExceptionType' for each column whose
  read threw on some returned rows; those cells are null.

- `excluded_fields`:

  object\[\]; never `NA`. Always empty: every field crosses, the
  fragments with the fragments flag.

- `error`:

  object; `NA` when always, for a single path: a library that cannot be
  read fails the call instead. {kind, type, message}; non-null only in a
  files\[\] entry under on-error skip.

- `fragments_included`:

  bool; never `NA`. Whether the rows are fragment ions (the fragments
  flag).

- `index_only_threshold`:

  int; in **bytes**; never `NA`. The threshold MslFileTypeHandler.Open
  applied: the one given, or mzLib's default, 1073741824 at the pin.

- `returned_count`:

  int; in **precursors**; never `NA`. Precursors returned in columns
  (the unit offset and limit count in); 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns: returned_count without
  fragments, one per fragment ion (or one for a precursor with none)
  with them; 0 when written to out.

- `offset`:

  int; in **precursors**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever precursors were left out by offset or
  limit.

- `records` (wire `columns`):

  table; `NA` when the table was written to out instead. Column name to
  per-row values, in column_names order.

- `output`:

  object; `NA` when out was not given. {path, format: "tsv", row_count}
  when out was given.

## Columns of `records`

- `precursor_index`:

  int; never `NA`. 0-based position of the precursor in the library
  (MslPrecursorIndexEntry.PrecursorIdx): the key MslLibrary.GetEntry
  takes, and the row order. Not m/z order.

- `full_sequence`:

  string; never `NA`. Modified sequence in mzLib bracket notation, e.g.
  PEPTM\[Common Variable:Oxidation on M\]IDE
  (MslLibraryEntry.FullSequence).

- `base_sequence`:

  string; never `NA`. The bare residues (MslLibraryEntry.BaseSequence).

- `precursor_mz`:

  float; in **m/z**; never `NA`. Precursor m/z, (M + zH)/z
  (MslLibraryEntry.PrecursorMz), stored as a 32-bit float.

- `charge_state`:

  int; in **charge**; never `NA`. Precursor charge
  (MslLibraryEntry.ChargeState).

- `retention_time`:

  float; in **iRT or min (per row)**; never `NA`. Indexed retention
  time, or a run-calibrated time in minutes, as retention_time_unit says
  (MslLibraryEntry.RetentionTime). iRT can be negative. 0 is passed
  through: mzLib documents no sentinel for it.

- `retention_time_unit`:

  string; never `NA`. iRT or minutes, from
  MslLibraryEntry.RtIsCalibrated (false: iRT; true: run-specific
  minutes). Per row because one library can hold both.

- `ion_mobility`:

  float; in **1/K0 (V s/cm2)**; `NA` when mzLib's 0.0, which
  MslLibraryEntry.IonMobility documents as not available; in
  absent_fields when no precursor has one. Ion mobility as 1/K0
  (MslLibraryEntry.IonMobility).

- `is_decoy`:

  bool; never `NA`. A decoy precursor, entrapment decoys included
  (MslLibraryEntry.IsDecoy).

- `is_entrapment`:

  bool; never `NA`. An entrapment precursor, target or decoy
  (MslLibraryEntry.IsEntrapment, mzLib \#1383). Independent of is_decoy.

- `decoy_contaminant_target`:

  string; never `NA`. T, D, ET or ED, from
  DecoyContaminantTargetLabel.For(is_decoy, false, is_entrapment). The
  .msl format has no contaminant flag, so C never appears.

- `protein_accession`:

  string; `NA` when mzLib's empty string: the precursor has no protein
  (MslLibraryEntry.ProteinAccession); in absent_fields when the library
  has no protein table. Accession of the source protein, as stored, e.g.
  P04637.

- `protein_name`:

  string; `NA` when mzLib's empty string: no protein, or no name for it
  (MslLibraryEntry.ProteinName); in absent_fields when the library has
  no protein table. Protein name, as stored.

- `gene_name`:

  string; `NA` when mzLib's empty string: no protein, or no gene for it
  (MslLibraryEntry.GeneName); in absent_fields when the library has no
  gene names. Gene symbol, as stored, e.g. TP53.

- `is_proteotypic`:

  bool; never `NA`. Expected to give a uniquely detectable signal
  (MslLibraryEntry.IsProteotypic), as the library's writer set it.

- `q_value`:

  float; in **fraction (0 to 1)**; `NA` when mzLib's NaN: no q-value is
  available, as for a purely predicted library (MslLibraryEntry.QValue).
  The library's q-value for the precursor; lower is better.

- `elution_group_id`:

  int; never `NA`. Shared by every precursor of one base sequence, so
  charge states and modified forms co-elute
  (MslLibraryEntry.ElutionGroupId); mzLib's writer assigns it.

- `molecule_type`:

  string; never `NA`. MslFormat.MoleculeType name: Peptide, Proteoform,
  Oligonucleotide or Glycopeptide.

- `source_type`:

  string; never `NA`. MslFormat.SourceType name
  (MslLibraryEntry.Source): Predicted, Empirical, or EmpiricalRefined
  (predicted intensities, measured retention time or ion mobility).

- `dissociation_type`:

  string; never `NA`. mzLib DissociationType name, e.g. HCD; Unknown
  when the library does not say (MslLibraryEntry.DissociationType).

- `nce`:

  int; in **NCE**; `NA` when mzLib's 0: unknown or not applicable
  (MslLibraryEntry.Nce). Normalized collision energy, a whole number:
  the format stores NCE x 10 and mzLib divides by 10 as an integer.

- `fragment_count`:

  int; in **fragment ions**; never `NA`. Fragment ions the precursor
  has, counted from the ions mzLib read; 0 for a precursor with none,
  whose one row has every fragment cell null. *Present only with
  `__fragments`.*

- `fragment_mz`:

  float; in **m/z**; `NA` when the precursor has no fragment ions
  (fragment_count 0). Fragment m/z (MslFragmentIon.Mz), a 32-bit float.
  Rows are in the stored order, which mzLib's writer sorts by m/z.
  *Present only with `__fragments`.*

- `fragment_intensity`:

  float; in **relative intensity (0 to 1)**; `NA` when the precursor has
  no fragment ions. Intensity relative to the precursor's most intense
  fragment, which mzLib's writer scales to 1 (MslFragmentIon.Intensity).
  Not comparable across precursors. *Present only with `__fragments`.*

- `product_type`:

  string; `NA` when the precursor has no fragment ions. mzLib
  ProductType name, e.g. b, y, c, zDot; for an internal ion, its
  N-terminal type (MslFragmentIon.ProductType). *Present only with
  `__fragments`.*

- `secondary_product_type`:

  string; `NA` when a terminal ion (MslFragmentIon.SecondaryProductType
  is null for every terminal ion), or no fragment ions. An internal
  ion's C-terminal ProductType, e.g. y for bIy. *Present only with
  `__fragments`.*

- `fragment_number`:

  int; `NA` when the precursor has no fragment ions. The ion's number in
  its series, e.g. 5 for b5; for an internal ion, its 0-based start
  residue (MslFragmentIon.FragmentNumber). *Present only with
  `__fragments`.*

- `secondary_fragment_number`:

  int; `NA` when a terminal ion, for which mzLib stores 0
  (MslFragmentIon.SecondaryFragmentNumber), or no fragment ions. An
  internal ion's end residue, exclusive. *Present only with
  `__fragments`.*

- `residue_position`:

  int; `NA` when mzLib's 0: always for a fragment with a custom neutral
  loss, whose stored position the format reuses to index its loss table
  (MslReader.ConvertFragment); or no fragment ions. Residue position the
  ion was annotated from, as stored (MslFragmentIon.ResiduePosition):
  for a library mzLib built from its own products, the
  Product.ResiduePosition convention, fragment_number for N-terminal
  ions and length - fragment_number + 1 for C-terminal ones (mzLib
  \#1408). *Present only with `__fragments`.*

- `fragment_charge`:

  int; in **charge**; `NA` when the precursor has no fragment ions.
  Fragment charge (MslFragmentIon.Charge). *Present only with
  `__fragments`.*

- `neutral_loss`:

  float; in **Da**; `NA` when the precursor has no fragment ions.
  Neutral-loss mass, negative for a loss (MslFragmentIon.NeutralLoss). 0
  means no loss, and crosses as 0, not null: it is the value of the
  shift, not a missing one. *Present only with `__fragments`.*

- `is_internal`:

  bool; `NA` when the precursor has no fragment ions. An internal ion
  (MslFragmentIon.IsInternalFragment). *Present only with
  `__fragments`.*

- `is_diagnostic`:

  bool; `NA` when the precursor has no fragment ions. A diagnostic ion,
  ProductType D (MslFragmentIon.IsDiagnosticIon). *Present only with
  `__fragments`.*

- `exclude_from_quant`:

  bool; `NA` when the precursor has no fragment ions. Marked to be left
  out of quantification, like DIA-NN's ExcludeFromAssay
  (MslFragmentIon.ExcludeFromQuant). *Present only with `__fragments`.*

- `annotation`:

  string; `NA` when the precursor has no fragment ions. mzLib's library
  rendering (MslFragmentIon.Annotation), e.g. b5, y3-18.0106, bIy\[2-5\]
  or y4^2+: the loss signed to four decimals and the charge appended,
  unlike Product.Annotation. *Present only with `__fragments`.*

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` comes
only from
[`library_read_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md).

- `mzlib_usage_error` (usage):

  path missing or blank; the file does not exist, or is a directory; its
  extension is not .msl (MslFileTypeHandler.IsMslFile)

- `mzlib_usage_error` (usage):

  fragments given a value; index-only-threshold, offset, limit or out
  given with no value; index-only-threshold not a whole number, or below
  0; offset or limit below 0 or not an integer; out resolving to the
  same file as path

- `mzlib_usage_error` (usage):

  out (with path or paths-stdin) whose extension is not .tsv, any case,
  including none (PYB-1); checked before any input is opened

- `mzlib_bridge_error` (correctness):

  mzLib cannot read the library: FormatException for a bad magic number,
  version, footer or precursor count, InvalidDataException for a CRC-32
  mismatch (MslReader); the error type is mzLib's exception name

- `mzlib_usage_error` (usage):

  paths-stdin with path, offset or limit; paths-stdin with a value; no
  paths on stdin; a repeated path; threads 0 or below -1; on-error not
  fail or skip; on-error skip with a single path; out naming an input

- `mzlib_usage_error` (usage):

  with on-error fail, an input that is missing or not a .msl: the usage
  error of that input, prefixed 'Input \<i\> (\<path\>)'

- `mzlib_bridge_error` (correctness):

  with on-error fail, mzLib failing to read an input: that input's own
  error type, message prefixed 'Input \<i\> (\<path\>)'

- `mzlib_usage_error` (usage):

  the answer is too large to return as one JSON document (PYB-3 part 1).
  Message: 'this read (N records) is too large to return as one
  document; use `--limit`/`--offset` or `--out`'. A library with
  fragments is the likeliest case: write it to out.

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Rows are in the library's stored order, by precursor_index, not by
  m/z: mzLib's own MslLibrary.GetAllEntries enumerates in ascending
  precursor m/z instead. Sort on precursor_mz for that order.

- retention_time is iRT for some precursors and run minutes for others:
  read retention_time_unit on each row (MslLibraryEntry.RtIsCalibrated).
  A library can hold both, so never plot or compare the column without
  it. There is no conversion: iRT becomes minutes only through a
  calibration against a run.

- The library is unfiltered: decoys (D, ED) and entrapment precursors
  (ET, ED) are rows. Filter on is_decoy, or on decoy_contaminant_target,
  before counting targets.

- A retention_time of 0 may mean the source had none:
  MslLibraryEntry.FromLibrarySpectrum, the path mzLib converts a text
  library through, writes RetentionTime ?? 0.0, and mzLib documents no
  sentinel for it, so 0 is passed through as a value.

- nce is a whole number: the format stores NCE x 10 and MslReader reads
  it back with integer division, so an NCE of 27.5 reads as 27.

- A library converted from a text library (.msp) has no NCE,
  dissociation type, q-value, ion mobility or protein: LibrarySpectrum
  has no field for them, and MslLibraryEntry.FromLibrarySpectrum fills
  mzLib's not-available defaults, so nce and q_value are null and
  dissociation_type is Unknown.

- precursor_mz, retention_time and ion_mobility are stored as 32-bit
  floats and widened to double by mzLib, so 465.74 reads back as
  465.739990234375. The extra digits are the float, not a measurement;
  round them for display.

- offset and limit count PRECURSORS. With fragments, each precursor
  gives one row per fragment ion.

- With fragments: residue_position is null for every fragment with a
  custom neutral loss, because the .msl format reuses that field to
  index its loss table and MslReader reads it back as 0.

- With fragments: fragment_intensity is relative within its precursor, 0
  to 1, because mzLib's writer scales each precursor's most intense
  fragment to 1 (MslWriter). It is not comparable across precursors.

- Index-only (at or above index_only_threshold) keeps memory to the
  precursors, but with fragments every fragment block is read twice from
  disk: once to count the rows, once to write them. mzLib keeps each
  precursor's stored fragment count private, so the bridge cannot size
  the table without reading them.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. For many files call
[`library_read_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md)
once rather than looping this function: one process, one start-up, and
the thread count stated on the wire.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.library.read`; many files:
  `pymzlib.library.read_many`

- Rust (mzLibRust): `mzlib::library::read_with` with
  `LibraryReadOptions`; many files: `mzlib::library::read_many`

- R (mzLibR): `library_read`; many files:
  [`library_read_many`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md)

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Text libraries (.msp, pDeep, ms2pip) could be read through the same
  table: mzLib's MslConverter.FromMspFile parses one into an in-memory
  MslLibrary. Left out of v1 because the conversion fills not-available
  defaults for every field LibrarySpectrum lacks, and because
  SpectralLibrary picks the text dialect from the file name. Revisit on
  demand.

- fragment_count is carried only with fragments: without them mzLib's
  index-only entries have an empty fragment list (MslLibrary.GetEntry(i,
  false), 'do not rely on MatchedFragmentIons being empty or full'), and
  the stored MslPrecursorRecord.FragmentCount is not public. A public
  count would give the precursor table this column at no I/O cost
  (PUSH-TO-MAINLAND).

- No decoy filter is offered: MslLibrary.GetAllEntries(includeDecoys:
  false) has one, but GetAllEntries yields m/z order with no index, so
  the bridge reads by index and the caller filters on is_decoy. A
  GetAllEntries overload with includeFragments, or one that yields the
  index, would let the bridge wrap it directly.

## See also

[`library_read_many`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md),
[`readers_read_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
for the runs a library is searched against

## Examples

``` r

lib <- library_read("synthetic_library.msl")
lib
#> <mzlibr_library_records> E:\CodeReview\wt-py-library-read\pkg\python\tests\fixtures\library\synthetic_library.msl: 6 precursors (4 target, 2 decoy), .msl format 3
#>   precursor m/z 406.242 to 952.493
#>   6 precursor(s) returned
#>   ! retention_time is in iRT and minutes: read retention_time_unit on each row
lib$records[, c("precursor_index", "full_sequence", "decoy_contaminant_target",
                "retention_time", "retention_time_unit")]
#>   precursor_index     full_sequence decoy_contaminant_target retention_time
#> 1               1          PEPTIDEK                        T           35.2
#> 2               2          KEDITPEP                        D           35.2
#> 3               3          PEPTLDEK                       ET           21.4
#> 4               4          KEDLTPEP                       ED           33.0
#> 5               5       ELVISLIVESK                        T           61.7
#> 6               6 MAHVRSGKHQIPAEALQ                        T           48.0
#>   retention_time_unit
#> 1                 iRT
#> 2                 iRT
#> 3             minutes
#> 4                 iRT
#> 5                 iRT
#> 6                 iRT
# A predicted precursor: mzLib's "not available" values are NA.
lib$records[5, c("ion_mobility", "q_value", "nce", "protein_accession")]
#>   ion_mobility q_value nce protein_accession
#> 5           NA      NA  NA              <NA>

# One precursor's fragment ions:
ions <- library_read("synthetic_library.msl", fragments = TRUE, limit = 1)
c(precursors = ions$returned_count, rows = ions$row_count)
#> precursors       rows 
#>          1          6 
ions$records[, c("annotation", "fragment_mz", "fragment_intensity", "neutral_loss")]
#>   annotation fragment_mz fragment_intensity neutral_loss
#> 1         y1    147.1128               0.40      0.00000
#> 2         b2    227.1026               0.30      0.00000
#> 3         y2    276.1554               1.00      0.00000
#> 4   bIy[2-5]    326.1710               0.10      0.00000
#> 5 y3-18.0106    373.2082               0.20    -18.01056
#> 6 y4-43.0058    461.2242               0.15    -43.00581
```
