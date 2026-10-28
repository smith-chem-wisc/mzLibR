# Read many spectral libraries into one long table, in one bridge call

Read an mzLib spectral library (.msl) as one row per precursor, or with
fragments one row per fragment ion, with mzLib's not-available sentinels
as nulls.

The many-files form of
[`library_read`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md):
a predicted library and the empirical one refined from it, say. ONE
bridge process reads every library, `threads` at a time, into ONE table,
in input order; each library's facts (counts, load mode, format version)
are its row of `files`.

## Usage

``` r
library_read_many(paths, fragments = FALSE, index_only_threshold = NULL, out = NULL,
  threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- paths:

  A character vector of `.msl` paths. Read in this order; a file's
  position is its `source_index`. A repeated path is refused rather than
  read twice.

- fragments:

  One row per fragment ion, in every library. `FALSE` by default.

- index_only_threshold:

  The size in **bytes** at and above which a library opens index-only,
  for every library. `NULL` (the default) leaves it to mzLib (1 GiB in
  the pinned build).

- out:

  Write the long table here as tab-separated text, one library at a
  time, in order, and return only a summary. A batch that stops on an
  error removes its partial table. The path must end in `.tsv` (any
  case).

- threads:

  Libraries read at once, a count of files: `1` (the default), more, or
  `-1` for one per core. The result is **byte-identical at any value**,
  so this trades memory for speed and never changes an answer. The
  default is 1 because a library loaded whole is held in memory while
  its rows are written.

- on_error:

  `"fail"` (the default) stops at the first library that cannot be read
  and raises, naming it; `"skip"` records the failure in that library's
  row of `files` and reads the rest.

- timeout:

  Seconds to allow for the whole batch, or `NULL` (the default) to wait.

## Value

An `mzlibr_library_batch` (an `mzlibr_batch`). `records` begins with
`source_index` (1-based, into `paths`) and `source_path`, then the
columns of
[`library_read`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md),
rows grouped by library in input order, each library's in
`precursor_index` order (1-based) - **whatever** `threads` **was**.
Their units are as there: `precursor_mz` and `fragment_mz` in m/z,
`charge_state` and `fragment_charge` a charge, `retention_time` in iRT
or minutes per row (read `retention_time_unit`), `ion_mobility` in 1/K0,
`q_value` a fraction, `nce` an NCE, `fragment_count` in fragment ions,
`fragment_intensity` a relative intensity (0 to 1) and `neutral_loss` in
Da. `file_count`, `read_count` and `failed_count` count files;
`record_count` counts precursors summed over the libraries read,
`returned_count` the precursors returned and `row_count` the rows (one
per fragment ion with `fragments`), both 0 with `out`.
`fragments_included` and `index_only_threshold` (in bytes) are as for
[`library_read`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md).

`files` has one row per input: the facts
[`library_read`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md)
reports at its top level, under the same names - `reader`,
`format_version`, `is_compressed`, `is_index_only`, `target_count`,
`decoy_count`, `proteoform_count`, `min_precursor_mz`,
`max_precursor_mz` (m/z), `record_count` (precursors), `caveats`,
`absent_fields` - and `error_kind`, `error_type` and `error_message`,
`NA` for a library that was read. \*\*For a library that could not be
read, the facts it could not establish are `NA`, not zero\*\*: nothing
was counted.

## Wraps

Wire verb `library read` with `--paths-stdin`. Generated from the
bridge's verb spec `library.read.yaml` (bridge commit `60b1d25968a1`) by
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

- `paths` (wire `--paths-stdin`):

  flag; default absent. Read MANY libraries: newline-delimited paths on
  stdin (blank lines ignored; a repeated path is a usage error).
  Exclusive with path; offset and limit are refused with it.

- `threads`:

  int; in **files**; default `1`; range `>= 1, or -1 for one per core`.
  Libraries read at once. The output is byte-identical at any value
  (libraries are emitted in input order); 1 is the default because a
  library loaded whole is held in memory while its rows are written.

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. fail: the first
  unreadable library stops the batch with its usage or correctness
  error, the message starting 'Input \<i\> (\<path\>)'. skip: the
  failure goes in that input's files\[i\].error and the rest are read.
  skip with a single path is a usage error.

- `out`:

  path; default absent; range `.tsv only, any case (PYB-1)`. Write the
  long table here, one library at a time in input order. A batch that
  stops on an error deletes its partial table. Must differ from every
  input. Must end in .tsv, any case (PYB-1); any other extension, or
  none, is a usage error before any input is read, and no extension is
  appended.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **files**; never `NA`. Paths given.

- `read_count`:

  int; in **files**; never `NA`. Libraries read.

- `failed_count`:

  int; in **files**; never `NA`. Libraries not read; non-zero only under
  on-error skip.

- `record_count`:

  int; in **precursors**; never `NA`. Precursors summed over the
  libraries read (BULK.md section 2).

- `returned_count`:

  int; in **precursors**; never `NA`. Precursors returned in columns,
  summed over the libraries; 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns (one per fragment ion
  with fragments); 0 when written to out.

- `on_error`:

  string; never `NA`. fail or skip, as requested. `--threads` is
  deliberately not echoed: the envelope is byte-identical at every
  thread count.

- `fragments_included`:

  bool; never `NA`. Whether the rows are fragment ions.

- `index_only_threshold`:

  int; in **bytes**; never `NA`. The threshold applied to every library:
  the one given, or mzLib's default.

- `column_names`:

  string\[\]; never `NA`. source_index, source_path, then the
  single-library columns.

- `records` (wire `columns`):

  table; `NA` when the table was written to out instead. One long table:
  rows grouped by input in input order, each library's rows in
  precursor_index order, whatever `--threads` is.

- `output`:

  object; `NA` when out was not given. {path, format: "tsv", row_count};
  written one library at a time.

- `files`:

  object\[\]; never `NA`. One per input, in input order: the per-file
  block (path, file_type, reader, format_version, is_compressed,
  is_index_only, target_count, decoy_count, proteoform_count,
  min_precursor_mz, max_precursor_mz, record_count, rows_not_read,
  retention_time_unit, caveats, column_names, absent_fields,
  failed_fields, excluded_fields, error). An unread input has error set
  and every fact it could not establish null or empty.

## Columns of `records`

The first two columns say which input each row came from -
`source_index` is 1-based in R, so it indexes `files` directly - and
rows are grouped by input, in input order, whatever `threads` is.

- `source_index`:

  int; never `NA`. 0-based position of the row's library in the stdin
  list.

- `source_path`:

  string; never `NA`. Absolute path of the row's library.

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
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.library.read`; many files:
  `pymzlib.library.read_many`

- Rust (mzLibRust): `mzlib::library::read_with` with
  `LibraryReadOptions`; many files: `mzlib::library::read_many`

- R (mzLibR): `library_read`; many files: `library_read_many`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

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

[`library_read`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md)

## Examples

``` r

batch <- library_read_many(c("synthetic_library.msl", "myPrositLib.msl"))
batch
#> <mzlibr_library_batch> library read
#>   2 inputs: 2 read, 0 failed (on_error = "fail")
#>   11 rows in `records`, 11 records
batch$files[, c("record_count", "target_count", "decoy_count", "format_version")]
#>   record_count target_count decoy_count format_version
#> 1            6            4           2              3
#> 2            5            5           0              3
batch$files$absent_fields[[2]]
#> [1] "ion_mobility"      "protein_accession" "protein_name"     
#> [4] "gene_name"        
table(batch$records$source_index)
#> 
#> 1 2 
#> 6 5 
```
