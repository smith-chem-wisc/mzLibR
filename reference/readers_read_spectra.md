# Read the scans of a spectra file: headers always, peaks on request

Read the scans of a spectra file: every scan's header always, and its
peak arrays only on request.

Seven file types offer the `spectra` view: `.mzML`, `.mgf`,
`_ms1.msalign`, `_ms2.msalign`, Thermo `.raw`, Bruker `.d` and timsTOF
`.d`.

## Usage

``` r
readers_read_spectra(path, limit = NULL, offset = 0, ms_order = NULL, peaks = FALSE,
  out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to a spectra file. A Bruker `.d` directory is also accepted.

- limit:

  Maximum scans to return. `NULL`, the default, returns all of them.

- offset:

  Scans to skip, applied **after** `ms_order`.

- ms_order:

  Keep only scans at this MS level - `1` for survey scans, `2` for
  fragment scans. `NULL`, the default, keeps every scan. Applied before
  `offset` and `limit`, so `ms_order = 2, limit = 10` means the first
  ten MS2 scans rather than the MS2 scans among the first ten.

- peaks:

  Include the `mz` and `intensity` arrays. `FALSE` by default, and worth
  leaving so unless you need them: a scan header is tens of bytes and
  its peak list is thousands, and a mid-size mzML holds tens of
  thousands of scans. `peak_count` still reports how many peaks each
  scan has.

- out:

  Write a tab-separated table here and return only a summary. With
  `peaks = TRUE` each cell holds a `;`-joined list.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely. A large `.raw`
  legitimately takes a while.

## Details

Retention times here **are** in minutes for every format - mzLib's
spectra readers convert at the boundary, unlike its result-file readers,
which pass the tool's own unit through untouched.

## Value

An `mzlibr_scan_records`. `scan_count` is the file's total in scans
**before** any `ms_order` filter, reported alongside `record_count` -
the scans that passed it - so a filter that matched nothing can never
look like an empty file. `returned_count` is the scans that came back,
starting `offset` scans in.

`records` has one row per scan. `retention_time` is in minutes for every
format; `injection_time` in ms; `total_ion_current` and
`selected_ion_intensity` in the instrument's intensity units;
`peak_count` in peaks; `compensation_voltage` in volts;
`selected_ion_charge_state_guess` a charge; and `scan_window_lower_mz`,
`scan_window_upper_mz`, `isolation_mz`, `isolation_width`,
`selected_ion_mz` and `selected_ion_monoisotopic_guess_mz` in m/z. With
`peaks = TRUE`, `mz` (m/z) and `intensity` (instrument units) are list
columns, one vector per scan. A precursor field is `NA` on an MS1 scan;
the generated sections below say what `NA` means for every column.

## Two of the seven need Windows

Bruker `.d` and timsTOF `.d` are read through vendor native libraries
(`baf2sql`, `timsdata`) and are **Windows-x64 only**. Thermo `.raw` uses
managed vendor assemblies and works everywhere. msalign files hold
**deconvolved neutral masses**, not raw m/z - do not re-deconvolve them.

## Wraps

Wire verb `readers read-spectra`. Generated from the bridge's verb spec
`readers.read-spectra.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`MsDataFileReader.GetDataFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/MsDataFileReader.cs)
  in `mzLib/Readers/MsDataFileReader.cs` at mzLib `23c2490e`

- [`MsDataFile.GetAllScansList`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/MassSpectrometry/MsDataFile.cs)
  in `mzLib/MassSpectrometry/MsDataFile.cs` at mzLib `23c2490e`

- [`SourceFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/MassSpectrometry/SourceFile.cs)
  in `mzLib/MassSpectrometry/SourceFile.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A spectra file: mzML, MGF, ms1/ms2 msalign, Thermo
  .raw, Bruker .d or timsTOF .d.

- `ms_order` (wire `--ms-order`):

  int; default absent; range `>= 1`. Keep only scans at this MS level.
  Applied BEFORE offset and limit, so ms-order 2 with limit 10 is the
  first ten MS2 scans. Given with no value is a usage error; absent
  means every level.

- `offset`:

  int; in **scans**; default `0`; range `>= 0`. Skip this many scans
  (after the ms-order filter).

- `limit`:

  int; in **scans**; default absent; range `>= 0`. Return at most this
  many scans. Absent means no limit.

- `peaks`:

  flag; default `FALSE`. Include the mz and intensity arrays. Off by
  default because peaks are roughly the size of the file; combine with
  limit, a scan filter, or out.

- `out`:

  path; default absent. Write the table to this path as TSV instead of
  returning columns. Must differ from path. Peak arrays are written as
  ';'-joined lists.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. Absolute path of the input.

- `file_type`:

  string; never `NA`. mzLib SupportedFileType name.

- `reader`:

  string; never `NA`. The mzLib MsDataFile class that read it.

- `scan_count`:

  int; in **scans**; never `NA`. The file's true total, before any
  filter, so a filter matching nothing never looks like an empty file.

- `ms_order`:

  int; `NA` when no ms-order filter was applied. The ms-order filter,
  echoed.

- `record_count`:

  int; in **scans**; never `NA`. Scans that passed the ms-order filter.

- `returned_count`:

  int; in **scans**; never `NA`. Scans in columns; 0 when written to
  out.

- `offset`:

  int; in **scans**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever scans were left out by offset or
  limit.

- `peaks_included`:

  bool; never `NA`. Whether mz and intensity columns are present.

- `retention_time_unit`:

  string; never `NA`. Always "minutes" for this verb: every mzLib
  MsDataFile reader converts at the boundary.

- `caveats`:

  string\[\]; never `NA`. Format-specific traps for THIS file; see
  caveats below.

- `column_names`:

  string\[\]; never `NA`. Column order.

- `records` (wire `columns`):

  table; `NA` when the table was written to out instead. Column name to
  per-scan values.

- `output`:

  object; `NA` when out was not given. {path, format: "tsv", row_count}
  when out was given.

- `error`:

  object; `NA` when always, for a single path: a file that cannot be
  read fails the call instead. {kind, type, message}; non-null only in a
  files\[\] entry under on-error skip.

## Columns of `records`

- `one_based_scan_number`:

  int; never `NA`. Scan number. For MGF it may be assigned in file
  order; see caveats.

- `ms_order`:

  int; never `NA`. MS level (1, 2, ...).

- `retention_time`:

  float; in **min**; never `NA`. Retention time.

- `polarity`:

  string; never `NA`. mzLib Polarity enum name.

- `mz_analyzer`:

  string; never `NA`. mzLib MZAnalyzerType enum name.

- `is_centroid`:

  bool; never `NA`. Whether the spectrum is centroided.

- `total_ion_current`:

  float; in **intensity (instrument units)**; never `NA`. Total ion
  current.

- `injection_time`:

  float; in **ms**; `NA` when the file did not record it. Ion injection
  time.

- `peak_count`:

  int; in **peaks**; never `NA`. Peaks in the scan; 0 when it has no
  spectrum.

- `scan_window_lower_mz`:

  float; in **m/z**; `NA` when no scan window recorded. Lower
  scan-window bound. DERIVED for MGF (first observed peak).

- `scan_window_upper_mz`:

  float; in **m/z**; `NA` when no scan window recorded. Upper
  scan-window bound. DERIVED for MGF (last observed peak).

- `scan_filter`:

  string; `NA` when the file did not record one. Vendor scan filter
  string.

- `native_id`:

  string; `NA` when the file did not record one. Vendor native ID.

- `scan_description`:

  string; `NA` when the file did not record one. Free-text scan
  description.

- `one_based_precursor_scan_number`:

  int; `NA` when MS1 scan, or the file does not record the survey scan
  (standard MGF). Survey scan this scan was selected from.

- `isolation_mz`:

  float; in **m/z**; `NA` when MS1 scan, or not recorded. Isolation
  window centre. May be calibration-adjusted.

- `isolation_width`:

  float; in **m/z**; `NA` when MS1 scan, or not recorded. Isolation
  window width.

- `selected_ion_mz`:

  float; in **m/z**; `NA` when MS1 scan, or not recorded. Selected
  precursor m/z. May be calibration-adjusted.

- `selected_ion_intensity`:

  float; in **intensity (instrument units)**; `NA` when MS1 scan, or not
  recorded. Selected precursor intensity.

- `selected_ion_charge_state_guess`:

  int; in **charge**; `NA` when MS1 scan, or no charge assigned.
  Instrument or file charge-state guess.

- `selected_ion_monoisotopic_guess_mz`:

  float; in **m/z**; `NA` when MS1 scan, or no monoisotopic guess.
  Monoisotopic precursor m/z guess.

- `dissociation_type`:

  string; `NA` when MS1 scan, or not recorded. mzLib DissociationType
  enum name.

- `hcd_energy`:

  string; `NA` when not HCD, or not recorded. Collision energy AS
  WRITTEN by the file (mzLib keeps it a string, e.g. "35" or
  "25,30,35"); not parsed.

- `compensation_voltage`:

  float; in **V**; `NA` when no FAIMS on this instrument. FAIMS
  compensation voltage.

- `mz`:

  float\[\]; in **m/z**; `NA` when scan has no spectrum. Peak m/z
  values. *Present only with `peaks`.*

- `intensity`:

  float\[\]; in **intensity (instrument units)**; `NA` when scan has no
  spectrum. Peak intensities, parallel to mz. *Present only with
  `peaks`.*

## On the wire but not projected yet

The bridge sends these, and this version of mzLibR does not return them
yet:

- `source`:

  arrives with the mzLib 1.0.592 port.

- `rows_not_read`:

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

  path missing or blank; file not found; the file is not one of the
  seven spectra types (the message names the views it does have)

- `mzlib_usage_error` (usage):

  ms-order, offset, limit or out given with no value; ms-order \< 1;
  offset or limit \< 0; non-integer value

- `mzlib_usage_error` (usage):

  out resolves to the same file as path

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

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Headers only unless peaks is set; peak_count tells you how many peaks,
  not what they are.

- Bruker .d and timsTOF .d need vendor native DLLs and read only on
  Windows x64. Thermo .raw is managed and reads everywhere.

- MGF: scan numbers may be file order; precursor scan is null unless the
  mzLib PRECURSORSCAN header is present; ms_order comes from MSLEVEL or
  is inferred from PEPMASS; the scan window is derived from observed
  peaks.

- msalign holds DECONVOLVED neutral monoisotopic masses, not raw m/z:
  the mz column is not comparable with an mzML or raw file's mz and must
  not be re-deconvolved.

- msalign: scan_window_lower_mz/\_upper_mz are in m/z while the mz
  column is neutral MASS, so they are not on the same axis. mzLib
  synthesises the window by converting each mass back to m/z at its
  reported charge (MsAlign.cs:526); msalign records no window of its
  own.

- timsTOF .d is ion-mobility-resolved and this view flattens it: the
  mobility dimension is collapsed into scans and no 1/K0 value is
  reported.

- source.acquisition_start_time from a Thermo .raw is the acquisition
  PC's local clock (is_utc false); ProteoWizard's mzML of the same run
  writes it as UTC assuming the CONVERTING machine's time zone, so the
  two can differ by the site's UTC offset
  (SourceFile.AcquisitionStartTime).

- Many files: `--paths-stdin` reads the list in ONE process, `--threads`
  at a time, into one long table with source_index and source_path
  first. The output does not depend on `--threads` (a tested property);
  the default of 1 is a memory choice, since every reader materialises
  the whole file.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_spectra`; many files:
  `pymzlib.readers.read_spectra_many`

- Rust (mzLibRust): `mzlib::readers::read_spectra_with` with
  `SpectraOptions`; many files: `mzlib::readers::read_spectra_many`

- R (mzLibR): `readers_read_spectra`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- injection_time unit 'ms' is from the mzML CV term (MS:1000927); not
  yet confirmed that every mzLib vendor reader (Thermo, Bruker, timsTOF,
  MGF) emits ms.

- source is read from MsDataFile.SourceFile after LoadAllStaticData;
  Bruker .d and timsTOF .d have not been checked for which of its fields
  they fill.

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

## See also

[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)

## Examples

``` r

scans <- readers_read_spectra("sliced_ethcd.mzML", limit = 3)
scans
#> <mzlibr_scan_records> E:\CodeReview\pymzlib-readers-wt\code\mzLib\mzLib\Test\DataFiles\sliced_ethcd.mzML (MzML, Mzml)
#>   6 scans in the file
#>   6 records in the file, 3 returned
#>   ! truncated - records were left behind
#>   retention_time_unit: minutes
#>   ! Peaks are not included. This is scan HEADERS only; peak_count reports how many peaks each scan has but not what they are. Pass peaks=true to include the mz and intensity arrays, and expect the payload to grow by roughly the size of the file.
#>   peaks not included - pass peaks = TRUE for the mz and intensity arrays
scans$records[, c("one_based_scan_number", "ms_order", "retention_time", "peak_count")]
#>   one_based_scan_number ms_order retention_time peak_count
#> 1                     1        1       38.92572        484
#> 2                     2        2       38.92606        409
#> 3                     3        2       38.93049        339
```
