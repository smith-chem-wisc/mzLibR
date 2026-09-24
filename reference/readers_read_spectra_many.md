# Read the scans of many spectra files, in one bridge call

Read the scans of a spectra file: every scan's header always, and its
peak arrays only on request.

The many-files form of
[`readers_read_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md):
every file read by one bridge process, `threads` at a time, into one
scan table - so a study's runs can be compared by instrument and
acquisition time from `files` without a loop.

## Usage

``` r
readers_read_spectra_many(paths, ms_order = NULL, peaks = FALSE, out = NULL,
  threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- paths:

  A character vector of paths to spectra files. Bruker `.d` directories
  are accepted.

- ms_order:

  Keep only scans at this MS level in every file, or `NULL` for every
  level.

- peaks:

  Include the `mz` and `intensity` arrays. `FALSE` by default.

- out:

  Write the long table here as tab-separated text and return only a
  summary.

- threads:

  Files read at once. `1`, the default, holds one whole file in memory
  at a time; `-1` uses one per core. The table does not depend on it.

- on_error:

  `"fail"`, the default, or `"skip"` to record an unreadable file in
  `files` and read the rest.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_read_batch`. `records` begins with `source_index` (1-based,
into `paths`) and `source_path`, then the columns of
[`readers_read_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md):
`retention_time` in minutes; `injection_time` in ms; `total_ion_current`
and `selected_ion_intensity` in the instrument's intensity units;
`peak_count` in peaks; `compensation_voltage` in volts;
`selected_ion_charge_state_guess` a charge; `scan_window_lower_mz`,
`scan_window_upper_mz`, `isolation_mz`, `isolation_width`,
`selected_ion_mz` and `selected_ion_monoisotopic_guess_mz` in m/z; and
with `peaks = TRUE` the list columns `mz` (m/z) and `intensity`. `files`
has one row per input: `scan_count` in scans, `record_count` in scans
after the `ms_order` filter, the run's `instrument_model`,
`instrument_model_accession`, `instrument_serial_number`,
`acquisition_start_time` and `acquisition_start_time_is_utc` (flattened
from each file's `source`; `NA` where the file does not record them),
`caveats` and the error columns. `file_count`, `read_count` and
`failed_count` count files; `record_count` and `returned_count` count
scans over the files read; `row_count` counts the rows of `records`.

## Wraps

Wire verb `readers read-spectra` with `--paths-stdin`. Generated from
the bridge's verb spec `readers.read-spectra.yaml` (bridge commit
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`MsDataFileReader.GetDataFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/MsDataFileReader.cs)
  in `mzLib/Readers/MsDataFileReader.cs` at mzLib `23c2490e`

- [`MsDataFile.GetAllScansList`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/MassSpectrometry/MsDataFile.cs)
  in `mzLib/MassSpectrometry/MsDataFile.cs` at mzLib `23c2490e`

- [`SourceFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/MassSpectrometry/SourceFile.cs)
  in `mzLib/MassSpectrometry/SourceFile.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `ms_order` (wire `--ms-order`):

  int; default absent; range `>= 1`. Keep only scans at this MS level.
  Applied BEFORE offset and limit, so ms-order 2 with limit 10 is the
  first ten MS2 scans. Given with no value is a usage error; absent
  means every level.

- `peaks`:

  flag; default `FALSE`. Include the mz and intensity arrays. Off by
  default because peaks are roughly the size of the file; combine with
  limit, a scan filter, or out.

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

  path; default absent. Write the long table here, one file at a time in
  input order, so memory holds at most threads files. A batch that stops
  on an error deletes its partial table. Must differ from every input.

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

  int; in **scans**; never `NA`. Scans summed over the files read
  (BULK.md section 2).

- `returned_count`:

  int; in **scans**; never `NA`. Scans returned in columns, summed over
  the files; 0 when written to out.

- `row_count`:

  int; in **rows**; never `NA`. Rows in columns (more than
  returned_count for a long view); 0 when written to out.

- `on_error`:

  string; never `NA`. fail or skip, as requested. `--threads` is
  deliberately not echoed: the envelope is byte-identical at every
  thread count.

- `ms_order`:

  int; `NA` when no ms-order filter was applied. The ms-order filter
  every file was read with.

- `peaks_included`:

  bool; never `NA`. Whether mz and intensity columns are present.

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
  block (path, file_type, reader, scan_count, source, record_count,
  rows_not_read, retention_time_unit, caveats, column_names,
  absent_fields, failed_fields, excluded_fields, error). An unread input
  has error set and every fact it could not establish null or empty.

## Columns of `records`

The first two columns say which input each row came from -
`source_index` is 1-based in R, so it indexes `files` directly - and
rows are grouped by input, in input order, whatever `threads` is.

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

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

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
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_spectra`; many files:
  `pymzlib.readers.read_spectra_many`

- Rust (mzLibRust): `mzlib::readers::read_spectra_with` with
  `SpectraOptions`; many files: `mzlib::readers::read_spectra_many`

- R (mzLibR): `readers_read_spectra`; many files:
  `readers_read_spectra_many`

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

[`readers_read_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)

## Examples

``` r

runs <- readers_read_spectra_many(
  c("sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"),
  threads = 2, on_error = "skip"
)
runs
#> <mzlibr_read_batch> readers read-spectra
#>   3 inputs: 2 read, 1 failed (on_error = "skip")
#>   8 rows in `records`, 8 records
#>   ! C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\no-such-run.mzML: File not found: 'C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\no-such-run.mzML'.
runs$files[, c("file_type", "scan_count", "instrument_model", "error_kind")]
#>   file_type scan_count instrument_model error_kind
#> 1      MzML          6  Orbitrap Fusion       <NA>
#> 2      MzML         NA             <NA>      usage
#> 3       Mgf          2             <NA>       <NA>
table(runs$records$source_index)
#> 
#> 1 3 
#> 6 2 
```
