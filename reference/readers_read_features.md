# Read deconvolved MS1 features, in the cross-format `ms1_features` view

Read deconvolved MS1 features through the cross-format ms1_features
view: m/z, charge, retention-time range, apex intensity and isotope
count.

Two file types offer it: TopFD/FLASHDeconv `_ms1.feature` and Dinosaur
`.feature.tsv`. A file without the view raises, with a message naming
the views it does have.

## Usage

``` r
readers_read_features(path, limit = NULL, offset = 0, out = NULL, timeout = NULL)
```

## Arguments

- path:

  Path to an `_ms1.feature` or Dinosaur `.feature.tsv`.

- limit:

  Maximum features to return. `NULL`, the default, returns all of them.

- offset:

  Features to skip.

- out:

  Write a tab-separated table here and return only a summary.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_feature_records`. `record_count` counts the features in the
whole file and `returned_count` the features returned, starting `offset`
features in.

`records` has one row per feature: `mz` in m/z; `charge`;
`retention_time_start` and `retention_time_end` in the unit
`retention_time_unit` names - `"unknown"` for `_ms1.feature`, see below;
`intensity`, the apex, in the instrument's intensity units; and
`number_of_isotopes`.

## One row is not one line of the file, for `_ms1.feature`

An `_ms1.feature` row is a deconvolved **neutral mass spanning a charge
range**, and mzLib expands it into one single-charge feature per charge
in that range. A hundred-feature file can read as a thousand rows.
**Dinosaur is one-for-one.** Either way
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
gives you the file's own rows.

`intensity` is the **apex** intensity, not the sum over the feature.
Both formats carry a summed intensity column too, and
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
has it.

## The retention-time unit is genuinely unknown for `_ms1.feature`

TopFD wrote seconds through v1.6.2 and minutes from v1.7.0 - \*within
the same file type\*, with nothing in the file to say which. mzLib
normalises neither, and its own deconvolution code resorts to a
heuristic (divide by 60 if the largest end time exceeds 500). This
package will not launder a guess into a stated fact, so
[`readers_retention_time_in_minutes`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
raises rather than converting. Dinosaur reports minutes and converts
without complaint.

## Wraps

Wire verb `readers read-features`. Generated from the bridge's verb spec
`readers.read-features.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`IMs1FeatureFile.GetMs1Features`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/BaseClasses/IMs1FeatureFile.cs)
  in `mzLib/Readers/BaseClasses/IMs1FeatureFile.cs` at mzLib `23c2490e`

- [`Ms1Feature.GetSingleChargeFeatures`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/ExternalResults/IndividualResultRecords/Ms1Feature.cs)
  in
  `mzLib/Readers/ExternalResults/IndividualResultRecords/Ms1Feature.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. A TopFD/FLASHDeconv \_ms1.feature or a Dinosaur
  .feature.tsv.

- `offset`:

  int; in **features**; default `0`; range `>= 0`. Skip this many
  features. A window, not a cursor: mzLib parses the whole file on every
  call, so paging re-reads it.

- `limit`:

  int; in **features**; default absent; range `>= 0`. Return at most
  this many features. Absent means no limit; there is deliberately no
  default cap.

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

  int; in **features**; never `NA`. Features in the whole file, before
  the window.

- `retention_time_unit`:

  string; never `NA`. minutes for Dinosaur; unknown for \_ms1.feature,
  whose unit changed at TopFD v1.7.0 without a format change.

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

  int; in **features**; never `NA`. Features returned in columns (the
  unit offset and limit count in); 0 when written to out.

- `offset`:

  int; in **features**; never `NA`. The offset applied, in features.

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

- `mz`:

  float; in **m/z**; never `NA`. Monoisotopic m/z of this single-charge
  feature.

- `charge`:

  int; in **charge**; never `NA`. Charge state.

- `retention_time_start`:

  float; in **retention_time_unit**; never `NA`. Feature start, in
  retention_time_unit.

- `retention_time_end`:

  float; in **retention_time_unit**; never `NA`. Feature end, in
  retention_time_unit.

- `intensity`:

  float; in **intensity (instrument units)**; `NA` when the file has no
  apex intensity (FLASHDeconv schema, then in absent_fields) or this
  feature has none. Apex intensity (Dinosaur intensityApex; TopFD
  Apex_intensity).

- `number_of_isotopes`:

  int; `NA` when \_ms1.feature: mzLib never sets it (then in
  absent_fields). Isotopes in the feature.

## On the wire but not projected yet

The bridge sends these, and this version of mzLibR does not return them
yet:

- `reader`:

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

  path missing or blank; file or directory not found; mzLib does not
  recognise the file type (the message points at readers formats)

- `mzlib_usage_error` (usage):

  offset, limit or out given with no value; offset or limit \< 0;
  non-integer value; out resolves to the same file as path

- `mzlib_usage_error` (usage):

  the file has no ms1_features view

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

- \_ms1.feature rows are expanded to one per charge in \[ChargeStateMin,
  ChargeStateMax\] (Ms1Feature.cs:84); Dinosaur is one row per file
  line.

- intensity is the per-charge APEX intensity (Ms1Feature.cs:86); a
  FLASHDeconv file has none and names intensity in absent_fields instead
  of reporting mzLib's zero.

- number_of_isotopes is absent for every \_ms1.feature
  (Ms1Feature.cs:91).

- Many files: `--paths-stdin` reads the list in ONE process, `--threads`
  at a time, into one long table with source_index and source_path
  first. The output does not depend on `--threads` (a tested property);
  the default of 1 is a memory choice, since every reader materialises
  the whole file.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.read_features`; many files:
  `pymzlib.readers.read_features_many`

- Rust (mzLibRust): `mzlib::readers::read_features_with` with
  `ReadOptions`; many files: `mzlib::readers::read_features_many`

- R (mzLibR): `readers_read_features`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- retention_time_start/\_end declare their unit as the envelope's
  retention_time_unit, not a fixed unit, because it genuinely varies per
  file.

## See also

[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md),
[`readers_retention_time_in_minutes`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)

## Examples

``` r

features <- readers_read_features("Ms1Feature_TopFDv1.6.2_ms1.feature", limit = 5)
features
#> <mzlibr_feature_records> E:\CodeReview\pymzlib-readers-wt\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature (Ms1Feature)
#>   25 records in the file, 5 returned
#>   ! truncated - records were left behind
#>   retention_time_unit: unknown
#>   ! One row here is one CHARGE STATE of one feature, not one row of the file. mzLib expands each deconvolved feature across [ChargeStateMin, ChargeStateMax] (Ms1Feature.cs:84), so record_count exceeds the file's line count, and a charge the tool never observed appears if the writer recorded a gapped charge range. Use read-records for the file's own rows.
#>   ! intensity is the per-charge APEX intensity (Ms1Feature.cs:86), not the summed intensity over the feature. The file's own Intensity column is a different number, and read-records has it.
#>   ! retention_time_start/_end are in UNKNOWN units for this format. TopFD wrote seconds through v1.6.2 and minutes from v1.7.0 without changing the file type, and mzLib normalises neither - its deconvolution parameters instead GUESS, dividing by 60 when the largest end time exceeds 500. Check the values against your gradient length before comparing them with anything.
#>   ! number_of_isotopes is null for every row of this format: the single-charge expansion mzLib builds never sets it (Ms1Feature.cs:91). Null means 'not reported', not 'no isotopes found'.
features$records
#>          mz charge retention_time_start retention_time_end intensity
#> 1 1548.9862      7              2372.27            2401.92 912795139
#> 2 1355.4889      8              2372.27            2401.92 912795139
#> 3 1204.9909      9              2372.27            2401.92 912795139
#> 4 1084.5925     10              2372.27            2401.92 912795139
#> 5  986.0848     11              2372.27            2401.92 912795139
#>   number_of_isotopes
#> 1                 NA
#> 2                 NA
#> 3                 NA
#> 4                 NA
#> 5                 NA
```
