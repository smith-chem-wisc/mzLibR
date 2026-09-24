# Every file type mzLib can recognise

List every file type mzLib recognises, with the extension it dispatches
on, the reader class that parses it, and the cross-format views it
offers.

Enumerated from mzLib itself rather than from a list maintained here, so
it reflects the installed version and cannot go stale.

## Usage

``` r
readers_formats(timeout = 60)
```

## Arguments

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

A data.frame with one row per format: `file_type`, `extension`,
`reader`, `is_quantifiable`, and a `views` list column.

## Views, and why most formats have none

It is tempting to read "31 formats" as "31 formats in one uniform
shape". They are not. The formats fall into disjoint families, and **14
of the 31 belong to none of them** — an empty `views` is a real and
common answer, meaning mzLib can parse the file but offers no
cross-format projection of it.

\- `"quantifiable"` — the cross-format record view, and the input
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
accepts. **Exactly four file types have it**: MetaMorpheus `psmtsv` and
`osmtsv`, `MsFraggerPsm`, and DIA-NN `DiaNnReport`. - `"ms1_features"` —
deconvolved MS1 features (TopFD `_ms1.feature`, Dinosaur). - `"spectra"`
— the file is spectra, not results. - `"spectral_match"` —
identifications that share no file-level interface.

Note also that `extension` is **not unique**: `BrukerD` and
`BrukerTimsTof` are both `.d`, told apart by what the directory holds,
and several formats share `.tsv`.

## Wraps

Wire verb `readers formats`. Generated from the bridge's verb spec
`readers.formats.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SupportedFileType`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Util/SupportedFileTypes.cs)
  in `mzLib/Readers/Util/SupportedFileTypes.cs` at mzLib `23c2490e`

- [`SupportedFileTypeExtensions.GetFileExtension`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Util/SupportedFileTypes.cs)
  in `mzLib/Readers/Util/SupportedFileTypes.cs` at mzLib `23c2490e`

- [`SupportedFileTypeExtensions.GetResultFileType`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Util/SupportedFileTypes.cs)
  in `mzLib/Readers/Util/SupportedFileTypes.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

This verb takes no parameters.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `format_count`:

  int; in **formats**; never `NA`. Number of entries in formats; equal
  to the number of SupportedFileType members. 36 at mzLib 1.0.592 (32 at
  1.0.591).

- `formats`:

  table; never `NA`. One entry per SupportedFileType, in enum
  declaration order. Enumerated from mzLib, not transcribed, so it
  cannot drift from what mzLib dispatches.

## Columns

- `file_type`:

  string; never `NA`. mzLib's SupportedFileType member name, e.g.
  "MsFraggerPsm", "psmtsv", "MzIdentMLGz".

- `extension`:

  string; `NA` when mzLib has no extension mapping for this member (a
  broken mzLib build); the listing reports it rather than failing. The
  suffix mzLib dispatches on, e.g. ".mzid.gz" or
  "QuantifiedProteinGroups.tsv". NOT unique: both Bruker types are ".d"
  and several share ".tsv"; DiaNnReport and PytheasResult dispatch on
  content, not on this value.

- `reader`:

  string; `NA` when mzLib has no reader mapping for this member (a
  broken mzLib build). The mzLib reader class name, e.g.
  "MzIdentMLResultFile".

- `views`:

  string\[\]; never `NA`. Cross-format views this type offers, from
  quantifiable, ms1_features, spectra, spectral_match. Empty is a real
  answer: 17 of 36 types offer none and are readable only through
  read-records.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  Never raised by this verb.

- `mzlib_bridge_error` (correctness):

  Never raised by this verb.

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Extensions are not unique and not always how dispatch works: identify
  a real path rather than matching on extension.

- An empty views list means mzLib parses the format into its own record
  type with no cross-format interface; read-records still reads it.

- At mzLib 1.0.592 the view families are: quantifiable 4, ms1_features
  2, spectral_match 6, spectra 7, none 17. These counts are the pin's,
  not the contract's: bindings should call formats rather than repeat
  them in prose.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.formats`

- Rust (mzLibRust): `mzlib::readers::formats`

- R (mzLibR): `readers_formats`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## See also

[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)

## Examples

``` r

formats <- readers_formats()
nrow(formats)
#> [1] 31
formats[formats$is_quantifiable, c("file_type", "extension", "reader")]
#>       file_type  extension           reader
#> 12       psmtsv    .psmtsv   PsmFromTsvFile
#> 13       osmtsv    .osmtsv   OsmFromTsvFile
#> 18 MsFraggerPsm    psm.tsv MsFraggerPsmFile
#> 30  DiaNnReport report.tsv  DiaNnReportFile

# An empty `views` is a real answer: mzLib reads the file but offers no cross-format view.
sum(lengths(formats$views) == 0)
#> [1] 14
```
