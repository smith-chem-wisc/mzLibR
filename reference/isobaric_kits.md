# The isobaric kits mzLib can name, with every channel's label and reporter-ion m/z

List the isobaric labelling kits mzLib can name, with every channel's
label and reporter-ion m/z.

Calls mzLib's `IsobaricMassTag.TryGetIsobaricMassTag` for each
`IsobaricMassTagType`, or for the one `kit` resolves to through
`IsobaricMassTag.TryGetTagType`. No file, no network.

## Usage

``` r
isobaric_kits(kit = NULL, timeout = 60)
```

## Arguments

- kit:

  One kit, matched by mzLib's whole-name rule: case-insensitive,
  optionally followed by `" on <motif>"` (`"TMT10"`, `"TMT6-plex"`,
  `"iTRAQ-4plex on K"`). Never a substring, so `"TMT10plex"` is refused
  rather than mistaken for a kit whose name it contains. `NULL`, the
  default, lists every kit.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**Nothing in the table is typed in.** mzLib holds only the channel
labels. Each m/z is a `DI HCD` diagnostic-ion line of the kit's
`Multiplex Label` modification in mzLib's embedded `TMT.txt`, plus one
proton, sorted ascending and paired with the labels by position - the
same m/z MetaMorpheus uses to quantify a TMT search.

## Value

An `mzlibr_isobaric_kits`: `kit`, the name you asked for exactly as
given (`NA` when every kit was listed); `kit_count` kits listed;
`record_count` channels listed over every kit; `absolute_tolerance`, the
half-width in Da of each channel's matching window; `kits`, a data.frame
of `kit` and `channel_count`, in mzLib's order; `caveats`; and
`records`, a data.frame with one row per (kit, channel): `kit`,
`channel_index` (1-based within the kit, ascending with the m/z),
`channel_label` as MetaMorpheus spells it, `reporter_ion_mz` (m/z,
charge 1, theoretical), and `mz_min` and `mz_max`, the matching window
in m/z.

## What the m/z are not

`reporter_ion_mz` is theoretical, at charge 1, from mzLib's `TMT.txt`:
not a calibrated or observed value. mzLib reads a reporter intensity as
the most intense peak between `mz_min` and `mz_max`. `TMT16` is the
lowest sixteen channels of TMTpro 18-plex, and `iTRAQ8` has no 120
channel (the phenylalanine immonium ion sits there); its eighth reagent
is 121.

## Wraps

Wire verb `isobaric kits`. Generated from the bridge's verb spec
`isobaric.kits.yaml` (bridge commit `e76157831b15`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`IsobaricMassTag.TryGetIsobaricMassTag`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/Modifications/IsobaricMassTag.cs)
  in `mzLib/Omics/Modifications/IsobaricMassTag.cs` at mzLib `0a808fec`

- [`IsobaricMassTag.TryGetTagType`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/Modifications/IsobaricMassTag.cs)
  in `mzLib/Omics/Modifications/IsobaricMassTag.cs` at mzLib `0a808fec`

- [`IsobaricMassTagType`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/Modifications/IsobaricMassTag.cs)
  in `mzLib/Omics/Modifications/IsobaricMassTag.cs` at mzLib `0a808fec`

## Parameters: units, ranges and defaults

- `kit`:

  string; default absent. One kit, resolved by mzLib's whole-name match,
  case-insensitive, optionally followed by " on \<motif\>": TMT10,
  tmt10, TMT6-plex, iTRAQ-4plex on K. Never a substring, so TMT10plex is
  unknown. An unknown or blank name is a usage error naming the known
  kits. Absent: every kit.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `kit`:

  string; `NA` when `--kit` was not given: every kit is listed. The
  `--kit` name exactly as given.

- `kit_count`:

  int; in **kits**; never `NA`. Kits listed.

- `record_count`:

  int; in **channels**; never `NA`. Rows in columns: the channels of
  every kit listed.

- `absolute_tolerance`:

  float; in **Da**; never `NA`. IsobaricMassTag.AbsoluteToleranceValue:
  the half-width of each channel's matching window.

- `kits`:

  object\[\]; never `NA`. One per kit in IsobaricMassTagType order:
  {kit, channel_count}.

- `column_names`:

  string\[\]; never `NA`. The table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. One row per (kit, channel), kits in enum order,
  channels in ascending m/z.

- `caveats`:

  string\[\]; never `NA`. Five static caveats.

## Columns of `records`

- `kit`:

  string; never `NA`. mzLib's IsobaricMassTagType member: TMT6, TMT10,
  TMT11, TMT16, TMT18, iTRAQ4, iTRAQ8, diLeu4, diLeu12 (at 0a808fec).
  TMT16 and TMT18 are TMTpro.

- `channel_index`:

  int; never `NA`. 0-based position within the kit, ascending with
  reporter_ion_mz.

- `channel_label`:

  string; never `NA`. The channel's name as MetaMorpheus spells it: 126,
  127N, 127C, ..., 115a.

- `reporter_ion_mz`:

  float; in **m/z**; never `NA`. Theoretical m/z at charge 1: a DI HCD
  diagnostic-ion mass of the kit's Multiplex Label modification in
  mzLib's TMT.txt, plus one proton.

- `mz_min`:

  float; in **m/z**; never `NA`. reporter_ion_mz minus
  absolute_tolerance.

- `mz_max`:

  float; in **m/z**; never `NA`. reporter_ion_mz plus
  absolute_tolerance.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--kit` blank, or not a name mzLib knows (the message lists the kits)

- `mzlib_bridge_error` (correctness):

  mzLib cannot build a kit from its embedded modifications (never at
  0a808fec)

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- reporter_ion_mz is theoretical, charge 1, derived from mzLib's
  TMT.txt; not calibrated or observed.

- mz_min and mz_max are mzLib's own reporter-ion matching window, the
  most intense peak inside it wins.

- TMT16 is the lowest sixteen channels of TMTpro 18-plex and has no
  modification entry of its own.

- iTRAQ8 has no 120 channel (the phenylalanine immonium ion sits at
  120.081); its eighth reagent is 121.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.isobaric.kits`

- Rust (mzLibRust): `mzlib::isobaric::kits`

- R (mzLibR): `isobaric_kits`

## Since

Wire protocol 1; pyMzLib 0.3.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- A second verb could read reporter intensities out of spectra through
  IsobaricMassTag.GetReporterIonIntensities; not built.

- Not yet ported: mzLibRust and mzLibR names above are proposals until
  their ports land (parity order).

## References

- [doi:10.1021/acs.jproteome.1c00168](https://doi.org/10.1021/acs.jproteome.1c00168):
  the TMTpro-18plex reagents paper (Li et al. 2021); mzLib cites it as
  the source of its 0.003 Da reporter-ion matching tolerance
  (IsobaricMassTag.AbsoluteToleranceValue, IsobaricMassTag.cs:41)

## See also

[`sdrf_design`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md),
for a label-free design from an SDRF.

## Examples

``` r

catalogue <- isobaric_kits()
catalogue$kits
#>       kit channel_count
#> 1    TMT6             6
#> 2   TMT10            10
#> 3   TMT11            11
#> 4   TMT16            16
#> 5   TMT18            18
#> 6  iTRAQ4             4
#> 7  iTRAQ8             8
#> 8  diLeu4             4
#> 9 diLeu12            12
catalogue$absolute_tolerance
#> [1] 0.003
catalogue$records$channel_label[catalogue$records$kit == "iTRAQ8"]
#> [1] "113" "114" "115" "116" "117" "118" "119" "121"

# TMTpro 18-plex, whose top channel sits at 135.15:
tmtpro <- isobaric_kits("TMT18")
utils::tail(tmtpro$records[, c("channel_index", "channel_label", "reporter_ion_mz")], 3)
#>    channel_index channel_label reporter_ion_mz
#> 16            16          134N        134.1482
#> 17            17          134C        134.1546
#> 18            18          135N        135.1516
```
