# Isobaric kits: TMT, TMTpro, iTRAQ and DiLeu channels

An isobaric experiment - TMT, TMTpro, iTRAQ, DiLeu - is read out of the
low-m/z **reporter ions**: one ion per channel, some only millidaltons
apart. Mislabel one channel, or shift its m/z, and every sample after it
is assigned to the wrong condition, with nothing in the result looking
wrong. So mzLibR keeps no table of its own.
[`isobaric_kits()`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
returns **mzLib’s**, the one MetaMorpheus quantifies TMT searches with:

| question | function | mzLib |
|----|----|----|
| every kit mzLib knows, with every channel | [`isobaric_kits()`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md) | `IsobaricMassTag`, `IsobaricMassTagType` (#1375) |
| one kit, by its name or a MetaMorpheus modification name | `isobaric_kits("TMT18")` | `IsobaricMassTag.TryGetTagType` |
| a channel’s reporter-ion m/z and matching window | `records$reporter_ion_mz`, `mz_min`, `mz_max` | `ReporterIonMzs`, `ReporterIonMzRanges` |

## Every kit

``` r

catalogue <- isobaric_kits()
catalogue
#> <mzlibr_isobaric_kits> 9 kit(s), 89 channels, matched within +/- 0.003 Da
#>   TMT6 (6), TMT10 (10), TMT11 (11), TMT16 (16), TMT18 (18), iTRAQ4 (4), iTRAQ8 (8), diLeu4 (4), diLeu12 (12)
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
```

`records` is one row per (kit, channel). `channel_index` is 1-based
within its kit, ascending with the m/z.

**`TMT16` and `TMT18` are TMTpro.** TMT16 is the lowest sixteen channels
of the 18-plex set and has no modification entry of its own, so its m/z
are TMT18’s first sixteen. **iTRAQ 8-plex has no 120 channel**: the
phenylalanine immonium ion sits at m/z 120.081, so the eighth reagent is
121.

``` r

catalogue$records$channel_label[catalogue$records$kit == "iTRAQ8"]
#> [1] "113" "114" "115" "116" "117" "118" "119" "121"
```

## One kit: TMTpro 18-plex

Ask for a kit the way MetaMorpheus names it. The match is on the **whole
name**, case-insensitive, optionally followed by `" on <motif>"`:
`"TMT18"`, `"tmt10"`, `"TMT6-plex"` and `"iTRAQ-4plex on K"` all work. A
substring never does, so `"TMT10plex"` is refused with a list of the
kits that exist, rather than quietly resolved to whichever kit’s name it
contains.

``` r

tmtpro <- isobaric_kits("TMT18")
head(tmtpro$records[, c("channel_index", "channel_label", "reporter_ion_mz", "mz_min", "mz_max")], 4)
#>   channel_index channel_label reporter_ion_mz   mz_min   mz_max
#> 1             1           126        126.1277 126.1247 126.1307
#> 2             2          127N        127.1248 127.1218 127.1278
#> 3             3          127C        127.1311 127.1281 127.1341
#> 4             4          128N        128.1281 128.1251 128.1311
```

## What the numbers are

Each `reporter_ion_mz` is a `DI HCD` diagnostic-ion line of the kit’s
`Multiplex Label` modification in mzLib’s embedded `TMT.txt`, plus one
proton: theoretical, at charge 1, never calibrated or observed. mzLib
reads a reporter intensity as the most intense peak between `mz_min` and
`mz_max`, `absolute_tolerance` Da either side:

``` r

tmtpro$absolute_tolerance
#> [1] 0.003
```

mzLib takes that tolerance from the TMTpro 18-plex paper (Li et al., J
Proteome Res 2021,
[doi:10.1021/acs.jproteome.1c00168](https://doi.org/10.1021/acs.jproteome.1c00168)),
as its `IsobaricMassTag.AbsoluteToleranceValue` says.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1021/acs.jproteome.1c00168](https://doi.org/10.1021/acs.jproteome.1c00168):
  the TMTpro-18plex reagents paper (Li et al. 2021); mzLib cites it as
  the source of its 0.003 Da reporter-ion matching tolerance
  (IsobaricMassTag.AbsoluteToleranceValue, IsobaricMassTag.cs:41).
