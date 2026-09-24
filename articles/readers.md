# Reading mass-spectrometry and search-result files

mzLib recognises dozens of file types: instrument data (mzML, Thermo
`.raw`, Bruker and timsTOF `.d`, MGF, msalign) and the output of a dozen
search and deconvolution tools. The `readers_` functions ask it what a
file is and read its records into data.frames.

The one idea to take away: **mzLib does not read every format into one
uniform shape.** The formats fall into a few disjoint *views*, and many
belong to none. This page shows how to tell which you have, and which
function reads it.

## What is this file?

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
```

A format’s `views` say which cross-format projection it offers. An empty
`views` is a real and common answer: mzLib can parse the file, but no
function other than
[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
reads it.

``` r

table(vapply(formats$views, function(v) if (length(v)) paste(v, collapse = "+") else "(none)",
             character(1)))
#> 
#>         (none)   ms1_features   quantifiable        spectra spectral_match 
#>             14              2              4              7              4
```

For one file,
[`readers_identify()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
resolves the type without parsing the contents:

``` r

readers_identify("PXD078927_msgf_1_1_0.mzid")
#> <mzlibr_file_info> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\DataFiles\PXD078927_msgf_1_1_0.mzid
#>   MzIdentML (.mzid), read by MzIdentMLResultFile
#>   views: spectral_match
```

## Five ways to read, and how to choose

| the file offers | read it with | columns |
|----|----|----|
| `spectra` | [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md) | one row per scan |
| `quantifiable` | [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md) | the uniform record view FlashLFQ takes |
| `ms1_features` | [`readers_read_features()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md) | deconvolved MS1 features |
| `spectral_match` | [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md) | identifications |
| anything | [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md) | the format’s own fields |

### Spectra

Scan headers always; peak arrays only when asked for with
`peaks = TRUE`, because they are roughly the size of the file.

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
scans$records[, c("one_based_scan_number", "ms_order", "retention_time", "peak_count",
                  "selected_ion_mz")]
#>   one_based_scan_number ms_order retention_time peak_count selected_ion_mz
#> 1                     1        1       38.92572        484              NA
#> 2                     2        2       38.92606        409        548.4539
#> 3                     3        2       38.93049        339        796.7652
```

`retention_time` is in minutes for every spectra format.
`selected_ion_mz` is `NA` on an MS1 scan, which has no precursor.
[`?readers_read_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
lists every column with its unit and what `NA` means in it.

### Search results in the uniform view

``` r

psms <- readers_read_results("FraggerPsm_FragPipev21.1_psm.tsv", limit = 2)
psms
#> <mzlibr_result_records> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv (MsFraggerPsm)
#>   5 records in the file, 2 returned
#>   ! truncated - records were left behind
#>   retention_time_unit: minutes
#>   ! is_decoy is null for this format: MSFragger's psm.tsv carries no target/decoy column, so mzLib cannot report decoy status (MsFraggerPsm.cs:231) and the field crosses as null. Null means 'unknown', not 'target' - do not filter this format on is_decoy == false.
#>   ! monoisotopic_mass is the THEORETICAL peptide mass (MsFraggerPsm.cs:233, CalculatedPeptideMass), not the observed precursor mass. The psmtsv formats report the theoretical mass here too, so the two are consistent - but neither is what the instrument measured.
#>   ! file_name is the full 'Spectrum File' path including its .pep.xml extension, whereas the psmtsv formats report a bare base name. The field is not a join key across formats.
psms$records[, c("base_sequence", "charge_state", "retention_time", "is_decoy")]
#>   base_sequence charge_state retention_time is_decoy
#> 1       KPVGAAK            2     0.03233000       NA
#> 2     KPAAAAGAK            2     0.04907667       NA
```

Read the `caveats` before trusting a column. Here they say that
MSFragger records no decoy label, so `is_decoy` is `NA` - unknown -
rather than a `FALSE` that would read as “target”. Convert retention
times with
[`readers_retention_time_in_minutes()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
rather than by hand: it uses the file’s own `retention_time_unit`, and
refuses when mzLib gives no basis for one.

``` r

readers_retention_time_in_minutes(psms)
#> [1] 0.03233000 0.04907667
```

### Features, matches, and everything else

``` r

features <- readers_read_features("Ms1Feature_TopFDv1.6.2_ms1.feature", limit = 5)
features$retention_time_unit
#> [1] "unknown"
try(readers_retention_time_in_minutes(features))
#> Error : Cannot convert retention time for 'Ms1Feature': mzLib gives no basis to say what unit it is in. Inspect the values against scan numbers before comparing them.
```

TopFD changed the unit of `_ms1.feature` retention times at v1.7.0
without changing the format, so the honest unit is `"unknown"`, and
conversion refuses rather than guessing.

``` r

matches <- readers_read_matches("Casanovo_5.0.0.mztab")
matches$records[, c("one_based_scan_number", "base_sequence", "modification_count")]
#>   one_based_scan_number   base_sequence modification_count
#> 1                     1      AGAHLQGGAK                  0
#> 2                     3       RGTGVENVK                  0
#> 3                     4 LDDPKEEDEEKEEGK                  0
#> 4                     5     HQGVMVGMGQK                  2
#> 5                     7         RQEFEMK                  1
```

[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
reads any recognised file into that format’s own fields, under mzLib’s
names, so a column can be looked up in the mzLib source:

``` r

prsms <- readers_read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv", limit = 3)
prsms$record_type
#> [1] "ToppicPrsm"
prsms$records[, c("one_based_scan_number", "base_sequence", "e_value")]
#>   one_based_scan_number base_sequence e_value
#> 1                   259        GYASDS  1e+300
#> 2                   270         THIGY  1e+300
#> 3                   354      SAECTKRF  1e+300
prsms$excluded_fields
#>                         field                      type
#> 1 alternative_identifications List<AlternativeToppicId>
#>                                                    reason
#> 1 a list of composite values has no faithful column shape
```

`excluded_fields` names what could not become a column, with the reason,
so a missing column is never mistaken for a field the format does not
have.

## Large files

Every call starts one bridge process and reads the whole file, so
`limit` and `offset` are a window, not a cursor: paging re-reads the
file each time. For a large file pass `out` to write a tab-separated
table and read it with
[`read.delim()`](https://rdrr.io/r/utils/read.table.html).
