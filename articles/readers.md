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
#> [1] 36
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
#>             17              2              4              7              6
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
| `spectral_match` | [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md) | identifications, with q-values and scores for mzIdentML |
| anything | [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md) | the format’s own fields |

MetaMorpheus and FlashLFQ quantification tables have readers of their
own, below.

### Spectra

Scan headers always; peak arrays only when asked for with
`peaks = TRUE`, because they are roughly the size of the file.

``` r

scans <- readers_read_spectra("sliced_ethcd.mzML", limit = 3)
scans
#> <mzlibr_scan_records> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\DataFiles\sliced_ethcd.mzML (MzML, Mzml)
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
#>   absent from this file (NA in every row): is_decoy
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

mzIdentML, which most search engines export, reads through the match
view with its q-value, rank and threshold - and with `scores = TRUE`,
one row per match and engine score:

``` r

matches <- readers_read_matches("PXD078927_msgf_1_1_0.mzid", limit = 3)
matches$records[, c("one_based_scan_number", "base_sequence", "q_value", "rank",
                    "pass_threshold")]
#>   one_based_scan_number  base_sequence   q_value rank pass_threshold
#> 1                 14316   HSNLNDATYQRT 0.0000000    1           TRUE
#> 2                 14316 KASAGQISVQPTFS 0.6588785    2           TRUE
#> 3                 14316   TRQYTADNLNSH 0.7383230    3           TRUE
scored <- readers_read_matches("PXD078927_msgf_1_1_0.mzid", limit = 1, scores = TRUE)
scored$records[, c("match_index", "score_name", "score_value")]
#>   match_index        score_name  score_value
#> 1           0    MS-GF:RawScore 1.150000e+02
#> 2           0 MS-GF:DeNovoScore 1.230000e+02
#> 3           0  MS-GF:SpecEValue 3.041116e-14
#> 4           0      MS-GF:EValue 8.393480e-12
#> 5           0      MS-GF:QValue 0.000000e+00
#> 6           0   MS-GF:PepQValue 0.000000e+00
#> 7           0      IsotopeError 1.000000e+00
```

Nothing in the match view is filtered: every candidate rank is a row,
and a column a format has no source for is `NA` in every row and named
in `absent_fields`. Casanovo is de novo, so it has no decoy label and no
q-value:

``` r

denovo <- readers_read_matches("Casanovo_5.0.0.mztab")
denovo$absent_fields
#> [1] "is_decoy"       "q_value"        "rank"           "pass_threshold"
```

[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
reads any recognised file into that format’s own fields, under mzLib’s
names, so a column can be looked up in the mzLib source:

``` r

prsms <- readers_read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv")
prsms$record_type
#> [1] "ToppicPrsm"
prsms$records[, c("one_based_scan_number", "base_sequence", "e_value")]
#>   one_based_scan_number base_sequence e_value
#> 1                   259        GYASDS  1e+300
#> 2                   270         THIGY  1e+300
#> 3                   354      SAECTKRF  1e+300
#> 4                   361      THASEKGY  1e+300
prsms$excluded_fields
#>                         field                      type
#> 1 alternative_identifications List<AlternativeToppicId>
#>                                                    reason verb
#> 1 a list of composite values has no faithful column shape   NA
```

`excluded_fields` names what could not become a column, with the reason,
so a missing column is never mistaken for a field the format does not
have.

## Quantification tables

mzLib 1.0.592 reads the tables MetaMorpheus and FlashLFQ write after
quantifying. Each is wide on disk, one column per sample, and arrives
long: one row per record per sample. `limit` and `offset` count records,
so one group’s samples are never split.

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

The table is unfiltered - decoys, contaminants and groups above 1% FDR
are all rows - and a blank intensity is `NA`, not zero.
[`readers_read_quantified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
reads the peptide table, where a FlashLFQ intensity of 0 is *not* a
measurement and `detection_type` says why; and
[`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
reads the per-site PTM occupancy, one row per group, sample, basis and
site:

``` r

sites <- readers_read_occupancy("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
head(sites$records[sites$records$basis == "count",
                   c("protein_group_name", "position", "modification", "numerator",
                     "denominator")])
#>    protein_group_name position                  modification numerator
#> 1              P68363      329              Deamidation on N         1
#> 3              P68363      329              Deamidation on N         1
#> 5              P68363      329              Deamidation on N         1
#> 7              P68363      329              Deamidation on N         2
#> 9              P68363      329              Deamidation on N         1
#> 11             P05141       52 N6,N6,N6-trimethyllysine on K         1
#>    denominator
#> 1            2
#> 3            2
#> 5            2
#> 7            4
#> 9            4
#> 11           1
```

## Many files in one call

Every reader has a `_many` form that hands the whole list to one bridge
process, which reads `threads` files at once and returns one long
table - the same table at any thread count. Never loop the one-file
function instead: each call starts a new process.

``` r

runs <- readers_read_spectra_many(c("sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"),
                                  threads = 2, on_error = "skip")
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

`records` starts with `source_index` - 1-based, so it indexes `files` -
and `source_path`. `files` has one row per input with everything a
one-file read reports about it, including the run’s instrument and
acquisition time for a spectra file; with `on_error = "skip"`, an input
that could not be read is a row with its error rather than a failed
call.

## Large files

Every call starts one bridge process and reads the whole file, so
`limit` and `offset` are a window, not a cursor: paging re-reads the
file each time. For a large file pass `out` to write a tab-separated
table and read it with
[`read.delim()`](https://rdrr.io/r/utils/read.table.html).
