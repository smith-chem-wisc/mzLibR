# Digesting and fragmenting peptidoforms

| question | function | mzLib |
|----|----|----|
| which peptidoforms does a protein give, with which fragments? | [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md) | `Loaders.LoadUniprot`, `Protein.Digest`, `PeptideWithSetModifications.Fragment` |
| which modified peptides, and how many distinct backbones? | [`digest_modified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_modified_peptides.md), [`digest_distinct_base_sequences()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_distinct_base_sequences.md) | the same digest |
| one ion series at a time | [`digest_fragments_by_series()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_fragments_by_series.md) | the same fragments |
| which annotated sites were not applied, and why? | [`census_excluded()`](https://smith-chem-wisc.github.io/mzLibR/reference/census_excluded.md), [`census_explain()`](https://smith-chem-wisc.github.io/mzLibR/reference/census_explain.md) | `ModificationCensus` |
| a peptide’s m/z at a charge | [`peptide_mz()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md) | none: the proton arithmetic is spelled out in [`?peptide_mz`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md) |
| a MetaMorpheus full sequence in Unimod accessions | [`peptidoform_convert()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_convert.md) | `SequenceConversionService.Convert`, Unimod serializer |
| the same sequence in ProForma | `peptidoform_convert(target = "ProForma")` | `SequenceConversionService.Convert`, ProForma serializer |
| which sequences mzLib could not convert, and why | `records$status`, `records$incompatible_items` | `ConversionWarnings` |
| which notations mzLib can read and write | `source_formats`, `target_formats` | `SequenceConversionService.AvailableSourceFormats`, `AvailableTargetFormats` |

[`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
fetches a UniProt entry, applies its annotated modifications, digests it
and computes fragment ions for every peptidoform. The recording used
here is the whole digest of human serum albumin, P02768, at mzLib’s
defaults: an argument left out is not sent, and mzLib chooses. The
digest says what was used.

``` r

digest <- peptidoform_fragments("P02768")
digest
#> <mzlibr_digest> P02768 ALBU_HUMAN (Homo sapiens)
#>   609 residues, trypsin, HCD, terminus Both
#>   min_length mzLib's default, missed cleavages mzLib's default, max 2 modifications
#>   310 peptidoforms over 202 distinct sequences
#>   10488 fragments: b=5244, y=5244
#>   ! 14 of 38 annotated modifications applied. See ?census_explain.
```

## Three data.frames that join on `peptide_index`

``` r

head(digest$peptides[, c("peptide_index", "base_sequence", "monoisotopic_mass",
                         "modification_count")])
#>   peptide_index         base_sequence monoisotopic_mass modification_count
#> 1             1     WVTFISLLFLFSSAYSR         2036.0771                  0
#> 2             2              DLGEENFK          950.4345                  0
#> 3             3 ALVLIAFAQYLQQCPFEDHVK         2432.2562                  0
#> 4             4            LVNEVTEFAK         1148.6077                  0
#> 5             5         TCVADESAENCDK         1383.5282                  0
#> 6             6         TCVADESAENCDK         1463.4946                  1
head(digest$fragments)
#>   peptide_index product_type fragment_number neutral_mass neutral_loss
#> 1             1            b               1     186.0793            0
#> 2             1            y               1     174.1117            0
#> 3             1            b               2     285.1477            0
#> 4             1            y               2     261.1437            0
#> 5             1            b               3     386.1954            0
#> 6             1            y               3     424.2070            0
#>   residue_position
#> 1                1
#> 2               17
#> 3                2
#> 4               16
#> 5                3
#> 6               15
digest$modifications
#>     peptide_index one_based_residue terminus                     id      mass
#> 1               6                 7     <NA>     Phosphoserine on S  79.96633
#> 2               8                 1     <NA>     Phosphoserine on S  79.96633
#> 3              11                 2     <NA>  Phosphothreonine on T  79.96633
#> 4              23                11     <NA>     Phosphoserine on S  79.96633
#> 5              33                 8     <NA>  Phosphothreonine on T  79.96633
#> 6              34                 6     <NA>  Phosphothreonine on T  79.96633
#> 7              35                 5     <NA>     Phosphoserine on S  79.96633
#> 8              36                 6     <NA>  Phosphothreonine on T  79.96633
#> 9              36                 8     <NA>  Phosphothreonine on T  79.96633
#> 10             37                 5     <NA>     Phosphoserine on S  79.96633
#> 11             37                 8     <NA>  Phosphothreonine on T  79.96633
#> 12             38                 5     <NA>     Phosphoserine on S  79.96633
#> 13             38                 6     <NA>  Phosphothreonine on T  79.96633
#> 14             42                 4     <NA>     Phosphoserine on S  79.96633
#> 15             44                19     <NA> N6-succinyllysine on K 100.01604
#> 16             46                 9     <NA>   N6-methyllysine on K  14.01565
#> 17             55                 5     <NA>     Phosphoserine on S  79.96633
#> 18             57                 1     <NA>     Phosphoserine on S  79.96633
#> 19             62                17     <NA>     Phosphoserine on S  79.96633
#> 20             64                14     <NA>     Phosphoserine on S  79.96633
#> 21             65                 7     <NA>     Phosphoserine on S  79.96633
#> 22             66                 7     <NA>     Phosphoserine on S  79.96633
#> 23             66                14     <NA>     Phosphoserine on S  79.96633
#> 24             68                 1     <NA>     Phosphoserine on S  79.96633
#> 25             70                10     <NA>  Phosphothreonine on T  79.96633
#> 26             72                 2     <NA>  Phosphothreonine on T  79.96633
#> 27             89                 8     <NA> N6-succinyllysine on K 100.01604
#> 28             91                 6     <NA> N6-succinyllysine on K 100.01604
#> 29            101                16     <NA>     Phosphoserine on S  79.96633
#> 30            103                11     <NA>     Phosphoserine on S  79.96633
#> 31            121                 9     <NA>  Phosphothreonine on T  79.96633
#> 32            122                 7     <NA>  Phosphothreonine on T  79.96633
#> 33            123                 6     <NA>     Phosphoserine on S  79.96633
#> 34            124                 7     <NA>  Phosphothreonine on T  79.96633
#> 35            124                 9     <NA>  Phosphothreonine on T  79.96633
#> 36            125                 6     <NA>     Phosphoserine on S  79.96633
#> 37            125                 9     <NA>  Phosphothreonine on T  79.96633
#> 38            126                 6     <NA>     Phosphoserine on S  79.96633
#> 39            126                 7     <NA>  Phosphothreonine on T  79.96633
#> 40            128                 8     <NA>  Phosphothreonine on T  79.96633
#> 41            129                 6     <NA>  Phosphothreonine on T  79.96633
#> 42            130                 5     <NA>     Phosphoserine on S  79.96633
#> 43            131                 6     <NA>  Phosphothreonine on T  79.96633
#> 44            131                 8     <NA>  Phosphothreonine on T  79.96633
#> 45            132                 5     <NA>     Phosphoserine on S  79.96633
#> 46            132                 8     <NA>  Phosphothreonine on T  79.96633
#> 47            133                 5     <NA>     Phosphoserine on S  79.96633
#> 48            133                 6     <NA>  Phosphothreonine on T  79.96633
#> 49            135                 8     <NA> N6-succinyllysine on K 100.01604
#> 50            137                 4     <NA> N6-succinyllysine on K 100.01604
#> 51            145                 5     <NA>     Phosphoserine on S  79.96633
#> 52            147                34     <NA> N6-succinyllysine on K 100.01604
#> 53            148                 4     <NA>     Phosphoserine on S  79.96633
#> 54            149                 4     <NA>     Phosphoserine on S  79.96633
#> 55            149                34     <NA> N6-succinyllysine on K 100.01604
#> 56            151                19     <NA> N6-succinyllysine on K 100.01604
#> 57            153                10     <NA>   N6-methyllysine on K  14.01565
#> 58            155                 9     <NA>   N6-methyllysine on K  14.01565
#> 59            160                 7     <NA> N6-succinyllysine on K 100.01604
#> 60            162                 4     <NA> N6-succinyllysine on K 100.01604
#> 61            171                 6     <NA>     Phosphoserine on S  79.96633
#> 62            173                 5     <NA>     Phosphoserine on S  79.96633
#> 63            175                 1     <NA>     Phosphoserine on S  79.96633
#> 64            179                38     <NA>     Phosphoserine on S  79.96633
#> 65            181                24     <NA>     Phosphoserine on S  79.96633
#> 66            182                17     <NA>     Phosphoserine on S  79.96633
#> 67            183                17     <NA>     Phosphoserine on S  79.96633
#> 68            183                24     <NA>     Phosphoserine on S  79.96633
#> 69            185                14     <NA>     Phosphoserine on S  79.96633
#> 70            186                 7     <NA>     Phosphoserine on S  79.96633
#> 71            187                 7     <NA>     Phosphoserine on S  79.96633
#> 72            187                14     <NA>     Phosphoserine on S  79.96633
#> 73            189                19     <NA>  Phosphothreonine on T  79.96633
#> 74            190                 1     <NA>     Phosphoserine on S  79.96633
#> 75            191                 1     <NA>     Phosphoserine on S  79.96633
#> 76            191                19     <NA>  Phosphothreonine on T  79.96633
#> 77            193                10     <NA>  Phosphothreonine on T  79.96633
#> 78            195                 2     <NA>  Phosphothreonine on T  79.96633
#> 79            213                10     <NA> N6-succinyllysine on K 100.01604
#> 80            215                 8     <NA> N6-succinyllysine on K 100.01604
#> 81            217                 6     <NA> N6-succinyllysine on K 100.01604
#> 82            226                33     <NA>     Phosphoserine on S  79.96633
#> 83            228                16     <NA>     Phosphoserine on S  79.96633
#> 84            230                11     <NA>     Phosphoserine on S  79.96633
#> 85            248                12     <NA>  Phosphothreonine on T  79.96633
#> 86            249                10     <NA>  Phosphothreonine on T  79.96633
#> 87            250                 9     <NA>     Phosphoserine on S  79.96633
#> 88            251                10     <NA>  Phosphothreonine on T  79.96633
#> 89            251                12     <NA>  Phosphothreonine on T  79.96633
#> 90            252                 9     <NA>     Phosphoserine on S  79.96633
#> 91            252                12     <NA>  Phosphothreonine on T  79.96633
#> 92            253                 9     <NA>     Phosphoserine on S  79.96633
#> 93            253                10     <NA>  Phosphothreonine on T  79.96633
#> 94            255                 9     <NA>  Phosphothreonine on T  79.96633
#> 95            256                 7     <NA>  Phosphothreonine on T  79.96633
#> 96            257                 6     <NA>     Phosphoserine on S  79.96633
#> 97            258                 7     <NA>  Phosphothreonine on T  79.96633
#> 98            258                 9     <NA>  Phosphothreonine on T  79.96633
#> 99            259                 6     <NA>     Phosphoserine on S  79.96633
#> 100           259                 9     <NA>  Phosphothreonine on T  79.96633
#> 101           260                 6     <NA>     Phosphoserine on S  79.96633
#> 102           260                 7     <NA>  Phosphothreonine on T  79.96633
#> 103           262                22     <NA> N6-succinyllysine on K 100.01604
#> 104           263                 8     <NA>  Phosphothreonine on T  79.96633
#> 105           264                 6     <NA>  Phosphothreonine on T  79.96633
#> 106           265                 5     <NA>     Phosphoserine on S  79.96633
#> 107           266                 8     <NA>  Phosphothreonine on T  79.96633
#> 108           266                22     <NA> N6-succinyllysine on K 100.01604
#> 109           267                 6     <NA>  Phosphothreonine on T  79.96633
#> 110           267                22     <NA> N6-succinyllysine on K 100.01604
#> 111           268                 6     <NA>  Phosphothreonine on T  79.96633
#> 112           268                 8     <NA>  Phosphothreonine on T  79.96633
#> 113           269                 5     <NA>     Phosphoserine on S  79.96633
#> 114           269                22     <NA> N6-succinyllysine on K 100.01604
#> 115           270                 5     <NA>     Phosphoserine on S  79.96633
#> 116           270                 8     <NA>  Phosphothreonine on T  79.96633
#> 117           271                 5     <NA>     Phosphoserine on S  79.96633
#> 118           271                 6     <NA>  Phosphothreonine on T  79.96633
#> 119           273                 8     <NA> N6-succinyllysine on K 100.01604
#> 120           275                 4     <NA> N6-succinyllysine on K 100.01604
#> 121           283                14     <NA>     Phosphoserine on S  79.96633
#> 122           285                35     <NA> N6-succinyllysine on K 100.01604
#> 123           286                 5     <NA>     Phosphoserine on S  79.96633
#> 124           287                 5     <NA>     Phosphoserine on S  79.96633
#> 125           287                35     <NA> N6-succinyllysine on K 100.01604
#> 126           289                34     <NA> N6-succinyllysine on K 100.01604
#> 127           290                 4     <NA>     Phosphoserine on S  79.96633
#> 128           291                 4     <NA>     Phosphoserine on S  79.96633
#> 129           291                34     <NA> N6-succinyllysine on K 100.01604
#> 130           293                19     <NA> N6-succinyllysine on K 100.01604
#> 131           295                13     <NA>   N6-methyllysine on K  14.01565
#> 132           297                10     <NA>   N6-methyllysine on K  14.01565
#> 133           299                 9     <NA>   N6-methyllysine on K  14.01565
#> 134           305                19     <NA> N6-succinyllysine on K 100.01604
#> 135           307                 7     <NA> N6-succinyllysine on K 100.01604
#> 136           309                 4     <NA> N6-succinyllysine on K 100.01604
#>     formal_charge
#> 1               0
#> 2               0
#> 3               0
#> 4               0
#> 5               0
#> 6               0
#> 7               0
#> 8               0
#> 9               0
#> 10              0
#> 11              0
#> 12              0
#> 13              0
#> 14              0
#> 15              0
#> 16              0
#> 17              0
#> 18              0
#> 19              0
#> 20              0
#> 21              0
#> 22              0
#> 23              0
#> 24              0
#> 25              0
#> 26              0
#> 27              0
#> 28              0
#> 29              0
#> 30              0
#> 31              0
#> 32              0
#> 33              0
#> 34              0
#> 35              0
#> 36              0
#> 37              0
#> 38              0
#> 39              0
#> 40              0
#> 41              0
#> 42              0
#> 43              0
#> 44              0
#> 45              0
#> 46              0
#> 47              0
#> 48              0
#> 49              0
#> 50              0
#> 51              0
#> 52              0
#> 53              0
#> 54              0
#> 55              0
#> 56              0
#> 57              0
#> 58              0
#> 59              0
#> 60              0
#> 61              0
#> 62              0
#> 63              0
#> 64              0
#> 65              0
#> 66              0
#> 67              0
#> 68              0
#> 69              0
#> 70              0
#> 71              0
#> 72              0
#> 73              0
#> 74              0
#> 75              0
#> 76              0
#> 77              0
#> 78              0
#> 79              0
#> 80              0
#> 81              0
#> 82              0
#> 83              0
#> 84              0
#> 85              0
#> 86              0
#> 87              0
#> 88              0
#> 89              0
#> 90              0
#> 91              0
#> 92              0
#> 93              0
#> 94              0
#> 95              0
#> 96              0
#> 97              0
#> 98              0
#> 99              0
#> 100             0
#> 101             0
#> 102             0
#> 103             0
#> 104             0
#> 105             0
#> 106             0
#> 107             0
#> 108             0
#> 109             0
#> 110             0
#> 111             0
#> 112             0
#> 113             0
#> 114             0
#> 115             0
#> 116             0
#> 117             0
#> 118             0
#> 119             0
#> 120             0
#> 121             0
#> 122             0
#> 123             0
#> 124             0
#> 125             0
#> 126             0
#> 127             0
#> 128             0
#> 129             0
#> 130             0
#> 131             0
#> 132             0
#> 133             0
#> 134             0
#> 135             0
#> 136             0
```

`monoisotopic_mass` and `neutral_mass` are **neutral** masses in
daltons, not m/z.
[`peptide_mz()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md)
converts at a charge:

``` r

head(peptide_mz(digest$peptides, charge = 2))
#> [1] 1019.0458  476.2245 1217.1354  575.3111  692.7714  732.7546
```

## Peptidoforms are not sequences

`peptides` has one row per sequence *and* modification placement, so its
row count is larger than the number of distinct sequences. Both are
legitimate answers to “how many peptides”, and they are not
interchangeable:

``` r

nrow(digest$peptides)
#> [1] 310
digest_distinct_base_sequences(digest)
#> [1] 202
nrow(digest_modified_peptides(digest))
#> [1] 108
head(digest_modified_peptides(digest)[, c("full_sequence", "modification_count")], 3)
#>                                  full_sequence modification_count
#> 6    TCVADES[UniProt:Phosphoserine on S]AENCDK                  1
#> 8        S[UniProt:Phosphoserine on S]LHTLFGDK                  1
#> 11 ET[UniProt:Phosphothreonine on T]YGEMADCCAK                  1
```

## Fragment series

HCD and CID produce `b` and `y` ions; ETD produces `c` and `z-dot`.
mzLib names no default dissociation type, so with none given the bridge
uses HCD until it does (mzLib#1488).

``` r

digest_fragments_by_series(digest)
#>   product_type    n
#> 1            b 5244
#> 2            y 5244
```

## Which modifications were applied, and which were not

mzLib loads only some UniProt feature types. The census says how many
annotations were applied, and why the rest were not:

``` r

census_explain(digest)
#> [1] "14 of 38 annotated modifications were applied, across 14 residue positions. Excluded by type: 24 x glycosylation site - mzLib loads only 'modified residue' and 'lipid moiety-binding region' annotations, so these were dropped on feature type alone. The exclusion is usually right: a glycation or glycosylation annotation describes a labile, heterogeneous adduct, so assigning it one exact mass and a clean fragment ladder would invent a species you cannot observe. But the reason is not reported, and the qualifier is not read - some annotations are marked 'in vitro' and some exist only in disease variants, which are different grounds for exclusion needing different judgements from you. Read the annotations on the UniProt entry before concluding anything about a specific site; this census can only tell you the count (smith-chem-wisc/mzLib#1112)."
head(census_excluded(digest))
#> [1] 24
```

## Before treating the list as complete

`max_isoforms` caps the peptidoforms generated per peptide, and the cap
truncates silently.
[`digest_truncated()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_truncated.md)
says whether it was reached:

``` r

digest_truncated(digest)
#> [1] FALSE
```

Two of mzLib’s defaults also decide which peptides exist at all: its
`min_length` drops everything shorter, and its protease, `"trypsin"`,
cuts before proline, where `"trypsin|P"` does not - the reverse of the
MaxQuant and Mascot spelling.
[`?peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
gives the numbers. To get what this page showed before mzLibR left the
defaults to mzLib, pass `protease = "trypsin|P", dissociation = "ETD"`.

## Converting full sequences to Unimod or ProForma

A MetaMorpheus result names each modification the way mzLib’s databases
do: `[UniProt:N-acetylserine on S]`, `[Common Variable:Oxidation on M]`.
Most other tools want a Unimod accession instead, such as `[UNIMOD:1]`.
[`peptidoform_convert()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_convert.md)
hands each sequence to mzLib’s `SequenceConversionService` and returns
one row per input, in order. mzLibR maps no modification itself; every
output is mzLib’s.

This example uses `BottomUpExample.psmtsv`, a MetaMorpheus search result
that ships with mzLib’s tests. Read its full sequences, then convert
them:

``` r

psms <- readers_read_results("BottomUpExample.psmtsv")
unimod <- peptidoform_convert(psms$records$full_sequence)
unimod
#> <mzlibr_sequence_conversions> 8 sequences, mzLib to Unimod (ReturnNull)
#>   8 converted, 0 with warnings, 0 failed
head(unimod$records[, c("status", "output")], 5)
#>      status                                 output
#> 1 converted            YPIEH[UNIMOD:34]GIVTNWDDMEK
#> 2 converted                     VAPEEHPVLLTEAPLNPK
#> 3 converted     AYHEQLSVAEITNAC[UNIMOD:4]FEPANQMVK
#> 4 converted         VEDMAELTC[UNIMOD:4]LNEASVLHNLK
#> 5 converted YPIEH[UNIMOD:34]GIVTNWDDM[UNIMOD:35]EK
unimod$converted_count == unimod$record_count
#> [1] TRUE
```

UniProt’s tele-methylhistidine became `UNIMOD:34` (Methyl).
Carbamidomethyl became `UNIMOD:4`, and oxidation became `UNIMOD:35`.
Unmodified sequences pass through unchanged.

Every result lists the notations mzLib has registered, so you never need
to guess a name:

``` r

unimod$source_formats
#> [1] "MassShift" "Modomics"  "ProForma"  "mzLib"
unimod$target_formats
#> [1] "Chronologer" "Essential"   "MassShift"   "ProForma"    "Unimod"     
#> [6] "mzLib"
```

### ProForma gives the same accessions

Ask for ProForma and this file gives the same accessions. The UniProt
modification is `UNIMOD:34` here too:

``` r

proforma <- peptidoform_convert(psms$records$full_sequence, target = "ProForma")
proforma$records$output[c(1, 5)]
#> [1] "YPIEH[UNIMOD:34]GIVTNWDDMEK"           
#> [2] "YPIEH[UNIMOD:34]GIVTNWDDM[UNIMOD:35]EK"
identical(proforma$records$output, unimod$records$output)
#> [1] TRUE
```

That is new with mzLib#1403. Before it, mzLib’s ProForma serializer
looked modifications up only in MetaMorpheus’s own list, which has no
UniProt entries, so it wrote `[UniProt:Tele-methylhistidine on H]` back
by name (mzLib#1401). It now falls back to every list mzLib loads, as
the Unimod serializer does. mzLibR never patched around the gap, so the
fix reached it, pyMzLib, mzLibRust and MetaMorpheus at once. The
`pro_forma` column that
[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
computes for a `.psmtsv` written without one comes from the same
serializer, so it gained the accessions too.

The two targets still differ on a modification that no list mzLib loads
knows. ProForma writes it back under its own name and the row says
`converted`; Unimod fails the row, or drops the modification, depending
on `mode`. In ProForma output, treat any bracket that is not a `UNIMOD:`
term as unresolved.

### When mzLib cannot convert a sequence

What happens to a modification the target cannot write depends on
`mode`, which is mzLib’s `SequenceConversionHandlingMode`:

| `mode` | the row | `output` |
|----|----|----|
| `"ReturnNull"` (default) | `failed` | `NA` |
| `"RemoveIncompatibleElements"` | `converted_with_warnings` | the sequence without that modification |
| `"UsePrimarySequence"` | `converted_with_warnings` | the sequence without that modification |
| `"ThrowException"` | no rows: the call raises an `mzlib_usage_error` naming the first such sequence | none |

`status` is mzLib’s own verdict: `converted` (an output, and nothing
recorded against it), `converted_with_warnings` (an output, but mzLib
noted something, such as a dropped modification or a skipped character)
or `failed` (no output).

These four sequences include one modification that has no Unimod
accession:

``` r

result <- peptidoform_convert(c(
  "[UniProt:N-acetylserine on S]SEQK",
  "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
  "PEPM[Common Variable:Oxidation on M]K",
  "PEPK[Made Up:Not a modification on K]R"
))
result$records$output
#> [1] "[UNIMOD:1]SEQK"   "PEPK[UNIMOD:36]R" "PEPM[UNIMOD:35]K" NA
bad <- result$records[result$records$status != "converted", ]
bad[, c("input", "status", "failure_reason")]
#>                                    input status failure_reason
#> 4 PEPK[Made Up:Not a modification on K]R failed           <NA>
bad$incompatible_items[[1]]
#> [1] "Made Up:Not a modification on K @3(K)"
```

`failure_reason` is `NA` here. Under `"ReturnNull"`, mzLib’s Unimod
serializer records the incompatible modification but no reason code, so
read `incompatible_items` to see what failed.

**Split ambiguous sequences first.** When MetaMorpheus cannot tell
candidates apart, it joins their full sequences with `|`. mzLib’s parser
does not refuse this (mzLib#1405): it skips each `|` with a warning and
joins the candidates into one sequence. The row is
`converted_with_warnings`, but the output is not a real peptide. Split
on `|` before you convert.

[`?peptidoform_convert`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_convert.md)
lists every parameter, field and caveat.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1002/pmic.200300744](https://doi.org/10.1002/pmic.200300744):
  Unimod, the modification database whose accessions the Unimod target
  writes (Creasy and Cottrell 2004, Proteomics).
- [doi:10.1021/acs.jproteome.1c00771](https://doi.org/10.1021/acs.jproteome.1c00771):
  ProForma 2.0, the ProForma target’s notation (LeDuc et al. 2022, J
  Proteome Res).
- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, the source of the entry, its sequence and its annotated
  modifications (UniProt Consortium 2025, NAR).
