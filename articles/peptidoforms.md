# Digesting and fragmenting peptidoforms

| question | function | mzLib |
|----|----|----|
| which peptidoforms does a protein give, with which fragments? | [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md) | `Loaders.LoadUniprot`, `Protein.Digest`, `PeptideWithSetModifications.Fragment` |
| which modified peptides, and how many distinct backbones? | [`digest_modified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_modified_peptides.md), [`digest_distinct_base_sequences()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_distinct_base_sequences.md) | the same digest |
| one ion series at a time | [`digest_fragments_by_series()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_fragments_by_series.md) | the same fragments |
| which annotated sites were not applied, and why? | [`census_excluded()`](https://smith-chem-wisc.github.io/mzLibR/reference/census_excluded.md), [`census_explain()`](https://smith-chem-wisc.github.io/mzLibR/reference/census_explain.md) | `ModificationCensus` |
| a peptide’s m/z at a charge | [`peptide_mz()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md) | none: the proton arithmetic is spelled out in [`?peptide_mz`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md) |

[`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
fetches a UniProt entry, applies its annotated modifications, digests it
and computes fragment ions for every peptidoform. The recording used
here is the whole digest of human serum albumin, P02768, at the
defaults.

``` r

digest <- peptidoform_fragments("P02768")
digest
#> <mzlibr_digest> P02768 ALBU_HUMAN (Homo sapiens)
#>   609 residues, trypsin|P, ETD, terminus Both
#>   min_length 7, 2 missed cleavages, max 2 modifications
#>   303 peptidoforms over 195 distinct sequences
#>   10728 fragments: c=5338, zDot=5390
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
#> 1             1            c               1     203.1059            0
#> 2             1         zDot               1     158.0930            0
#> 3             1            c               2     302.1743            0
#> 4             1         zDot               2     245.1250            0
#> 5             1            c               3     403.2220            0
#> 6             1         zDot               3     408.1883            0
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
#> 5              34                 8     <NA>  Phosphothreonine on T  79.96633
#> 6              35                 6     <NA>  Phosphothreonine on T  79.96633
#> 7              36                 5     <NA>     Phosphoserine on S  79.96633
#> 8              37                 6     <NA>  Phosphothreonine on T  79.96633
#> 9              37                 8     <NA>  Phosphothreonine on T  79.96633
#> 10             38                 5     <NA>     Phosphoserine on S  79.96633
#> 11             38                 8     <NA>  Phosphothreonine on T  79.96633
#> 12             39                 5     <NA>     Phosphoserine on S  79.96633
#> 13             39                 6     <NA>  Phosphothreonine on T  79.96633
#> 14             43                 5     <NA>     Phosphoserine on S  79.96633
#> 15             45                19     <NA> N6-succinyllysine on K 100.01604
#> 16             47                 9     <NA>   N6-methyllysine on K  14.01565
#> 17             56                 5     <NA>     Phosphoserine on S  79.96633
#> 18             58                 1     <NA>     Phosphoserine on S  79.96633
#> 19             63                17     <NA>     Phosphoserine on S  79.96633
#> 20             65                14     <NA>     Phosphoserine on S  79.96633
#> 21             66                 7     <NA>     Phosphoserine on S  79.96633
#> 22             67                 7     <NA>     Phosphoserine on S  79.96633
#> 23             67                14     <NA>     Phosphoserine on S  79.96633
#> 24             69                 1     <NA>     Phosphoserine on S  79.96633
#> 25             71                10     <NA>  Phosphothreonine on T  79.96633
#> 26             73                 2     <NA>  Phosphothreonine on T  79.96633
#> 27             89                 8     <NA> N6-succinyllysine on K 100.01604
#> 28             91                 6     <NA> N6-succinyllysine on K 100.01604
#> 29            101                16     <NA>     Phosphoserine on S  79.96633
#> 30            103                11     <NA>     Phosphoserine on S  79.96633
#> 31            119                 9     <NA>  Phosphothreonine on T  79.96633
#> 32            120                 7     <NA>  Phosphothreonine on T  79.96633
#> 33            121                 6     <NA>     Phosphoserine on S  79.96633
#> 34            122                 7     <NA>  Phosphothreonine on T  79.96633
#> 35            122                 9     <NA>  Phosphothreonine on T  79.96633
#> 36            123                 6     <NA>     Phosphoserine on S  79.96633
#> 37            123                 9     <NA>  Phosphothreonine on T  79.96633
#> 38            124                 6     <NA>     Phosphoserine on S  79.96633
#> 39            124                 7     <NA>  Phosphothreonine on T  79.96633
#> 40            126                 8     <NA>  Phosphothreonine on T  79.96633
#> 41            127                 6     <NA>  Phosphothreonine on T  79.96633
#> 42            128                 5     <NA>     Phosphoserine on S  79.96633
#> 43            129                 6     <NA>  Phosphothreonine on T  79.96633
#> 44            129                 8     <NA>  Phosphothreonine on T  79.96633
#> 45            130                 5     <NA>     Phosphoserine on S  79.96633
#> 46            130                 8     <NA>  Phosphothreonine on T  79.96633
#> 47            131                 5     <NA>     Phosphoserine on S  79.96633
#> 48            131                 6     <NA>  Phosphothreonine on T  79.96633
#> 49            133                 8     <NA> N6-succinyllysine on K 100.01604
#> 50            135                 4     <NA> N6-succinyllysine on K 100.01604
#> 51            142                14     <NA>     Phosphoserine on S  79.96633
#> 52            144                35     <NA> N6-succinyllysine on K 100.01604
#> 53            145                 5     <NA>     Phosphoserine on S  79.96633
#> 54            146                 5     <NA>     Phosphoserine on S  79.96633
#> 55            146                35     <NA> N6-succinyllysine on K 100.01604
#> 56            148                19     <NA> N6-succinyllysine on K 100.01604
#> 57            150                10     <NA>   N6-methyllysine on K  14.01565
#> 58            152                 9     <NA>   N6-methyllysine on K  14.01565
#> 59            158                 7     <NA> N6-succinyllysine on K 100.01604
#> 60            160                 4     <NA> N6-succinyllysine on K 100.01604
#> 61            169                 6     <NA>     Phosphoserine on S  79.96633
#> 62            171                 5     <NA>     Phosphoserine on S  79.96633
#> 63            173                 1     <NA>     Phosphoserine on S  79.96633
#> 64            177                38     <NA>     Phosphoserine on S  79.96633
#> 65            179                24     <NA>     Phosphoserine on S  79.96633
#> 66            180                17     <NA>     Phosphoserine on S  79.96633
#> 67            181                17     <NA>     Phosphoserine on S  79.96633
#> 68            181                24     <NA>     Phosphoserine on S  79.96633
#> 69            183                14     <NA>     Phosphoserine on S  79.96633
#> 70            184                 7     <NA>     Phosphoserine on S  79.96633
#> 71            185                 7     <NA>     Phosphoserine on S  79.96633
#> 72            185                14     <NA>     Phosphoserine on S  79.96633
#> 73            187                19     <NA>  Phosphothreonine on T  79.96633
#> 74            188                 1     <NA>     Phosphoserine on S  79.96633
#> 75            189                 1     <NA>     Phosphoserine on S  79.96633
#> 76            189                19     <NA>  Phosphothreonine on T  79.96633
#> 77            191                10     <NA>  Phosphothreonine on T  79.96633
#> 78            193                 2     <NA>  Phosphothreonine on T  79.96633
#> 79            210                10     <NA> N6-succinyllysine on K 100.01604
#> 80            212                 8     <NA> N6-succinyllysine on K 100.01604
#> 81            214                 6     <NA> N6-succinyllysine on K 100.01604
#> 82            223                33     <NA>     Phosphoserine on S  79.96633
#> 83            225                16     <NA>     Phosphoserine on S  79.96633
#> 84            227                11     <NA>     Phosphoserine on S  79.96633
#> 85            243                12     <NA>  Phosphothreonine on T  79.96633
#> 86            244                10     <NA>  Phosphothreonine on T  79.96633
#> 87            245                 9     <NA>     Phosphoserine on S  79.96633
#> 88            246                10     <NA>  Phosphothreonine on T  79.96633
#> 89            246                12     <NA>  Phosphothreonine on T  79.96633
#> 90            247                 9     <NA>     Phosphoserine on S  79.96633
#> 91            247                12     <NA>  Phosphothreonine on T  79.96633
#> 92            248                 9     <NA>     Phosphoserine on S  79.96633
#> 93            248                10     <NA>  Phosphothreonine on T  79.96633
#> 94            250                 9     <NA>  Phosphothreonine on T  79.96633
#> 95            251                 7     <NA>  Phosphothreonine on T  79.96633
#> 96            252                 6     <NA>     Phosphoserine on S  79.96633
#> 97            253                 7     <NA>  Phosphothreonine on T  79.96633
#> 98            253                 9     <NA>  Phosphothreonine on T  79.96633
#> 99            254                 6     <NA>     Phosphoserine on S  79.96633
#> 100           254                 9     <NA>  Phosphothreonine on T  79.96633
#> 101           255                 6     <NA>     Phosphoserine on S  79.96633
#> 102           255                 7     <NA>  Phosphothreonine on T  79.96633
#> 103           257                22     <NA> N6-succinyllysine on K 100.01604
#> 104           258                 8     <NA>  Phosphothreonine on T  79.96633
#> 105           259                 6     <NA>  Phosphothreonine on T  79.96633
#> 106           260                 5     <NA>     Phosphoserine on S  79.96633
#> 107           261                 8     <NA>  Phosphothreonine on T  79.96633
#> 108           261                22     <NA> N6-succinyllysine on K 100.01604
#> 109           262                 6     <NA>  Phosphothreonine on T  79.96633
#> 110           262                22     <NA> N6-succinyllysine on K 100.01604
#> 111           263                 6     <NA>  Phosphothreonine on T  79.96633
#> 112           263                 8     <NA>  Phosphothreonine on T  79.96633
#> 113           264                 5     <NA>     Phosphoserine on S  79.96633
#> 114           264                22     <NA> N6-succinyllysine on K 100.01604
#> 115           265                 5     <NA>     Phosphoserine on S  79.96633
#> 116           265                 8     <NA>  Phosphothreonine on T  79.96633
#> 117           266                 5     <NA>     Phosphoserine on S  79.96633
#> 118           266                 6     <NA>  Phosphothreonine on T  79.96633
#> 119           268                 8     <NA> N6-succinyllysine on K 100.01604
#> 120           270                 4     <NA> N6-succinyllysine on K 100.01604
#> 121           277                17     <NA>     Phosphoserine on S  79.96633
#> 122           279                44     <NA> N6-succinyllysine on K 100.01604
#> 123           280                14     <NA>     Phosphoserine on S  79.96633
#> 124           281                14     <NA>     Phosphoserine on S  79.96633
#> 125           281                44     <NA> N6-succinyllysine on K 100.01604
#> 126           283                35     <NA> N6-succinyllysine on K 100.01604
#> 127           284                 5     <NA>     Phosphoserine on S  79.96633
#> 128           285                 5     <NA>     Phosphoserine on S  79.96633
#> 129           285                35     <NA> N6-succinyllysine on K 100.01604
#> 130           287                19     <NA> N6-succinyllysine on K 100.01604
#> 131           289                13     <NA>   N6-methyllysine on K  14.01565
#> 132           291                10     <NA>   N6-methyllysine on K  14.01565
#> 133           293                 9     <NA>   N6-methyllysine on K  14.01565
#> 134           298                19     <NA> N6-succinyllysine on K 100.01604
#> 135           300                 7     <NA> N6-succinyllysine on K 100.01604
#> 136           302                 4     <NA> N6-succinyllysine on K 100.01604
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
#> [1] 303
digest_distinct_base_sequences(digest)
#> [1] 195
nrow(digest_modified_peptides(digest))
#> [1] 108
head(digest_modified_peptides(digest)[, c("full_sequence", "modification_count")], 3)
#>                                  full_sequence modification_count
#> 6    TCVADES[UniProt:Phosphoserine on S]AENCDK                  1
#> 8        S[UniProt:Phosphoserine on S]LHTLFGDK                  1
#> 11 ET[UniProt:Phosphothreonine on T]YGEMADCCAK                  1
```

## Fragment series

ETD, the default, produces `c` and `z-dot` ions; HCD and CID produce `b`
and `y`.

``` r

digest_fragments_by_series(digest)
#>   product_type    n
#> 1            c 5338
#> 2         zDot 5390
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

Two defaults also decide which peptides exist at all: `min_length = 7`
drops everything shorter, and mzLib’s `"trypsin|P"` *applies* the
proline rule - the reverse of the MaxQuant and Mascot spelling.
[`?peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
gives the numbers.

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
