# Digesting and fragmenting peptidoforms

[`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
fetches a UniProt entry, applies its annotated modifications, digests it
and computes fragment ions for every peptidoform. The recording used
here is human serum albumin, P02768, with at most one modification per
peptide.

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
digest
#> <mzlibr_digest> P02768 ALBU_HUMAN (Homo sapiens)
#>   609 residues, trypsin|P, ETD, terminus Both
#>   min_length 7, 2 missed cleavages, max 1 modifications
#>   2 peptidoforms over 2 distinct sequences
#>   127 fragments: c=64, zDot=63
#>   ! 14 of 38 annotated modifications applied. See ?census_explain.
```

## Three data.frames that join on `peptide_index`

``` r

head(digest$peptides[, c("peptide_index", "base_sequence", "monoisotopic_mass",
                         "modification_count")])
#>   peptide_index                       base_sequence monoisotopic_mass
#> 1             1     ALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK          3562.853
#> 2             2 RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK          4136.902
#>   modification_count
#> 1                  0
#> 2                  1
head(digest$fragments)
#>   peptide_index product_type fragment_number neutral_mass neutral_loss
#> 1             1            c               1     88.06366            0
#> 2             1         zDot               1    130.08680            0
#> 3             1            c               2    201.14773            0
#> 4             1         zDot               2    201.12392            0
#> 5             1            c               3    300.21614            0
#> 6             1         zDot               3    348.19233            0
#>   residue_position
#> 1                1
#> 2               31
#> 3                2
#> 4               30
#> 5                3
#> 6               29
digest$modifications
#>   peptide_index one_based_residue terminus                     id    mass
#> 1             2                35     <NA> N6-succinyllysine on K 100.016
#>   formal_charge
#> 1             0
```

`monoisotopic_mass` and `neutral_mass` are **neutral** masses in
daltons, not m/z.
[`peptide_mz()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md)
converts at a charge:

``` r

peptide_mz(digest$peptides, charge = 2)
#> [1] 1782.434 2069.458
```

## Peptidoforms are not sequences

`peptides` has one row per sequence *and* modification placement, so its
row count is larger than the number of distinct sequences. Both are
legitimate answers to “how many peptides”, and they are not
interchangeable:

``` r

nrow(digest$peptides)
#> [1] 2
digest_distinct_base_sequences(digest)
#> [1] 2
digest_modified_peptides(digest)[, c("full_sequence", "modification_count")]
#>                                                         full_sequence
#> 2 RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK[UniProt:N6-succinyllysine on K]
#>   modification_count
#> 2                  1
```

## Fragment series

ETD, the default, produces `c` and `z-dot` ions; HCD and CID produce `b`
and `y`.

``` r

digest_fragments_by_series(digest)
#>   product_type  n
#> 1            c 64
#> 2         zDot 63
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
