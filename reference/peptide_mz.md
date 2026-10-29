# The m/z of intact peptides at a given charge

Two conventions are handled explicitly, because getting either wrong is
invisible in the answer.

## Usage

``` r
peptide_mz(peptides, charge)
```

## Arguments

- peptides:

  Rows of a `digest$peptides` data.frame.

- charge:

  Total charge, a positive whole number.

## Details

**The proton mass, not the hydrogen atom.** 1.00727646677 against
1.007825 - a difference of 0.55 mDa, or 1.1 ppm at m/z 500, which on an
Orbitrap is a match versus a miss.

**Permanently charged modifications are not double-counted.**
Trimethylation of a lysine epsilon-amine gives a quaternary ammonium,
which carries its charge with no proton. UniProt records the delta as
43.054227 (C3H7 \*minus an electron\*); mzLib loads it one proton
lighter, at 42.046951 (C3H6), so adding `charge` protons to
`monoisotopic_mass` gives the right m/z. Adding a full complement to
UniProt's mass would put a 2+ trimethylated peptide half a Thomson high.

`fixed_charges` is deprecated and always 0 from the pyMzLib 0.7.0
bridge, so the formula below adds `charge` protons. It still subtracts
`fixed_charges`, which keeps the m/z right for a digest from an older
bridge, whose `monoisotopic_mass` carried UniProt's charged mass: the
m/z is the same from either.

Fragments carry `neutral_mass`, charged modifications included at their
corrected mass, so a fragment's m/z at charge z is
`(neutral_mass + z * 1.00727646677) / z`.

## Value

A numeric vector of m/z, one per row.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
peptide_mz(digest$peptides, charge = 2)
#> [1] 1782.434 2069.458
```
