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

**Fixed charges are not double-counted.** A peptide whose modification
leaves a permanently charged residue - trimethyl-lysine gives a
quaternary ammonium, and UniProt records the delta as 43.054227, which
is C3H7 \*minus an electron\* - already carries that charge in its
`monoisotopic_mass`. Only `charge - fixed_charges` protons are added.
Adding a full complement would put a 2+ trimethylated peptide half a
Thomson high, on the most important histone modification there is.

A peptide with a fixed charge is observable at that charge with no
protonation at all, which is why `charge` may not be below it.

Note that \*\*fragments carry `neutral_mass` and deliberately have no
m/z\*\*: a c or z ion carries only the fixed charges within its own
span, and per-fragment charge accounting does not exist on this wire.

## Value

A numeric vector of m/z, one per row.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
peptide_mz(digest$peptides, charge = 2)
#> [1] 1782.434 2069.458
```
