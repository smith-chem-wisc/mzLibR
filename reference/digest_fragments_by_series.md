# Fragment ions per product type

**Prefer this to `nrow(digest$fragments)` whenever the ion series
matter, which for ETD is always.** A bare total folds in the one extra
full-length `z-dot` per peptide, which is not a backbone-cleavage ion
between two residues.

## Usage

``` r
digest_fragments_by_series(digest)
```

## Arguments

- digest:

  An
  [`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  result.

## Details

The `zDot` series runs `1..length`, not `1..length-1`. The extra ion
numbered `length` is the whole peptide minus NH2 - the N-Ca cleavage at
residue 1 - which is **correct and deliberate**, not a defect. It is
absent when the peptide starts with proline.

Note also that `zDot` counts come in **below** `length` per peptide
rather than above it, because z-dot ions are suppressed N-terminal to
proline while the complementary `c` ions are not. On albumin that is 138
proline sites, 138 suppressed `z-dot` and **0** suppressed `c`
(smith-chem-wisc/mzLib#1110).

## Value

A data.frame of `product_type` and `n`, ordered by product type.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
digest_fragments_by_series(digest)
#>   product_type  n
#> 1            c 64
#> 2         zDot 63
```
