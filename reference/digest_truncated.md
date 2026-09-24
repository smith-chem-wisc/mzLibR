# Whether a digest hit the isoform cap, and is therefore incomplete

A short answer and a truncated answer look identical from the outside.
Check this before treating a peptide list as exhaustive.

## Usage

``` r
digest_truncated(digest)
```

## Arguments

- digest:

  An
  [`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  result.

## Value

`TRUE` if any peptide hit `max_isoforms`.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
digest_truncated(digest)
#> [1] FALSE
```
