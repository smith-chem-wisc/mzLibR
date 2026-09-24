# How many distinct base sequences a digest produced

`nrow(digest$peptides)` counts **peptidoforms** - one per
sequence-and-modification-placement. This counts distinct sequences. On
albumin at two modifications the two are **303** and **195**.

## Usage

``` r
digest_distinct_base_sequences(digest)
```

## Arguments

- digest:

  An
  [`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  result.

## Value

A single count.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
nrow(digest$peptides)                  # peptidoforms
#> [1] 2
digest_distinct_base_sequences(digest)  # distinct sequences
#> [1] 2
```
