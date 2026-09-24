# Annotated features that could not be used

`annotated - applied`. See
[`census_explain`](https://smith-chem-wisc.github.io/mzLibR/reference/census_explain.md)
for why they were excluded, and why that exclusion is correct.

## Usage

``` r
census_excluded(digest)
```

## Arguments

- digest:

  An
  [`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  result, or its `census`.

## Value

A single count.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
census_excluded(digest)
#> [1] 24
```
