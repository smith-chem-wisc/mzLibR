# Only the peptidoforms carrying at least one modification

Only the peptidoforms carrying at least one modification

## Usage

``` r
digest_modified_peptides(digest)
```

## Arguments

- digest:

  An
  [`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  result.

## Value

The rows of `digest$peptides` with `modification_count > 0`.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
digest_modified_peptides(digest)[, c("full_sequence", "modification_count")]
#>                                                         full_sequence
#> 2 RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK[UniProt:N6-succinyllysine on K]
#>   modification_count
#> 2                  1
```
