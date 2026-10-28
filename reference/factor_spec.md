# A condition factor of a model, as mzLib's FactorSpec

Names a sample-sheet column for
[`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
to treat-code against a reference level.

## Usage

``` r
factor_spec(name, reference_level = NULL)
```

## Arguments

- name:

  The sample sheet column holding each sample's level.

- reference_level:

  The level every other level is compared against, or `NULL` to compare
  every pair of levels (mzLib then warns, because the comparisons are
  many).

## Value

An `mzlibr_factor_spec`: a list of `name` and `reference_level` (`NA`
when none).

## See also

[`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md),
[`covariate_spec`](https://smith-chem-wisc.github.io/mzLibR/reference/covariate_spec.md)

## Examples

``` r

factor_spec("malat_load", "500ng")
#> <mzlibr_factor_spec> malat_load against 500ng
factor_spec("condition")
#> <mzlibr_factor_spec> condition (no reference: every pair of levels)
```
