# A numeric covariate of a model, as mzLib's CovariateSpec

Names a sample-sheet column for
[`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
to enter as `(value - centre) / scale`. Leave a field `NULL` and nothing
is sent for it, so mzLib's own default applies.

## Usage

``` r
covariate_spec(name, unit, scale = NULL, scaled_unit = NULL, centre = NULL)
```

## Arguments

- name:

  The sample sheet column holding each sample's value, in `unit`.

- unit:

  The stored unit, e.g. `"years"`.

- scale:

  Stored units per scaled unit, e.g. `10` years per decade. `NULL`:
  mzLib's 1.

- scaled_unit:

  The unit a slope is per, e.g. `"decade"`. `NULL`: `scale` stored units
  (`"years"` at scale 1, `"10 years"` at scale 10).

- centre:

  The value subtracted first, in `unit`. `NULL`: the stratum's mean.

## Value

An `mzlibr_covariate_spec`: a list of `name`, `unit`, `scale`,
`scaled_unit` and `centre`, each of the last three `NA` when not given
(mzLib's default applies).

## See also

[`diff_design`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md),
[`factor_spec`](https://smith-chem-wisc.github.io/mzLibR/reference/factor_spec.md)

## Examples

``` r

covariate_spec("age", "years", scale = 10, scaled_unit = "decade", centre = 50)
#> <mzlibr_covariate_spec> age in years, scale 10 per decade, centred at 50
covariate_spec("bmi", "kg/m2")
#> <mzlibr_covariate_spec> bmi in kg/m2, centred at the stratum's mean
```
