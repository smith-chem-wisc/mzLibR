# Retention times in minutes, whatever unit the format wrote

The conversion you would otherwise write by hand, using
`retention_time_unit`.

## Usage

``` r
readers_retention_time_in_minutes(records, column = NULL)
```

## Arguments

- records:

  A
  [`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md),
  [`readers_read_features`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md),
  [`readers_read_spectra`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  or
  [`readers_read_quantified_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
  result, or the `_many` form of one, whose rows are converted by each
  file's own unit.

- column:

  Which retention-time column to convert. Defaults to `retention_time`
  for a record or scan table and to `retention_time_start` for a feature
  table, which has two. Name the other explicitly —
  `"retention_time_end"` — to convert it.

## Details

\*\*Raises when the unit is `"unknown"` rather than guessing.\*\* A
silently unconverted time axis is the specific mistake this module
exists to prevent: TopPIC still writes seconds while MetaMorpheus and
MSFragger write minutes (mzLib normalises MSFragger since PR \#1116, but
not TopPIC), and a 60x error in a retention-time comparison looks like a
chromatography problem rather than a units problem.

## Value

A numeric vector of retention times in minutes.

## When it raises

Whenever the unit is `"unknown"`, which is not hypothetical: it is the
honest answer for TopFD `_ms1.feature`, where seconds became minutes at
v1.7.0 without the file type changing. Guessing there is what mzLib's
own deconvolution code does, and what this deliberately does not.

## Examples

``` r

psms <- readers_read_results("FraggerPsm_FragPipev21.1_psm.tsv", limit = 2)
psms$retention_time_unit
#> [1] "minutes"
readers_retention_time_in_minutes(psms)
#> [1] 0.03233000 0.04907667

# TopFD's _ms1.feature gives no basis for a unit, so conversion refuses rather than guess.
features <- readers_read_features("Ms1Feature_TopFDv1.6.2_ms1.feature", limit = 5)
try(readers_retention_time_in_minutes(features))
#> Error : Cannot convert retention time for 'Ms1Feature': mzLib gives no basis to say what unit it is in. Inspect the values against scan numbers before comparing them.
```
