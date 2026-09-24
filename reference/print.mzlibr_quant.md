# Print a quantification result

A compact summary that names the two things most likely to be misread:
how many match-between-runs transfers the peptide roll-up is hiding, and
how many protein intensities are `NA` rather than `0`.

## Usage

``` r
# S3 method for class 'mzlibr_quant'
print(x, ...)
```

## Arguments

- x:

  A
  [`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  result.

- ...:

  Ignored.

## Value

`x`, invisibly.
