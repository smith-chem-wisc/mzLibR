# Every cell under a column, one list element per row

The accessor for a multi-cardinality column such as
`comment[modification parameters]`, which legitimately repeats - up to
eight times in one corpus file. Positions a row is too short to reach
are skipped, so each element holds only cells that exist.

## Usage

``` r
sdrf_all(doc, column)
```

## Arguments

- doc:

  An
  [`sdrf_read`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  or
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  result.

- column:

  A column name.

## Value

A list with one character vector per returned row.

## Examples

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
sdrf_all(doc, "comment[modification parameters]")[[1]]
#> [1] "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"           
#> [2] "NT=Oxidation;AC=UNIMOD:35;TA=M,W,H;MT=Variable"         
#> [3] "NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=Variable"    
#> [4] "NT=Acetyl;AC=UNIMOD:1;TA=K;MT=Variable"                 
#> [5] "NT=Ammonia-loss;AC=UNIMOD:385;TA=Q,C;MT=Variable"       
#> [6] "NT=Glu->pyro-Glu;AC=UNIMOD:27;TA=E;MT=Variable"         
#> [7] "NT=Deamidated;AC=UNIMOD:7;PP=Protein N-term;MT=Variable"
#> [8] "NT=Phospho;AC=UNIMOD:21;TA=S,T,Y,A;MT=Variable"         
```
