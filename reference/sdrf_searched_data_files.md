# The searched file each row of a restricted SDRF was joined to

The `comment[searched data file]` cell of each row: the searched file's
name, with no folder.

## Usage

``` r
sdrf_searched_data_files(doc)
```

## Arguments

- doc:

  An
  [`sdrf_restrict`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)
  result, or an SDRF read back from one it wrote.

## Value

A character vector, one element per row; `NA` for a row the document
does not reach, or every row when the document has no such column.

## See also

[`sdrf_restrict`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)

## Examples

``` r

pool1 <- "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
pool5 <- "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
searched <- c(paste0("D:\\search\\", pool1, "1-calib.mzML"), paste0(pool1, "2-calib.mzML"),
  paste0(pool1, "2.raw"), "sample_A.mzML", "not_in_sdrf.mzML")
names(searched) <- c("", "", "", paste0(pool5, "1.raw"), "")
sdrf_searched_data_files(sdrf_restrict("PXD008841.subset.sdrf.tsv", searched, limit = 2))
#> [1] "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
#> [2] "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
```
