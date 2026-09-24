# The published locations of each file, as controlled-vocabulary terms

One long data.frame rather than a nested list, so it pipes into
anything.

## Usage

``` r
pride_locations(files)
```

## Arguments

- files:

  A
  [`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  data.frame.

## Details

You rarely need this: `https_url` already carries the fetchable URL,
resolved by mzLib's `TryGetHttpsDownloadUrl`, which **searches** the
locations rather than taking the first. That distinction matters - **the
order is not stable.** In PXD000001 the mztab lists FTP first while the
MGF lists **Aspera** first, so code that took `locations[[1]]` would get
an unfetchable `prd_ascp@fasp.ebi.ac.uk:...` address for some files and
a working one for others, in the same project. Do not re-implement the
search.

## Value

A data.frame with `file_name`, `accession`, `name` and `value`, one row
per location.

## Examples

``` r

files <- pride_list_files("PXD000001")
head(pride_locations(files))
#>                                                     file_name     accession
#> 1                  PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz PRIDE:0000469
#> 2                  PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz PRIDE:0000468
#> 3                    PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz PRIDE:0000468
#> 4                    PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz PRIDE:0000469
#> 5 TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML PRIDE:0000469
#> 6 TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML PRIDE:0000468
#>              name
#> 1    FTP Protocol
#> 2 Aspera Protocol
#> 3 Aspera Protocol
#> 4    FTP Protocol
#> 5    FTP Protocol
#> 6 Aspera Protocol
#>                                                                                                                        value
#> 1        ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz
#> 2          prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2012/03/PXD000001/generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz
#> 3            prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2012/03/PXD000001/generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz
#> 4          ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz
#> 5 ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
#> 6   prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2012/03/PXD000001/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
```
