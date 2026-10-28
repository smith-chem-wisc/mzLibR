# One species pair's counts, in either order

The row of a snapshot's `pairs` for species `a` and `b`, which may be
given either way round: side A is always the species that sorts first.

## Usage

``` r
proteins_orthology_pair(snapshot, a, b)
```

## Arguments

- snapshot:

  A
  [`proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)
  result.

- a,b:

  Two species of the snapshot, Ensembl production names; the same one
  twice for its paralogs.

## Value

A one-row data.frame: `species_a`, `species_b`, and `row_count`,
`ortholog_count`, `paralog_count` and `homoeolog_count`, in rows.

## See also

[`proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)

## Examples

``` r

gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
          mus_musculus = "Mus_musculus.GRCm39.116.gtf",
          rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
snap <- proteins_build_orthology("116", gtfs = gtfs,
  dumps = paste0("homologies/", names(gtfs), "/Compara.116.protein_default.homologies.tsv.gz"),
  gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
  pair_status = c("homo_sapiens", "mus_musculus"))
proteins_orthology_pair(snap, "mus_musculus", "mus_musculus")$paralog_count
#> [1] 1
```
