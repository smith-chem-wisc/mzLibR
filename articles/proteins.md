# Protein databases: organisms, genes and unique peptides

A search result names proteins by accession. Three questions it cannot
answer on its own are answered by mzLib from the protein database you
searched:

| question | function |
|----|----|
| what is this accession: organism, taxon, genes, GO terms? | [`proteins_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md) |
| which Ensembl gene is it, reproducibly? | [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md) |
| does this peptide identify one protein? | [`proteins_classify_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md) |

All three take one database or several - UniProt XML or FASTA,
optionally gzipped - and read several in **one** bridge call, `threads`
at a time. Mark contaminant databases with `contaminants`: it changes
answers.

## What is this accession?

``` r

db <- proteins_read(c("human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"),
                    contaminants = "contaminants.fasta",
                    tables = c("proteins", "go_terms", "ensembl_genes"))
db
#> <mzlibr_protein_database> 4 of 4 database(s) read, 8 proteins
#>   tables: proteins, go_terms, ensembl_genes
#>   go_terms: 47 rows
#>   ensembl_genes: 6 rows
#>   3 database(s) cannot say: go_terms, ensembl_genes, ensembl_gene_ids (see files$absent_fields)
db$proteins[, c("accession", "organism", "ncbi_taxonomy_id", "primary_gene_name",
                "length", "monoisotopic_mass", "is_contaminant")]
#>   accession     organism ncbi_taxonomy_id primary_gene_name length
#> 1    P04406 Homo sapiens             9606             GAPDH    335
#> 2    O43653 Homo sapiens             9606              PSCA    114
#> 3    P02768 Homo sapiens             9606               ALB    609
#> 4    Q13409 Homo sapiens             9606           DYNC1I2    638
#> 5  Q13409-2 Homo sapiens             9606           DYNC1I2    632
#> 6  Q13409-3 Homo sapiens             9606           DYNC1I2    612
#> 7    Q9Z0X1 Mus musculus            10090             Aifm1    612
#> 8    P02769   Bos taurus             9913               ALB    607
#>   monoisotopic_mass is_contaminant
#> 1          36030.40          FALSE
#> 2          11950.87          FALSE
#> 3          69321.50          FALSE
#> 4          71412.11          FALSE
#> 5          70600.75          FALSE
#> 6          68383.72          FALSE
#> 7          66723.84          FALSE
#> 8          69248.44           TRUE
```

`length` is in residues and `monoisotopic_mass` in daltons, of the
unmodified sequence as written plus one water - the precursor, not the
mature protein. `ncbi_taxonomy_id` is kept a string: it is an
identifier, not a quantity.

GO terms and Ensembl gene links are long tables, one row per protein and
term or transcript, returned only when asked for:

``` r

head(db$go_terms[, c("accession", "go_id", "aspect", "term_name")])
#>   accession      go_id            aspect
#> 1    P04406 GO:0005737 CellularComponent
#> 2    P04406 GO:0005829 CellularComponent
#> 3    P04406 GO:0070062 CellularComponent
#> 4    P04406 GO:0097452 CellularComponent
#> 5    P04406 GO:0043231 CellularComponent
#> 6    P04406 GO:0005811 CellularComponent
#>                                  term_name
#> 1                                cytoplasm
#> 2                                  cytosol
#> 3                    extracellular exosome
#> 4                             GAIT complex
#> 5 intracellular membrane-bounded organelle
#> 6                            lipid droplet
db$ensembl_genes[, c("accession", "gene_id", "versioned_gene_id", "gene_version")]
#>   accession         gene_id  versioned_gene_id gene_version
#> 1    P04406 ENSG00000111640 ENSG00000111640.15           15
#> 2    P04406 ENSG00000111640 ENSG00000111640.15           15
#> 3    P04406 ENSG00000111640 ENSG00000111640.15           15
#> 4    P04406 ENSG00000111640 ENSG00000111640.15           15
#> 5    P04406 ENSG00000111640 ENSG00000111640.15           15
#> 6    O43653 ENSG00000167653    ENSG00000167653           NA
```

### A FASTA cannot say

A FASTA header carries organism, taxon and gene name, and nothing else.
So a FASTA contributes no GO or Ensembl rows - and says so, rather than
leaving the tables silently short:

``` r

db$files[, c("file_type", "contaminant", "protein_count")]
#>    file_type contaminant protein_count
#> 1 UniProtXml       FALSE             2
#> 2      Fasta       FALSE             4
#> 3      Fasta       FALSE             1
#> 4      Fasta        TRUE             1
db$files$absent_fields
#> [[1]]
#> character(0)
#> 
#> [[2]]
#> [1] "go_terms"         "ensembl_genes"    "ensembl_gene_ids"
#> 
#> [[3]]
#> [1] "go_terms"         "ensembl_genes"    "ensembl_gene_ids"
#> 
#> [[4]]
#> [1] "go_terms"         "ensembl_genes"    "ensembl_gene_ids"
```

Empty there means “this format is silent”, never “this protein has
none”. Read the UniProt XML of the same proteome for those.

`source_index` in every table is 1-based, so it indexes `files`
directly:

``` r

table(basename(db$files$path[db$proteins$source_index]))
#> 
#> contaminants.fasta  human_extra.fasta   human_subset.xml  mouse_aifm1.fasta 
#>                  1                  4                  2                  1
```

## Which gene is it, reproducibly?

[`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
resolves each protein to stable Ensembl gene ids, **counted against a
gene set you supply**: an Ensembl GTF for one release. Nothing is
downloaded or defaulted, because the release is part of the answer.

``` r

genes <- proteins_resolve_genes(c("human_subset.xml", "human_extra.fasta"),
                                contaminants = "contaminants.fasta",
                                gtf = "Homo_sapiens.GRCh38.116.gtf",
                                xref = "Homo_sapiens.GRCh38.116.uniprot.tsv")
genes$gene_set[c("source_file_name", "release", "genome_build")]
#> $source_file_name
#> [1] "Homo_sapiens.GRCh38.116.gtf"
#> 
#> $release
#> [1] "116"
#> 
#> $genome_build
#> [1] "GRCh38.p14"
genes$outcome_counts
#>               resolved             multi_gene       off_primary_only 
#>                      1                      0                      1 
#>          not_in_source unrecognized_accession contaminant_not_mapped 
#>                      4                      0                      1
genes$resolutions[, c("accession", "outcome", "n_genes", "gene_id", "gene_symbol",
                      "off_primary_genes", "ensembl_xref_agrees")]
#>   accession                outcome n_genes         gene_id gene_symbol
#> 1    P04406               resolved       1 ENSG00000111640       GAPDH
#> 2    O43653       off_primary_only       0            <NA>        <NA>
#> 3    P02768          not_in_source       0            <NA>        <NA>
#> 4    Q13409          not_in_source       0            <NA>        <NA>
#> 5  Q13409-2          not_in_source       0            <NA>        <NA>
#> 6  Q13409-3          not_in_source       0            <NA>        <NA>
#> 7    P02769 contaminant_not_mapped       0            <NA>        <NA>
#>   off_primary_genes ensembl_xref_agrees
#> 1                 0                TRUE
#> 2                 1                  NA
#> 3                 0                  NA
#> 4                 0                  NA
#> 5                 0                  NA
#> 6                 0                  NA
#> 7                 0                  NA
```

Every protein gets exactly one outcome. A `multi_gene` protein has one
row per gene, never a pick. `off_primary_only` means its genes exist
only on ALT haplotypes or patches outside the gene set - use Ensembl’s
primary-assembly GTF. `ensembl_xref_agrees` is `NA` when there is
nothing to compare, which is unknown, not false. Keep
`genes$gene_set$sha256` with your results: it pins the exact gene set.

## Does this peptide identify one protein?

``` r

calls <- proteins_classify_peptides(
  c("VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"),
  c("human_subset.xml", "human_extra.fasta"),
  contaminants = "contaminants.fasta")
calls$peptides[, c("peptide", "sharing", "accession_count")]
#>           peptide           sharing accession_count
#> 1        VGVNGFGR            Unique               1
#> 2  LVLNGNPLTLFQER            Unique               1
#> 3 ALSEQINIFFDYSGR  SharedWithinGene               3
#> 4         YLYEIAR SharedAcrossGenes               2
#> 5        AEFVEVTK            Unique               1
#> 6        PEPTIDEK     NotInDatabase               0
calls$sharing_counts
#>     NotInDatabase            Unique  SharedWithinGene SharedAcrossGenes 
#>                 1                 3                 1                 1
```

- `LVLNGNPLTLFQER` is GAPDH’s `LVINGNPITIFQER`: **I and L are one
  residue**, because a mass spectrometer cannot tell them apart.
- `ALSEQINIFFDYSGR` is in three isoforms of one gene: `SharedWithinGene`
  supports the gene, not an isoform.
- `YLYEIAR` is in human and bovine albumin, which share no gene. Leave
  the contaminant database out and it would read as `Unique`.

``` r

calls$peptides$accessions[[4]]
#> [1] "P02768" "P02769"
calls$peptides$shared_gene_keys[[3]]
#> [1] "entry:Q13409"              "gene:Homo sapiens:DYNC1I2"
```

Containment, not digestion: a protein contains a peptide if its sequence
contains it anywhere, whatever the protease. That is deliberately
conservative, so a peptide called `Unique` cannot be explained by
another entry at a site the search’s cleavage rules skipped.
