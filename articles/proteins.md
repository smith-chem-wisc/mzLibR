# Protein databases: organisms, genes, unique peptides and GO

A search result names proteins by accession. The questions it cannot
answer on its own are answered by mzLib from the protein database you
searched:

| question | function | mzLib |
|----|----|----|
| what is this accession: organism, taxon, genes, GO terms? | [`proteins_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md) | `ProteinDbLoader` |
| which Ensembl gene is it, reproducibly? | [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md) | `EnsemblGeneResolver` |
| does this peptide identify one protein? | [`proteins_classify_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md) | `PeptideUniquenessClassifier` |
| which GO terms does each protein group carry, every member kept? | [`proteins_annotate_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md) | `GoGroupAnnotator` (#1353), `ToGoAnnotationGroups` (#1366) |
| where do I get a go.obo, and keep it pinned? | [`proteins_update_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_update_go.md) | `Loaders.UpdateGeneOntology` (#1353) |

The first three take one database or several - UniProt XML or FASTA,
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

To look up only some accessions, pass them as `accessions`. The match is
**exact**, so an isoform suffix that the database does not use is
reported, not guessed at:

``` r

found <- proteins_read(c("human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"),
                       contaminants = "contaminants.fasta",
                       accessions = c("P04406", "Q9Z0X1", "P02769", "P04406-1"))
setNames(found$proteins$ncbi_taxonomy_id, found$proteins$accession)
#>  P04406  Q9Z0X1  P02769 
#>  "9606" "10090"  "9913"
found$accessions_not_found
#> [1] "P04406-1"
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

## Which GO terms does each protein group carry?

A search reports *protein groups*, and a group can have several members:
two histone H2A.Z variants no peptide could tell apart are one group,
`P0C0S5|Q71UI9`. Taking the first accession throws away half the
evidence, and it is not even a choice MetaMorpheus made: it sorts
members and never picks a leading protein. And UniProt records the most
specific terms, so a protein annotated to *mitochondrial inner membrane*
is not, as written, “in the mitochondrion” at all.

[`proteins_annotate_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)
answers both, with mzLib’s `GoGroupAnnotator` reading the protein-group
table MetaMorpheus wrote: **one row per (group, term)** that *any*
member holds, directly or through an ancestor (`is_a` and `part_of`),
and every row says which members carry it. The data here is real:
mzLib’s copy of a MetaMorpheus 1.1.11 protein-group table from
PXD036557, the UniProt entries it names, and GO release 2026-07-26,
trimmed to the terms those proteins reach.

### First, get a go.obo - once, and keep it

Terms are added, renamed and moved between GO releases, so a GO result
means something only relative to one release.
[`proteins_annotate_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)
never downloads anything; fetching is its own function, called on
purpose:

``` r

update <- proteins_update_go("go.obo")
update$go[c("release", "term_count")]
#> $release
#> [1] "releases/2026-07-26"
#> 
#> $term_count
#> [1] 48340
```

Call it again later and a newer release replaces the file, while the old
one is kept beside it as `go.obo.<timestamp>`, so a result you already
published can still be reproduced.

### Annotate the groups

``` r

go <- proteins_annotate_go("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
                           go_obo = "go-pxd036557.obo", category_map = "organelle_map.tsv")
go
#> <mzlibr_go_annotations> 5 groups, 563 (group, term) rows
#>   go.obo: go-pxd036557.obo, releases/2026-07-26, sha256 c1cdfd098c5c
#>   groups at q <= 0.01: annotated 4, no_go_terms 0, no_entry 0, contaminant 1
#>   categories: organelle v1, 30 rows
c(table_rows = go$table_row_count, decoys = go$decoy_group_count, groups = go$group_count,
  rows = go$row_count)
#> table_rows     decoys     groups       rows 
#>          6          1          5        563
```

The decoy group is skipped, because decoys carry no GO. Albumin, a
contaminant, gets one row with no term, and its `annotation_status` says
why:

``` r

go$annotations[go$annotations$protein_group == "P02769",
               c("protein_group", "go_id", "annotation_status", "n_with", "n_members")]
#>     protein_group go_id annotation_status n_with n_members
#> 358        P02769  <NA>       contaminant      0         1
```

Both histones are annotated to the nucleosome, each on its own evidence:

``` r

histones <- go$annotations[go$annotations$protein_group == "P0C0S5|Q71UI9", ]
nucleosome <- histones[histones$go_id %in% "GO:0000786", ]
nucleosome[, c("go_name", "aspect", "n_with", "n_members", "propagated", "inherited")]
#>        go_name             aspect n_with n_members propagated inherited
#> 253 nucleosome cellular_component      2         2      FALSE     FALSE
nucleosome$accession_used[[1]]
#> [1] "P0C0S5" "Q71UI9"
nucleosome$evidence_by_member[[1]]
#> $P0C0S5
#> [1] "ECO:0000353"
#> 
#> $Q71UI9
#> [1] "ECO:0000353"
```

### Union, consensus, direct only: filters, not modes

The table is the **union** over a group’s members. Everything narrower
is one comparison on the rows:

``` r

c(union = nrow(histones), consensus = sum(histones$n_with == histones$n_members),
  direct = sum(!histones$propagated))
#>     union consensus    direct 
#>       106        57        18
histones$accession_used[histones$go_name %in% "euchromatin"]
#> [[1]]
#> [1] "P0C0S5"
```

*Euchromatin* is one histone’s term, not the other’s. For an enrichment
tool that propagates on its own, keep the rows with `propagated` false.

### In your own vocabulary

`category_map` maps terms to the categories you care about, through
anchors: a term belongs to a category when the anchor is the term or one
of its ancestors. mzLib ships no vocabulary.

``` r

head(go$categories$records)
#>        go_id      category       subcategory
#> 1 GO:0000785       nucleus nucleus:chromatin
#> 2 GO:0000786       nucleus nucleus:chromatin
#> 3 GO:0000791       nucleus nucleus:chromatin
#> 4 GO:0000792       nucleus nucleus:chromatin
#> 5 GO:0005576 extracellular              <NA>
#> 6 GO:0005634       nucleus              <NA>
```

### Pin what you used

Keep `go$go$sha256` and `go$annotation_database$sha256` with your
results. The header’s counters count **groups at q \<= 0.01**, not rows;
the rows themselves are never filtered on q.

``` r

go$go[c("release", "sha256")]
#> $release
#> [1] "releases/2026-07-26"
#> 
#> $sha256
#> [1] "c1cdfd098c5cf395ce4bbefdea017fd698fdf22d135b28623e74f6f714e1c5fa"
go$header[c("counter_q_value_max", "status_annotated", "status_contaminant")]
#> counter_q_value_max    status_annotated  status_contaminant 
#>              "0.01"                 "4"                 "1"
```

When the UniProt XML is newer than the go.obo it cites terms the release
lacks, and the call fails naming every one. Fetch a newer go.obo, or
pass `skip_unknown_go_ids = TRUE` to drop those ids and have them listed
in `unresolved_go_ids`.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl releases, GTFs and cross-references (Ensembl 2026, NAR).
- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, whose XML dbReferences carry the gene links resolved (UniProt
  Consortium 2025, NAR).
- [doi:10.1093/genetics/iyad031](https://doi.org/10.1093/genetics/iyad031):
  The Gene Ontology knowledgebase in 2023 (GO Consortium, Genetics 224,
  iyad031), the ontology and annotations used.
- [doi:10.1038/75556](https://doi.org/10.1038/75556): Gene Ontology:
  tool for the unification of biology (Ashburner et al., Nat Genet 25,
  25-29, 2000), the is_a / part_of graph propagated over.
