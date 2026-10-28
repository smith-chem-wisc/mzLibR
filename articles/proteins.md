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
| where is residue 120 of this protein in its isoform, or in the bovine protein? | [`proteins_align()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md), [`proteins_map_source_to_target()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md) | `PairwiseAligner` (#1428) |
| which entry, isoform and variants does `P04637-2_S70N` or `NP_000537.3_R72P` name? | [`proteins_parse_accessions()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_parse_accessions.md) | `VariantApplication.ParseAccession` (#1382) |
| which mouse genes are orthologous to this human gene, and if none, why not? | [`proteins_build_orthology()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md) | `OrthologySnapshot`, `OrthologySnapshotWriter` (#1381) |

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

## Map a residue from one sequence onto another

A position means something only in the sequence it was counted in. A
phosphosite reported on the canonical isoform, a variant called against
last year’s UniProt release, a cysteine you know in the human protein:
before you look for it in an isoform, an ortholog or a newer entry, you
need the residue that *corresponds* to it there. Counting the same
number along the other sequence is wrong as soon as there is one
insertion or deletion before it.

[`proteins_align()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
aligns pairs of sequences with mzLib’s `PairwiseAligner` (#1428) and
maps every residue to its counterpart in the other sequence, or to `NA`
when it has none. Here the canonical DYNC1I2 sequence (Q13409) is
aligned with its isoforms 2B and 2C, read from the database with
`sequences = TRUE`:

``` r

seqs <- proteins_read("human_extra.fasta", contaminants = "contaminants.fasta",
                      accessions = c("Q13409", "Q13409-2", "Q13409-3", "P02768", "P02769"),
                      sequences = TRUE)
seq <- setNames(seqs$proteins$sequence, seqs$proteins$accession)
isoforms <- proteins_align(list(c(seq[["Q13409"]], seq[["Q13409-2"]]),
                                c(seq[["Q13409"]], seq[["Q13409-3"]])))
two <- proteins_alignments(isoforms)
p <- c(75, 76, 82, 100, 120, 133)
data.frame(canonical = p,
           isoform_2B = proteins_map_source_to_target(two[[1]], p),
           isoform_2C = proteins_map_source_to_target(two[[2]], p))
#>   canonical isoform_2B isoform_2C
#> 1        75         75         75
#> 2        76         NA         NA
#> 3        82         76         76
#> 4       100         94         94
#> 5       120        114         NA
#> 6       133        127        107
```

Reading it row by row:

- **Up to residue 75 the numbers agree.** Both isoforms begin with the
  same 75 residues.
- **Residues 76 to 81 have no counterpart in either isoform**: `NA`,
  never the nearest residue.
- **After the gap the numbers are shifted**: canonical residue 100 is
  residue 94 of each isoform.
- **Isoform 2C lacks 20 more** (canonical 113 to 132), so canonical 120
  has no counterpart in 2C, and canonical 133 is 2C’s residue 107.

[`proteins_map_target_to_source()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_target_to_source.md)
goes the other way. The aligned rows hold `-` wherever a residue faces a
gap, so printing them one above the other shows where the gaps fall:

``` r

writeLines(substr(c(isoforms$records$aligned_source[[2]], isoforms$records$aligned_target[[2]]), 61, 140))
#> EALLQSMGLTPESPIVFSEYWVPPPMSPSSKSVSTPSEAGSQDSGDGAVGSRTLHWDTDPSVLQLHSDSDLGRGPIKLGM
#> EALLQSMGLTPESPI------VPPPMSPSSKSVSTPSEAGSQDSGDGAVGSR--------------------RGPIKLGM
isoforms$records[, c("pair_index", "identical", "substituted", "source_residues_facing_gap")]
#>   pair_index identical substituted source_residues_facing_gap
#> 1          1       632           0                          6
#> 2          2       612           0                         26
isoforms$aligner_id
#> [1] "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"
```

`identical` and `substituted` count the paired residues that are the
same letter and that differ; `source_residues_facing_gap` and
`target_residues_facing_gap` count the residues with no counterpart.
`score` is the alignment’s BLOSUM62 score in half-bits (see the traps).
`aligner_id` names the method and every parameter, so store it beside
any position you map.

### Orthologs: the counterpart can be a different amino acid

Across species most residues still pair, but many pair with a
*different* residue. Human serum albumin’s free thiol, Cys34 of the
mature protein, is residue 58 of the UniProt sequence, which counts the
24-residue signal peptide and propeptide. In bovine albumin it is a
cysteine at 58 too:

``` r

albumin <- proteins_align(list(c(seq[["P02768"]], seq[["P02769"]])))
albumin$records[, c("identical", "substituted", "source_residues_facing_gap")]
#>   identical substituted source_residues_facing_gap
#> 1       465         142                          2
map <- albumin$residue_map
map[map$sequence == "source" & map$position %in% c(26, 58),
    c("position", "residue", "mapped_position", "mapped_residue")]
#>    position residue mapped_position mapped_residue
#> 26       26       A              26              T
#> 58       58       C              58              C
```

Residue 26 pairs, but it is A in human and T in bovine albumin. **Always
check the residue at the mapped position** before you carry a
modification across: a phosphosite on S is not a site on A. The residue
map puts both letters side by side.

### The options, and whose defaults they are

Every option you leave out is **mzLib’s to set**. mzLibR sends nothing
for it, and the result says what the aligner used. The isoform alignment
above passed no options, so these are mzLib’s own defaults in this
build:

``` r

str(unclass(isoforms)[c("matrix", "gap_open", "gap_extend", "end_gaps", "max_cells")])
#> List of 5
#>  $ matrix    : chr "BLOSUM62"
#>  $ gap_open  : num 11
#>  $ gap_extend: num 1
#>  $ end_gaps  : chr "free"
#>  $ max_cells : num 2e+08
```

| argument | meaning |
|----|----|
| `gap_open` | cost of opening a gap, in BLOSUM62’s half-bits; a gap of *k* residues costs `gap_open + k * gap_extend` (BLAST’s convention; mzLib’s defaults are BLASTP’s) |
| `gap_extend` | cost of each residue in a gap, the first included |
| `end_gaps` | `"free"`: an overhang at either end costs nothing (EMBOSS needle’s default), so a signal peptide or a longer tail is not penalised. `"charged"`: every gap costs |
| `max_cells` | the largest table mzLib will build, `(nchar(source) + 1) * (nchar(target) + 1)` cells of one byte; a larger pair is refused rather than allowed to exhaust memory |
| `residue_map` | `FALSE` leaves out the residue map, for a large batch where you want only the scores and counts |

``` r

charged <- proteins_align(list(c("MLLLPEPTIDE", "PEPTIDE")), end_gaps = "charged")
charged$records[, c("aligned_target", "score")]
#>   aligned_target score
#> 1    ----PEPTIDE    24
charged$aligner_id
#> [1] "global-affine;BLOSUM62;open=11;extend=1;end-gaps=charged"
```

With free end gaps the same pair scores 39: the four-residue overhang
costs `11 + 4 * 1` only when it is charged.

### The traps

- **`score` is not an identity, a bit score or an E-value.** It is a raw
  sum of BLOSUM62 scores in half-bits, less the gaps, and it grows with
  length. Compare scores only between alignments with the same
  `aligner_id`. mzLibR computes no percent identity, because there are
  several definitions: divide `identical` by the length you mean.
- **It is a global alignment.** Every residue of both sequences is
  placed, paired or facing a gap. With free end gaps a short sequence is
  placed on the best-scoring stretch of a long one, even when it has no
  real match there.
- **One optimal alignment, chosen the same way every time.** Isoform 2B
  above is a real case: canonical V76 and V82 are both V, either could
  pair with the isoform’s V76, and mzLib pairs V82. Another aligner may
  give the equally good other answer.
- **Pass bare one-letter sequences.** mzLib refuses only `-` and `.`, as
  an `mzlib_usage_error`. Any other character BLOSUM62 does not name -
  selenocysteine `U`, pyrrolysine `O`, a digit, a space, a modification
  bracket - is scored as an unknown residue `X`, silently.
- **Long pairs need a larger `max_cells`.** Two titin isoforms need
  about 1.2 billion cells, past the default; mzLib accepts up to
  2,147,483,591.
- **Many pairs are one call.**
  [`proteins_align()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
  takes the whole list at once, in one bridge process, so a few thousand
  pairs cost one start-up; never loop it.

## Read an accession into entry, isoform and variants

The accessions a search reports are not all plain entries. UniProt
isoforms carry `-2`; RefSeq accessions carry a `.3` version; and mzLib
names a protein with applied sequence variants after its entry,
`P04637_S70N_A80T` or `NP_000537.3_R72P`. To join results to anything
you keep per entry - genes, orthologs, GO terms - you need the entry
back, without losing the isoform or the variants. Cutting at the first
`_` turns `NP_000537` into `NP`.

[`proteins_parse_accessions()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_parse_accessions.md)
reads each accession with mzLib’s own grammars
(`VariantApplication.ParseAccession`, \#1382):

``` r

parsed <- proteins_parse_accessions(c("P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P",
  "P04637_S70N_A80T", "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"))
parsed$records[, c("input", "namespace", "entry_accession", "isoform", "version", "applied_variants")]
#>                 input    namespace     entry_accession isoform version
#> 1              P04637      uniprot              P04637      NA      NA
#> 2            P04637-2      uniprot              P04637       2      NA
#> 3            Q13409-2      uniprot              Q13409       2      NA
#> 4    NP_000537.3_R72P       refseq           NP_000537      NA       3
#> 5    P04637_S70N_A80T      uniprot              P04637      NA      NA
#> 6            P04637_2 unrecognized            P04637_2      NA      NA
#> 7        DECOY_P02768 unrecognized        DECOY_P02768      NA      NA
#> 8 sp|P04637|P53_HUMAN unrecognized sp|P04637|P53_HUMAN      NA      NA
#> 9       P04637|P02768 unrecognized       P04637|P02768      NA      NA
#>   applied_variants
#> 1             <NA>
#> 2             <NA>
#> 3             <NA>
#> 4             R72P
#> 5        S70N_A80T
#> 6             <NA>
#> 7             <NA>
#> 8             <NA>
#> 9             <NA>
parsed$namespace_counts
#> unrecognized      uniprot       refseq 
#>            4            4            1
```

- **`entry_accession`** is the UniProt entry or the unversioned RefSeq
  accession. It answers “same entry?”, **not** “same sequence?”:
  `P04637` and `P04637-2` share it, but they are different sequences.
  Join per-entry data on it; join sequence-level data on `input`.
- **`isoform`** and **`version`** are the `-n` and `.n` numbers, `NA`
  when there is none.
- **`parent_accession`** is the entry the variants were applied to,
  isoform or version kept: `NP_000537.3` for `NP_000537.3_R72P`. Look
  the parent sequence up under it.
- **`applied_variants`** is the variant suffix exactly as written, `NA`
  when there is none. It is not checked against any sequence.
- **`namespace`** is `uniprot`, `refseq` or `unrecognized`, spelled as
  [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
  spells it, so the two tables join.

**It parses and never repairs.** Everything outside the two grammars
comes back `unrecognized`, with `entry_accession` equal to the input, so
a data problem stays visible instead of being fixed silently.
**`P04637_2` is not a variant of P04637**: it is mzLib’s loader counter
for a *different* entry whose accession collided with P04637’s. Decoy,
contaminant and entrapment prefixes and FASTA headers are left as they
are. A protein group (`P04637|P02768`) is not split: split groups on `|`
first, if that is what you mean. Lower case and spaces count, because
mzLibR sends every accession exactly as you give it. Only RefSeq
*protein* accessions (`NP_`, `XP_`, `YP_`, `WP_`, `AP_`) are recognized,
so a transcript (`NM_`) or an Ensembl id is unrecognized too.

## Orthologs across species: an Ensembl Compara snapshot

You have a mouse result and a human question: which human gene is this
mouse gene’s ortholog? Ensembl Compara answers that, per species pair,
in bulk files with three traps:

1.  **The dumps partition.** Each species has its own homology dump, and
    Compara puts each of a pair’s relationships in *one* of the pair’s
    two dumps, arbitrarily. In release 116 every human-mouse orthology
    is in the **mouse** dump and none is in the human one. Reading the
    dump named after your species gives a complete-looking, badly
    incomplete answer.
2.  **Orthology is not transitive.** Human A ~ mouse B and mouse B ~ rat
    C does not make A ~ C. Chaining pairs into groups invents tuples
    Compara never asserted.
3.  **“No ortholog” means three things.** The gene’s tree may hold a
    gene of the other species with no ortholog called, hold none of that
    species, or the gene may be in no tree at all.

[`proteins_build_orthology()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)
builds mzLib’s orthology snapshot (`OrthologySnapshot`, \#1381), which
makes the first trap structural - **mzLib refuses a call without every
species’ dump** - and keeps the other two visible. This example is
mzLib’s own test world: four human genes, three mouse and two rat, with
every human-mouse row in the mouse dump. The species are named by
Ensembl production name, the names of `gtfs`:

``` r

gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
          mus_musculus = "Mus_musculus.GRCm39.116.gtf",
          rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
dumps <- paste0("homologies/", names(gtfs), "/Compara.116.protein_default.homologies.tsv.gz")
snap <- proteins_build_orthology("116", gtfs = gtfs, dumps = dumps,
  gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
  pair_status = c("homo_sapiens", "mus_musculus"))   # every human gene, in mouse
snap$pair_status$records[, c("gene_id", "gene_biotype", "status")]
#>           gene_id   gene_biotype                    status
#> 1 ENSG00000000001 protein_coding              has_ortholog
#> 2 ENSG00000000002 protein_coding              has_ortholog
#> 3 ENSG00000000003 protein_coding tree_lacks_target_species
#> 4 ENSG00000000004 protein_coding           not_in_any_tree
substr(snap$snapshot_id, 1, 12)
#> [1] "5487ebcabe48"
```

Record `snap$snapshot_id` with any result you draw from the snapshot. It
is mzLib’s sha256 over the format, the release, the species and every
input’s sha256, so the same inputs always give the same id, and a
different id means different inputs.

### Every species’ dump, always

The result shows where the rows came from. The human dump holds one row
(human-rat); all three human-mouse rows came from the mouse dump:

``` r

snap$dumps[, c("species", "rows_read", "rows_kept")]
#>             species rows_read rows_kept
#> 1      homo_sapiens         1         1
#> 2      mus_musculus         5         5
#> 3 rattus_norvegicus         1         1
snap$pairs
#>           species_a         species_b row_count ortholog_count paralog_count
#> 1      homo_sapiens      homo_sapiens         0              0             0
#> 2      homo_sapiens      mus_musculus         3              3             0
#> 3      homo_sapiens rattus_norvegicus         1              1             0
#> 4      mus_musculus      mus_musculus         1              0             1
#> 5      mus_musculus rattus_norvegicus         2              2             0
#> 6 rattus_norvegicus rattus_norvegicus         0              0             0
#>   homoeolog_count
#> 1               0
#> 2               0
#> 3               0
#> 4               0
#> 5               0
#> 6               0
proteins_orthology_pair(snap, "mus_musculus", "homo_sapiens")$row_count   # either order
#> [1] 3
```

Leave the mouse dump out and mzLib refuses: *“No homology dump for
mus_musculus. Each pair’s relationships are split between its two
species’ dumps, so every species’ dump is required.”* That is an
`mzlib_usage_error`. A pair with nothing in it, such as the human
paralogs here, still has its row.

### What to download, once per release

Nothing is downloaded for you. For release 116, from
`https://ftp.ensembl.org/pub/release-116/`:

| input | path under the release | size (human, 116) | argument |
|----|----|----|----|
| each species’ primary-assembly GTF | `gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` | 141 MB | `gtfs = c(species = path)` |
| each species’ homology dump | `tsv/ensembl-compara/homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz` | 109 MB | `dumps` |
| the gene trees, one file for all species | `compara/vertebrates.GeneTree_content.default.e116.txt.gz` | 79 MB | `gene_trees` |

**Take the dump in each species’ own folder**: every species’ dump has
the same file name, and one folder up, `homologies/` also holds a 10 GB
file of that name, which is not one species’ dump. **Keep Ensembl’s file
names**: mzLib holds the gene trees and every dump to `release` only
when their names carry it, and a caveat says which were not checked. A
compact gene table made with mzLib’s `EnsemblGeneSetWriter` for
[`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
works too, in `gene_sets` instead of `gtfs`; it carries the GTF’s
provenance, so the snapshot and its id are the same.

### What each status means

`pair_status = c(from, to)` gives every gene of `from`’s gene set
exactly one status in `to`. It says how far Compara’s inference got, not
what biology holds:

| status | what it means | read it as |
|----|----|----|
| `has_ortholog` | at least one ortholog in `to` | an answer; the Parquet has the rows |
| `no_edge_in_shared_tree` | the gene’s tree holds a `to` gene, and Compara called no ortholog to it | the nearest the data comes to “no ortholog exists” |
| `tree_lacks_target_species` | the gene’s tree holds no `to` gene | nothing was compared |
| `not_in_any_tree` | the gene is in the gene set but in no gene tree | no call was attempted |

`snap$pair_status$status_counts` counts genes per status, every status
present. **A status is not symmetric**: ask for
`c("mus_musculus", "homo_sapiens")` separately. **Every biotype is a
row**, so choose your denominator by filtering on `gene_biotype`.

### Write the snapshot, then query it

The relationships themselves go to disk. Give `out` a directory that
does not exist yet:

``` r

written <- proteins_build_orthology("116", gtfs = gtfs, dumps = dumps,
  gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
  restrict_gene_trees = TRUE, out = "compara-116")
identical(written$snapshot_id, snap$snapshot_id)
#> [1] TRUE
written$written$files[, c("path", "rows", "bytes")]
#>                                                  path rows bytes
#> 1                          genes/homo_sapiens.parquet    4  1444
#> 2                          genes/mus_musculus.parquet    3  1459
#> 3                     genes/rattus_norvegicus.parquet    2  1499
#> 4                        members/homo_sapiens.parquet    3  1554
#> 5                        members/mus_musculus.parquet    3  1609
#> 6                   members/rattus_norvegicus.parquet    2  1639
#> 7            pairs/homo_sapiens__homo_sapiens.parquet    0   653
#> 8            pairs/homo_sapiens__mus_musculus.parquet    3  4221
#> 9       pairs/homo_sapiens__rattus_norvegicus.parquet    1  4141
#> 10           pairs/mus_musculus__mus_musculus.parquet    1  4195
#> 11      pairs/mus_musculus__rattus_norvegicus.parquet    2  4272
#> 12 pairs/rattus_norvegicus__rattus_norvegicus.parquet    0   653
#> 13                                          views.sql   NA  3297
```

The directory holds `genes/SPECIES.parquet` (the gene set),
`members/SPECIES.parquet` (gene to Compara tree), `pairs/A__B.parquet`
(every relationship between A and B, A sorting first; A == B holds the
paralogs), `views.sql` and `manifest.json` (every input and file with
its sha256). **A directory that is not empty is refused**: a snapshot is
written once, and next release’s goes in a new directory to compare with
this one. `restrict_gene_trees = TRUE` keeps only your species’ genes
from the gene-tree file; without it, mzLib’s default keeps every gene in
it - 4.2 million for the vertebrates of release 116 - and a caveat says
how many. The snapshot is the same either way, as the equal ids show.

**mzLib writes Parquet and ships no reader, and mzLibR reads none** (it
has no dependencies). The `arrow` package reads the tables, and `duckdb`
runs the macros mzLib writes in `views.sql` (`orthologs`, `pair_status`,
`species_set3`):

``` r

# Not run: needs the arrow package, which mzLibR does not depend on.
pairs <- arrow::read_parquet("compara-116/pairs/homo_sapiens__mus_musculus.parquet")
pairs[pairs$relationship_class == "ortholog",
      c("gene_a", "gene_b", "relationship_type", "identity_a", "identity_b")]
```

`pairs/A__B.parquet` has one row per relationship, ordered by Compara’s
`homology_id`, with Compara’s values verbatim: `relationship_type`
(`ortholog_one2one`, `ortholog_one2many`, …), `relationship_class`, the
genes and proteins of each side, `identity_a` and `identity_b` (percent,
each over its **own** protein, so they differ), `goc_score`,
`wga_coverage`, `dn`, `ds`, `is_high_confidence` and `source_dump`. **A
NULL stays null, never 0 or false**, and arrow reads it as `NA`.
[`?proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)
lists every column of the three tables.

### The traps

- **Never chain pairs.** Use `views.sql`’s `species_set3` (or mzLib’s
  `OrthologySnapshot.SpeciesSet` for more species), which requires every
  two members of a tuple to be orthologs.
- **One-to-many is several rows.** Nothing is picked: a human gene with
  two mouse orthologs has two rows, both `ortholog_one2many`. Decide per
  question.
- **The gene set and the dumps must be one release.** A homology naming
  a gene its gene set lacks is refused as an `mzlib_bridge_error`
  (`InvalidDataException`): that is what a GTF and dumps from different
  releases look like.
- **`source_dump` cannot tell species apart.** Ensembl gives every
  species’ dump the same file name, so the column holds that one name
  everywhere, and a caveat says so. `manifest.json`’s `inputs` list each
  dump with its species and sha256.
- **Join with
  [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
  on the gene set.** Its `gene_set_sha256` column and
  `snap$gene_sets$sha256` are the same value when both were counted
  against the same GTF. Check that before joining a resolution to a
  snapshot.
- **It is slow and large.** Every input is read whole. For human, mouse
  and rat at release 116, mzLib reports about 20 s to build and 40 s
  more to write 13 files (9.0 MB).

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1101/gr.073585.107](https://doi.org/10.1101/gr.073585.107):
  Ensembl Compara gene trees, the duplication-aware phylogenies its
  orthology and paralogy calls are read from (Vilella et al. 2009,
  Genome Res 19:327).
- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl releases, their GTFs and the Compara homology dumps (Ensembl
  2026, NAR).
- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, whose XML dbReferences carry the gene links resolved (UniProt
  Consortium 2025, NAR).
- [doi:10.1016/0022-2836(70)90057-4](https://doi.org/10.1016/0022-2836%2870%2990057-4):
  Global alignment by dynamic programming (Needleman and Wunsch 1970, J
  Mol Biol 48:443).
- [doi:10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836%2882%2990398-9):
  Affine gap costs in the same time and memory order, the form mzLib
  implements (Gotoh 1982, J Mol Biol 162:705).
- [doi:10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915):
  BLOSUM62, the substitution matrix and its half-bit units (Henikoff and
  Henikoff 1992, PNAS 89:10915).
- [doi:10.1093/genetics/iyad031](https://doi.org/10.1093/genetics/iyad031):
  The Gene Ontology knowledgebase in 2023 (GO Consortium, Genetics 224,
  iyad031), the ontology and annotations used.
- [doi:10.1038/75556](https://doi.org/10.1038/75556): Gene Ontology:
  tool for the unification of biology (Ashburner et al., Nat Genet 25,
  25-29, 2000), the is_a / part_of graph propagated over.
- [doi:10.1093/nar/gkv1189](https://doi.org/10.1093/nar/gkv1189):
  RefSeq, whose protein accession prefixes and .n versions mzLib matches
  (O’Leary et al. 2016, NAR).
