# Resolve every protein to stable Ensembl gene ids, against a gene set you pin

Resolve every protein in the given databases to stable Ensembl gene ids,
counted against a caller-supplied, release-pinned gene set, with one
outcome per protein and the hash of every input.

Wraps mzLib's `EnsemblGeneResolver`. The gene links come from the
database itself (UniProt's Ensembl `dbReference`s), and each is
**counted against the gene set you pass**: a gene in the set counts; a
gene outside it (an ALT haplotype, a patch) is dropped from `n_genes`
but counted in `off_primary_genes`, and a protein with only such genes
is `off_primary_only`. Every protein gets exactly one `outcome`:
`resolved`, `multi_gene`, `off_primary_only`, `not_in_source`,
`unrecognized_accession` or `contaminant_not_mapped`.

## Usage

``` r
proteins_resolve_genes(databases, gtf = NULL, gene_set = NULL, xref = NULL,
  contaminants = NULL, threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- databases:

  The search database or databases, UniProt XML for real answers. A
  FASTA carries no gene links, so every FASTA protein is
  `not_in_source`; a caveat says so.

- gtf:

  An Ensembl GTF, plain or `.gz`. Only its `gene` rows are read. Give
  exactly one of `gtf` and `gene_set`.

- gene_set:

  Instead of `gtf`, a compact gene table written by mzLib's
  `EnsemblGeneSetWriter` (about 0.5 MB for human, against a 141 MB GTF).
  It carries the GTF's provenance, so rows are keyed exactly as against
  the GTF.

- xref:

  Optional: Ensembl's `Species.Assembly.Release.uniprot.tsv.gz`, for a
  second opinion. Each gene row then says whether Ensembl agrees, and a
  gene only Ensembl links gets its own row (`source == "ensembl_xref"`).
  Without it, `ensembl_xref_agrees` is `NA` - unknown, not false.

- contaminants:

  Databases to load as contaminants. Their proteins are
  `contaminant_not_mapped`: never mapped, and never silently dropped.

- threads:

  Databases read (and hashed) at once. Default 1; `-1` means every core.
  Same answer at any value.

- on_error:

  `"fail"` (the default) or `"skip"`, as for
  [`proteins_read`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md).

- timeout:

  Seconds to allow, or `NULL` (the default) to wait: a gzipped GTF takes
  tens of seconds.

## Details

**You supply the gene set, and you should pin it.** Gene ids, versions,
symbols and the very set of genes on the primary assembly change between
Ensembl releases, so a resolution means something only relative to one
release. Use Ensembl's **primary-assembly** GTF
(`Species.Assembly.Release.gtf.gz`, not the `chr_patch_hapl_scaff` one),
keep Ensembl's file name so the release is recorded, and keep
`gene_set$sha256` with your results. Nothing is downloaded or defaulted
here.

## Value

An `mzlibr_gene_resolutions`. `file_count`, `read_count` and
`failed_count` count database files; `protein_count` counts the proteins
resolved; `record_count` counts the rows of `resolutions` - at least one
per protein, one per gene for `multi_gene`, plus any `ensembl_xref`
rows. `outcome_counts` is a named vector: outcome to number of proteins,
every outcome present, zeros included. `gene_set` is a list
(`source_file_name`, `sha256`, `release`, `genome_build`,
`genebuild_last_updated`, `gene_count`), and `xref` one like it, or
`NULL` without `xref`. `files` has one row per input, with its
`search_database_sha256` (of the decompressed database).

`resolutions` is a data.frame whose columns after `source_index`
(1-based) and `source_path` are mzLib's own `GeneResolutionTsv` schema.
`n_genes` counts genes on the gene set's assembly; `off_primary_genes`
counts linked genes outside it. `gene_id`, `versioned_gene_id` and
`gene_biotype` are `NA` when the row has no gene; `isoform` is `NA`
without a `-N` suffix.

## One row per gene, never a pick

A `multi_gene` protein has one row per gene, and a protein with no gene
one outcome row. No cell joins several genes, and nothing chooses among
them for you.

## Wraps

Wire verb `genes resolve`. Generated from the bridge's verb spec
`genes.resolve.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`EnsemblGeneResolver.Resolve`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneResolver.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneResolver.cs` at
  mzLib `23c2490e`

- [`GeneResolutionTsv.Schema`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneResolver.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneResolver.cs` at
  mzLib `23c2490e`

- [`EnsemblGeneSet.LoadGtf`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSet.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSet.cs` at
  mzLib `23c2490e`

- [`EnsemblGeneSetReader.Load`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSetReader.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSetReader.cs`
  at mzLib `23c2490e`

- [`EnsemblXrefTable.Load`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblXrefTable.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblXrefTable.cs` at
  mzLib `23c2490e`

- [`ProteinDbLoader.LoadProteinXML`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

## Parameters: units, ranges and defaults

- `databases` (wire `--path`):

  path; default absent. One search database (UniProt XML for real
  answers; a FASTA carries no gene links). Exactly one of path and
  paths-stdin.

- `contaminants` (wire `--contaminant`):

  flag; default `FALSE`. Load the `--path` database as contaminants:
  every protein is contaminant_not_mapped.

- wire `--paths-stdin`:

  flag; default `FALSE`. Databases on stdin, one per line, optionally
  followed by a tab and 'contaminant' or 'target' (BULK.md §1). *Not an
  argument here: chosen automatically when more than one database is
  given.*

- `gtf`:

  path; default absent. An Ensembl GTF, plain or .gz; only gene rows are
  read. Use the primary-assembly GTF (Species.Assembly.Release.gtf.gz)
  and keep Ensembl's file name, which carries the release. Exactly one
  of gtf and gene-set; there is no default.

- `gene_set` (wire `--gene-set`):

  path; default absent. Instead of gtf: a compact gene table from
  mzLib's EnsemblGeneSetWriter. It carries the GTF's provenance, so rows
  are keyed as against the GTF (tested byte-identical).

- `xref`:

  path; default absent. Optional Ensembl
  Species.Assembly.Release.uniprot.tsv(.gz): per-row agreement, and a
  row per gene only Ensembl links (source ensembl_xref).

- `threads`:

  int; default `1`; range `>= 1, or -1 for every core`. Databases loaded
  (and hashed) at once. Output identical at any value.

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. As proteins read.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **files**; never `NA`. Databases given.

- `read_count`:

  int; in **files**; never `NA`. Databases read.

- `failed_count`:

  int; in **files**; never `NA`. Databases that failed (on-error skip
  only).

- `protein_count`:

  int; in **proteins**; never `NA`. Proteins resolved.

- `record_count`:

  int; in **rows**; never `NA`. Rows in columns: at least one per
  protein; one per gene for multi_gene; plus ensembl_xref rows.

- `gene_set`:

  object; never `NA`. {source_file_name, sha256 (of the file's bytes as
  read: compressed for .gz), release (null: not in the file name),
  genome_build (null: no \#!genome-build header), genebuild_last_updated
  (null: no header), gene_count}.

- `xref`:

  object; `NA` when no `--xref` was given. {source_file_name, sha256,
  release (nullable), accession_count}.

- `outcome_counts`:

  map\<string,int\>; in **proteins**; never `NA`. Outcome name -\>
  proteins with it; all six keys always present.

- `caveats`:

  string\[\]; never `NA`. FASTA proteins (all not_in_source), a GTF with
  no release or genome build, no xref (agreement unknown, not false).

- `files`:

  object\[\]; never `NA`. As proteins read; search_database_sha256 is
  set: lower-case hex sha256 of the DECOMPRESSED database, so .xml and
  .xml.gz agree.

- `column_names`:

  string\[\]; never `NA`. source_index, source_path, then mzLib's
  GeneResolutionTsv.Schema headers in its order (a bridge test holds
  them equal).

- `resolutions` (wire `columns`):

  table; never `NA`. Column name to per-row values.

## Columns of `resolutions`

- `source_index`:

  int; never `NA`. Input position.

- `source_path`:

  string; never `NA`. Input path.

- `accession`:

  string; never `NA`. Verbatim, as the search reported it (a variant
  suffix kept).

- `entry_accession`:

  string; never `NA`. UniProt entry or unversioned RefSeq accession;
  verbatim when unrecognized.

- `isoform`:

  int; `NA` when no -N isoform suffix. UniProt isoform number.

- `namespace`:

  string; never `NA`. uniprot \| refseq \| unrecognized.

- `outcome`:

  string; never `NA`. resolved \| multi_gene \| off_primary_only \|
  not_in_source \| unrecognized_accession \| contaminant_not_mapped
  (GeneResolutionTsv.OutcomeName).

- `n_genes`:

  int; in **genes**; never `NA`. Genes the search database links on the
  gene set's assembly.

- `gene_id`:

  string; `NA` when the outcome has no gene on the assembly. Stable
  Ensembl gene id.

- `versioned_gene_id`:

  string; `NA` when no gene on the row. As the database wrote it (e.g.
  ENSG00000111640.15); for an ensembl_xref row, built from the GTF's
  version, or null when the GTF has none.

- `gene_symbol`:

  string; `NA` when no gene on the row, or the GTF gives it no
  gene_name. The release's gene_name: a display label, never a key.

- `gene_biotype`:

  string; `NA` when no gene on the row. The release's gene_biotype
  ('unknown' when the GTF omits it).

- `off_primary_genes`:

  int; in **genes**; never `NA`. Linked genes outside the gene set (ALT
  haplotypes, patches, scaffolds): dropped from n_genes, counted here.

- `uniprot_gene_name`:

  string; `NA` when the entry names no primary gene. The database's
  primary gene name, as a label.

- `source`:

  string; never `NA`. search_database_dbreference, or ensembl_xref for a
  gene only Ensembl links (it carries the search database's outcome and
  never changes it).

- `search_database_sha256`:

  string; never `NA`. sha256 of the decompressed database.

- `gene_set_release`:

  string; `NA` when the GTF's file name carries no release. Ensembl
  release.

- `gene_set_sha256`:

  string; never `NA`. sha256 of the GTF as read.

- `ensembl_xref_agrees`:

  bool; `NA` when no gene on the row, or no `--xref`: unknown, not
  false. Whether Ensembl's xref links this accession to this gene.

- `ensembl_xref_info_type`:

  string; `NA` when Ensembl's xref does not link the pair, or no
  `--xref`. Ensembl's strongest evidence: DIRECT, SEQUENCE_MATCH,
  INFERRED_PAIR.

- `ensembl_xref_sha256`:

  string; `NA` when no `--xref` was given. sha256 of the xref file.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  neither or both of gtf and gene-set; a missing GTF, gene set or xref;
  any database-selection error proteins read raises

- `mzlib_bridge_error` (correctness):

  a GTF gene row without gene_id (InvalidDataException); an xref whose
  header is not Ensembl's nine columns; a compact gene table mzLib's
  reader refuses; a database mzLib cannot parse

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- The caller supplies and pins the gene set. Use the primary-assembly
  GTF, not chr_patch_hapl_scaff: counting ALT haplotypes as genes made
  human release 116 look 6.99% multi-gene when it is 0.36% (mzLib
  EnsemblGeneSet remarks).

- A FASTA carries no Ensembl links, so every FASTA protein is
  not_in_source regardless of biology; the result caveats count them.

- Contaminants are never mapped (contaminant_not_mapped), even when they
  would resolve.

- multi_gene is one row per gene, never a pick.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.resolve_genes`

- Rust (mzLibRust): `mzlib::proteins::resolve_genes_with` with
  `GeneResolveOptions`

- R (mzLibR): `proteins_resolve_genes`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Rust and R spellings are intended, not yet implemented. The wire
  module is 'genes' but every binding projects it in its proteins module
  (BULK.md §6); confirm.

- No verb writes the compact gene set (EnsemblGeneSetWriter); callers
  need mzLib or a future 'genes gene-set' verb to make one.

- The example gene set is a hand-built three-gene GTF in Ensembl's
  layout; O43653 is off_primary_only only because PSCA is left out of
  it.

## References

- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl releases, GTFs and cross-references (Ensembl 2026, NAR)

- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, whose XML dbReferences carry the gene links resolved (UniProt
  Consortium 2025, NAR)

## See also

[`proteins_read`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md),
[`proteins_classify_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md)

## Examples

``` r

genes <- proteins_resolve_genes(c("human_subset.xml", "human_extra.fasta"),
  contaminants = "contaminants.fasta",
  gtf = "Homo_sapiens.GRCh38.116.gtf", xref = "Homo_sapiens.GRCh38.116.uniprot.tsv")
genes
#> <mzlibr_gene_resolutions> 7 proteins, 7 rows
#>   gene set: Homo_sapiens.GRCh38.116.gtf, Ensembl release 116
#>   outcomes: resolved 1, off_primary_only 1, not_in_source 4, contaminant_not_mapped 1
#>   ! 4 protein(s) came from a FASTA, which carries no Ensembl links, so each is not_in_source whatever Ensembl knows about it. Resolve against the UniProt XML of the same proteome, or pass --xref to add Ensembl's own links as ensembl_xref rows.
genes$gene_set[c("release", "genome_build", "sha256")]
#> $release
#> [1] "116"
#> 
#> $genome_build
#> [1] "GRCh38.p14"
#> 
#> $sha256
#> [1] "cfe0494115dddf652f4658ebd9e89c0cb77e007d464e67c8d99a2a13cd004fd6"
#> 
genes$outcome_counts
#>               resolved             multi_gene       off_primary_only 
#>                      1                      0                      1 
#>          not_in_source unrecognized_accession contaminant_not_mapped 
#>                      4                      0                      1 
genes$resolutions[, c("accession", "outcome", "n_genes", "gene_id", "gene_symbol")]
#>   accession                outcome n_genes         gene_id gene_symbol
#> 1    P04406               resolved       1 ENSG00000111640       GAPDH
#> 2    O43653       off_primary_only       0            <NA>        <NA>
#> 3    P02768          not_in_source       0            <NA>        <NA>
#> 4    Q13409          not_in_source       0            <NA>        <NA>
#> 5  Q13409-2          not_in_source       0            <NA>        <NA>
#> 6  Q13409-3          not_in_source       0            <NA>        <NA>
#> 7    P02769 contaminant_not_mapped       0            <NA>        <NA>
```
