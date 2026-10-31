# Map residues across species through Ensembl Compara's gene-tree alignments

Map a residue of one protein to the residue of another protein in the
same column of Ensembl Compara's gene-tree alignment, or say why there
is none, for a batch of queries.

Answers **which mouse residue is human p53's Ser15?** Ensembl Compara
has already aligned every protein of each gene tree, across all its
species, to build the tree. A residue's column in that alignment is the
published answer. This wraps mzLib's `ComparaGeneTreeAlignment` (mzLib
\#1465): `Load` reads Compara's peptide alignment of every gene tree,
keeping the canonical proteins of your species' genes, and `MapResidue`
answers each query: which residue of the target protein sits in the same
alignment column as residue `position` of the source protein. Local
files only; nothing is downloaded.

## Usage

``` r
proteins_map_residues(queries, alignment, gene_trees, gtfs = NULL, gene_sets = NULL,
  proteins = FALSE, timeout = NULL)
```

## Arguments

- queries:

  The `(source protein, position, target protein)` queries, one result
  row each, in this order, duplicates included: a data.frame with three
  columns, `source_protein_id`, `position` and `target_protein_id`, or a
  list of `list(source_protein_id, position, target_protein_id)`. The
  ids are Ensembl protein stable ids; the position is a whole number,
  the 1-based residue in the source protein. They travel on stdin after
  the species lines.

- alignment:

  Compara's peptide alignments,
  `Compara.<release>.protein_<collection>.aa.fasta.gz` (908 MB for
  `protein_default` in release 116), from
  `emf/ensembl-compara/homologies/`. Required.

- gene_trees:

  Compara's gene-tree file of the same release,
  `<collection>.GeneTree_content.default.e<release>.txt.gz`. It names
  each gene's tree and its canonical protein. Required.

- gtfs:

  A named character vector from species to its Ensembl primary-assembly
  GTF, as for
  [`proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md):
  `c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf.gz")`. Only the genes
  in these are kept, and so only their canonical proteins. `NULL` (the
  default) with no `gene_sets`: mzLib's default keeps every gene in the
  gene-tree file - millions for the vertebrates - and a caveat says so.
  The answers for your species' proteins are the same.

- gene_sets:

  Instead of a GTF, a compact gene table mzLib's `EnsemblGeneSetWriter`
  made from one, by species in the same way. A species may be in `gtfs`
  or `gene_sets`, not both.

- proteins:

  Also return the kept proteins with their ungapped sequences, as
  `proteins`. `FALSE` (the default): `proteins` is `NULL`.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait: for human, mouse
  and rat at release 116 the whole call takes about a minute, nearly all
  of it reading the alignment.

## Details

Every query gets exactly one outcome, mzLib's own name for it:
`"Aligned"` (with the column, the target position and its residue),
`"GapInTarget"` (no target residue in that column - never the nearest
one), `"DifferentAlignments"` (the proteins are in different gene
trees), or `"SourceNotInAlignment"` / `"TargetNotInAlignment"` (that
protein was not kept: only each gene's canonical protein is aligned).

**Positions are Ensembl's.** To start from a UniProt position, align the
UniProt sequence to the canonical Ensembl protein with
[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
first, and align the target's Ensembl protein to its UniProt sequence
after. `proteins = TRUE` returns the Ensembl sequences for that.

**A position outside the source protein fails the whole call.** mzLib's
`MapResidue` throws for it rather than returning an outcome, and mzLibR
passes positions through unchecked, so keep them within each source
protein: with `proteins = TRUE`, a sequence's length is its last valid
position.

## Value

An `mzlibr_residue_mappings`.

`records` is a data.frame with one row per query, in the order given:
`source_protein_id`, `position` (the 1-based residue in the source, as
given) and `target_protein_id`, as given; `outcome`, one of the five
outcomes above; `column`, the source residue's 1-based column in the
gene tree's alignment, `NA` when the two proteins share no alignment;
`target_position`, the 1-based position of the target's residue in that
column, in the target's ungapped sequence, `NA` unless `"Aligned"`; and
`target_residue`, the target's residue there, one letter, `NA` unless
`"Aligned"`. It may be a different amino acid from the source's: a
column pairs homologous positions, not identical ones.

`row_count` counts the queries, one row each (queries). `outcome_counts`
is a named vector, outcome to the queries with it (queries), every
outcome present in mzLib's order, zeros included. `column_names` lists
the table's columns.

`alignment` is a list describing the alignment file: `source_file_name`;
`sha256`, of the file's bytes as read (the compressed bytes for a
`.gz`); `release` and `collection`, read from the file name (`NA` when
the name carries none); `alignment_count`, the alignments in the file,
one per gene tree, kept proteins or not (alignments); `protein_count`,
the canonical proteins kept (proteins); and
`restricted_to_gene_tree_content`, the sha256 of the gene-tree file
whose canonical proteins were kept, the same as `gene_trees$sha256`.

`gene_trees` is a list, as
[`proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)
reports it: `source_file_name`, `sha256`, `release`, `collection`,
`member_count`, the genes kept from the file (genes), and
`restricted_to_gene_sets`, the gene sets' sha256s, or `NULL` when no
gene set was given and every gene was kept.

`gene_sets` is a data.frame, one row per species given (none without
`gtfs` or `gene_sets`): `species`, `source_file_name`, `sha256`,
`release`, `genome_build` and `genebuild_last_updated` (`NA` when the
GTF's name or headers carry none), and `gene_count`, in genes, every
biotype.

`proteins_included` says whether `proteins = TRUE` was given. `proteins`
is `NULL` without it; with it, a list: `row_count` (proteins kept),
`column_names`, and `records`, a data.frame with one row per kept
protein in order of protein id: `protein_id`, the Ensembl protein stable
id; `alignment_index`, which alignment holds it (**1-based in R**, in
file order; two proteins share a column only when they share this);
`tree_id`, its gene's Compara gene tree; and `sequence`, its residues
with the alignment's gaps removed, the sequence every position counts
in.
`setNames(m$proteins$records$sequence, m$proteins$records$protein_id)`
gives the sequences by id, for
[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md).
`caveats` names this call's traps: gene trees read unrestricted, files
named for different releases or for none, queries naming a protein the
alignment did not keep.

## What each outcome means

`"Aligned"`: the target has a residue in the source residue's column;
`column`, `target_position` and `target_residue` are set.
`"GapInTarget"`: the target has a gap in that column, so no residue
corresponds; `column` is set, the target's are `NA`.
`"DifferentAlignments"`: the two proteins are in different gene trees,
so this alignment says nothing about them. It does not say they are
unrelated. `"SourceNotInAlignment"` and `"TargetNotInAlignment"`: that
protein was not kept - not in the file, not its gene's canonical
protein, or outside your gene sets. Another isoform is never mapped
through the canonical one.

## The traps

**A position past the end of the source protein fails the whole call**
with mzLib's message (for example "ENSP00000269305 has 393 residues."),
as an `mzlib_usage_error`, and no row is returned, not even for the
queries that were fine. A source protein that was not kept is an
outcome, `"SourceNotInAlignment"`, checked before the position. **Use
Ensembl ids and Ensembl positions**: a UniProt accession is
`"SourceNotInAlignment"`, never an error. **Use the alignment and the
gene trees of one release**: mzLib does not compare them, and a gene's
canonical protein and tree can change between releases. Keep Ensembl's
file names, so the result can report the release each one carries.
**Sharing a gene tree is not orthology**: a tree holds paralogs too, so
ask
[`proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)
which genes are orthologs, then map their residues here. \*\*Check
`target_residue`\*\* before carrying a modification across: a
phosphosite on S is not a site on A.

## Wraps

Wire verb `genes map-residues`. Generated from the bridge's verb spec
`genes.map-residues.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`ComparaGeneTreeAlignment.Load`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs)
  in
  `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs`
  at mzLib `cc67b7a2`

- [`ComparaGeneTreeAlignment.MapResidue`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs)
  in
  `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs`
  at mzLib `cc67b7a2`

- [`ComparaGeneTreeAlignment.Proteins`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs)
  in
  `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs`
  at mzLib `cc67b7a2`

- [`ComparaColumnOutcome`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs)
  in
  `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeAlignment.cs`
  at mzLib `cc67b7a2`

- [`ComparaGeneTreeContent.Load`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeContent.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeContent.cs`
  at mzLib `cc67b7a2`

- [`EnsemblGeneSet.LoadGtf`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSet.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSet.cs` at
  mzLib `cc67b7a2`

- [`EnsemblGeneSetReader.Load`](https://github.com/smith-chem-wisc/mzLib/blob/cc67b7a29/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSetReader.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSetReader.cs`
  at mzLib `cc67b7a2`

## Parameters: units, ranges and defaults

- `alignment`:

  path; required. Compara's peptide alignments of its gene trees,
  emf/ensembl-compara/homologies/Compara.RELEASE.protein_COLLECTION.aa.fasta.gz
  (plain or .gz): one gapped multiple alignment per tree, as FASTA,
  separated by a // line. Read with ComparaGeneTreeAlignment.Load, which
  keeps only the canonical proteins of the genes kept from gene-trees,
  and checks every row's characters and width whether it is kept or not.
  908 MB compressed for protein_default in release 116.

- `gene_trees` (wire `--gene-trees`):

  path; required. Compara's gene-tree membership dump,
  COLLECTION.GeneTree_content.default.eRELEASE.txt.gz (plain or .gz),
  read with ComparaGeneTreeContent.Load. It names each gene's tree and
  its canonical protein (the row flagged Y), which is the only protein
  of the gene the alignment holds.

- wire `--species`:

  string\[\]; default absent; range
  `>= 0 lines; each exactly three tab-separated fields`. Read from stdin
  (not an option), first: one line per species, the species' Ensembl
  production name (homo_sapiens), a tab, gtf or gene-set, a tab, and the
  path, exactly as genes orthology reads them. gtf is an Ensembl
  primary-assembly GTF, read with EnsemblGeneSet.LoadGtf; gene-set is
  the compact table mzLib's EnsemblGeneSetWriter makes from one, read
  with EnsemblGeneSetReader.Load. Their only use is
  ComparaGeneTreeContent.Load's restrictTo: only the genes in them are
  kept, and so only their canonical proteins. None (the stdin starts
  with –): mzLib's default, restrictTo null, keeps every gene in the
  gene-tree file and every canonical protein in the alignment, millions
  for the vertebrate collection, and a caveat says so. A species listed
  twice, a kind other than gtf or gene-set, or a missing file is a usage
  error. *Not an argument here: rendered from gtfs and gene_sets, two
  named character vectors from Ensembl production name to path, one
  stdin line per species (none: mzLib keeps every gene).*

- `queries`:

  string\[\]; required; range
  `>= 1 line; each exactly three tab-separated fields`. Read from stdin
  (not an option), after a line holding only –: one query per line, the
  source protein's Ensembl stable id, a tab, a 1-based residue position
  in it, a tab, and the target protein's id (ENSP00000269305, 15,
  ENSMUSP00000104298). One row per query, in input order, duplicates
  included. The position is parsed as an integer and passed to
  MapResidue unchecked; the bridge does not compare it with the
  protein's length.

- `proteins`:

  flag; default `FALSE`. Also return the kept proteins
  (ComparaGeneTreeAlignment.Proteins: id, alignment, tree and ungapped
  sequence), so a protein's Ensembl sequence can be aligned to a UniProt
  sequence with proteins align without a second file. Absent: proteins
  is null.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `alignment`:

  object; never `NA`. The alignment file and what was kept from it
  (typed under result.tables.alignment).

- `gene_trees`:

  object; never `NA`. The gene-tree file and what was kept from it, the
  same fields as genes orthology's gene_trees (typed under
  result.tables.gene_trees).

- `gene_sets`:

  map\<string,object\>; never `NA`. Species to its gene set, in ordinal
  order of species, each with genes resolve's gene_set fields (typed
  under result.tables.gene_sets). Empty when no species lines were
  given.

- `row_count`:

  int; in **queries**; never `NA`. Rows: the queries, one each.

- `outcome_counts`:

  map\<string,int\>; in **queries**; never `NA`. Outcome to queries with
  it, every ComparaColumnOutcome in mzLib's declaration order, zeros
  included.

- `column_names`:

  string\[\]; never `NA`. The mapping table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per query, in input order.

- `proteins_included`:

  bool; never `NA`. Whether `--proteins` was given.

- `proteins`:

  object; `NA` when `--proteins` was not given. {row_count,
  column_names, columns}: one row per kept protein, in ordinal order of
  protein id (typed under result.tables.proteins).

- `caveats`:

  string\[\]; never `NA`. This call's traps: gene trees read
  unrestricted (every gene kept), an alignment and gene-tree file named
  for different releases or for none, queries naming a protein the
  alignment did not keep.

## Columns of `records`

- `source_protein_id`:

  string; never `NA`. The query's source protein id, as given.

- `position`:

  int; never `NA`. The query's 1-based residue position in the source
  protein, as given.

- `target_protein_id`:

  string; never `NA`. The query's target protein id, as given.

- `outcome`:

  string; never `NA`. mzLib's ComparaColumnOutcome name, verbatim:
  Aligned (the target has a residue in the source residue's column),
  GapInTarget (the target has a gap there: no residue corresponds, never
  the nearest one), DifferentAlignments (the two proteins are in
  different gene trees, so the alignment says nothing about them),
  SourceNotInAlignment or TargetNotInAlignment (that protein was not
  kept: not in the file, not a canonical protein, or outside the gene
  sets). Exactly one per query.

- `column`:

  int; `NA` when the two proteins do not share an alignment
  (DifferentAlignments, SourceNotInAlignment or TargetNotInAlignment).
  The source residue's 1-based column in the gene tree's alignment. Set
  for Aligned and GapInTarget.

- `target_position`:

  int; `NA` when the outcome is not Aligned. The 1-based position of the
  target's residue in that column, in the target's ungapped sequence.

- `target_residue`:

  string; `NA` when the outcome is not Aligned. The target's residue
  there, one letter as the alignment writes it. It may differ from the
  source's residue: a column pairs homologous positions, not identical
  ones.

## Fields of `alignment`

a list.

- `alignment$source_file_name` (wire `source_file_name`):

  string; never `NA`. The file's name
  (ComparaGeneTreeAlignment.SourceFileName).

- `alignment$sha256` (wire `sha256`):

  string; never `NA`. Lower-case hex sha256 of the file's bytes as read,
  the compressed bytes for a .gz (SourceSha256).

- `alignment$release` (wire `release`):

  string; `NA` when the file name does not start Compara.RELEASE.. The
  Ensembl release parsed from the file name (Release).

- `alignment$collection` (wire `collection`):

  string; `NA` when the file name does not match
  Compara.RELEASE.COLLECTION.aa.fasta. The collection parsed from the
  file name, e.g. protein_default (Collection).

- `alignment$alignment_count` (wire `alignment_count`):

  int; in **alignments**; never `NA`. Alignments (gene trees) in the
  file, kept proteins or not (AlignmentCount).

- `alignment$protein_count` (wire `protein_count`):

  int; in **proteins**; never `NA`. Proteins kept (Count): the canonical
  proteins of the genes kept from gene-trees that the file holds.

- `alignment$restricted_to_gene_tree_content` (wire
  `restricted_to_gene_tree_content`):

  string; never `NA`. The sha256 of the gene-tree file whose canonical
  proteins were kept (RestrictedToGeneTreeContent); the same as
  gene_trees.sha256.

## Fields of `gene_trees`

a list, as proteins_build_orthology() reports it.

- `gene_trees$source_file_name` (wire `source_file_name`):

  string; never `NA`. The file's name.

- `gene_trees$sha256` (wire `sha256`):

  string; never `NA`. Lower-case hex sha256 of the file's bytes as read
  (compressed for a .gz).

- `gene_trees$release` (wire `release`):

  string; `NA` when the file name carries no .eRELEASE.txt suffix. The
  release parsed from the file name.

- `gene_trees$collection` (wire `collection`):

  string; `NA` when the file name does not start
  COLLECTION.GeneTree_content.. The Compara collection, e.g.
  vertebrates.

- `gene_trees$member_count` (wire `member_count`):

  int; in **genes**; never `NA`. Genes kept from the file
  (ComparaGeneTreeContent.Count): those in the gene sets, or every gene
  when none was given.

- `gene_trees$restricted_to_gene_sets` (wire `restricted_to_gene_sets`):

  string\[\]; `NA` when no species lines were given: every gene in the
  file was kept. The sha256s of the gene sets the members were
  restricted to, in ordinal order (RestrictedToGeneSets).

## Fields of `gene_sets`

a data.frame, one row per species, `species` first.

- `source_file_name`:

  string; never `NA`. The GTF's file name; for a gene-set table, still
  the GTF's.

- `sha256`:

  string; never `NA`. Lower-case hex sha256 of the GTF's bytes as read;
  the value genes resolve and genes orthology report for the same GTF.

- `release`:

  string; `NA` when the GTF's file name carries no release. The
  annotation release from the GTF's name.

- `genome_build`:

  string; `NA` when the GTF has no \#!genome-build header. e.g.
  GRCh38.p14.

- `genebuild_last_updated`:

  string; `NA` when the GTF has no \#!genebuild-last-updated header. The
  header's value.

- `gene_count`:

  int; in **genes**; never `NA`. Genes in the set, every biotype.

## Fields of `proteins$records`

a data.frame, one row per kept protein, inside the proteins list.

- `protein_id`:

  string; never `NA`. The Ensembl protein stable id, as the alignment's
  header gives it (ComparaAlignedProtein.ProteinId).

- `alignment_index`:

  int; never `NA`. Which alignment holds it, 0-based in file order
  (AlignmentIndex). Two proteins share a column only when they share
  this.

- `tree_id`:

  string; never `NA`. Its gene's Compara gene tree, e.g.
  ENSGT00950000183153, from gene-trees (TreeId).

- `sequence`:

  string; never `NA`. The protein's residues with the alignment's gaps
  removed (Sequence): the sequence positions count in.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--alignment` or `--gene-trees` missing or given without a value; a
  value after `--proteins`

- `mzlib_usage_error` (usage):

  no line holding only – on stdin, or no queries after it; a species
  line that is not species, tab, gtf or gene-set, tab, path; a species
  listed twice; a query that is not three tab-separated fields, or whose
  position is not an integer

- `mzlib_usage_error` (usage):

  a GTF, gene-set table, gene-trees file or alignment that does not
  exist (checked before any input is read)

- `mzlib_usage_error` (usage):

  mzLib's MapResidue refused a query's position as outside the source
  protein (ArgumentOutOfRangeException): 0, negative, or past its last
  residue. The message names the query and quotes mzLib, which gives the
  protein's length. The whole call fails on the first such query, in
  input order, and no row is returned. A source protein that is not kept
  is an outcome, SourceNotInAlignment, checked before the position.

- `mzlib_bridge_error` (correctness):

  a file mzLib's reader refuses (InvalidDataException, naming file and
  line): in the alignment, a sequence line before any header, an empty
  header, a character that is not a letter, \* or -, rows of different
  widths, an empty alignment, a kept protein twice, or an alignment
  whose kept proteins belong to two gene trees; in the gene-tree file, a
  row without four cells, a flag other than Y or N, or a kept gene in
  two trees or without exactly one canonical protein; a GTF gene row
  without gene_id

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- A column pairs positions that Compara's alignment puts together, which
  is homology, not identity: the target residue may be a different amino
  acid, so check target_residue before carrying a modification across (a
  phosphosite on S is not a site on A).

- Only each gene's canonical protein is in the alignment (the gene-tree
  row flagged Y). A query naming another isoform's Ensembl id is
  SourceNotInAlignment or TargetNotInAlignment, never mapped through the
  canonical one. To start from a UniProt position, align the UniProt
  sequence to the canonical Ensembl protein with proteins align first,
  then map, then align the target's canonical protein to its UniProt
  sequence: three steps, each of which can meet a gap.

- A gap is GapInTarget, never the nearest residue: no residue of the
  target sits in that column.

- DifferentAlignments means the two proteins are in different gene
  trees, so this alignment says nothing about them. It does not say the
  genes are not homologous, nor that they are orthologs when they share
  one: a tree holds paralogs as well, so read orthology from genes
  orthology.

- mzLib's MapResidue throws for a position outside the source protein
  instead of returning an outcome, so one such query fails the whole
  call, and the bridge does not check lengths for it. Keep positions
  inside each source protein; `--proteins` gives each kept protein's
  sequence.

- mzLib does not compare the alignment's release with the gene-tree
  file's. Use the two files of one release: a gene's canonical protein
  and its tree can change between releases. The envelope reports the
  release each file name carries, and a caveat says when they differ or
  one is absent.

- Every input is read whole, and the alignment keeps only the canonical
  proteins of the genes kept, so give the species' gene sets: with none,
  every gene in the gene-tree file is kept (millions for the vertebrate
  collection), and so every canonical protein of the alignment.

## Performance

ComparaGeneTreeAlignment.Load reads the whole alignment file once,
checking every row, and is the cost; mapping a query afterwards is a
binary search. For human, mouse and rat at release 116 (64,114 proteins
kept from 54,308 alignments), mzLib \#1465 reports about 33 s to load,
and this bridge measured 55 s for the whole call on 2026-10-28, reading
the three GTFs, the gene-tree file and the alignment. Give all of a
batch's queries to one call. Neither load takes a degree of parallelism,
so there is no `--threads`. Every call starts one bridge process, which
costs a .NET start-up before any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.map_residues`

- Rust (mzLibRust): `mzlib::proteins::map_residues_with` with
  `MapResiduesOptions`

- R (mzLibR): `proteins_map_residues`

## Since

Wire protocol 1; pyMzLib 0.7.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The outcome is mzLib's ComparaColumnOutcome name verbatim (Aligned,
  GapInTarget, ...), as classify-peptides projects PeptideSharing. genes
  orthology's statuses are snake_case because mzLib names them
  (OrthologySnapshot.StatusName). PUSH-TO-MAINLAND: a
  ComparaGeneTreeAlignment.OutcomeName helper in the same style; when it
  lands the wire names change to its spelling, a breaking change to call
  out.

- MapResidue throws ArgumentOutOfRangeException for a position outside
  the source, which fails a whole batch. PUSH-TO-MAINLAND: a
  non-throwing outcome (PositionOutsideSource) so one bad query is one
  row. Until then the bridge projects the throw and does not check
  lengths itself.

- mzLib compares no release between the alignment and the gene-tree
  content (OrthologySnapshot.Build does, for the dumps).
  PUSH-TO-MAINLAND: ComparaGeneTreeAlignment.Load could refuse a
  gene-tree content named for another release, as Build refuses one; the
  bridge only reports it as a caveat.

- result.tables (alignment, gene_trees, gene_sets, proteins) is not
  validated by tools/check_verbs.py, which checks only envelope_fields,
  columns and bulk. Bridge tests pin each nested object's keys
  (GenesMapResiduesTests).

- Not yet ported to mzLibRust or mzLibR; the binding names above are
  proposals.

## References

- [doi:10.1101/gr.073585.107](https://doi.org/10.1101/gr.073585.107):
  Ensembl Compara gene trees, each built on a multiple alignment of its
  member proteins, the alignment this verb reads columns from (Vilella
  et al. 2009, Genome Res 19:327)

- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl releases, their GTFs and the Compara gene-tree files (Ensembl
  2026, NAR)

## See also

[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
for the UniProt legs,
[`proteins_build_orthology`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)

## Examples

``` r

# A real cut of Ensembl release 116: the TP53 and GAPDH gene trees' alignments, with the human,
# mouse and rat rows exactly as Ensembl publishes them.
gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
          mus_musculus = "Mus_musculus.GRCm39.116.gtf",
          rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
site <- proteins_map_residues(list(list("ENSP00000269305", 15, "ENSMUSP00000104298")),
  alignment = "Compara.116.protein_default.aa.fasta.gz",
  gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz", gtfs = gtfs)
site
#> <mzlibr_residue_mappings> 1 queries through Compara.116.protein_default.aa.fasta.gz (Ensembl 116, 5 proteins kept)
#>   Aligned 1
# Human p53's Ser15 is mouse p53's Ser18:
site$records[, c("outcome", "column", "target_position", "target_residue")]
#>   outcome column target_position target_residue
#> 1 Aligned    476              18              S
c(site$alignment$release, site$gene_trees$collection)
#> [1] "116"         "vertebrates"
```
