# Build an Ensembl Compara orthology snapshot for your species, and write it as Parquet

Build an Ensembl Compara orthology snapshot for a set of species from
their gene sets, the gene trees and every species' homology dump, and
optionally write it as Parquet and say, gene by gene, how far Compara's
inference got in another species.

Answers **which mouse genes are orthologous to this human gene, and if
none, why not?** Wraps mzLib's `OrthologySnapshot.Build` and
`OrthologySnapshotWriter.WriteAsync` (mzLib \#1381): which genes of one
species are orthologous to which genes of another, kept exactly as
Compara states it. One-to-many and many-to-many orthologies are several
rows, paralogs are kept, and nothing is picked, merged or collapsed.
Local files only; nothing is downloaded.

## Usage

``` r
proteins_build_orthology(release, gene_trees, dumps, gtfs = NULL, gene_sets = NULL,
  restrict_gene_trees = FALSE, out = NULL, pair_status = NULL, timeout = NULL)
```

## Arguments

- release:

  The Ensembl release of the Compara inputs, e.g. `"116"` (a whole
  number is sent as its digits). Required: mzLib holds the gene trees
  and every dump to it when their file names carry a release. A GTF's
  own number is not compared with it (Ensembl publishes some species
  with their own annotation release).

- gene_trees:

  Compara's gene-tree file,
  `<collection>.GeneTree_content.default.e<release>.txt.gz`. It tells a
  gene whose tree holds no gene of the other species from one whose tree
  does but got no ortholog call.

- dumps:

  One Compara homology dump per species
  (`homologies/<species>/Compara.<release>.protein_default.homologies.tsv.gz`),
  in any order: each dump's species is read from its rows. Sent on stdin
  after the species lines.

- gtfs:

  A named character vector from species to its Ensembl primary-assembly
  GTF, plain or `.gz`:
  `c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf.gz")`. The names are
  Ensembl production names, as the dumps' `species` column spells them,
  and they are the snapshot's species. Sent on stdin, one line per
  species.

- gene_sets:

  Instead of a GTF, a compact gene table mzLib's `EnsemblGeneSetWriter`
  made from one, by species in the same way (it keeps the GTF's
  provenance, so the snapshot is the same). A species may be in `gtfs`
  or `gene_sets`, not both.

- restrict_gene_trees:

  Keep only the gene-tree members in your gene sets. `FALSE` (the
  default) is mzLib's default and keeps every gene in the file - 4.2
  million for the vertebrates of release 116. The snapshot is the same
  either way; only memory differs.

- out:

  A directory to write the snapshot to: Parquet tables (`genes/`,
  `members/`, `pairs/`), `views.sql` and `manifest.json`. It must not
  exist yet or be empty: mzLib writes a snapshot once and never updates
  one in place. `NULL` (the default): nothing is written.

- pair_status:

  `c(from_species, to_species)`: return every gene of `from_species`
  with its status in `to_species`. Two different species of the
  snapshot, in that order: a status is not symmetric. `NULL` (the
  default): not computed.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait: for human, mouse
  and rat at release 116, mzLib reports about 20 s to build and 40 s
  more to write.

## Details

**Every species' homology dump is required, and mzLib refuses a call
without one.** Compara puts each of a pair's relationships in one of the
pair's two dumps, arbitrarily, and in practice they partition: in
release 116 every human-mouse orthology is in the \*mouse\* dump and
none is in the human one. Reading the dump named after your species
gives a complete-looking, badly incomplete answer.

**Orthology is pairwise and not transitive.** A~B and B~C do not make
A~C, so never chain pairs into groups. The written `views.sql` has a
`species_set3` macro that forms a tuple only when every two members are
orthologs.

The relationships are in the Parquet files `out` writes, which mzLib
cannot read back and mzLibR does not read (it has no dependencies):
query them with the `arrow` or `duckdb` packages, as the proteins
article shows. Without `out`, the call checks that the inputs fit
together and, with `pair_status`, says gene by gene how far Compara's
inference got in one other species.

## Value

An `mzlibr_orthology_snapshot`. The relationships themselves are in the
Parquet files `out` writes, not here: this says what went in, what came
out, and, with `pair_status`, every gene's status in one other species.

`release` is the release you gave; `species` the snapshot's species, in
sorted order; `format` is `"ensembl-orthology-snapshot"` and
`format_version` its version, as `manifest.json` names them.
`snapshot_id` is mzLib's id for the snapshot, a sha256 over the format,
the release, the species and every input's sha256: the same inputs give
the same id, written or not. Record it with any result drawn from the
snapshot. `identical_duplicates_dropped` counts homology rows found in
both of a pair's dumps, identical, and counted once, in rows: 0 in
release 116, whose dumps partition.

`gene_sets` is a data.frame, one row per species: `species`,
`source_file_name`, `sha256` (the value
[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
reports as `gene_set_sha256`, so a resolution and a snapshot can be
checked against one gene set), `release`, `genome_build` and
`genebuild_last_updated` (`NA` when the GTF's name or headers carry
none), and `gene_count`, in genes, every biotype.

`gene_trees` is a list: `source_file_name`, `sha256`, `release` (`NA`
when the name carries none, so it was not checked), `collection` (e.g.
`"vertebrates"`), `member_count`, the genes held from the file, and
`restricted_to_gene_sets`, the gene sets' sha256s the members were
restricted to, or `NULL` when every gene in the file was kept.

`dumps` is a data.frame, one row per species, the species read from the
dump's rows: `source_file_name` (the same for every species, so it alone
does not say whose dump it is), `sha256`, `release`, `rows_read`, every
data row read and checked, in rows, and `rows_kept`, the rows whose two
species are both snapshot species, in rows.

`pairs` is a data.frame, one row per species pair `(a, b)` with `a`
sorting first, `(s, s)` included: `species_a`, `species_b`, and the
relationships between them in rows - `row_count`, one row each, never
merged or picked; `ortholog_count`, always 0 for `(s, s)`;
`paralog_count`, always 0 between two species; and `homoeolog_count`.
See
[`proteins_orthology_pair`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_orthology_pair.md).

`written` is `NULL` without `out`. With it, a list: `directory`,
absolute; `manifest`, the path of `manifest.json`, which lists every
input and file with its sha256; and `files`, a data.frame of every file
written but the manifest, in order of path: `path` (relative to
`directory`, with forward slashes), `rows` (rows in the Parquet table,
`NA` for `views.sql`, 0 for a pair with no relationship, which is still
written), `bytes` (the file's size in bytes) and `sha256`.

`pair_status` is `NULL` without `pair_status`. With it, a list:
`from_species`, `to_species`, `row_count` (every gene of `from_species`,
in genes), `status_counts` (a named vector, status to genes with it,
every status present in mzLib's order, zeros included), `column_names`,
and `records`, a data.frame with one row per gene in order of gene id:
`gene_id`, `gene_name` (`NA` when the GTF gives none: a label, never a
key), `gene_biotype` and `status`. `caveats` names this call's traps.

## What each status means

It says how far Compara's inference got, not what biology holds.
`has_ortholog`: at least one ortholog in `to_species`.
`no_edge_in_shared_tree`: the gene's tree holds a `to_species` gene, and
Compara called no ortholog to it - the nearest the data comes to "no
ortholog exists". `tree_lacks_target_species`: the gene's tree holds no
`to_species` gene, so nothing was compared. `not_in_any_tree`: the gene
is in the gene set but in no gene tree, so no call was attempted.
**Every biotype is a row**, so choose your denominator by filtering on
`gene_biotype`. One status's genes are
`ps$records$gene_id[ps$records$status == "not_in_any_tree"]`, with
`ps <- snap$pair_status`.

## The traps

**The gene set and the dumps must be one release.** A homology naming a
gene its gene set lacks is refused as an `mzlib_bridge_error`
(`InvalidDataException`): that is what a GTF and dumps from different
releases look like. Use the primary-assembly GTF, as for
[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md).
**Keep Ensembl's file names**: mzLib holds the gene trees and the dumps
to `release` only when their names carry it, and a caveat says which
were not checked. **A directory that is not empty is refused** as `out`:
a snapshot is written once, and next release's goes in a new directory.

## Wraps

Wire verb `genes orthology`. Generated from the bridge's verb spec
`genes.orthology.yaml` (bridge commit `60b1d25968a1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`OrthologySnapshot.Build`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshot.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshot.cs` at
  mzLib `d10017ce`

- [`OrthologySnapshot.PairStatus`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshot.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshot.cs` at
  mzLib `d10017ce`

- [`OrthologySnapshot.StatusName`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshot.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshot.cs` at
  mzLib `d10017ce`

- [`OrthologySnapshotWriter.WriteAsync`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshotWriter.cs)
  in
  `mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshotWriter.cs`
  at mzLib `d10017ce`

- [`OrthologySnapshotWriter.SnapshotId`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshotWriter.cs)
  in
  `mzLib/UsefulProteomicsDatabases/Ensembl/OrthologySnapshotWriter.cs`
  at mzLib `d10017ce`

- [`ComparaHomologyDump.Load`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaHomologyDump.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaHomologyDump.cs` at
  mzLib `d10017ce`

- [`ComparaGeneTreeContent.Load`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeContent.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/ComparaGeneTreeContent.cs`
  at mzLib `d10017ce`

- [`EnsemblGeneSet.LoadGtf`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSet.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSet.cs` at
  mzLib `d10017ce`

- [`EnsemblGeneSetReader.Load`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSetReader.cs)
  in `mzLib/UsefulProteomicsDatabases/Ensembl/EnsemblGeneSetReader.cs`
  at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `release`:

  string; required; range `an Ensembl release number, e.g. 116`. The
  Ensembl release of the Compara inputs. Required, never defaulted: the
  release is part of the answer. mzLib holds the gene trees and every
  dump to it, but only when their file names carry a release
  (Compara.116.protein_default.homologies.tsv.gz,
  vertebrates.GeneTree_content.default.e116.txt.gz); a renamed file is
  not checked, and the result's caveats say which. A GTF's own number is
  not compared with it: Ensembl publishes some species with their own
  annotation release (Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz in
  release 116), so mzLib holds a gene set to the dumps by content
  instead.

- wire `--species`:

  string\[\]; required; range
  `>= 1 line; each exactly three tab-separated fields`. Read from stdin
  (not an option), first: one line per species, the species' Ensembl
  production name (homo_sapiens, as the dumps' species column spells
  it), a tab, gtf or gene-set, a tab, and the path. gtf is an Ensembl
  primary-assembly GTF, plain or .gz, read with EnsemblGeneSet.LoadGtf;
  gene-set is the compact table mzLib's EnsemblGeneSetWriter makes from
  one, read with EnsemblGeneSetReader.Load, which keeps the GTF's
  provenance. These lines are the snapshot's species. A species listed
  twice, a kind other than gtf or gene-set, or a missing file is a usage
  error. Blank lines are skipped. *Not an argument here: rendered from
  gtfs and gene_sets, two named character vectors from Ensembl
  production name to path, one stdin line per species.*

- `dumps`:

  string\[\]; required; range `exactly one per species`. Read from stdin
  (not an option), after a line holding only –: one Compara homology
  dump path per line
  (homologies/SPECIES/Compara.RELEASE.protein_default.homologies.tsv.gz,
  plain or .gz), each read with ComparaHomologyDump.Load for every
  species given. The dump's genome is read from its rows, not its name.
  mzLib requires EVERY species' dump: Compara puts each of a pair's
  relationships in one of the pair's two dumps, arbitrarily, so a
  missing dump would silently lose them. A missing dump, two of one
  species, or one of a species not given is mzLib's refusal, a usage
  error.

- `gene_trees` (wire `--gene-trees`):

  path; required. Compara's gene-tree membership dump,
  COLLECTION.GeneTree_content.default.eRELEASE.txt.gz (plain or .gz),
  read with ComparaGeneTreeContent.Load. It is what tells a gene whose
  tree holds no gene of the other species from one whose tree does but
  got no ortholog call.

- `restrict_gene_trees` (wire `--restrict-gene-trees`):

  flag; default `FALSE`. Keep only the gene-tree members in the given
  gene sets (ComparaGeneTreeContent.Load with restrictTo, the gene
  sets). Absent: mzLib's default, restrictTo null, which keeps every
  gene in the file, 4.2 million for the vertebrate collection of release
  116 (mzLib's ComparaGeneTreeContent remarks). The snapshot, its id and
  every written byte are the same either way; only memory,
  gene_trees.member_count and gene_trees.restricted_to_gene_sets differ.
  mzLib checks the one-tree and one-canonical-protein rules only on the
  genes kept.

- `out`:

  path; default absent; range
  `a directory that does not exist or is empty`. Write the snapshot to
  this directory with OrthologySnapshotWriter.WriteAsync: Parquet tables
  (zstd) genes/SPECIES.parquet, members/SPECIES.parquet and
  pairs/A\_\_B.parquet, plus views.sql (DuckDB table macros) and
  manifest.json. The directory and its subfolders are created. A
  DIRECTORY, not a .tsv file: the PYB-1 extension rule is for options
  that write one table. mzLib refuses a directory that is not empty,
  because a snapshot is written once and compared with a later one,
  never updated in place; the bridge applies that same refusal, and
  refuses a path that is a file, before any input is read. Absent:
  nothing is written.

- `pair_status` (wire `--pair-status`):

  string; default absent; range
  `FROM,TO: two different species of the snapshot`. Return pair_status:
  every gene of FROM's gene set, in ordinal order of gene id, with its
  status in TO (OrthologySnapshot.PairStatus). Two production names
  separated by a comma, e.g. homo_sapiens,mus_musculus. Both must be
  species given on stdin, and different; the bridge checks that before
  any input is read, with the rule mzLib's PairStatus applies. One
  direction per call: status is not symmetric. Absent: pair_status is
  null.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `release`:

  string; never `NA`. The `--release` given (OrthologySnapshot.Release).

- `species`:

  string\[\]; never `NA`. The snapshot's species, Ensembl production
  names, in ordinal order (OrthologySnapshot.Species).

- `format`:

  string; never `NA`. ensembl-orthology-snapshot
  (OrthologySnapshotWriter.FormatName): what manifest.json's format
  says.

- `format_version`:

  int; never `NA`. 1 at d10017ce
  (OrthologySnapshotWriter.FormatVersion).

- `snapshot_id`:

  string; never `NA`. OrthologySnapshotWriter.SnapshotId: lower-case hex
  sha256 over the format, the release, the species and every input's
  sha256. The same inputs give the same id, written or not; with `--out`
  it is manifest.json's snapshot_id. Record it with any result drawn
  from the snapshot.

- `identical_duplicates_dropped`:

  int; in **rows**; never `NA`. Homology rows met in both of a pair's
  dumps, identical in every value, and counted once
  (OrthologySnapshot.IdenticalDuplicatesDropped). 0 in release 116,
  whose dumps partition; counted so a change is visible.

- `gene_sets`:

  map\<string,object\>; never `NA`. Species to its gene set, in species
  order; each value has genes resolve's gene_set fields (typed under
  result.tables.gene_sets).

- `gene_trees`:

  object; never `NA`. The gene-tree file and what was kept from it
  (typed under result.tables.gene_trees).

- `dumps`:

  map\<string,object\>; never `NA`. Species (the dump's genome, read
  from its rows) to its homology dump, in species order (typed under
  result.tables.dumps).

- `pairs`:

  object\[\]; never `NA`. One entry per species pair (a, b) with a \<= b
  in ordinal order, (s, s) included, in OrthologySnapshot.Pairs order:
  how many relationships of each class it holds (typed under
  result.tables.pairs).

- `written`:

  object; `NA` when no `--out` was given; nothing was written. The
  written snapshot: its directory and every file but manifest.json,
  which sits beside them (typed under result.tables.written).

- `pair_status`:

  object; `NA` when no `--pair-status` was given. {from_species,
  to_species, row_count, status_counts, column_names, columns}: one row
  per gene of from_species (typed under result.tables.pair_status and
  result.tables.pair_status_columns).

- `caveats`:

  string\[\]; never `NA`. This call's traps: an input whose file name
  carries no release (not checked against `--release`), dumps sharing
  one file name (source_dump cannot tell them apart), a GTF named for
  another release, identical duplicate rows, unrestricted gene trees
  holding genes outside the gene sets.

## Fields of `gene_sets`

- `source_file_name`:

  string; never `NA`. The GTF's file name; for a gene-set table, still
  the GTF's (EnsemblGeneSet.SourceFileName).

- `sha256`:

  string; never `NA`. Lower-case hex sha256 of the GTF's bytes as read,
  the compressed bytes for a .gz (EnsemblGeneSet.SourceSha256).
  manifest.json's gene_set_sha256 and genes resolve's gene_set_sha256
  carry the same value, so a resolution and a snapshot can be checked
  against one gene set.

- `release`:

  string; `NA` when the GTF's file name carries no release. The
  annotation release from the GTF's name; it may differ from `--release`
  (see the release param).

- `genome_build`:

  string; `NA` when the GTF has no \#!genome-build header. e.g.
  GRCh38.p14.

- `genebuild_last_updated`:

  string; `NA` when the GTF has no \#!genebuild-last-updated header. The
  header's value.

- `gene_count`:

  int; in **genes**; never `NA`. Genes in the set, every biotype.

## Fields of `gene_trees`

- `gene_trees$source_file_name` (wire `source_file_name`):

  string; never `NA`. The file's name.

- `gene_trees$sha256` (wire `sha256`):

  string; never `NA`. Lower-case hex sha256 of the file's bytes as read
  (compressed for a .gz).

- `gene_trees$release` (wire `release`):

  string; `NA` when the file name carries no .eRELEASE.txt suffix, so
  the release was not checked. The release parsed from the file name.

- `gene_trees$collection` (wire `collection`):

  string; `NA` when the file name does not start
  COLLECTION.GeneTree_content.. The Compara collection, e.g.
  vertebrates; manifest.json's collection.

- `gene_trees$member_count` (wire `member_count`):

  int; in **genes**; never `NA`. Genes held from the file
  (ComparaGeneTreeContent.Count): every gene in it, or with
  `--restrict-gene-trees` only those in the gene sets.

- `gene_trees$restricted_to_gene_sets` (wire `restricted_to_gene_sets`):

  string\[\]; `NA` when no `--restrict-gene-trees`: every gene in the
  file was kept. The sha256s of the gene sets the members were
  restricted to, in ordinal order
  (ComparaGeneTreeContent.RestrictedToGeneSets).

## Fields of `dumps`

- `source_file_name`:

  string; never `NA`. The dump's file name. Ensembl gives every species'
  dump the same name, so this alone does not say which dump it is: the
  key does.

- `sha256`:

  string; never `NA`. Lower-case hex sha256 of the file's bytes as read
  (compressed for a .gz).

- `release`:

  string; `NA` when the file name does not start Compara.RELEASE., so
  the release was not checked. The release parsed from the file name.

- `rows_read`:

  int; in **rows**; never `NA`. Every data row read and checked, kept or
  not (ComparaHomologyDump.RowsRead).

- `rows_kept`:

  int; in **rows**; never `NA`. Rows whose two species are both snapshot
  species (ComparaHomologyDump.Rows). The rest name a species not given.

## Fields of `pairs`

- `species_a`:

  string; never `NA`. The species that sorts first; side A of every row
  in pairs/A\_\_B.parquet.

- `species_b`:

  string; never `NA`. The other species; the same as species_a for the
  paralogs of one species.

- `row_count`:

  int; in **rows**; never `NA`. Relationships between the two
  (OrthologySnapshot.Pair), one row each, never merged or picked; the
  rows of pairs/A\_\_B.parquet.

- `ortholog_count`:

  int; in **rows**; never `NA`. Of those, ortholog_one2one,
  ortholog_one2many or ortholog_many2many (relationship_class ortholog).
  Always 0 for (s, s).

- `paralog_count`:

  int; in **rows**; never `NA`. within_species_paralog, other_paralog or
  gene_split. Always 0 between two species.

- `homoeolog_count`:

  int; in **rows**; never `NA`. homoeolog_one2one, homoeolog_one2many or
  homoeolog_many2many: between subgenomes of one polyploid species.

## Fields of `written`

- `written$directory` (wire `directory`):

  string; never `NA`. The `--out` directory, absolute. manifest.json is
  in it.

- `written$files` (wire `files`):

  object\[\]; never `NA`. Every file written but manifest.json, in
  ordinal order of path (OrthologySnapshotWriter's
  SnapshotManifest.Files): {path, rows, bytes, sha256} (typed under
  result.tables.written_files).

## Fields of `written$files`

a data.frame, one row per file written but manifest.json.

- `path`:

  string; never `NA`. Relative to directory, with forward slashes, e.g.
  pairs/homo_sapiens\_\_mus_musculus.parquet.

- `rows`:

  int; in **rows**; `NA` when the file is not a table (views.sql). Rows
  in the Parquet table; 0 for a pair with no relationship, which is
  still written.

- `bytes`:

  int; in **bytes**; never `NA`. The file's size.

- `sha256`:

  string; never `NA`. Lower-case hex sha256 of the file; manifest.json
  lists the same.

## Fields of `pair_status`

- `pair_status$from_species` (wire `from_species`):

  string; never `NA`. FROM: whose genes the rows are.

- `pair_status$to_species` (wire `to_species`):

  string; never `NA`. TO: the species each status is about.

- `pair_status$row_count` (wire `row_count`):

  int; in **genes**; never `NA`. Every gene of from_species' gene set,
  all biotypes.

- `pair_status$status_counts` (wire `status_counts`):

  map\<string,int\>; in **genes**; never `NA`. Status name to genes with
  it, every OrthologyStatus in mzLib's order, zeros included.
  not_in_gene_set is always 0 here: every row is a gene of the gene set.

- `pair_status$column_names` (wire `column_names`):

  string\[\]; never `NA`. gene_id, gene_name, gene_biotype, status.

- `pair_status$records` (wire `columns`):

  table; never `NA`. Column name to one value per gene.

## Fields of `pair_status$records`

a data.frame, one row per gene.

- `gene_id`:

  string; never `NA`. Stable Ensembl gene id, unversioned as Ensembl
  GTFs write it.

- `gene_name`:

  string; `NA` when the GTF gives the gene no gene_name (common for
  novel lncRNAs). The release's gene_name: a display label, never a key.

- `gene_biotype`:

  string; never `NA`. The release's gene_biotype, unknown when the GTF
  omits it. Filter on it to choose the denominator: all biotypes are
  rows, and protein_coding alone is not the same set.

- `status`:

  string; never `NA`. OrthologySnapshot.StatusName: has_ortholog (at
  least one ortholog in to_species), no_edge_in_shared_tree (the gene's
  tree holds a to_species gene and Compara called no ortholog to it: the
  closest the data comes to no ortholog exists),
  tree_lacks_target_species (the tree holds no to_species gene, so
  nothing was compared), or not_in_any_tree (in the gene set but in no
  gene tree: no call was attempted).

## Fields of `parquet_genes`

the columns of genes/SPECIES.parquet in a written snapshot; mzLibR does
not read Parquet.

- `gene_id`:

  string; never `NA`. genes/SPECIES.parquet, one row per gene of the
  gene set.

- `gene_version`:

  int; `NA` when the GTF gives no gene_version (absent is not 0). The
  gene_version attribute.

- `gene_biotype`:

  string; never `NA`. unknown when the GTF omits it.

- `gene_name`:

  string; `NA` when the GTF gives no gene_name. A display label.

- `seq_region`:

  string; never `NA`. The GTF's first column: chromosome or scaffold.

- `species`:

  string; never `NA`. The production name.

## Fields of `parquet_members`

the columns of members/SPECIES.parquet in a written snapshot.

- `gene_id`:

  string; never `NA`. members/SPECIES.parquet, one row per gene of the
  gene set that is in a gene tree.

- `species`:

  string; never `NA`. The production name.

- `group_id`:

  string; never `NA`. compara-RELEASE:TREE_ID, so a tree id is never
  compared across releases.

- `source_group_id`:

  string; never `NA`. Compara's gene-tree stable id, e.g.
  ENSGT00940000153582.

- `canonical_protein_id`:

  string; never `NA`. The protein Compara used for the gene (the row
  flagged Y).

## Fields of `parquet_pairs`

the columns of pairs/A\_\_B.parquet in a written snapshot.

- `homology_id`:

  string; never `NA`. pairs/A\_\_B.parquet, one row per relationship,
  ordered by homology id. Compara's id, unique within a release.

- `relationship_type`:

  string; never `NA`. Compara's homology_type verbatim, e.g.
  ortholog_one2many.

- `relationship_class`:

  string; never `NA`. ortholog, paralog or homoeolog
  (ComparaHomologyDump.ClassName).

- `species_a`:

  string; never `NA`. The species that sorts first.

- `species_b`:

  string; never `NA`. The other species.

- `gene_a`:

  string; never `NA`. Side A's gene; every value moves with its gene
  when a row is reoriented.

- `protein_a`:

  string; never `NA`. Side A's protein in the alignment.

- `identity_a`:

  float; in **percent**; `NA` when Compara gave NULL. Percentage of A's
  protein identical to B's in the alignment, 0 to 100.

- `gene_b`:

  string; never `NA`. Side B's gene.

- `protein_b`:

  string; never `NA`. Side B's protein.

- `identity_b`:

  float; in **percent**; `NA` when Compara gave NULL. Percentage of B's
  protein identical to A's, 0 to 100. Not identity_a: each is over its
  own protein's length.

- `dn`:

  float; in **nonsynonymous substitutions per nonsynonymous site**; `NA`
  when Compara gave NULL (it computes dN for few comparisons). Compara's
  dn verbatim.

- `ds`:

  float; in **synonymous substitutions per synonymous site**; `NA` when
  Compara gave NULL. Compara's ds verbatim.

- `goc_score`:

  int; in **score, 0 to 100**; `NA` when Compara gave NULL (not computed
  for this pair). Gene order conservation score.

- `wga_coverage`:

  float; in **percent**; `NA` when Compara gave NULL (no whole-genome
  alignment for this pair). Whole-genome alignment coverage, 0 to 100.

- `is_high_confidence`:

  bool; `NA` when Compara gave NULL: it declines to say, which is not
  false. Compara's high-confidence flag.

- `source_dump`:

  string; never `NA`. The FILE NAME of the dump the row was read from.
  Which of a pair's two dumps holds a row is arbitrary, so this is
  provenance, not meaning; and as Ensembl names every dump the same, it
  cannot say which species' dump.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  `--release` or `--gene-trees` missing or given without a value;
  `--out` or `--pair-status` given without a value; a value after
  `--restrict-gene-trees`

- `mzlib_usage_error` (usage):

  no species lines on stdin, no line holding only –, or no dump paths
  after it; a species line that is not species, tab, gtf or gene-set,
  tab, path; a species listed twice

- `mzlib_usage_error` (usage):

  a GTF, gene-set table, dump or gene-trees file that does not exist;
  `--out` naming a file or a directory that is not empty
  (OrthologySnapshotWriter.WriteAsync's refusal, applied before any
  input is read)

- `mzlib_usage_error` (usage):

  `--pair-status` not two comma-separated names, naming a species not
  given on stdin, or one species twice (the rule of
  OrthologySnapshot.PairStatus, applied before any input is read)

- `mzlib_usage_error` (usage):

  mzLib's OrthologySnapshot.Build refused the inputs as not fitting
  together (ArgumentException, with mzLib's message): the gene trees or
  a dump name a release other than `--release`; a dump is of a species
  not given, two dumps are of one species, or a species has no dump

- `mzlib_bridge_error` (correctness):

  a file mzLib's reader refuses (InvalidDataException, naming file and
  line): a GTF gene row without gene_id; a table that is not a gene-set
  table; a dump with other columns, a short row, a value that is not a
  number, an unknown homology_type, a paralog between two species, an
  ortholog within one, or a second genome; a gene-tree row without four
  cells or with a flag other than Y or N, or a kept gene in two trees or
  without exactly one canonical protein

- `mzlib_bridge_error` (correctness):

  OrthologySnapshot.Build found data that breaks the snapshot's rules
  (InvalidDataException): one gene in two species' gene sets, a homology
  naming a gene its species' gene set lacks (the gene set and the dumps
  are from different releases), or one homology id with two different
  rows

- `mzlib_bridge_error` (correctness):

  the snapshot could not be written (IOException from the disk)

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Every species' homology dump is required, and mzLib refuses a call
  without one. Compara's README warns that each genome's dump holds an
  arbitrary subset of its relationships, and in practice they partition:
  in release 116 every human-mouse orthology is in the mouse dump and
  none is in the human one, so the dump named after one species gives a
  complete-looking, badly incomplete answer.

- Orthology is asserted per species pair and is not transitive: A~B and
  B~C do not make A~C. Never chain pairs into groups. The written
  views.sql's species_set3 macro, and mzLib's
  OrthologySnapshot.SpeciesSet, form a tuple only when every two members
  are orthologs.

- Nothing is picked, merged or collapsed: a one-to-many or many-to-many
  orthology is several rows, and a gene's paralogs are the pair (s, s).

- A status says how far Compara's inference got, not what biology holds.
  no_edge_in_shared_tree is the nearest the data comes to no ortholog
  exists; tree_lacks_target_species and not_in_any_tree mean nothing was
  compared.

- Every biotype is in pair_status and in genes/SPECIES.parquet. Which
  genes form a denominator is the caller's choice: filter on
  gene_biotype. mzLib's own human-to-mouse check counted protein-coding
  plus IG and TR genes, which is not the same set as protein_coding
  alone.

- Values are Compara's, verbatim: a NULL in a dump stays null in the
  Parquet (identity, dn, ds, goc_score, wga_coverage,
  is_high_confidence), never 0 or false. identity_a and identity_b
  differ because each is over its own protein.

- The gene set and the dumps must be from one release: mzLib refuses a
  homology whose gene is not in its species' gene set. Use the
  primary-assembly GTF (Species.Assembly.Release.gtf.gz), as for genes
  resolve.

- A snapshot is written once. The same inputs give the same bytes and
  the same snapshot_id, so a later release is a new directory compared
  with the old one; a directory that is not empty is refused, never
  overwritten.

- mzLib ships a writer and no reader. Query the Parquet with DuckDB
  (load views.sql for the orthologs, pair_status and species_set3
  macros), pyarrow, polars or R's arrow.

- Ensembl names every species' dump
  Compara.RELEASE.protein_default.homologies.tsv.gz in its own folder,
  and the pairs tables' source_dump column records only the file name,
  so it cannot tell two species' dumps apart. manifest.json's inputs
  list each dump with its species and sha256.

## Performance

Everything is read into memory, and the inputs are large: human release
116's GTF is 141 MB compressed (mzLib's EnsemblGeneSet remarks), and the
vertebrate gene-tree file holds 4.2 million genes unless
`--restrict-gene-trees` keeps only the given species'. For human, mouse
and rat at release 116, mzLib \#1381 reports the build from the three
dumps at about 20 s and the written snapshot (13 files, 9.0 MB) at about
40 s more. No call here takes a degree of parallelism, so there is no
`--threads`. Every call starts one bridge process, which costs a .NET
start-up before any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.build_orthology`

- Rust (mzLibRust): `mzlib::proteins::build_orthology_with` with
  `OrthologyOptions`

- R (mzLibR): `proteins_build_orthology`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- result.tables (gene_sets, gene_trees, dumps, pairs, written,
  pair_status and the Parquet tables) is not validated by
  tools/check_verbs.py, which checks only envelope_fields, columns and
  bulk. Bridge tests pin each nested object's keys
  (GenesOrthologyTests).

- The Parquet tables' units are mzLib's doc comments (ComparaHomology:
  identity percent, goc_score and wga_coverage 0 to 100). dn and ds are
  written as dN and dS conventionally are, substitutions per site;
  neither mzLib nor this bridge has checked Compara's own definition.

- Only one pair-status direction per call, and no SpeciesSet tuples or
  ortholog rows on the wire: the written Parquet and views.sql answer
  those. Whether a binding needs them without a Parquet reader (pyMzLib
  has no third-party dependency) is open.

- source_dump records the dump's file name, which Ensembl gives every
  species' dump (PUSH-TO-MAINLAND: mzLib could record the dump's genome,
  or the path under homologies/, beside it).

- Not yet ported to mzLibRust or mzLibR; the binding names above are
  proposals.

## References

- [doi:10.1101/gr.073585.107](https://doi.org/10.1101/gr.073585.107):
  Ensembl Compara gene trees, the duplication-aware phylogenies its
  orthology and paralogy calls are read from (Vilella et al. 2009,
  Genome Res 19:327)

- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl releases, their GTFs and the Compara homology dumps (Ensembl
  2026, NAR)

## See also

[`proteins_orthology_pair`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_orthology_pair.md),
[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)

## Examples

``` r

# mzLib's own test world: four human genes, three mouse and two rat, with every human-mouse
# row in the mouse dump, as in release 116.
gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
          mus_musculus = "Mus_musculus.GRCm39.116.gtf",
          rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
dumps <- paste0("homologies/", names(gtfs), "/Compara.116.protein_default.homologies.tsv.gz")
snap <- proteins_build_orthology("116", gtfs = gtfs, dumps = dumps,
  gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
  pair_status = c("homo_sapiens", "mus_musculus"))
snap
#> <mzlibr_orthology_snapshot> Ensembl 116, 3 species: homo_sapiens, mus_musculus, rattus_norvegicus
#>   snapshot_id 5487ebcabe48... (record it with any result)
#>   orthologs: homo_sapiens~mus_musculus 3, homo_sapiens~rattus_norvegicus 1, mus_musculus~rattus_norvegicus 2
#>   homo_sapiens in mus_musculus: has_ortholog 2, tree_lacks_target_species 1, not_in_any_tree 1
#>   ! The dumps of homo_sapiens, mus_musculus, rattus_norvegicus share the file name Compara.116.protein_default.homologies.tsv.gz, and the pairs tables' source_dump column holds only the file name, so it cannot say which species' dump a row came from. manifest.json's inputs list each dump with its species and sha256.
#>   ! The gene trees were read unrestricted: 9 genes were held in memory, 8 of them in the given gene sets. --restrict-gene-trees holds only those and gives the same snapshot.
snap$pair_status$records[, c("gene_id", "gene_biotype", "status")]
#>           gene_id   gene_biotype                    status
#> 1 ENSG00000000001 protein_coding              has_ortholog
#> 2 ENSG00000000002 protein_coding              has_ortholog
#> 3 ENSG00000000003 protein_coding tree_lacks_target_species
#> 4 ENSG00000000004 protein_coding           not_in_any_tree
snap$dumps[, c("species", "rows_kept")]
#>             species rows_kept
#> 1      homo_sapiens         1
#> 2      mus_musculus         5
#> 3 rattus_norvegicus         1
proteins_orthology_pair(snap, "mus_musculus", "homo_sapiens")$row_count
#> [1] 3

# With out = the snapshot is written, and its id is the same:
written <- proteins_build_orthology("116", gtfs = gtfs, dumps = dumps,
  gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
  restrict_gene_trees = TRUE, out = "compara-116")
identical(written$snapshot_id, snap$snapshot_id)
#> [1] TRUE
head(written$written$files[, c("path", "rows")], 7)
#>                                       path rows
#> 1               genes/homo_sapiens.parquet    4
#> 2               genes/mus_musculus.parquet    3
#> 3          genes/rattus_norvegicus.parquet    2
#> 4             members/homo_sapiens.parquet    3
#> 5             members/mus_musculus.parquet    3
#> 6        members/rattus_norvegicus.parquet    2
#> 7 pairs/homo_sapiens__homo_sapiens.parquet    0
```
