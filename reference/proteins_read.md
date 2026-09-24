# Read protein databases: one row per protein, with GO terms and Ensembl genes on request

Read protein databases (UniProt XML or FASTA) into one row per protein,
with GO terms and Ensembl gene links as long tables on request.

Loads each database with mzLib's `ProteinDbLoader`, exactly as a search
would but with no decoys generated, and answers "what is this
accession?": organism, NCBI taxon, gene names, length and mass, plus -
when you ask for them - its Gene Ontology terms and Ensembl gene links
as long tables.

## Usage

``` r
proteins_read(databases, tables = "proteins", accessions = NULL, contaminants = NULL,
  sequences = FALSE, threads = 1, on_error = "fail", timeout = NULL)
```

## Arguments

- databases:

  A UniProt XML (`.xml`) or FASTA (`.fasta`, `.fa`, `.faa`, `.fas`),
  each optionally gzipped, or a character vector of them. Read in the
  order given, in one bridge call.

- tables:

  Which tables to return, from `"proteins"`, `"go_terms"` and
  `"ensembl_genes"`. Default `"proteins"`. Ask for the other two when
  you want them: a whole human proteome has about twenty GO rows per
  protein.

- accessions:

  Keep only proteins whose accession is one of these - an **exact**
  match, so `"P04406-1"` does not find `"P04406"`. Applies to every
  table; misses are listed in `accessions_not_found`. Sent on stdin, so
  tens of thousands are fine. `NULL` keeps every protein.

- contaminants:

  Databases to load as contaminants; `is_contaminant` is `TRUE` on their
  rows. They come after `databases` in `files`.

- sequences:

  Add a `sequence` column to the `proteins` table.

- threads:

  Databases read at once. Default 1; `-1` means every core. A resource
  choice only: the result is identical at any value.

- on_error:

  `"fail"` (the default) raises on the first unreadable database, in
  input order; `"skip"` records the failure in that database's row of
  `files` and reads the rest.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait: a proteome XML
  takes a while.

## Details

**A FASTA knows no GO terms and no Ensembl genes.** Its headers carry
organism (`OS=`), taxon (`OX=`) and gene name (`GN=`) and nothing else,
so a FASTA contributes no `go_terms` or `ensembl_genes` rows. That is
the format's silence, not a biological absence, and `files` says so in
its `absent_fields` column. Read the UniProt XML of the same proteome
for those.

## Value

An `mzlibr_protein_database`. `file_count`, `read_count` and
`failed_count` count database files; `record_count` counts the proteins
that passed the accession filter; `accession_filter_count` counts the
distinct accessions in the filter, and is `NA` with no filter;
`accessions_not_found` is `NULL` with no filter. `tables` names what was
returned, `caveats` the traps in the result as a whole, and `files` has
one row per input (its `source_index`, `contaminant`, `protein_count`,
`decoy_count`, `record_count`, `caveats`, `absent_fields` and any
`error_*`).

`proteins` is a data.frame with one row per protein, or `NULL` when it
was not requested: `length` in residues and `monoisotopic_mass` in Da
(of the unmodified sequence as written, plus one water; `NA` when the
sequence holds a letter with no defined mass). `gene_names` and
`ensembl_gene_ids` are list columns. `go_terms` and `ensembl_genes` are
data.frames, or `NULL` when not requested. \*\*`source_index` is
1-based\*\* in all of them, so `db$files[db$proteins$source_index, ]`
finds each row's database.

## Masses are of the precursor as written

`monoisotopic_mass` is the unmodified sequence plus one water - the
initiator methionine, signal peptide and propeptide included - so it is
not the mass of the mature protein.

## No decoys, no variants

A database is read as written. A decoy already in the file (an accession
starting `DECOY`) is kept and flagged `is_decoy`; none are generated. A
UniProt XML that records a **genotype** still has those variants applied
by mzLib, which renames the accession (`P38936_C117Y`); the file's
`caveats` say how many.

## Wraps

Wire verb `proteins read`. Generated from the bridge's verb spec
`proteins.read.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`ProteinDbLoader.LoadProteinXML`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

- [`ProteinDbLoader.LoadProteinFasta`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

- [`Protein.NcbiTaxonomyId`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/Protein.cs)
  in `mzLib/Proteomics/Protein/Protein.cs` at mzLib `23c2490e`

- [`Protein.GoTerms`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/GoTerm.cs)
  in `mzLib/Proteomics/Protein/GoTerm.cs` at mzLib `23c2490e`

- [`Protein.EnsemblGeneReferences`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/EnsemblGeneReference.cs)
  in `mzLib/Proteomics/Protein/EnsemblGeneReference.cs` at mzLib
  `23c2490e`

- [`Peptide.MonoisotopicMass`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/AminoAcidPolymer/AminoAcidPolymer.cs)
  in `mzLib/Proteomics/AminoAcidPolymer/AminoAcidPolymer.cs` at mzLib
  `23c2490e`

## Parameters: units, ranges and defaults

- `databases` (wire `--path`):

  path; default absent. One database: UniProt XML (.xml) or FASTA
  (.fasta, .fa, .faa, .fas), each optionally .gz. Exactly one of path
  and paths-stdin.

- `contaminants` (wire `--contaminant`):

  flag; default `FALSE`. Load the `--path` database as contaminants
  (mzLib isContaminant). Usage error with paths-stdin.

- wire `--paths-stdin`:

  flag; default `FALSE`. Databases on stdin, one per line, each
  optionally followed by a tab and 'contaminant' or 'target'. Blank
  lines ignored; a repeated path is a usage error (BULK.md §1). *Not an
  argument here: chosen automatically when more than one database is
  given.*

- `accessions` (wire `--accessions-stdin`):

  flag; default `FALSE`. Keep only proteins whose accession exactly
  equals one read from stdin (one per line). With paths-stdin, stdin
  holds the paths, then a line holding only '–', then the accessions.
  Applies to every table.

- `tables`:

  string\[\]; default `proteins`; range
  `subset of proteins, go_terms, ensembl_genes`. Comma list of the
  tables to return. A table not asked for is null.

- `sequences`:

  flag; default `FALSE`. Add a sequence column to the proteins table.

- `threads`:

  int; default `1`; range `>= 1, or -1 for every core`. Databases loaded
  at once. Output is byte-identical at any value (tested at 1 and 4).
  mzLib's own maxThreads is pinned to 1 per load, so the degree is owned
  here alone.

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. fail raises the failure
  of the LOWEST failing input index; skip records it in files\[i\].error
  and reads the rest.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **files**; never `NA`. Databases given.

- `read_count`:

  int; in **files**; never `NA`. Databases read.

- `failed_count`:

  int; in **files**; never `NA`. Databases that failed (on-error skip
  only).

- `record_count`:

  int; in **proteins**; never `NA`. Proteins that passed the accession
  filter (rows of the proteins table when requested).

- `tables`:

  string\[\]; never `NA`. The tables returned, in the order proteins,
  go_terms, ensembl_genes.

- `accession_filter_count`:

  int; in **accessions**; `NA` when no accession filter was given.
  Distinct accessions in the filter.

- `accessions_not_found`:

  string\[\]; `NA` when no accession filter was given. Filter accessions
  no database contained, in the order given. Matching is exact: P04406-1
  does not find P04406.

- `caveats`:

  string\[\]; never `NA`. Whole-result traps, e.g. proteins with no
  computable mass.

- `files`:

  object\[\]; never `NA`. One per input in input order (BULK.md §3):
  source_index, path (absolute), file_type (UniProtXml \| Fasta \| null
  if it failed first), reader, contaminant, protein_count, decoy_count,
  record_count, search_database_sha256 (null here), caveats,
  absent_fields, error (null or {kind: usage\|correctness, type,
  message}). Present with `--path` too, so one shape serves both.

- `column_names`:

  string\[\]; `NA` when the proteins table was not requested. Column
  order of the proteins table.

- `proteins` (wire `columns`):

  table; `NA` when the proteins table was not requested. Column name to
  per-protein values.

- `go_terms`:

  object; `NA` when go_terms was not requested. {row_count,
  column_names, columns}: one row per (protein, GO id), with columns
  source_index, source_path, accession, go_id, aspect, term_name,
  evidence_codes, projects (typed under result.tables.go_terms in this
  spec).

- `ensembl_genes`:

  object; `NA` when ensembl_genes was not requested. {row_count,
  column_names, columns}: one row per (protein, Ensembl transcript),
  with columns source_index, source_path, accession, transcript_id,
  protein_id, gene_id, versioned_gene_id, gene_version (typed under
  result.tables.ensembl_genes in this spec).

## Columns of `proteins`

- `source_index`:

  int; never `NA`. Position of the input in files (paths first, in
  order).

- `source_path`:

  string; never `NA`. Absolute path of that input.

- `accession`:

  string; never `NA`. As the database wrote it. A FASTA with a repeated
  accession gets \_2, \_3 appended by mzLib; an XML entry with an
  applied genotype variant gets a variant suffix (see caveats).

- `name`:

  string; `NA` when the database records no entry name. Entry name, e.g.
  G3P_HUMAN.

- `full_name`:

  string; `NA` when the database records no protein name. Recommended
  protein name.

- `organism`:

  string; `NA` when the database records no organism. Scientific name
  (the XML organism element, FASTA OS=).

- `ncbi_taxonomy_id`:

  string; `NA` when the entry records no taxonomy id (a FASTA header
  without OX=). NCBI taxon, e.g. "9606". A string because it is an
  identifier.

- `primary_gene_name`:

  string; `NA` when the entry names no primary gene. Primary gene name
  (the XML gene name of type primary, FASTA GN=).

- `gene_names`:

  string\[\]; never `NA`. Every gene name the entry gives, in file order
  (primary, synonyms, ORF, locus). May be empty; may repeat when the
  file repeats.

- `length`:

  int; in **residues**; never `NA`. Sequence length.

- `monoisotopic_mass`:

  float; in **Da**; `NA` when the sequence holds a letter with no
  defined residue mass (X, B, Z, J). Monoisotopic mass of the unmodified
  sequence as written plus one water: the precursor, initiator Met and
  signal peptide included.

- `is_contaminant`:

  bool; never `NA`. Loaded from a database marked contaminant.

- `is_decoy`:

  bool; never `NA`. The file marks it a decoy (accession prefix DECOY).
  None are generated.

- `is_entrapment`:

  bool; never `NA`. mzLib's entrapment flag (FASTA accession containing
  'Random', case-insensitive).

- `ensembl_gene_ids`:

  string\[\]; never `NA`. Distinct stable Ensembl gene ids
  (Protein.EnsemblGeneIds), ordinal order. Empty for a FASTA: see
  absent_fields.

- `sequence`:

  string; never `NA`. The base sequence. *Present only with
  `sequences`.*

## Fields of `go_terms`

- `source_index`:

  int; never `NA`. Input position.

- `source_path`:

  string; never `NA`. Input path.

- `accession`:

  string; never `NA`. The protein.

- `go_id`:

  string; never `NA`. GO accession, e.g. GO:0005737. Unique per protein.

- `aspect`:

  string; never `NA`. mzLib GoAspect name: BiologicalProcess,
  CellularComponent, MolecularFunction, or Unknown when UniProt gave no
  C:/F:/P: prefix (never inferred).

- `term_name`:

  string; `NA` when the reference carried no term property. Term name
  with the aspect prefix removed.

- `evidence_codes`:

  string\[\]; never `NA`. ECO ids unioned over UniProt's repeats of this
  GO id, sorted ordinal. Empty is normal.

- `projects`:

  string\[\]; never `NA`. Annotating projects (UniProtKB, MGI, HPA,
  ...), unioned, sorted.

## Fields of `ensembl_genes`

- `source_index`:

  int; never `NA`. Input position.

- `source_path`:

  string; never `NA`. Input path.

- `accession`:

  string; never `NA`. The protein.

- `transcript_id`:

  string; `NA` when the reference named no transcript. Ensembl
  transcript as UniProt wrote it, versioned.

- `protein_id`:

  string; `NA` when UniProt gave no Ensembl protein id. Ensembl protein
  id as written.

- `gene_id`:

  string; never `NA`. Stable gene id, e.g. ENSG00000111640: the one to
  join on.

- `versioned_gene_id`:

  string; never `NA`. As UniProt wrote it, e.g. ENSG00000111640.15;
  equal to gene_id when unversioned.

- `gene_version`:

  int; `NA` when the id carried no numeric version (absent is not 0).
  The gene version.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  neither or both of path and paths-stdin; contaminant with paths-stdin;
  no paths on stdin; a blank or repeated path; a role other than
  target/contaminant; accessions-stdin with no accessions; paths-stdin
  with accessions-stdin but no '–' line

- `mzlib_usage_error` (usage):

  tables names something else; threads 0 or \< -1; on-error not
  fail/skip; an option given without a value

- `mzlib_usage_error` (usage):

  a database not named .xml or a FASTA extension; a missing database
  (under skip: recorded in files\[i\].error with kind usage)

- `mzlib_bridge_error` (correctness):

  mzLib cannot parse a database (malformed XML); under skip, recorded
  with kind correctness

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- A FASTA carries no GO or Ensembl cross-references: its
  files\[i\].absent_fields is \[go_terms, ensembl_genes,
  ensembl_gene_ids\] and it contributes no rows to either table. Empty
  there means the format is silent, not that no annotation exists.

- Decoys are not generated (DecoyType.None), but decoys already in the
  file are kept with is_decoy true.

- Sequence variants are not expanded (maxHeterozygousVariants 0), except
  that mzLib still applies variants a VCF-annotated XML records a
  genotype for; such proteins carry a variant accession suffix, and the
  file's caveats count them.

- GO terms are the terms as annotated, not propagated up the GO graph.

- Parallel XML loads share ProteinDbLoader's static modification state;
  it affects only modification resolution, which this verb does not
  report, and output was byte-identical at threads 1 and 4.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.read`

- Rust (mzLibRust): `mzlib::proteins::read_with` with
  `ProteinReadOptions`

- R (mzLibR): `proteins_read`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Rust (mzlib::proteins::read_with + ProteinReadOptions, Table enum) and
  R (proteins_read) spellings are intended, not yet implemented; the
  checker cannot verify them.

- files\[\] is always present, with `--path` too, instead of BULK.md
  §3's flattened single-file block: a new verb has no pre-existing
  single-file shape to preserve, and one shape is simpler for every
  binding. Confirm or amend BULK.md.

- Two lists sharing stdin, separated by a line holding only '–', is new
  wire grammar (also used by proteins classify-peptides). It belongs in
  CONTRACT.md or BULK.md if adopted.

- result.tables (go_terms, ensembl_genes) is not validated by
  tools/check_verbs.py, which checks only 'columns'; the pyMzLib tests
  pin them.

- No `--out` TSV mode yet (BULK.md §2 mentions one): multi-table output
  is JSON only, since three tables would need three files.

- proteins_read_filtered.json (an accession-filtered read) is a second
  recording used by the pyMzLib tests; it is not listed as an example
  because the doctest replay bridge cannot tell two recordings of one
  verb apart by stdin.

## References

- [doi:10.1038/75556](https://doi.org/10.1038/75556): The Gene Ontology
  and its three aspects (Ashburner et al. 2000, Nature Genetics)

- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, the source of the XML's organism, taxonomy, GO and Ensembl
  cross-references (UniProt Consortium 2025, NAR)

- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl, whose transcript and gene ids the ensembl_genes table carries
  (Ensembl 2026, NAR)

## See also

[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md),
[`proteins_classify_peptides`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md)

## Examples

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
db$proteins[, c("accession", "organism", "ncbi_taxonomy_id", "length", "monoisotopic_mass")]
#>   accession     organism ncbi_taxonomy_id length monoisotopic_mass
#> 1    P04406 Homo sapiens             9606    335          36030.40
#> 2    O43653 Homo sapiens             9606    114          11950.87
#> 3    P02768 Homo sapiens             9606    609          69321.50
#> 4    Q13409 Homo sapiens             9606    638          71412.11
#> 5  Q13409-2 Homo sapiens             9606    632          70600.75
#> 6  Q13409-3 Homo sapiens             9606    612          68383.72
#> 7    Q9Z0X1 Mus musculus            10090    612          66723.84
#> 8    P02769   Bos taurus             9913    607          69248.44
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
# A FASTA says what it cannot know, rather than leaving the tables silently empty:
db$files[, c("path", "file_type", "contaminant")]
#>                                   path  file_type contaminant
#> 1   fixtures/proteins/human_subset.xml UniProtXml       FALSE
#> 2  fixtures/proteins/human_extra.fasta      Fasta       FALSE
#> 3  fixtures/proteins/mouse_aifm1.fasta      Fasta       FALSE
#> 4 fixtures/proteins/contaminants.fasta      Fasta        TRUE
db$files$absent_fields[[2]]
#> [1] "go_terms"         "ensembl_genes"    "ensembl_gene_ids"
```
