# Classify peptides by how widely they are shared across protein databases, with I = L

Classify peptides as unique to one sequence, shared within a gene,
shared across genes, or not in the databases, treating I and L as the
same residue.

Wraps mzLib's `PeptideUniquenessClassifier.Classify`. Each peptide is
one of:

## Usage

``` r
proteins_classify_peptides(peptides, databases, contaminants = NULL, threads = 1,
  timeout = NULL)
```

## Arguments

- peptides:

  A character vector of unmodified base sequences, upper case
  (`"PEPTIDEK"`), one result row each in the same order, duplicates
  included. Sent on stdin. Strip modifications first: mzLib refuses a
  sequence with a bracketed modification, or in lower case
  (`"peptidek"`), rather than guess.

- databases:

  The search database or databases: UniProt XML or FASTA, optionally
  gzipped.

- contaminants:

  Contaminant databases, searched alongside: they are real sequences in
  the search space, so a peptide a contaminant shares with a target is
  shared.

- threads:

  Databases read at once (the classification itself is one pass).
  Default 1; `-1` means every core. Same answer at any value.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait.

## Details

`"Unique"` - every protein containing it has the **same sequence**
(identical entries under two accessions count as one; both are listed).
`"SharedWithinGene"` - several distinct sequences, all sharing a gene:
isoforms, so the peptide supports the gene, not an isoform.
`"SharedAcrossGenes"` - sequences with no gene in common.
`"NotInDatabase"` - no target protein contains it.

A protein contains a peptide when its sequence contains it **anywhere,
whatever the protease** - deliberately conservative, so a peptide called
`"Unique"` cannot be explained by another entry at a site the search's
cleavage rules happened to skip. **I and L are the same residue**
throughout, because a mass spectrometer cannot tell them apart.

## Value

An `mzlibr_peptide_classification`. `file_count` counts the database
files searched; `peptide_count` counts the peptides (rows);
`target_protein_count` counts the proteins searched, contaminants
included; `decoy_proteins_ignored` counts the decoy proteins, never
searched. `sharing_counts` is a named vector: class to number of
peptides, every class present. `i_and_l_equivalent` is always `TRUE`,
carried so the rule travels with the result.

`peptides` is a data.frame with one row per peptide, in the order given:
`peptide` (exactly as given, not I/L-folded), `sharing`,
`accession_count` (the target proteins containing it), and the list
columns `accessions` and `shared_gene_keys` (the gene keys common to all
of them, such as `"ensembl:ENSG00000111640"`; empty for
`"NotInDatabase"` and `"SharedAcrossGenes"`).

## Why there is no on_error

Every database is part of one search space. Skipping one that failed to
load would report the peptides it contains as unique, so any failure
raises.

## Wraps

Wire verb `proteins classify-peptides`. Generated from the bridge's verb
spec `proteins.classify-peptides.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PeptideUniquenessClassifier.Classify`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/PeptideUniqueness.cs)
  in `mzLib/Proteomics/Protein/PeptideUniqueness.cs` at mzLib `23c2490e`

- [`PeptideUniquenessClassifier.DefaultGeneKeys`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/PeptideUniqueness.cs)
  in `mzLib/Proteomics/Protein/PeptideUniqueness.cs` at mzLib `23c2490e`

- [`ProteinDbLoader.LoadProteinXML`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

- [`ProteinDbLoader.LoadProteinFasta`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

## Parameters: units, ranges and defaults

- `databases` (wire `--path`):

  path; default absent. One database. Exactly one of path and
  paths-stdin. With path, all of stdin is peptides.

- `contaminants` (wire `--contaminant`):

  flag; default `FALSE`. Load the `--path` database as contaminants
  (they are still searched).

- wire `--paths-stdin`:

  flag; default `FALSE`. stdin holds database paths (each optionally tab
  'contaminant'), then a line holding only '–', then the peptides. *Not
  an argument here: chosen automatically when more than one database is
  given.*

- `threads`:

  int; default `1`; range `>= 1, or -1 for every core`. Databases loaded
  at once. mzLib's classifier has no thread knob of its own; it is one
  pass. Output identical at any value.

- wire `--on-error`:

  string; default `fail`; range `fail`. Only 'fail' is accepted: every
  database is one search space, and skipping one would report the
  peptides it holds as unique. *Not an argument here: the only accepted
  value is fail, the default; there is no skip to choose.*

- `peptides`:

  string\[\]; required; range `non-empty, letters A-Z only`. Read from
  stdin (not an option). Unmodified base sequences, upper case, one per
  line; blank lines ignored. One result per peptide, in order,
  duplicates included.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **files**; never `NA`. Databases searched.

- `peptide_count`:

  int; in **peptides**; never `NA`. Rows (input peptides).

- `target_protein_count`:

  int; in **proteins**; never `NA`. Non-decoy proteins searched,
  contaminants included.

- `decoy_proteins_ignored`:

  int; in **proteins**; never `NA`. Decoy entries present in the files
  and never searched.

- `i_and_l_equivalent`:

  bool; never `NA`. Always true: I and L are one residue for matching.
  On the wire so the rule travels with the result.

- `sharing_counts`:

  map\<string,int\>; in **peptides**; never `NA`. Class -\> peptides;
  all four keys always present.

- `files`:

  object\[\]; never `NA`. As proteins read; record_count is the target
  proteins searched from that file.

- `column_names`:

  string\[\]; never `NA`. Column order.

- `peptides` (wire `columns`):

  table; never `NA`. Column name to per-peptide values.

## Columns of `peptides`

- `peptide`:

  string; never `NA`. Exactly as given (not I/L-folded).

- `sharing`:

  string; never `NA`. mzLib PeptideSharing name: NotInDatabase \| Unique
  (every containing protein has one sequence) \| SharedWithinGene
  (several sequences sharing a gene key) \| SharedAcrossGenes.

- `accession_count`:

  int; in **proteins**; never `NA`. Target proteins containing the
  peptide.

- `accessions`:

  string\[\]; never `NA`. Those proteins, distinct, ordinal order; empty
  for NotInDatabase.

- `shared_gene_keys`:

  string\[\]; never `NA`. Gene keys common to all of them
  (ensembl:ENSG..., gene:ORGANISM:PRIMARY_NAME, entry:UNIPROT_ENTRY),
  ordinal order. Empty for NotInDatabase and SharedAcrossGenes.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  no peptides; a peptide that is not non-empty A-Z (mzLib's
  ArgumentException, reclassified; the message names the peptide) -
  lower case and bracketed modifications included

- `mzlib_usage_error` (usage):

  on-error other than fail; any database-selection error proteins read
  raises; paths-stdin without the '–' line

- `mzlib_bridge_error` (correctness):

  a database mzLib cannot parse

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Containment, not digestion: a protein contains a peptide if its
  sequence contains it anywhere, whatever the protease. Deliberately
  conservative.

- Identical sequences under different accessions are one sequence: the
  peptide is Unique and lists both accessions.

- Contaminants are searched; leaving the contaminant database out can
  turn a SharedAcrossGenes peptide (YLYEIAR, human and bovine albumin)
  into Unique.

- A FASTA protein has no ensembl: key, so gene sharing for it rests on
  organism + primary gene name and the UniProt entry.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.classify_peptides`

- Rust (mzLibRust): `mzlib::proteins::classify_peptides_with` with
  `ClassifyOptions`

- R (mzLibR): `proteins_classify_peptides`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Rust and R spellings are intended, not yet implemented.

- The gene-key function is mzLib's DefaultGeneKeys and is not selectable
  on the wire.

- The peptides param travels on stdin, not as an option; the spec format
  has no stdin field yet, so its doc says so.

## References

- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, the entries and gene names the gene keys come from (UniProt
  Consortium 2025, NAR)

- [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239):
  Ensembl gene ids, the first gene key (Ensembl 2026, NAR)

## See also

[`proteins_read`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md),
[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)

## Examples

``` r

peptides <- c("VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK")
calls <- proteins_classify_peptides(peptides, c("human_subset.xml", "human_extra.fasta"),
  contaminants = "contaminants.fasta")
calls
#> <mzlibr_peptide_classification> 6 peptides against 7 proteins in 3 database(s)
#>   NotInDatabase 1, Unique 3, SharedWithinGene 1, SharedAcrossGenes 1
#>   I and L are one residue for matching
calls$peptides[, c("peptide", "sharing", "accession_count")]
#>           peptide           sharing accession_count
#> 1        VGVNGFGR            Unique               1
#> 2  LVLNGNPLTLFQER            Unique               1
#> 3 ALSEQINIFFDYSGR  SharedWithinGene               3
#> 4         YLYEIAR SharedAcrossGenes               2
#> 5        AEFVEVTK            Unique               1
#> 6        PEPTIDEK     NotInDatabase               0
# YLYEIAR is in human and bovine albumin, which share no gene:
calls$peptides$accessions[[4]]
#> [1] "P02768" "P02769"
```
