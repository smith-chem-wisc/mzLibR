# Annotate MetaMorpheus protein groups with Gene Ontology terms, keeping every member

Annotate a stored MetaMorpheus protein-group table with Gene Ontology
terms: one row per (group, term) that any member holds, directly or
through an ancestor, with no member privileged.

Wraps mzLib's `GoGroupAnnotator` over a stored protein-group table
(`ProteinGroupFromTsv.ToGoAnnotationGroups`). For each group it returns
one row per GO term that **any** member holds - directly, or by
propagation up `is_a` and `part_of` - and each row names the members
that carry it. Nothing is collapsed: MetaMorpheus picks no leading
protein, so the table is the union, and the views people usually want
are filters on its rows.

## Usage

``` r
proteins_annotate_go(groups, database, go_obo, category_map = NULL,
  skip_unknown_go_ids = FALSE, out = NULL, categories_out = NULL, limit = NULL,
  offset = 0, timeout = NULL)
```

## Arguments

- groups:

  A MetaMorpheus protein-group table: `AllQuantifiedProteinGroups.tsv`,
  `AllProteinGroups.tsv` (a search without quantification), or one
  file's `<file>_ProteinGroups.tsv`.

- database:

  The UniProt XML (`.xml` or `.xml.gz`) whose GO terms annotate the
  members: the proteome the search used. A FASTA is refused, because it
  carries no GO and every group would read `"no_go_terms"` whatever the
  proteins are.

- go_obo:

  A go.obo file. It must exist; nothing is fetched.

- category_map:

  Optional: your own term-to-category map in mzLib's format
  (`#!category_map_format 1`, `#!map_name`, `#!map_version`, then
  `category`, `subcategory` and `anchor_go_id` columns). mzLib ships no
  vocabulary. Adds `categories`.

- skip_unknown_go_ids:

  `FALSE` (the default) fails when the database cites a GO id the
  ontology release lacks - usually a UniProt release newer than the
  go.obo - and names every missing id. `TRUE` drops each such id and
  lists it in `unresolved_go_ids`.

- out:

  Write the **whole** table here with mzLib's own `GoAnnotationTsv`
  writer, provenance header included, whatever `limit` and `offset` are.
  Must end in `.tsv`. For a large run pair it with `limit = 0`, so only
  the summary comes back.

- categories_out:

  Write the category table here with mzLib's `GoCategoryTsv` writer.
  Must end in `.tsv`, and needs `category_map`.

- limit:

  At most this many rows (rows, not groups) in `annotations`. `NULL`,
  the default, returns them all. Never shortens `out`.

- offset:

  Rows to skip before the window. Default 0.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait: a whole proteome
  XML takes a while.

## Details

**Every non-decoy group gets at least one row.** A group with no term
gets a single row whose `annotation_status` says why: `"no_go_terms"`
(its members have none), `"no_entry"` (a member is not in the database -
annotate against the database the search used), or `"contaminant"`.
Decoy groups get none.

**Pin the ontology.** Terms and their ancestors change between GO
releases, so a result means something only relative to one go.obo. This
reads the file you name and never downloads one: fetch a release on
purpose with
[`proteins_update_go`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_update_go.md),
keep the file, and record `go$sha256` with your results.

## Value

An `mzlibr_go_annotations`:

\- `annotations`, a data.frame with one row per (group, term) in the
returned window, under mzLib's own `GoAnnotationTsv` column names and
order (`column_names`). `accession_used`, `accession_direct`,
`accession_inherited`, `evidence` and `entrapment_members` are list
columns of character vectors; `evidence_by_member` is a list column of
named lists, member to its own evidence codes. `go_id`, `go_name`,
`aspect`, `inherited` and `propagated` are `NA` on a term-less row.
`n_members` and `n_with` count proteins (members); `n_with` is 0 on a
term-less row. `q_value` is the group's q-value, a fraction (0 to 1),
never filtered on. - `group_count` groups annotated (every non-decoy
group, each once); `table_row_count` rows in the table read, decoys
included; `decoy_group_count` decoy groups skipped. - `row_count` rows
in the whole table; `returned_count` rows in `annotations`; `offset`
rows skipped; `truncated` is `TRUE` when `limit` or `offset` left rows
out. - Provenance: `groups_file` and `groups_file_sha256`;
`annotation_database` (`path`, `file_type`, `reader`, `protein_count`,
`sha256` of the decompressed database, `caveats`); `go`
(`source_file_name`, `sha256`, `release`, `term_count`); and `header`,
mzLib's `#!` header as a named character vector. **Its counters count
groups at q \<= 0.01, not rows.** - `skip_unknown_go_ids`,
`unresolved_go_ids`, `caveats`. - `written` (`path`, `row_count`) and
`categories_written` (`path`), or `NULL` when not asked. - `categories`,
or `NULL` without `category_map`: `map_name`, `map_version`,
`source_file_name`, `sha256`, `anchor_count`, `row_count`, and
`records`, one row per (term, category, subcategory). `subcategory` is
`NA` where the term reaches only the category's own anchor; a term under
no anchor has no row. Join to `annotations` on `go_id`.

## The views are filters, not options

\- **Consensus** - every member carries the term:
`n_with == n_members`. - **Direct annotations only**: `!propagated`. -
**No borrowing**: `!inherited`. An isoform (`P04406-2`) or sequence
variant (`P04406_A20T`) absent from the database takes its entry's terms
and is listed in `accession_inherited`; an isoform can differ from its
entry in cellular component, which is why the row says so.

## Wraps

Wire verb `proteins annotate-go`. Generated from the bridge's verb spec
`proteins.annotate-go.yaml` (bridge commit `e76157831b15`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`GoGroupAnnotator.AnnotateAll`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GoGroupAnnotator.cs)
  in `mzLib/UsefulProteomicsDatabases/GeneOntology/GoGroupAnnotator.cs`
  at mzLib `0a808fec`

- [`GeneOntologyGraph.Load`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GeneOntologyGraph.cs)
  in `mzLib/UsefulProteomicsDatabases/GeneOntology/GeneOntologyGraph.cs`
  at mzLib `0a808fec`

- [`ProteinGroupFromTsvExtensions.ToGoAnnotationGroups`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsvExtensions.cs)
  in
  `mzLib/Readers/InternalResults/IndividualResultRecords/ProteinGroupFromTsvExtensions.cs`
  at mzLib `0a808fec`

- [`GoAnnotationTsv.Schema`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GoAnnotationTsv.cs)
  in `mzLib/UsefulProteomicsDatabases/GeneOntology/GoAnnotationTsv.cs`
  at mzLib `0a808fec`

- [`GoAnnotationTsv.Write`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GoAnnotationTsv.cs)
  in `mzLib/UsefulProteomicsDatabases/GeneOntology/GoAnnotationTsv.cs`
  at mzLib `0a808fec`

- [`GoCategoryResolver.Categorize`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GoCategoryResolver.cs)
  in
  `mzLib/UsefulProteomicsDatabases/GeneOntology/GoCategoryResolver.cs`
  at mzLib `0a808fec`

- [`GoCategoryTsv.Write`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GoCategoryTsv.cs)
  in `mzLib/UsefulProteomicsDatabases/GeneOntology/GoCategoryTsv.cs` at
  mzLib `0a808fec`

- [`ProteinDbLoader.LoadProteinXML`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `0a808fec`

## Parameters: units, ranges and defaults

- `groups`:

  path; required. A MetaMorpheus protein-group table:
  AllQuantifiedProteinGroups.tsv, AllProteinGroups.tsv (no
  quantification, mzLib \#1365) or a file's \<file\>\_ProteinGroups.tsv.
  Read by ProteinGroupFromTsvFile.

- `database`:

  path; required. The UniProt XML (.xml or .xml.gz) whose GO terms
  annotate the members: the proteome the search used. Exactly one,
  because every row carries that database's sha256. A FASTA is a usage
  error: it carries no GO.

- `go_obo` (wire `--go-obo`):

  path; required. A go.obo file, which must exist. Nothing is downloaded
  (see proteins update-go), so the release is always one the caller
  chose to keep.

- `skip_unknown_go_ids` (wire `--skip-unknown-go-ids`):

  flag; default `FALSE`. Drop each GO id the database cites that the
  release lacks, listing it in unresolved_go_ids, instead of failing on
  the first run that meets one.

- `category_map` (wire `--category-map`):

  path; default absent. A consumer term-to-category map in mzLib's
  format (#!category_map_format 1, \#!map_name, \#!map_version; columns
  category, subcategory, anchor_go_id). Adds the categories table. mzLib
  ships no map.

- `out`:

  path; default absent; range `.tsv only, any case (PYB-1)`. Write the
  WHOLE table here with GoAnnotationTsv.Write, header included, whatever
  offset and limit are. Written atomically through \<out\>.partial. Must
  not name an input; its folder must exist.

- `categories_out` (wire `--categories-out`):

  path; default absent; range `.tsv only, any case (PYB-1)`. Write the
  category table here with GoCategoryTsv.Write. Needs category-map.

- `offset`:

  int; in **rows**; default `0`; range `>= 0`. Skip this many rows in
  the payload.

- `limit`:

  int; in **rows**; default absent; range `>= 0`. Return at most this
  many rows. Absent means every row. 0 returns only the summary, the
  usual pairing with out for a large run.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `groups_file`:

  string; never `NA`. Absolute path of the protein-group table.

- `groups_file_sha256`:

  string; never `NA`. Lower-case hex sha256 of the table's bytes; also
  the header's source_file_sha256.

- `table_row_count`:

  int; in **rows**; never `NA`. Rows in the table, decoys included.

- `decoy_group_count`:

  int; in **groups**; never `NA`. Decoy rows (label contains D, so an
  entrapment decoy ED too), skipped: decoys carry no GO.

- `group_count`:

  int; in **groups**; never `NA`. Groups annotated: every non-decoy
  group, contaminants included, each once even when MetaMorpheus wrote
  it twice (mzLib \#1366).

- `annotation_database`:

  object; never `NA`. {path, file_type (UniProtXml), reader,
  protein_count, sha256 (of the DECOMPRESSED bytes), caveats}.

- `go`:

  object; never `NA`. {source_file_name, sha256, release (the
  data-version; null when the file has none), term_count (obsolete terms
  included)}.

- `skip_unknown_go_ids`:

  bool; never `NA`. Whether skip-unknown-go-ids was given.

- `unresolved_go_ids`:

  string\[\]; never `NA`. GO ids the database cites that the release
  lacks, in ordinal order. Always empty without skip-unknown-go-ids (the
  call fails instead).

- `header`:

  map\<string,string\>; never `NA`. The \#!key value header
  GoAnnotationTsv.Write produced, read back: go_annotation_format,
  mzlib_version, mzlib_release, go_release, go_obo_sha256,
  annotation_db_sha256, source_file_sha256, unresolved_go_ids (only when
  some were dropped), counter_q_value_max, n_multi_member_groups and
  status\_\<status\>. The counters count GROUPS at q \<=
  counter_q_value_max, not rows.

- `row_count`:

  int; in **rows**; never `NA`. Rows in the whole table.

- `returned_count`:

  int; in **rows**; never `NA`. Rows in columns.

- `offset`:

  int; in **rows**; never `NA`. The offset applied.

- `truncated`:

  bool; never `NA`. True whenever rows were left out of columns by limit
  or offset.

- `caveats`:

  string\[\]; never `NA`. Ids dropped from the database, a go.obo with
  no data-version, and the database's own load caveats.

- `written`:

  object; `NA` when out was not given. {path, row_count}; row_count is
  the whole table.

- `categories`:

  object; `NA` when category-map was not given. {map_name, map_version,
  source_file_name, sha256, anchor_count, row_count, column_names
  \[go_id, category, subcategory\], columns}: GoCategoryTsv's rows read
  back. subcategory is "category:subcategory", or null where the term
  reaches only the category's own anchor. A term under no anchor has no
  row.

- `categories_written`:

  object; `NA` when categories-out was not given. {path}.

- `column_names`:

  string\[\]; never `NA`. Exactly GoAnnotationTsv.Schema's headers, in
  its order (a bridge test holds them equal).

- `annotations` (wire `columns`):

  table; never `NA`. Column name to per-row values, for the returned
  window.

## Columns of `annotations`

- `protein_group`:

  string; never `NA`. The group as MetaMorpheus named it, members
  \|-joined.

- `accession_used`:

  string\[\]; never `NA`. Members carrying the term, directly or by
  propagation, ordinal order. Empty on a term-less row.

- `accession_direct`:

  string\[\]; never `NA`. Of those, members annotated to this exact
  term.

- `accession_inherited`:

  string\[\]; never `NA`. Of those, members whose terms were borrowed: a
  UniProt isoform or mzLib sequence variant takes its entry's terms when
  it has none of its own.

- `go_id`:

  string; `NA` when a term-less row: the group has no term
  (annotation_status says why). The term's primary id (an alternative id
  in the database resolves to it).

- `go_name`:

  string; `NA` when a term-less row. The term's name in this release.

- `aspect`:

  string; `NA` when a term-less row. biological_process \|
  cellular_component \| molecular_function \| unknown, spelled by
  mzLib's own schema column.

- `evidence`:

  string\[\]; never `NA`. Evidence codes (ECO) pooled over the carrying
  members, ordinal order.

- `evidence_by_member`:

  map\<string,string\[\]\>; never `NA`. Member to its own evidence codes
  for this term; empty on a term-less row. The file spells it
  member=code,code;member=code.

- `inherited`:

  bool; `NA` when a term-less row. True when every carrying member's
  terms were borrowed (accession_inherited covers accession_used).

- `propagated`:

  bool; `NA` when a term-less row. True when no member is annotated to
  this exact term: it is implied by a more specific one, over is_a or
  part_of.

- `n_members`:

  int; in **proteins**; never `NA`. Members in the group.

- `n_with`:

  int; in **proteins**; never `NA`. Members carrying the term; 0 on a
  term-less row. n_with == n_members is the consensus view.

- `entrapment_members`:

  string\[\]; never `NA`. The group's entrapment members
  (ProteinDbLoader.IsEntrapmentAccession, or marked so by the database),
  on every row of the group.

- `annotation_status`:

  string; never `NA`. annotated \| no_go_terms \| no_entry \|
  contaminant (GoAnnotationTsv.StatusName). A row with a term is always
  annotated.

- `q_value`:

  float; in **fraction (0 to 1)**; never `NA`. The group's q-value from
  the table. Never filtered on.

- `go_release`:

  string; `NA` when the go.obo has no data-version header. The ontology
  release, e.g. releases/2026-07-26.

- `go_obo_sha256`:

  string; never `NA`. sha256 of the go.obo.

- `annotation_db_sha256`:

  string; never `NA`. sha256 of the decompressed annotation database.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  a missing groups table, database, go.obo or category map; a FASTA
  database; out or categories-out not .tsv, naming an input, naming each
  other, or in a folder that does not exist; categories-out without
  category-map; a negative offset or limit. Every one is raised before
  any file is read.

- `mzlib_bridge_error` (correctness):

  the database cites a GO id the release lacks, without
  skip-unknown-go-ids (InvalidDataException naming every id); a go.obo
  mzLib cannot parse, including a truncated one whose edges name absent
  terms (InvalidDataException with the line); a category map it refuses;
  a group written twice with different members (InvalidDataException); a
  table or database it cannot parse

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- The table is the UNION over a group's members, never a pick:
  MetaMorpheus chooses no leading protein. Consensus (n_with ==
  n_members), direct-only (propagated false) and no-borrowing (inherited
  false) are filters on the rows.

- Every non-decoy group gets at least one row; a group with no term gets
  one term-less row whose annotation_status says why. Decoys get none.

- The header counters count groups at q \<= 0.01, not rows; the rows are
  unfiltered and every group is present.

- Pin the ontology: terms and ancestors change between releases. This
  verb reads only an existing go.obo and never downloads one; record
  go.sha256 with the results.

- mzlib_release in the header is 'none' and mzlib_version is
  1.0.0+\<commit\>, because the bridge builds mzLib from source at the
  pinned commit; the commit names the release.

- An isoform or variant member absent from the database (or present
  without GO of its own) borrows its entry's terms and is listed in
  accession_inherited; an isoform can differ from its entry in cellular
  component, which is why the row says so.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.annotate_go`

- Rust (mzLibRust): `mzlib::proteins::annotate_go_with` with
  `GoAnnotateOptions`

- R (mzLibR): `proteins_annotate_go`

## Since

Wire protocol 1; pyMzLib 0.3.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Only one annotation database: a search over a target XML plus a
  separate contaminant database annotates contaminant members as
  no_entry unless their entries are in the XML. Whether GoGroupAnnotator
  should take several databases (and stamp several sha256s) is mzLib's
  call.

- No bulk (`--paths-stdin`) form yet. One table per call; a cohort of
  many tables against one proteome reloads the proteome each time.

- PYB-3 (a document too large to return is a usage error) is not built
  in the bridge yet; until it is, use out with limit 0 for a large
  table.

## References

- [doi:10.1093/genetics/iyad031](https://doi.org/10.1093/genetics/iyad031):
  The Gene Ontology knowledgebase in 2023 (GO Consortium, Genetics 224,
  iyad031), the ontology and annotations used

- [doi:10.1038/75556](https://doi.org/10.1038/75556): Gene Ontology:
  tool for the unification of biology (Ashburner et al., Nat Genet 25,
  25-29, 2000), the is_a / part_of graph propagated over

- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, whose XML GO dbReferences are the direct annotations (UniProt
  Consortium 2025, NAR)

## See also

[`proteins_update_go`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_update_go.md),
[`proteins_read`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)
for a database's own GO terms

## Examples

``` r

go <- proteins_annotate_go("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
  go_obo = "go-pxd036557.obo", category_map = "organelle_map.tsv")
go
#> <mzlibr_go_annotations> 5 groups, 563 (group, term) rows
#>   go.obo: go-pxd036557.obo, releases/2026-07-26, sha256 c1cdfd098c5c
#>   groups at q <= 0.01: annotated 4, no_go_terms 0, no_entry 0, contaminant 1
#>   categories: organelle v1, 30 rows
go$header[c("status_annotated", "status_contaminant")]
#>   status_annotated status_contaminant 
#>                "4"                "1" 
histones <- go$annotations[go$annotations$protein_group == "P0C0S5|Q71UI9", ]
c(rows = nrow(histones), consensus = sum(histones$n_with == histones$n_members))
#>      rows consensus 
#>       106        57 
head(go$categories$records)
#>        go_id      category       subcategory
#> 1 GO:0000785       nucleus nucleus:chromatin
#> 2 GO:0000786       nucleus nucleus:chromatin
#> 3 GO:0000791       nucleus nucleus:chromatin
#> 4 GO:0000792       nucleus nucleus:chromatin
#> 5 GO:0005576 extracellular              <NA>
#> 6 GO:0005634       nucleus              <NA>
```
