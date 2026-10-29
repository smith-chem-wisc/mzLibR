# Protein databases: proteins read, genes resolve, proteins classify-peptides. See
# scripts/spec-facts.R for what each list means.

R_TABLE["proteins read"] <- "proteins"
R_TABLE["genes resolve"] <- "resolutions"
R_TABLE["proteins classify-peptides"] <- "peptides"

# The protein verbs take one database or a character vector of them as their first argument, and
# contaminant databases as a separate vector, so the wire's --path / --paths-stdin / --contaminant
# choice is made from the arguments' shape rather than spelled by the caller - as in pyMzLib.
PROTEINS_DATABASE_PARAMS <- list(
  "param.path" = as_r("databases", "one database or several; several travel as --paths-stdin"),
  "param.contaminant" = as_r("contaminants", "contaminant databases are their own vector, tagged on stdin or sent as --contaminant"),
  "param.paths-stdin" = not_here("chosen automatically when more than one database is given")
)

# Every wire `source_index` is 0-based; R's is 1-based, so it indexes `files` directly.
PROTEINS_SOURCE_INDEX <- as_r("source_index", "1-based in R, so it indexes `files` directly")

R_DEVIATIONS[["proteins read"]] <- c(PROTEINS_DATABASE_PARAMS, list(
  "param.accessions-stdin" = as_r("accessions", "the accession vector itself; giving it sets the flag and fills stdin"),
  "field.columns" = as_r("proteins", "the proteins table is a data.frame, so it is `proteins`"),
  "field.source_index" = PROTEINS_SOURCE_INDEX,
  "field.go_terms.source_index" = PROTEINS_SOURCE_INDEX,
  "field.ensembl_genes.source_index" = PROTEINS_SOURCE_INDEX
))
R_DEVIATIONS[["genes resolve"]] <- c(PROTEINS_DATABASE_PARAMS, list(
  "field.columns" = as_r("resolutions", "the resolution table is a data.frame, so it is `resolutions`"),
  "field.source_index" = PROTEINS_SOURCE_INDEX
))
R_DEVIATIONS[["proteins classify-peptides"]] <- c(PROTEINS_DATABASE_PARAMS, list(
  "param.on-error" = not_here("the only accepted value is fail, the default; there is no skip to choose"),
  "field.columns" = as_r("peptides", "the classification table is a data.frame, one row per peptide, so it is `peptides`")
))

PARENT_MAP[c("proteins_read", "proteins_resolve_genes", "proteins_classify_peptides")] <- c(
  "proteins.read", "proteins.resolve_genes", "proteins.classify_peptides"
)
PARENT_OMISSIONS[c(
  "ProteinDatabase.taxonomy", "ProteinDatabase.organisms", "PeptideClassification.sharing_of"
)] <- c(
  "two columns of `proteins`: `setNames(db$proteins$ncbi_taxonomy_id, db$proteins$accession)`",
  "two columns of `proteins`: `setNames(db$proteins$organism, db$proteins$accession)`",
  "two columns of `peptides`: `setNames(calls$peptides$sharing, calls$peptides$peptide)`"
)

FIELD_CHECKS[["proteins read"]] <- list("proteins_read_human.json", function(d) mz$proteins_parse_database(d))
FIELD_CHECKS[["genes resolve"]] <- list("genes_resolve_human.json", function(d) mz$proteins_parse_resolutions(d))
FIELD_CHECKS[["proteins classify-peptides"]] <- list("proteins_classify_peptides.json", function(d) mz$proteins_parse_classification(d))

# ---------------------------------------------------------------- annotate-go, update-go

R_TABLE["proteins annotate-go"] <- "annotations"
R_DEVIATIONS[["proteins annotate-go"]] <- list(
  "field.columns" = as_r("annotations", "the annotation table is a data.frame, one row per (group, term), so it is `annotations`")
)

PARENT_MAP[c("proteins_annotate_go", "proteins_update_go")] <- c("proteins.annotate_go", "proteins.update_go")

FIELD_CHECKS[["proteins annotate-go"]] <- list("proteins_annotate_go_pxd036557.json", function(d) mz$proteins_parse_go_annotations(d))
FIELD_CHECKS[["proteins update-go"]] <- list("proteins_update_go.json", function(d) mz$proteins_parse_go_update(d))

# Real-data recordings the articles replay (pyMzLib #72), beyond the spec's examples.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "proteins read", fixture = "proteins_read_accessions.json")

# ---------------------------------------------------------------- align, parse-accession (pyMzLib 0.6.0)

R_TABLE["proteins align"] <- "records"
PROTEINS_PAIR_INDEX <- as_r("pair_index", "1-based in R, so it indexes `pairs` directly")
R_DEVIATIONS[["proteins align"]] <- list(
  "param.no-residue-map" = as_r("residue_map", "stated positively: residue_map = FALSE sends --no-residue-map"),
  "field.columns" = RECORDS,
  "field.pair_index" = PROTEINS_PAIR_INDEX,
  "field.residue_map.pair_index" = PROTEINS_PAIR_INDEX
)
R_TABLE["proteins parse-accession"] <- "records"
R_DEVIATIONS[["proteins parse-accession"]] <- list("field.columns" = RECORDS)

PARENT_MAP[c(
  "proteins_align", "proteins_alignments", "proteins_map_source_to_target",
  "proteins_map_target_to_source", "proteins_parse_accessions"
)] <- c(
  "proteins.align_many", "Alignments.alignments", "Alignment.map_source_to_target",
  "Alignment.map_target_to_source", "proteins.parse_accessions"
)
PARENT_OMISSIONS[c("proteins.align", "ParsedAccessions.entry_of")] <- c(
  "proteins_align(list(c(source, target))): one pair is a list of one, as the spec spells the R call",
  "two columns of `records`: `setNames(parsed$records$entry_accession, parsed$records$input)`"
)

FIELD_CHECKS[["proteins align"]] <- list("proteins_align_many.json", function(d) mz$proteins_parse_alignments(d))
FIELD_CHECKS[["proteins parse-accession"]] <- list("proteins_parse_accession.json", function(d) mz$proteins_parse_accessions_result(d))

# The proteins article reads the DYNC1I2 isoforms and both albumins with their sequences, which its
# alignment section aligns.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "proteins read", fixture = "proteins_read_sequences.json")

# ---------------------------------------------------------------- genes orthology (pyMzLib 0.6.0)

# A list in R rather than a data.frame: one gene-tree file, one written snapshot, one pair status.
orthology_in <- function(element, why) {
  function(field) as_r(paste0(element, "$", field), why)
}
orthology_trees <- orthology_in("gene_trees", "an element of the gene_trees list")
orthology_written <- orthology_in("written", "an element of the written list, NULL without out")
orthology_status <- orthology_in("pair_status", "an element of the pair_status list, NULL without pair_status")
ORTHOLOGY_FILES <- not_here("a column of the written$files data.frame")
ORTHOLOGY_STATUS_ROWS <- not_here("a column of the pair_status$records data.frame")
ORTHOLOGY_PARQUET <- not_here("a column of the Parquet files out = writes, which mzLibR does not read")

R_DEVIATIONS[["genes orthology"]] <- c(
  list(
    "param.species" = not_here(paste(
      "rendered from gtfs and gene_sets, two named character vectors from Ensembl production name",
      "to path, one stdin line per species"
    )),
    "field.gene_sets" = as_r("gene_sets", "a data.frame, one row per species, `species` first"),
    "field.dumps" = as_r("dumps", "a data.frame, one row per species, `species` first"),
    "field.pairs" = as_r("pairs", "a data.frame, one row per species pair"),
    "field.gene_trees.source_file_name" = orthology_trees("source_file_name"),
    "field.gene_trees.sha256" = orthology_trees("sha256"),
    "field.gene_trees.release" = orthology_trees("release"),
    "field.gene_trees.collection" = orthology_trees("collection"),
    "field.gene_trees.member_count" = orthology_trees("member_count"),
    "field.gene_trees.restricted_to_gene_sets" = orthology_trees("restricted_to_gene_sets"),
    "field.written.directory" = orthology_written("directory"),
    "field.written.files" = orthology_written("files"),
    "field.pair_status.from_species" = orthology_status("from_species"),
    "field.pair_status.to_species" = orthology_status("to_species"),
    "field.pair_status.row_count" = orthology_status("row_count"),
    "field.pair_status.status_counts" = orthology_status("status_counts"),
    "field.pair_status.column_names" = orthology_status("column_names"),
    "field.pair_status.columns" = as_r("pair_status$records", "the gene table is a data.frame, so it is `records`"),
    "table.written_files" = as_r("written$files", "a data.frame, one row per file written but manifest.json"),
    "table.pair_status_columns" = as_r("pair_status$records", "a data.frame, one row per gene"),
    "table.parquet_genes" = not_here("the columns of genes/SPECIES.parquet in a written snapshot; mzLibR does not read Parquet"),
    "table.parquet_members" = not_here("the columns of members/SPECIES.parquet in a written snapshot"),
    "table.parquet_pairs" = not_here("the columns of pairs/A__B.parquet in a written snapshot")
  ),
  stats::setNames(rep(list(ORTHOLOGY_FILES), 4L), paste0("field.written_files.", c("path", "rows", "bytes", "sha256"))),
  stats::setNames(rep(list(ORTHOLOGY_STATUS_ROWS), 4L), paste0("field.pair_status_columns.", c("gene_id", "gene_name", "gene_biotype", "status"))),
  stats::setNames(rep(list(ORTHOLOGY_PARQUET), 6L), paste0("field.parquet_genes.", c(
    "gene_id", "gene_version", "gene_biotype", "gene_name", "seq_region", "species"
  ))),
  stats::setNames(rep(list(ORTHOLOGY_PARQUET), 5L), paste0("field.parquet_members.", c(
    "gene_id", "species", "group_id", "source_group_id", "canonical_protein_id"
  ))),
  stats::setNames(rep(list(ORTHOLOGY_PARQUET), 17L), paste0("field.parquet_pairs.", c(
    "homology_id", "relationship_type", "relationship_class", "species_a", "species_b", "gene_a",
    "protein_a", "identity_a", "gene_b", "protein_b", "identity_b", "dn", "ds", "goc_score",
    "wga_coverage", "is_high_confidence", "source_dump"
  )))
)
rm(orthology_in, orthology_trees, orthology_written, orthology_status)

PARENT_MAP[c("proteins_build_orthology", "proteins_orthology_pair")] <- c(
  "proteins.build_orthology", "OrthologySnapshot.pair"
)
PARENT_OMISSIONS[c("PairStatus.genes", "WrittenSnapshot.manifest")] <- c(
  "a row filter on `records`: `ps$records$gene_id[ps$records$status == status]`",
  "an element, `written$manifest`, not a function"
)

# The pair-status recording has no written snapshot and the written one no pair status, so the
# check of the spec's fields reads the one and takes `written` from the other.
FIELD_CHECKS[["genes orthology"]] <- list("genes_orthology_pair_status.json", function(d) {
  snapshot <- mz$proteins_parse_orthology(d)
  snapshot$written <- mz$proteins_parse_orthology(read_fixture("genes_orthology_written.json"))$written
  snapshot
})

# Not yet ported: the pyMzLib 0.7.0 bridge's genes map-residues (mzLib #1465), a later mzLibR port.
PARENT_OMISSIONS[c("proteins.map_residues", "AlignedProteins.sequences")] <- c(
  "not ported yet: wire `genes map-residues` (mzLib #1465), new in the pyMzLib 0.7.0 bridge",
  "not ported yet: part of `proteins.map_residues`'s result"
)
