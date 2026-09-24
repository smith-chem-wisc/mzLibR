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
