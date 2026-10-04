# Peptidoforms. See scripts/spec-facts.R for what each list means.

R_TABLE["peptidoform fragments"] <- "peptides"

R_DEVIATIONS[["peptidoform fragments"]] <- list(
  "param.no-modifications" = as_r("modifications", "stated positively: modifications = FALSE sends --no-modifications"),
  "param.max-mods" = as_r("max_modifications", "spelled out"),
  "field.annotated_modification_sites" = as_r("census$sites", "gathered into the census; see census_explain()"),
  "field.annotated_modifications_loaded" = as_r("census$applied", "gathered into the census; see census_explain()"),
  "field.uniprot_annotated_features" = as_r("census$annotated", "gathered into the census; see census_explain()"),
  "field.unresolved_modifications" = as_r("census$unresolved", "gathered into the census; see census_excluded()"),
  "field.uniprot_features_by_type" = as_r("census$by_type", "gathered into the census, as a data.frame"),
  "field.peptide_count" = not_here("nrow() of `peptides`"),
  "field.modifications" = as_r("modifications", "its own long data.frame, joined to `peptides` on peptide_index"),
  "field.fragments" = as_r("fragments", "its own long data.frame, joined to `peptides` on peptide_index")
)

PARENT_MAP[c(
  "peptidoform_fragments", "digest_truncated", "digest_modified_peptides", "census_explain",
  "census_excluded", "peptide_mz"
)] <- c(
  "peptidoform.fragments", "Digest.truncated", "Digest.modified_peptides",
  "ModificationCensus.explain", "ModificationCensus.excluded", "Peptide.mz"
)
PARENT_ADDITIONS["digest_distinct_base_sequences"] <-
  "from mzLibRust (Digest::distinct_base_sequences); pyMzLib has no equivalent"
PARENT_ADDITIONS["digest_fragments_by_series"] <-
  "from mzLibRust (Digest::fragments_by_series); pyMzLib has only fragment_count"
PARENT_OMISSIONS[c("Peptide.is_modified", "Digest.fragment_count")] <- c(
  "a column comparison, `modification_count > 0`",
  "`nrow(digest$fragments)`; not promoted to a function because a bare total folds in the spurious ETD y series — see ?digest_fragments_by_series"
)

FIELD_CHECKS[["peptidoform fragments"]] <- list("peptidoform_P02768_small.json", function(d) mz$peptidoform_parse(d))

# Real-data recordings the articles replay (pyMzLib #72), beyond the spec's examples.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "peptidoform fragments", fixture = "peptidoform_P02768.json")

# peptidoform convert: pyMzLib #75.
R_TABLE["peptidoform convert"] <- "records"
R_DEVIATIONS[["peptidoform convert"]] <- list(
  "param.from" = as_r("source", "pyMzLib's name, the pair to target; mzLib calls it sourceFormat"),
  "param.to" = as_r("target", "the pair to source; mzLib calls it targetFormat"),
  "field.columns" = RECORDS
)
PARENT_MAP["peptidoform_convert"] <- "peptidoform.convert"
PARENT_OMISSIONS[c(
  "SequenceConversions.outputs", "SequenceConversions.sequences",
  "SequenceConversions.not_converted", "ConvertedSequence.ok"
)] <- c(
  "`records$output`",
  "the rows are `records`, a data.frame",
  "`records[records$status != \"converted\", ]`",
  "`records$status == \"converted\"`"
)
FIELD_CHECKS[["peptidoform convert"]] <- list("peptidoform_convert_unimod.json", function(d) mz$peptidoform_parse_conversions(d))

# The peptidoforms article reads MetaMorpheus's BottomUpExample.psmtsv (mzLib's test data) before
# converting its full sequences (pyMzLib #75).
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "readers read-results", fixture = "readers_results_psmtsv.json")
