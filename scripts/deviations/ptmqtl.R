# PTM site occupancy and PTM-QTL: ptmqtl occupancy, pairs and trait-effects. See
# scripts/spec-facts.R for what each list means.

R_TABLE[c("ptmqtl occupancy", "ptmqtl pairs", "ptmqtl trait-effects")] <- "records"
R_DEVIATIONS[["ptmqtl occupancy"]] <- list("field.columns" = RECORDS)
R_DEVIATIONS[["ptmqtl pairs"]] <- list("field.columns" = RECORDS)
R_DEVIATIONS[["ptmqtl trait-effects"]] <- list("field.columns" = RECORDS)

PARENT_MAP[c("ptmqtl_occupancy", "ptmqtl_pairs", "ptmqtl_trait_effects")] <- c(
  "ptmqtl.occupancy", "ptmqtl.pairs", "ptmqtl.trait_effects"
)

FIELD_CHECKS[["ptmqtl occupancy"]] <- list("ptmqtl_occupancy_fractions.json", function(d) mz$ptmqtl_parse_occupancy(d))
FIELD_CHECKS[["ptmqtl pairs"]] <- list("ptmqtl_pairs_batch_run_level.json", function(d) mz$ptmqtl_parse_pairs(d))
FIELD_CHECKS[["ptmqtl trait-effects"]] <- list("ptmqtl_trait_effects.json", function(d) mz$ptmqtl_parse_trait_effects(d))
