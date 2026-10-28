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

# Not yet ported: three options the pyMzLib 0.7.0 bridge added (mzLib #1464, #1468). The fields that
# report them (include_bounds, min_untied_share, run_covariate_*, scale, and the pair columns
# bounds_included and run_covariate_removed) are projected; the options are a later port.
PTMQTL_NOT_YET <- "not ported yet: the pyMzLib 0.7.0 bridge's option, a later mzLibR port"
R_DEVIATIONS[["ptmqtl pairs"]][["param.include-bounds"]] <- not_here(PTMQTL_NOT_YET)
R_DEVIATIONS[["ptmqtl pairs"]][["param.run-covariate"]] <- not_here(PTMQTL_NOT_YET)
R_DEVIATIONS[["ptmqtl trait-effects"]][["param.scale"]] <- not_here(paste0(
  PTMQTL_NOT_YET, "; mzLib's default, natural_log_odds, always applies, and the result's scale says so"
))
PARENT_PARAM_OMISSIONS[c(
  "ptmqtl_pairs.include_bounds", "ptmqtl_pairs.run_covariate", "ptmqtl_trait_effects.scale"
)] <- c(
  "not ported yet (mzLib #1468, wire `--include-bounds`); `include_bounds` in the result reports mzLib's default",
  "not ported yet (mzLib #1468, wire `--run-covariate`); `run_covariate_file` in the result is `NA`",
  "not ported yet (mzLib #1464, wire `--scale`); `scale` in the result reports mzLib's default"
)
