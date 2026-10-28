# Differential design and observations: diff design, diff observations. See scripts/spec-facts.R
# for what each list means.

R_TABLE[c("diff design", "diff observations")] <- "records"
DIFF_NORMALIZATION <- function(field) {
  as_r(paste0("normalization$", field), "the normalization is one object, a named list `normalization`, not a table")
}
R_DEVIATIONS[["diff design"]] <- list(
  "param.stdin" = not_here(paste(
    "the model is the factors, covariates, interactions, strata and batch arguments, sent one stdin",
    "line per term"
  )),
  "field.columns" = RECORDS,
  "table.strata" = as_r("strata", "a list of `mzlibr_design_stratum`, named by stratum key")
)
R_DEVIATIONS[["diff observations"]] <- list(
  "param.stdin" = as_r("runs", "a data.frame with file and sample columns, and optionally fraction and technical_replicate, sent one run per stdin line"),
  "field.columns" = RECORDS,
  "field.normalization.setting" = DIFF_NORMALIZATION("setting"),
  "field.normalization.reference_set_size" = DIFF_NORMALIZATION("reference_set_size"),
  "field.normalization.minimum_reference_set_size" = DIFF_NORMALIZATION("minimum_reference_set_size"),
  "field.normalization.shifts" = DIFF_NORMALIZATION("shifts"),
  "field.normalization.warnings" = DIFF_NORMALIZATION("warnings")
)

PARENT_MAP[c("diff_design", "diff_observations")] <- c("diff.design", "diff.observations")
PARENT_ADDITIONS[c("factor_spec", "covariate_spec")] <- c(
  "pyMzLib's `diff.FactorSpec(name, reference_level)`, a dataclass rather than a function",
  "pyMzLib's `diff.CovariateSpec(name, unit, scale, scaled_unit, centre)`, a dataclass rather than a function"
)
PARENT_OMISSIONS["Design.stratum"] <- "`strata[[key]]`: the strata are a list named by stratum key"

FIELD_CHECKS[["diff design"]] <- list("diff_design_malat.json", function(d) mz$diff_parse_design(d))
FIELD_CHECKS[["diff observations"]] <- list("diff_observations_mm_mbr_kept.json", function(d) mz$diff_parse_observations(d))

# The stats article's contrasts chain, beyond the spec's examples (pyMzLib's conftest REPLAY_EXTRA):
# the dilution series' design with no reference level, every pair of loads a contrast.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "diff design", fixture = "diff_design_malat_pairs.json")
