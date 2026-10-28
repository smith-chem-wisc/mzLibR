# Differential abundance: stats fit, adjust, meta, combine-p, correlate, fit-logistic, fit-mixed
# and entrapment-fdp. See scripts/spec-facts.R for what each list means.

R_TABLE[c(
  "stats fit", "stats adjust", "stats meta", "stats combine-p", "stats correlate", "stats fit-logistic",
  "stats fit-mixed", "stats entrapment-fdp"
)] <- "records"
STATS_PRIOR <- function(field) as_r(paste0("prior$", field), "the prior is one object, a named list `prior`, not a table")
STATS_CONTRASTS <- function(field) {
  as_r(paste0("contrasts$", field), "the contrasts are one object, an `mzlibr_moderated_contrasts`, beside the coefficient table")
}
STATS_COEFFICIENTS <- as_r("coefficients", "a character vector, sent one name per stdin line")
R_DEVIATIONS[["stats fit"]] <- list(
  "param.stdin" = STATS_COEFFICIENTS,
  "field.columns" = RECORDS,
  "field.prior.estimator" = STATS_PRIOR("estimator"),
  "field.prior.df" = STATS_PRIOR("df"),
  "field.prior.df_infinite" = STATS_PRIOR("df_infinite"),
  "field.prior.trended" = STATS_PRIOR("trended"),
  "field.prior.spline_basis_count" = STATS_PRIOR("spline_basis_count"),
  "field.prior.scale" = STATS_PRIOR("scale"),
  "field.contrasts.confidence_level" = STATS_CONTRASTS("confidence_level"),
  "field.contrasts.tested_contrasts" = STATS_CONTRASTS("tested_contrasts"),
  "field.contrasts.weights" = STATS_CONTRASTS("weights"),
  "field.contrasts.row_count" = STATS_CONTRASTS("row_count"),
  "field.contrasts.column_names" = STATS_CONTRASTS("column_names"),
  "field.contrasts.columns" = as_r("contrasts$records", "the contrasts table is a data.frame, so it is `contrasts$records`")
)
R_DEVIATIONS[["stats adjust"]] <- list(
  "param.stdin" = as_r("p_values", "a numeric vector, sent one per stdin line; NA crosses as an empty line"),
  "field.columns" = RECORDS
)
R_DEVIATIONS[["stats meta"]] <- list(
  "param.stdin" = as_r("studies", "a data.frame with feature, estimate and standard_error columns, sent one study per stdin line"),
  "field.columns" = RECORDS
)
R_DEVIATIONS[["stats combine-p"]] <- list(
  "param.stdin" = as_r("tests", "a data.frame with feature and p_value columns, and optionally weight and effect, sent one test per stdin line"),
  "field.columns" = RECORDS
)
R_DEVIATIONS[["stats correlate"]] <- list("field.columns" = RECORDS)
R_DEVIATIONS[["stats fit-logistic"]] <- list("param.stdin" = STATS_COEFFICIENTS, "field.columns" = RECORDS)
R_DEVIATIONS[["stats fit-mixed"]] <- list("param.stdin" = STATS_COEFFICIENTS, "field.columns" = RECORDS)
R_DEVIATIONS[["stats entrapment-fdp"]] <- list("field.columns" = RECORDS)

PARENT_MAP[c(
  "stats_fit", "stats_adjust", "stats_meta", "stats_combine_p", "stats_correlate", "stats_fit_logistic",
  "stats_fit_mixed", "stats_entrapment_fdp"
)] <- c(
  "stats.fit", "stats.adjust", "stats.meta", "stats.combine_p", "stats.correlate", "stats.fit_logistic",
  "stats.fit_mixed", "stats.entrapment_fdp"
)
PARENT_OMISSIONS[c(
  "ModeratedFit.rows", "ModeratedContrasts.rows", "Adjusted.bh_adjusted", "CombinedPValues.p_values",
  "EntrapmentFdpEstimates.at"
)] <- c(
  "`records[records$coefficient == name, ]`",
  "`contrasts$records[contrasts$records$contrast == id, ]`",
  "`records$bh_adjusted`",
  "`setNames(records$p_value, records$feature)`",
  "`records[records$q_value_threshold == threshold, ]`, zero rows where nothing was discovered"
)

# The contrasts recording, so the field check sees the contrasts object as well as the (empty)
# coefficient table.
FIELD_CHECKS[["stats fit"]] <- list("stats_fit_malat_contrasts.json", function(d) mz$stats_parse_fit(d))
FIELD_CHECKS[["stats adjust"]] <- list("stats_adjust.json", function(d) mz$stats_parse_adjusted(d))
FIELD_CHECKS[["stats meta"]] <- list("stats_meta_metafor.json", function(d) mz$stats_parse_meta(d))
FIELD_CHECKS[["stats combine-p"]] <- list("stats_combine_p_signed.json", function(d) mz$stats_parse_combined(d))
FIELD_CHECKS[["stats correlate"]] <- list("stats_correlate_limma_age_given_sex.json", function(d) mz$stats_parse_correlations(d))
FIELD_CHECKS[["stats fit-logistic"]] <- list("stats_fit_logistic_glm.json", function(d) mz$stats_parse_logistic(d))
FIELD_CHECKS[["stats fit-mixed"]] <- list("stats_fit_mixed_lmm_reml.json", function(d) mz$stats_parse_mixed(d))
FIELD_CHECKS[["stats entrapment-fdp"]] <- list("stats_entrapment_fdp_reference.json", function(d) mz$stats_parse_entrapment(d))

# The diff article's recordings, beyond the spec's examples (pyMzLib's conftest REPLAY_EXTRA): stats
# fit on the design file diff design wrote for the MALAT1 dilution series, without and with its
# contrasts file.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "stats fit", fixture = "stats_fit_malat_diff_design.json")
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "stats fit", fixture = "stats_fit_malat_diff_design_contrasts.json")

# Not yet ported: the pyMzLib 0.7.0 bridge's latent-scale mixed model and its detection link
# (stats fit-latent-mixed, stats fit-link; mzLib #1469), a later mzLibR port.
PARENT_OMISSIONS[c("stats.fit_latent_mixed", "stats.fit_link")] <- c(
  "not ported yet: wire `stats fit-latent-mixed` (mzLib #1469), new in the pyMzLib 0.7.0 bridge",
  "not ported yet: wire `stats fit-link` (mzLib #1469), new in the pyMzLib 0.7.0 bridge"
)
