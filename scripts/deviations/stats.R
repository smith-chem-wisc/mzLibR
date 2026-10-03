# Differential abundance: stats fit, stats adjust, stats meta. See scripts/spec-facts.R for what
# each list means.

R_TABLE[c("stats fit", "stats adjust", "stats meta")] <- "records"
STATS_PRIOR <- function(field) as_r(paste0("prior$", field), "the prior is one object, a named list `prior`, not a table")
R_DEVIATIONS[["stats fit"]] <- list(
  "param.stdin" = as_r("coefficients", "a character vector, sent one name per stdin line"),
  "field.columns" = RECORDS,
  "field.prior.df" = STATS_PRIOR("df"),
  "field.prior.df_infinite" = STATS_PRIOR("df_infinite"),
  "field.prior.trended" = STATS_PRIOR("trended"),
  "field.prior.spline_basis_count" = STATS_PRIOR("spline_basis_count"),
  "field.prior.scale" = STATS_PRIOR("scale")
)
R_DEVIATIONS[["stats adjust"]] <- list(
  "param.stdin" = as_r("p_values", "a numeric vector, sent one per stdin line; NA crosses as an empty line"),
  "field.columns" = RECORDS
)
R_DEVIATIONS[["stats meta"]] <- list(
  "param.stdin" = as_r("studies", "a data.frame with feature, estimate and standard_error columns, sent one study per stdin line"),
  "field.columns" = RECORDS
)

PARENT_MAP[c("stats_fit", "stats_adjust", "stats_meta")] <- c("stats.fit", "stats.adjust", "stats.meta")
PARENT_OMISSIONS[c("ModeratedFit.rows", "Adjusted.bh_adjusted")] <- c(
  "`records[records$coefficient == name, ]`",
  "`records$bh_adjusted`"
)

FIELD_CHECKS[["stats fit"]] <- list("stats_fit_limma.json", function(d) mz$stats_parse_fit(d))
FIELD_CHECKS[["stats adjust"]] <- list("stats_adjust.json", function(d) mz$stats_parse_adjusted(d))
FIELD_CHECKS[["stats meta"]] <- list("stats_meta_metafor.json", function(d) mz$stats_parse_meta(d))
