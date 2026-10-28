# PTM site occupancy and PTM-QTL statistics, computed by mzLib's Quantification.PtmQtl.
#
#   ptmqtl_occupancy()      SiteOccupancyCalculator.Calculate: every site's intensity occupancy in
#                           every run, each cell with its state
#   ptmqtl_pairs()          PtmPairEngine.Physical / .CoVarying: modifications on one molecule, or
#                           rising and falling together across runs
#   ptmqtl_trait_effects()  SiteTraitEffects.Fit: per site, a trait's effect on logit occupancy,
#                           with a random intercept per biological replicate
#
# Ported from pyMzLib's `ptmqtl.py` (mzLib #1369, #1430, #1446). Every number is mzLib's: this file
# names the files and options and projects each result into a data.frame, `records`.

# The first pyMzLib whose bridge has the ptmqtl verbs (their specs' since.pymzlib).
PTMQTL_SINCE <- "0.6.0"

# Every value of the occupancy `state` column, as mzLib's OccupancyState is written on the wire.
# Only "quantified" cells have a `fraction`.
PTMQTL_OCCUPANCY_STATES <- c("quantified", "floor", "count_only", "not_detected")

# What ptmqtl_pairs() takes as `type`. There is no default; you choose.
PTMQTL_PAIR_TYPES <- c("physical", "covarying", "both")

# The named rules ptmqtl_pairs() takes as `run_level_group`. For any other grouping, map each
# modification to a group in a file and pass it as `run_level_groups`.
PTMQTL_RUN_LEVEL_GROUPS <- c("category", "modification")

# Every value of the ptmqtl_trait_effects() `status` column. Only "fitted" rows carry statistics.
PTMQTL_TRAIT_STATUSES <- c(
  "fitted", "below_support", "too_few_observations", "rank_deficient", "too_few_groups",
  "variance_ratio_at_limit", "not_converged"
)

# ---------------------------------------------------------------- argument shaping

# A character vector of names as the wire's one tab-separated option.
ptmqtl_names <- function(value, what) {
  if (!is.character(value) || length(value) == 0L) {
    stop(mzlib_usage_error(paste0(what, " takes a non-empty character vector of names; got ",
      paste(deparse(value), collapse = " "), ".")))
  }
  for (name in value) {
    if (is.na(name) || !nzchar(trimws(name))) {
      stop(mzlib_usage_error(paste0("Every name in ", what, " must be a non-empty string; got a blank or NA.")))
    }
    if (grepl("[\t\r\n]", name)) {
      stop(mzlib_usage_error(paste0("A name in ", what, " contains a tab or a line break: '", name, "'.")))
    }
  }
  paste(unname(value), collapse = "\t")
}

ptmqtl_filter_args <- function(modifications, categories) {
  if (!is.null(modifications) && !is.null(categories)) {
    stop(mzlib_usage_error(
      "Give modifications or categories, not both: mzLib takes one filter on the modification name."
    ))
  }
  if (!is.null(modifications)) {
    return(c("--modifications", ptmqtl_names(modifications, "modifications")))
  }
  if (!is.null(categories)) {
    return(c("--categories", ptmqtl_names(categories, "categories")))
  }
  character(0)
}

ptmqtl_build_occupancy_args <- function(observations, samples, modifications, categories, out) {
  args <- c("ptmqtl", "occupancy", "--observations", readers_normalise_path(observations, "observations"))
  if (!is.null(samples)) {
    args <- c(args, "--samples", readers_normalise_path(samples, "samples"))
  }
  args <- c(args, ptmqtl_filter_args(modifications, categories))
  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  args
}

ptmqtl_build_pairs_args <- function(observations, type, samples, modifications, categories,
                                    min_quantified_fraction, keep_ceiling, remove_run_level,
                                    run_level_group, run_level_groups) {
  if (!is.character(type) || length(type) != 1L || !type %in% PTMQTL_PAIR_TYPES) {
    stop(mzlib_usage_error(paste0(
      "type must be one of ", paste(PTMQTL_PAIR_TYPES, collapse = ", "), "; got ",
      paste(deparse(type), collapse = " "), ". mzLib has no default between its two methods."
    )))
  }
  args <- c("ptmqtl", "pairs", "--observations", readers_normalise_path(observations, "observations"), "--type", type)
  if (!is.null(samples)) {
    args <- c(args, "--samples", readers_normalise_path(samples, "samples"))
  }
  args <- c(args, ptmqtl_filter_args(modifications, categories))
  if (!is.null(min_quantified_fraction)) {
    args <- c(args, "--min-quantified-fraction", stats_number(min_quantified_fraction, "min_quantified_fraction"))
  }
  # A flag left FALSE sends nothing, so mzLib's default applies.
  if (stats_flag(keep_ceiling, "keep_ceiling")) {
    args <- c(args, "--keep-ceiling")
  }
  if (stats_flag(remove_run_level, "remove_run_level")) {
    args <- c(args, "--remove-run-level")
  }
  if (!is.null(run_level_group)) {
    if (!is.character(run_level_group) || length(run_level_group) != 1L || !run_level_group %in% PTMQTL_RUN_LEVEL_GROUPS) {
      stop(mzlib_usage_error(paste0(
        "run_level_group must be one of ", paste(PTMQTL_RUN_LEVEL_GROUPS, collapse = ", "), "; got ",
        paste(deparse(run_level_group), collapse = " "),
        ". For any other grouping, pass a modification-to-group file as run_level_groups."
      )))
    }
    args <- c(args, "--run-level-group", run_level_group)
  }
  if (!is.null(run_level_groups)) {
    args <- c(args, "--run-level-groups", readers_normalise_path(run_level_groups, "run_level_groups"))
  }
  args
}

ptmqtl_build_trait_args <- function(occupancy, runs, min_runs, min_replicates, logit_epsilon, keep_ceiling) {
  args <- c(
    "ptmqtl", "trait-effects",
    "--occupancy", readers_normalise_path(occupancy, "occupancy"),
    "--runs", readers_normalise_path(runs, "runs")
  )
  if (!is.null(min_runs)) {
    args <- c(args, "--min-runs", stats_whole(min_runs, "min_runs"))
  }
  if (!is.null(min_replicates)) {
    args <- c(args, "--min-replicates", stats_whole(min_replicates, "min_replicates"))
  }
  if (!is.null(logit_epsilon)) {
    args <- c(args, "--logit-epsilon", stats_number(logit_epsilon, "logit_epsilon"))
  }
  if (stats_flag(keep_ceiling, "keep_ceiling")) {
    args <- c(args, "--keep-ceiling")
  }
  args
}

# ---------------------------------------------------------------- parsing

ptmqtl_optional_count <- function(data, name) {
  value <- data[[name]]
  if (wire_null(value)) NA_real_ else as.numeric(value)
}

ptmqtl_optional_flag <- function(data, name) {
  value <- data[[name]]
  if (wire_null(value)) NA else isTRUE(value)
}

# A filter list: NULL when no filter was given, which differs from an empty filter.
ptmqtl_optional_names <- function(value) {
  if (wire_null(value)) NULL else wire_strings(value)
}

ptmqtl_parse_occupancy <- function(data) {
  written <- data[["written"]]
  structure(
    list(
      observations_file = wire_field(data, "observations_file", "character", NA_character_),
      samples_file = wire_field(data, "samples_file", "character", NA_character_),
      observation_count = stats_count(data, "observation_count"),
      run_count = stats_count(data, "run_count"),
      sample_count = ptmqtl_optional_count(data, "sample_count"),
      run_names = wire_strings(data[["run_names"]]),
      site_count = stats_count(data, "site_count"),
      state_counts = proteins_parse_counts(data[["state_counts"]]),
      modifications = ptmqtl_optional_names(data[["modifications"]]),
      categories = ptmqtl_optional_names(data[["categories"]]),
      ignored_columns = wire_strings(data[["ignored_columns"]]),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data),
      written = if (is.list(written)) {
        list(
          path = wire_field(written, "path", "character", NA_character_),
          format = wire_field(written, "format", "character", NA_character_),
          row_count = stats_count(written, "row_count")
        )
      } else {
        NULL
      },
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_site_occupancy"
  )
}

ptmqtl_parse_pairs <- function(data) {
  structure(
    list(
      observations_file = wire_field(data, "observations_file", "character", NA_character_),
      samples_file = wire_field(data, "samples_file", "character", NA_character_),
      type = wire_field(data, "type", "character", NA_character_),
      observation_count = stats_count(data, "observation_count"),
      run_count = stats_count(data, "run_count"),
      sample_count = ptmqtl_optional_count(data, "sample_count"),
      run_names = wire_strings(data[["run_names"]]),
      covarying_run_count = ptmqtl_optional_count(data, "covarying_run_count"),
      min_quantified_fraction = ptmqtl_optional_count(data, "min_quantified_fraction"),
      exclude_ceiling = ptmqtl_optional_flag(data, "exclude_ceiling"),
      remove_run_level = ptmqtl_optional_flag(data, "remove_run_level"),
      run_level_group = wire_field(data, "run_level_group", "character", NA_character_),
      run_level_groups_file = wire_field(data, "run_level_groups_file", "character", NA_character_),
      min_runs_for_run_level = ptmqtl_optional_count(data, "min_runs_for_run_level"),
      min_sites_for_run_level = ptmqtl_optional_count(data, "min_sites_for_run_level"),
      modifications = ptmqtl_optional_names(data[["modifications"]]),
      categories = ptmqtl_optional_names(data[["categories"]]),
      ignored_columns = wire_strings(data[["ignored_columns"]]),
      physical_count = ptmqtl_optional_count(data, "physical_count"),
      covarying_count = ptmqtl_optional_count(data, "covarying_count"),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_ptm_pairs"
  )
}

ptmqtl_parse_trait_effects <- function(data) {
  structure(
    list(
      occupancy_file = wire_field(data, "occupancy_file", "character", NA_character_),
      runs_file = wire_field(data, "runs_file", "character", NA_character_),
      trait_name = wire_field(data, "trait_name", "character", NA_character_),
      covariate_names = wire_strings(data[["covariate_names"]]),
      cell_count = stats_count(data, "cell_count"),
      run_count = stats_count(data, "run_count"),
      trait_run_count = stats_count(data, "trait_run_count"),
      occupancy_runs_not_used = wire_strings(data[["occupancy_runs_not_used"]]),
      min_runs = stats_count(data, "min_runs"),
      min_replicates = stats_count(data, "min_replicates"),
      logit_epsilon = stats_count(data, "logit_epsilon"),
      exclude_ceiling = isTRUE(data[["exclude_ceiling"]]),
      status_counts = proteins_parse_counts(data[["status_counts"]]),
      ignored_columns = wire_strings(data[["ignored_columns"]]),
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data, c(covariate_effects = "numeric")),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_trait_effects"
  )
}

# ---------------------------------------------------------------- the public surface

#' Compute every site's occupancy in every run, with mzLib's SiteOccupancyCalculator
#'
#' For a site (protein, position, modification) in a run, occupancy is the summed intensity of the
#' quantified peptidoforms that carry the modification there, over the summed intensity of every
#' quantified peptidoform that covers the position - the estimator MetaMorpheus writes as
#' `IntensityOccupancy_`. Each cell keeps its state, so a censored cell is never mistaken for a 0: a
#' `"floor"` (the modified form was identified but not quantified) is low but not zero, and
#' `"not_detected"` is unknown. Neither has a `fraction`.
#'
#' mzLib always leaves out the categories `Common Variable` and `Common Fixed`, peptide-terminal
#' modifications and every C-terminal modification. An N-terminal modification is a site only on a
#' peptide that starts at residue 1, or at 2 after the initiator Met is removed.
#'
#' @param observations A TSV with a header row naming its columns, one row per (run, peptidoform,
#'   protein) identified by MS/MS: `run`, `full_sequence` (MetaMorpheus notation),
#'   `protein_accession`, `start_residue` and `end_residue` (1-based, in the protein), `intensity`
#'   (instrument units; blank or `NA` for identified but not quantified), and optionally
#'   `previous_residue` (one letter). A shared peptide is one row per protein. Never one row per
#'   PSM, and never a match-between-runs transfer. mzLib has no reader that builds this table; the
#'   ptmqtl article shows how to assemble it.
#' @param samples A TSV with a header row, then `run`, `sample`, to sum each sample's fractions
#'   before the ratio; the `run` column of the result then holds the sample. Map only fractions
#'   together, never technical replicates. Leave it `NULL` for per-run occupancy.
#' @param modifications Count only these modifications as sites: exact names,
#'   `Category:IdWithMotif`, as a character vector. Not with `categories`. `NULL`: no filter.
#' @param categories Count only modifications of these categories (the text before the first `:`),
#'   as a character vector, e.g. `"Common Biological"`. Not with `modifications`.
#' @param out Write the table to this `.tsv` instead of returning it - the file
#'   [ptmqtl_trait_effects()] reads. Then `records` has no rows and `written` says where.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_site_occupancy`:
#'
#'   - `records`, a data.frame with one row per (site, run) where the site was identified or
#'     covered, by protein, then run, then position and modification: `site_key`,
#'     `accession:residue position:modification`, the key [ptmqtl_pairs()] and
#'     [ptmqtl_trait_effects()] use; `protein_accession`; `position`, the 1-based residue (a protein
#'     N-terminal modification at the residue it sits on); `residue`; `modification`,
#'     `Category:IdWithMotif`; `category`, the text before the first `:`; `run`, the run or with
#'     `samples` the sample; `state`, one of `"quantified"`, `"floor"`, `"count_only"` or
#'     `"not_detected"`; `modified_intensity` (intensity, instrument units), the summed intensity of
#'     the quantified forms carrying the modification here, 0 unless quantified; `covering_intensity`
#'     (intensity, instrument units), of every quantified form covering the position, 0 for
#'     `"count_only"`; `fraction` (fraction, 0 to 1), modified over covering intensity, `NA` unless
#'     `state` is `"quantified"` - censored, never 0; and `unmodified_quantified`, where `FALSE`
#'     with `fraction` 1 is a **ceiling**: the site was seen only modified.
#'   - `observations_file` and `samples_file` (`NA` without `samples`), as absolute paths;
#'     `observation_count` rows read; `run_count` runs; `sample_count` samples, `NA` without
#'     `samples`; `run_names`; `site_count` sites; `state_counts`, cells per state;
#'     `modifications` and `categories`, the filter given, `NULL` for none; `ignored_columns`,
#'     columns that were not read (a misspelled `previous_residue` shows up here); `row_count` rows,
#'     returned or written; `column_names`; `written` (`path`, `format`, `row_count`), or `NULL`.
#'   - `caveats`: what the numbers are and are not, plus any that apply to this input.
#'
#' @section Positions:
#'
#' `position` is 1-based, unlike [readers_read_occupancy()], whose N-terminus is 0. A peptide shared
#' by several proteins counts in full for each; MetaMorpheus gives it a proportional share per
#' protein group, so its `IntensityOccupancy_` can differ.
#'
#' @seealso [ptmqtl_pairs()], [ptmqtl_trait_effects()], [readers_read_occupancy()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own hand-computed cases, one protein each:
#' cases <- system.file("extdata", "ptmqtl", "cases.tsv", package = "mzLibR")
#' occ <- ptmqtl_occupancy(cases)
#' occ
#' occ$records[occ$records$protein_accession == "STATES", c("run", "state", "fraction")]
#' # RATIO is 30 / (30 + 70 + 100): a missed-cleavage form covering the site counts below the line.
#' occ$records[occ$records$protein_accession == "RATIO", c("modification", "fraction")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec ptmqtl.occupancy
#' @export
ptmqtl_occupancy <- function(observations, samples = NULL, modifications = NULL, categories = NULL,
                             out = NULL, timeout = 600) {
  args <- ptmqtl_build_occupancy_args(observations, samples, modifications, categories, out)
  bridge_require_verb("ptmqtl occupancy", PTMQTL_SINCE)
  ptmqtl_parse_occupancy(bridge_invoke(args, timeout = timeout))
}

#' Find PTM-PTM pairs within one dataset, with mzLib's PtmPairEngine
#'
#' Two kinds of evidence, which you choose between with `type`:
#'
#' - **physical** (`PtmPairEngine.Physical`): both modifications were identified on one
#'   peptidoform, so on one molecule. `statistic` is the median co-occupancy over the runs where both
#'   were quantified. There is no test.
#' - **covarying** (`PtmPairEngine.CoVarying`): the two sites' occupancies rise and fall together (or
#'   oppositely) across runs. `statistic` is Spearman's rho on the runs where both are quantified,
#'   with a Benjamini-Hochberg `q` per family. Floors, undetected cells and, by default, ceilings are
#'   left out: two unrelated sites both read 1 in low-load runs, which would correlate them through
#'   detection alone.
#'
#' Occupancy is computed inside by the same mzLib call [ptmqtl_occupancy()] makes.
#'
#' @param observations [ptmqtl_occupancy()]'s observations table, with the same rules.
#' @param type `"physical"`, `"covarying"` or `"both"` (physical rows first). Required: mzLib has no
#'   default between its two methods.
#' @param samples [ptmqtl_occupancy()]'s run-to-sample table; pairs are then across samples.
#' @param modifications As [ptmqtl_occupancy()]: only these modifications are sites.
#' @param categories As [ptmqtl_occupancy()]: only these categories are sites.
#' @param min_quantified_fraction Co-varying only. A site enters pairs only if quantified in at least
#'   this fraction of the runs (fraction of the runs, above 0 and at most 1), rounded up. `NULL`:
#'   mzLib's default (0.7 at mzLib `d10017ce`), which `min_quantified_fraction` in the result
#'   reports.
#' @param keep_ceiling Co-varying only. Keep runs where a site was seen only modified. `FALSE` sends
#'   nothing, so mzLib's default (leave them out) applies.
#' @param remove_run_level Co-varying only. Correlate given each run's level - the median occupancy
#'   of the other tested sites of the pair's group - to remove a batch or loading shift. **Batch
#'   only**: where the level follows the design, it removes the effect you study.
#' @param run_level_group With `remove_run_level`: `"category"` or `"modification"`, for one level
#'   per group of sites. `NULL`: one level for every site, mzLib's default. Not with
#'   `run_level_groups`.
#' @param run_level_groups With `remove_run_level`: a TSV with a header row, then `modification`,
#'   `group`, for any grouping, such as chemistry. Every modification among the sites needs a group.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_ptm_pairs`:
#'
#'   - `records`, a data.frame with one row per pair of sites, physical pairs first: `result_type`,
#'     `"physical"` (both on one peptidoform) or `"covarying"` (across runs); `site_a` and `site_b`,
#'     the two sites' keys, `site_a` first by ordinal order; `same_protein`; `overlapping`, one
#'     peptide covers both (such a co-varying pair is correlated by construction); `statistic` (physical:
#'     fraction (0 to 1), the median co-occupancy; covarying: unitless correlation (-1 to 1),
#'     Spearman's rho), `NA` when physical pairs were never quantified together, or not estimable;
#'     `p_value`, a fraction (0 to 1), two-sided Spearman, `NA` for physical (no test) or not
#'     estimable; `q`, a fraction (0 to 1), Benjamini-Hochberg within `fdr_family`, `NA` when not
#'     tested (physical, overlapping, not estimable), not a target-decoy q-value; `fdr_family`, `P`
#'     or `A` then `intra` or `inter` protein, e.g. `"A:inter"`; `n` (runs), physical: runs
#'     identified together, covarying: runs where both were quantified; `statistic_n` (runs), the
#'     runs `statistic` was computed on; `spearman_method`, `"exact"`, `"asymptotic"` or
#'     `"not_estimable"`, `NA` for physical; and `run_level_removed`.
#'   - `observations_file`, `samples_file` and `run_level_groups_file` (`NA` when not given), as
#'     absolute paths; `type`; `observation_count` rows read; `run_count` runs; `sample_count`
#'     samples, `NA` without `samples`; `run_names`; and, `NA` for `type = "physical"`:
#'     `covarying_run_count` (runs), the runs where any site was identified or covered, the
#'     denominator of `min_quantified_fraction` (fraction of the runs, 0 to 1), the guard used;
#'     `exclude_ceiling` and `remove_run_level`; `run_level_group`, `"none"`, `"category"`,
#'     `"modification"` or `"file"`, `NA` when the run level was not removed;
#'     `min_runs_for_run_level` (runs) and `min_sites_for_run_level` (sites), mzLib's thresholds,
#'     `NA` when not removed; `modifications` and `categories`, `NULL` for no filter;
#'     `ignored_columns`; `physical_count` pairs (`NA` for `"covarying"`) and `covarying_count`
#'     pairs (`NA` for `"physical"`); `row_count` rows; `column_names`.
#'   - `caveats`: what each pair type means, and the run level when it was removed, plus any that
#'     apply to this input.
#'
#' @seealso [ptmqtl_occupancy()], [stats_correlate()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # Two modifications on one molecule, mzLib's PhysicalPairsCarryCoOccupancy case:
#' cases <- system.file("extdata", "ptmqtl", "cases.tsv", package = "mzLibR")
#' found <- ptmqtl_pairs(cases, type = "physical")
#' found$records[, c("site_a", "statistic", "n", "statistic_n", "p_value")]
#'
#' # Twelve runs in two batches, every site 0.2 higher in the second. The batch alone correlates
#' # two sites; given each run's level, they are unrelated:
#' batch <- system.file("extdata", "ptmqtl", "batch.tsv", package = "mzLibR")
#' plain <- ptmqtl_pairs(batch, type = "covarying")
#' given <- ptmqtl_pairs(batch, type = "covarying", remove_run_level = TRUE)
#' pick <- function(r) r$statistic[startsWith(r$site_a, "PAIR_S") & startsWith(r$site_b, "PAIR_T")]
#' round(c(plain = pick(plain$records), given = pick(given$records)), 2)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec ptmqtl.pairs
#' @export
ptmqtl_pairs <- function(observations, type, samples = NULL, modifications = NULL, categories = NULL,
                         min_quantified_fraction = NULL, keep_ceiling = FALSE, remove_run_level = FALSE,
                         run_level_group = NULL, run_level_groups = NULL, timeout = 600) {
  if (missing(type)) {
    stop(mzlib_usage_error(paste0(
      "type is required: one of ", paste(PTMQTL_PAIR_TYPES, collapse = ", "),
      ". mzLib has no default between its two methods."
    )))
  }
  args <- ptmqtl_build_pairs_args(observations, type, samples, modifications, categories,
    min_quantified_fraction, keep_ceiling, remove_run_level, run_level_group, run_level_groups
  )
  bridge_require_verb("ptmqtl pairs", PTMQTL_SINCE)
  ptmqtl_parse_pairs(bridge_invoke(args, timeout = timeout))
}

#' Fit, per site, the effect of a per-run trait on logit occupancy, with mzLib's SiteTraitEffects
#'
#' The model, per site, is logit(occupancy) = b0 + b_trait x trait + covariates + a random intercept
#' per biological replicate + error, fitted by REML (mzLib's `MixedModel`). The random intercept is
#' the point: repeated runs of one animal or donor are not counted as independent, so a trait that is
#' constant within each replicate is tested at the between-replicate degrees of freedom.
#'
#' Only quantified cells are responses. Floors, count-only and undetected cells are left out, never
#' entered as 0, so this answers "does the modified share change?", not "is the modification
#' detected more often?". Ceilings are left out by default.
#'
#' @param occupancy The table [ptmqtl_occupancy()] writes with `out` (columns found by header name).
#'   A `fraction` cell is used as written, as mzLib does for a stored fraction.
#' @param runs A TSV with a header row, then one row per run: the run, its biological replicate, the
#'   trait (its header names it; any unit), then any covariates (each header names one). A run with
#'   a blank trait is not used.
#' @param min_runs Fewest runs with a response for a site to be fitted (runs). `NULL`: mzLib's
#'   default (6 at mzLib `d10017ce`), which `min_runs` in the result reports.
#' @param min_replicates Fewest distinct replicates among them (replicates). `NULL`: mzLib's default
#'   (4 at mzLib `d10017ce`).
#' @param logit_epsilon Fractions are clamped to `[epsilon, 1 - epsilon]` before the logit (fraction,
#'   above 0 and below 0.5). `NULL`: mzLib's default (5e-5 at mzLib `d10017ce`).
#' @param keep_ceiling Use runs where the site was seen only modified. `FALSE` sends nothing, so
#'   mzLib's default (leave them out) applies.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_trait_effects`:
#'
#'   - `records`, a data.frame with one row per site with a usable response, by site key:
#'     `site_key`, `protein_accession`, `position`, `residue` and `modification`, as
#'     [ptmqtl_occupancy()] writes them; `status`, one of `"fitted"`, `"below_support"`,
#'     `"too_few_observations"`, `"rank_deficient"`, `"too_few_groups"`, `"variance_ratio_at_limit"`
#'     or `"not_converged"`; `effect` (logit occupancy per trait unit), the change in
#'     logit(occupancy) per unit of the trait, so `exp(effect)` is the odds ratio, `NA` when not
#'     fitted; `standard_error` (logit occupancy per trait unit); `df` (degrees of freedom), by the
#'     containment rule; `p_value`, a fraction (0 to 1), two-sided, for the trait; `q`, a fraction
#'     (0 to 1), Benjamini-Hochberg over this call's fitted sites; `runs`, runs contributing a
#'     response; `replicates`, distinct replicates among them; `median_fraction` (fraction, 0 to 1),
#'     median occupancy over those runs; `covariate_effects` (logit occupancy per covariate unit), a
#'     list column, one per covariate in `covariate_names` order, each `NA` when not fitted;
#'     `replicate_variance` (squared logit occupancy), between replicates; and `residual_variance`
#'     (squared logit occupancy), run to run within a replicate. Statistics are `NA` where `status`
#'     is not `"fitted"`.
#'   - `occupancy_file` and `runs_file`, as absolute paths; `trait_name`, the runs table's third
#'     header; `covariate_names`; `cell_count` occupancy rows (cells) read; `run_count` runs in the
#'     runs table; `trait_run_count`, runs with a trait value (runs); `occupancy_runs_not_used`;
#'     `min_runs` (runs) and `min_replicates` (replicates), the thresholds used, yours or mzLib's
#'     default; `logit_epsilon` (fraction, 0 to 1), the clamp used; `exclude_ceiling`;
#'     `status_counts`, sites (sites) per status; `ignored_columns`; `row_count` rows;
#'     `column_names`.
#'   - `caveats`: the model and what it leaves out, plus any that apply to this input.
#'
#' @section What it is not:
#'
#' It says nothing about whether *detecting* the modification changes with the trait: only
#' quantified cells enter. `q` adjusts over this call's fitted sites only - not across datasets, and
#' not a target-decoy q-value.
#'
#' @seealso [ptmqtl_occupancy()], [stats_fit_mixed()] for the same model on any responses
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's TraitDesign: six replicates at three doses, two injections each, with
#' # logit(occupancy) = -2 + 0.5 x dose + a replicate offset. trait_runs.tsv gives each run's
#' # replicate, dose and an injection covariate; trait_occupancy.tsv is ptmqtl_occupancy(out = )'s.
#' data <- system.file("extdata", "ptmqtl", package = "mzLibR")
#' fit <- ptmqtl_trait_effects(file.path(data, "trait_occupancy.tsv"),
#'   file.path(data, "trait_runs.tsv"))
#' fit
#' fit$records[, c("protein_accession", "status", "effect", "df", "runs", "replicates")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec ptmqtl.trait-effects
#' @export
ptmqtl_trait_effects <- function(occupancy, runs, min_runs = NULL, min_replicates = NULL,
                                 logit_epsilon = NULL, keep_ceiling = FALSE, timeout = 600) {
  args <- ptmqtl_build_trait_args(occupancy, runs, min_runs, min_replicates, logit_epsilon, keep_ceiling)
  bridge_require_verb("ptmqtl trait-effects", PTMQTL_SINCE)
  ptmqtl_parse_trait_effects(bridge_invoke(args, timeout = timeout))
}

# ---------------------------------------------------------------- printing

#' Print site occupancy
#'
#' @param x A [ptmqtl_occupancy()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_site_occupancy <- function(x, ...) {
  cat("<mzlibr_site_occupancy> ", format(x$site_count), " sites in ", format(length(x$run_names)),
    if (is.na(x$sample_count)) " runs, " else " samples, ", format(x$row_count), " cells",
    if (!is.null(x$written)) paste0(", written to ", x$written$path) else "", "\n",
    sep = ""
  )
  if (length(x$state_counts) > 0L) {
    cat("  ", paste0(names(x$state_counts), " ", x$state_counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Print PTM pairs
#'
#' @param x A [ptmqtl_pairs()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_ptm_pairs <- function(x, ...) {
  cat("<mzlibr_ptm_pairs> ", x$type, ": ", format(x$row_count), " pair(s) over ",
    format(length(x$run_names)), " runs",
    if (isTRUE(x$remove_run_level)) paste0(", given the run level (", x$run_level_group, ")") else "", "\n",
    sep = ""
  )
  invisible(x)
}

#' Print trait effects
#'
#' @param x A [ptmqtl_trait_effects()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_trait_effects <- function(x, ...) {
  cat("<mzlibr_trait_effects> ", format(x$row_count), " sites against ", x$trait_name,
    if (length(x$covariate_names) > 0L) paste0(", given ", paste(x$covariate_names, collapse = ", ")) else "",
    " (min ", format(x$min_runs), " runs, ", format(x$min_replicates), " replicates)\n",
    sep = ""
  )
  if (length(x$status_counts) > 0L) {
    cat("  ", paste0(names(x$status_counts), " ", x$status_counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
