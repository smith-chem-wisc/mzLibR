# ptmqtl_occupancy(), ptmqtl_pairs() and ptmqtl_trait_effects().
#
# Ported from pyMzLib's test_ptmqtl.py and test_ptmqtl_live.py. The payloads are pyMzLib's
# recordings, verbatim; the inputs are mzLib's hand-computed PtmQtl cases, shipped in
# extdata/ptmqtl/. Three layers: what is sent and parsed, and refusals before the bridge; replayed
# end to end; and LIVE, the bridge today returns what was recorded.

ptmqtl_path <- function(name) system.file("extdata", "ptmqtl", name, package = "mzLibR")

# ---------------------------------------------------------------- occupancy

test_that("occupancy sends only what was given", {
  expect_identical(mz$ptmqtl_build_occupancy_args("cases.tsv", NULL, NULL, NULL, NULL),
    c("ptmqtl", "occupancy", "--observations", "cases.tsv"))
  expect_identical(mz$ptmqtl_build_occupancy_args("cases.tsv", "s.tsv", NULL, c("Common Biological", "UniProt"), "o.tsv"),
    c("ptmqtl", "occupancy", "--observations", "cases.tsv", "--samples", "s.tsv",
      "--categories", "Common Biological\tUniProt", "--out", "o.tsv"))
  expect_identical(utils::tail(mz$ptmqtl_build_occupancy_args("cases.tsv", NULL, "Common Biological:Phosphorylation on S", NULL, NULL), 2L),
    c("--modifications", "Common Biological:Phosphorylation on S"))
})

test_that("the cases recording parses; a censored cell has no fraction, never a 0", {
  occ <- mz$ptmqtl_parse_occupancy(recorded_data("ptmqtl_occupancy_cases.json"))
  expect_true(inherits(occ, "mzlibr_site_occupancy"))
  expect_identical(occ$column_names, c("site_key", "protein_accession", "position", "residue", "modification", "category",
    "run", "state", "modified_intensity", "covering_intensity", "fraction", "unmodified_quantified"))
  expect_identical(sort(names(occ$state_counts)), sort(mz$PTMQTL_OCCUPANCY_STATES))
  expect_identical(c(occ$row_count, nrow(occ$records), sum(occ$state_counts)), c(14, 14, 14))
  expect_identical(c(occ$site_count, occ$state_counts[["quantified"]], occ$state_counts[["count_only"]]), c(9, 9, 3))
  expect_true(is.na(occ$samples_file) && is.na(occ$sample_count) && is.null(occ$written))
  expect_true(is.null(occ$modifications) && is.null(occ$categories))
  expect_identical(occ$ignored_columns, character(0))
  r <- occ$records
  ratio <- r[r$site_key == "RATIO:S13:Common Biological:Phosphorylation on S" & r$run == "r1", ]
  expect_true(abs(ratio$fraction - 30 / 200) < 1e-15)
  expect_identical(c(ratio$modified_intensity, ratio$covering_intensity), c(30, 200))
  expect_identical(is.na(r$fraction), r$state != "quantified")
  expect_false(any(r$protein_accession %in% c("EXCLUDED", "PEPTIDE_NTERM", "NTERM_AFTER_K")), info = "mzLib's exclusions")
  expect_identical(r$position[r$site_key == "NTERM_AFTER_MET:S2:UniProt:N-acetylserine on S" & r$run == "r1"], 2)
  states <- r[r$protein_accession == "STATES", ]
  expect_identical(states$state, c("count_only", "floor", "not_detected", "quantified"))
  expect_identical(states$fraction, c(NA, NA, NA, 0.1))
  ceiling <- r[r$protein_accession == "CEILING", ]
  expect_identical(c(ceiling$fraction, ceiling$unmodified_quantified), c(1, 0))
  expect_true(any(startsWith(occ$caveats, "position is 1-based")))
})

test_that("with samples the fractions are summed before the ratio", {
  occ <- mz$ptmqtl_parse_occupancy(recorded_data("ptmqtl_occupancy_fractions.json"))
  expect_identical(occ$run_names, c("A", "B"))
  expect_identical(c(occ$run_count, occ$sample_count), c(3, 2))
  expect_true(endsWith(occ$samples_file, "fractions_samples.tsv"))
  expect_true(any(startsWith(occ$caveats, "With --samples")))
  expect_identical(occ$records$state, c("floor", "quantified", "not_detected"))
  expect_identical(occ$records$fraction, c(NA, 0.2, NA))
})

test_that("occupancy written to out has no rows, and says where it went", {
  occ <- mz$ptmqtl_parse_occupancy(recorded_data("ptmqtl_occupancy_trait_out.json"))
  expect_identical(nrow(occ$records), 0L)
  expect_identical(names(occ$records), occ$column_names)
  expect_identical(c(occ$written$format, occ$written$row_count, occ$row_count), c("tsv", "24", "24"))
  expect_true(endsWith(occ$written$path, "trait_occupancy.tsv"))
})

test_that("bad occupancy arguments are refused before the bridge", {
  build <- function(...) {
    given <- list(observations = "cases.tsv", samples = NULL, modifications = NULL, categories = NULL, out = NULL)
    given[names(list(...))] <- list(...)
    do.call(mz$ptmqtl_build_occupancy_args, given)
  }
  expect_error(build(modifications = "a", categories = "b"), class = "mzlib_usage_error", contains = "not both")
  expect_error(build(categories = character(0)), class = "mzlib_usage_error", contains = "non-empty character vector")
  expect_error(build(categories = c("ok", "")), class = "mzlib_usage_error", contains = "non-empty string")
  expect_error(build(modifications = "a\tb"), class = "mzlib_usage_error", contains = "a tab or a line break")
  expect_error(build(samples = 3), class = "mzlib_usage_error", contains = "samples must be")
  expect_error(build(out = 1), class = "mzlib_usage_error", contains = "out must be")
})

# ---------------------------------------------------------------- pairs

test_that("pairs sends only what was given", {
  build <- function(...) {
    given <- list(observations = "batch.tsv", type = "covarying", samples = NULL, modifications = NULL, categories = NULL,
      min_quantified_fraction = NULL, keep_ceiling = FALSE, remove_run_level = FALSE, run_level_group = NULL,
      run_level_groups = NULL)
    given[names(list(...))] <- list(...)
    do.call(mz$ptmqtl_build_pairs_args, given)
  }
  expect_identical(build(), c("ptmqtl", "pairs", "--observations", "batch.tsv", "--type", "covarying"))
  every <- build(type = "both", samples = "s.tsv", min_quantified_fraction = 0.5, keep_ceiling = TRUE,
    remove_run_level = TRUE, run_level_group = "category")
  expect_identical(every, c("ptmqtl", "pairs", "--observations", "batch.tsv", "--type", "both", "--samples", "s.tsv",
    "--min-quantified-fraction", "0.5", "--keep-ceiling", "--remove-run-level", "--run-level-group", "category"))
  expect_identical(utils::tail(build(remove_run_level = TRUE, run_level_groups = "groups.tsv"), 2L),
    c("--run-level-groups", "groups.tsv"))
  expect_error(build(type = "all"), class = "mzlib_usage_error", contains = "type must be one of")
  expect_error(build(type = NULL), class = "mzlib_usage_error", contains = "type must be one of")
  expect_error(ptmqtl_pairs("batch.tsv"), class = "mzlib_usage_error", contains = "type is required")
  expect_error(build(min_quantified_fraction = "0.7"), class = "mzlib_usage_error", contains = "must be a number")
  expect_error(build(min_quantified_fraction = TRUE), class = "mzlib_usage_error", contains = "must be a number")
  expect_error(build(keep_ceiling = "yes"), class = "mzlib_usage_error", contains = "keep_ceiling must be TRUE or FALSE")
  expect_error(build(remove_run_level = 1), class = "mzlib_usage_error", contains = "remove_run_level must be TRUE or FALSE")
  expect_error(build(remove_run_level = TRUE, run_level_group = "chemistry"), class = "mzlib_usage_error",
    contains = "run_level_group must be one of")
  expect_error(build(categories = "a", modifications = "b"), class = "mzlib_usage_error", contains = "not both")
})

test_that("the physical recording parses; physical pairs carry no test and echo no co-varying settings", {
  found <- mz$ptmqtl_parse_pairs(recorded_data("ptmqtl_pairs_physical.json"))
  expect_true(inherits(found, "mzlibr_ptm_pairs"))
  expect_identical(c(found$type, found$physical_count), c("physical", "1"))
  expect_true(is.na(found$covarying_count))
  expect_true(is.na(found$min_quantified_fraction) && is.na(found$exclude_ceiling))
  expect_true(is.na(found$remove_run_level) && is.na(found$covarying_run_count))
  pair <- found$records
  expect_identical(nrow(pair), 1L)
  expect_identical(pair$site_a, "TWO_ON_ONE:K14:Common Biological:Acetylation on K")
  expect_identical(c(pair$result_type, pair$fdr_family), c("physical", "P:intra"))
  expect_true(pair$overlapping)
  expect_true(abs(pair$statistic - 0.25) < 1e-15)
  expect_identical(c(pair$n, pair$statistic_n), c(2, 1))
  expect_true(is.na(pair$p_value) && is.na(pair$q) && is.na(pair$spearman_method))
})

test_that("the co-varying recordings report mzLib's defaults, and the run level when removed", {
  plain <- mz$ptmqtl_parse_pairs(recorded_data("ptmqtl_pairs_batch.json"))
  expect_identical(plain$min_quantified_fraction, 0.7)
  expect_true(plain$exclude_ceiling)
  expect_false(plain$remove_run_level)
  expect_true(is.na(plain$run_level_group) && is.na(plain$physical_count) && is.na(plain$min_runs_for_run_level))
  expect_identical(plain$covarying_run_count, 12)
  given <- mz$ptmqtl_parse_pairs(recorded_data("ptmqtl_pairs_batch_run_level.json"))
  expect_true(given$remove_run_level)
  expect_identical(given$run_level_group, "none")
  expect_identical(c(given$min_runs_for_run_level, given$min_sites_for_run_level), c(10, 3))
  expect_true(all(given$records$run_level_removed))
  not_estimable <- given$records[given$records$spearman_method == "not_estimable", ]
  expect_true(nrow(not_estimable) > 0L && all(is.na(not_estimable$statistic) & is.na(not_estimable$q)))
  expect_true(any(startsWith(given$caveats, "With --remove-run-level")))
  pick <- function(r) r$statistic[startsWith(r$site_a, "PAIR_S") & startsWith(r$site_b, "PAIR_T")]
  expect_identical(sprintf("%.2f", c(pick(plain$records), pick(given$records))), c("0.76", "0.03"))
})

# ---------------------------------------------------------------- trait_effects

test_that("trait_effects sends only what was given", {
  expect_identical(mz$ptmqtl_build_trait_args("occ.tsv", "runs.tsv", NULL, NULL, NULL, FALSE),
    c("ptmqtl", "trait-effects", "--occupancy", "occ.tsv", "--runs", "runs.tsv"))
  every <- mz$ptmqtl_build_trait_args("occ.tsv", "runs.tsv", 8, 5, 1e-6, TRUE)
  expect_identical(every[7:11], c("--min-runs", "8", "--min-replicates", "5", "--logit-epsilon"))
  # The wire text itself, not as.numeric() of it: R reads 17-digit text through long double, which is only
  # double on Apple silicon, so the round trip lands a bit off there while mzLib (C#) parses it exactly.
  expect_identical(every[[12L]], "9.9999999999999995e-07")
  expect_identical(every[[13L]], "--keep-ceiling")
  expect_error(mz$ptmqtl_build_trait_args("occ.tsv", "runs.tsv", 6.5, NULL, NULL, FALSE), class = "mzlib_usage_error",
    contains = "min_runs must be a whole number")
  expect_error(mz$ptmqtl_build_trait_args("occ.tsv", "runs.tsv", NULL, TRUE, NULL, FALSE), class = "mzlib_usage_error",
    contains = "min_replicates must be a whole number")
  expect_error(mz$ptmqtl_build_trait_args("occ.tsv", "runs.tsv", NULL, NULL, "small", FALSE), class = "mzlib_usage_error",
    contains = "logit_epsilon must be a number")
  expect_error(mz$ptmqtl_build_trait_args("occ.tsv", "runs.tsv", NULL, NULL, NULL, NA), class = "mzlib_usage_error",
    contains = "keep_ceiling must be TRUE or FALSE")
  expect_error(ptmqtl_trait_effects("occ.tsv", NULL), class = "mzlib_usage_error", contains = "runs must be")
})

test_that("the trait recording recovers the dose effect, with mzLib's support thresholds", {
  fit <- mz$ptmqtl_parse_trait_effects(recorded_data("ptmqtl_trait_effects.json"))
  expect_true(inherits(fit, "mzlibr_trait_effects"))
  expect_identical(c(fit$min_runs, fit$min_replicates, fit$logit_epsilon), c(6, 4, 5e-05))
  expect_true(fit$exclude_ceiling)
  expect_identical(names(fit$status_counts), mz$PTMQTL_TRAIT_STATUSES)
  expect_identical(c(fit$trait_name, fit$covariate_names), c("dose", "injection"))
  expect_identical(c(fit$trait_run_count, fit$run_count), c(12, 12))
  expect_identical(fit$occupancy_runs_not_used, character(0))
  trait <- fit$records[fit$records$protein_accession == "TRAIT_SITE", ]
  expect_identical(trait$status, "fitted")
  expect_true(abs(trait$effect - 0.5) < 0.1)
  expect_identical(c(trait$df, trait$runs, trait$replicates), c(4, 12, 6))
  expect_true(is.list(fit$records$covariate_effects))
  expect_identical(length(trait$covariate_effects[[1L]]), 1L)
  expect_true(abs(trait$covariate_effects[[1L]] - 0.02) < 0.005)
  censored <- fit$records[fit$records$protein_accession == "CENSORED_SITE", ]
  expect_identical(c(sprintf("%.2f", censored$effect), censored$df, censored$runs, censored$replicates), c("0.54", "3", "10", "5"))
})

test_that("a site below support has NA statistics, and NA covariate effects", {
  fit <- mz$ptmqtl_parse_trait_effects(recorded_data("ptmqtl_trait_effects_strict.json"))
  expect_identical(fit$min_replicates, 6)
  expect_identical(fit$status_counts[["below_support"]], 1)
  censored <- fit$records[fit$records$protein_accession == "CENSORED_SITE", ]
  expect_identical(censored$status, "below_support")
  expect_true(is.na(censored$effect) && is.na(censored$p_value) && is.na(censored$q))
  expect_identical(censored$covariate_effects[[1L]], NA_real_)
  expect_true(!is.na(censored$median_fraction))
  expect_identical(censored$replicates, 5)
})

# ---------------------------------------------------------------- end to end, replayed

test_that("the three verbs replay end to end, a flag or option choosing the recording", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(ptmqtl_occupancy(ptmqtl_path("cases.tsv"))$site_count, 9)
  expect_identical(ptmqtl_occupancy(ptmqtl_path("fractions.tsv"), samples = ptmqtl_path("fractions_samples.tsv"))$sample_count, 2)
  written <- ptmqtl_occupancy(ptmqtl_path("trait.tsv"), out = file.path(tempdir(), "trait_occupancy.tsv"))
  expect_identical(written$written$row_count, 24)
  expect_identical(ptmqtl_pairs(ptmqtl_path("cases.tsv"), type = "physical")$physical_count, 1)
  expect_false(ptmqtl_pairs(ptmqtl_path("batch.tsv"), type = "covarying")$remove_run_level)
  expect_true(ptmqtl_pairs(ptmqtl_path("batch.tsv"), type = "covarying", remove_run_level = TRUE)$remove_run_level)
  default <- ptmqtl_trait_effects(ptmqtl_path("trait_occupancy.tsv"), ptmqtl_path("trait_runs.tsv"))
  strict <- ptmqtl_trait_effects(ptmqtl_path("trait_occupancy.tsv"), ptmqtl_path("trait_runs.tsv"), min_replicates = 6)
  expect_identical(c(default$min_replicates, strict$min_replicates), c(4, 6))
})

# ---------------------------------------------------------------- against a real mzLib

ptmqtl_live <- function() {
  skip_unless_bridge_has("ptmqtl occupancy")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
}

test_that("LIVE: the ptmqtl recordings are what the bridge returns today", {
  ptmqtl_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  calls <- list(
    ptmqtl_occupancy_cases.json = function() ptmqtl_occupancy(ptmqtl_path("cases.tsv")),
    ptmqtl_occupancy_fractions.json = function() ptmqtl_occupancy(ptmqtl_path("fractions.tsv"), samples = ptmqtl_path("fractions_samples.tsv")),
    ptmqtl_pairs_physical.json = function() ptmqtl_pairs(ptmqtl_path("cases.tsv"), type = "physical"),
    ptmqtl_pairs_batch.json = function() ptmqtl_pairs(ptmqtl_path("batch.tsv"), type = "covarying"),
    ptmqtl_pairs_batch_run_level.json = function() ptmqtl_pairs(ptmqtl_path("batch.tsv"), type = "covarying", remove_run_level = TRUE),
    ptmqtl_trait_effects.json = function() ptmqtl_trait_effects(ptmqtl_path("trait_occupancy.tsv"), ptmqtl_path("trait_runs.tsv")),
    ptmqtl_trait_effects_strict.json = function() {
      ptmqtl_trait_effects(ptmqtl_path("trait_occupancy.tsv"), ptmqtl_path("trait_runs.tsv"), min_replicates = 6)
    }
  )
  parsers <- list(ptmqtl_occupancy = mz$ptmqtl_parse_occupancy, ptmqtl_pairs = mz$ptmqtl_parse_pairs,
    ptmqtl_trait = mz$ptmqtl_parse_trait_effects)
  for (name in names(calls)) {
    live <- calls[[name]]()
    recorded <- parsers[[which(startsWith(name, names(parsers)))]](recorded_data(name))
    expect_identical(records_disagree(live$records, recorded$records), "", info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
    for (counts in c("state_counts", "status_counts")) {
      # Keys from mzLib enums; a repin that grows one must not slip past the recording.
      if (!is.null(recorded[[counts]])) expect_identical(live[[counts]], recorded[[counts]], info = name)
    }
  }
})

test_that("LIVE: occupancy(out = ) writes the table trait_effects() reads", {
  ptmqtl_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  out <- file.path(tempfile("nested"), "occupancy.tsv")
  occ <- ptmqtl_occupancy(ptmqtl_path("trait.tsv"), out = out)
  expect_identical(nrow(occ$records), 0L)
  expect_identical(occ$written$row_count, occ$row_count)
  # Compared line by line: a Windows checkout may carry the shipped copy with CRLF endings.
  expect_identical(readLines(out, warn = FALSE), readLines(ptmqtl_path("trait_occupancy.tsv"), warn = FALSE))
  fit <- ptmqtl_trait_effects(out, ptmqtl_path("trait_runs.tsv"))
  recorded <- mz$ptmqtl_parse_trait_effects(recorded_data("ptmqtl_trait_effects.json"))
  expect_identical(records_disagree(fit$records, recorded$records), "")
})

test_that("LIVE: the filters reach mzLib", {
  ptmqtl_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  every <- ptmqtl_occupancy(ptmqtl_path("cases.tsv"))
  biological <- ptmqtl_occupancy(ptmqtl_path("cases.tsv"), categories = "Common Biological")
  expect_identical(unique(biological$records$category), "Common Biological")
  expect_identical(biological$categories, "Common Biological")
  expect_true(biological$row_count < every$row_count)
  phospho <- ptmqtl_occupancy(ptmqtl_path("cases.tsv"), modifications = "Common Biological:Phosphorylation on S")
  expect_identical(unique(phospho$records$modification), "Common Biological:Phosphorylation on S")
})

test_that("LIVE: refusals cross as usage errors", {
  ptmqtl_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  expect_error(ptmqtl_pairs(ptmqtl_path("batch.tsv"), type = "physical", remove_run_level = TRUE),
    class = "mzlib_usage_error", contains = "shapes co-varying pairs only")
  expect_error(ptmqtl_pairs(ptmqtl_path("batch.tsv"), type = "covarying", run_level_group = "category"),
    class = "mzlib_usage_error", contains = "give --remove-run-level as well")
  condition <- expect_error(ptmqtl_occupancy(ptmqtl_path("cases.tsv"), out = tempfile(fileext = ".csv")),
    class = "mzlib_usage_error")
  expect_true(grepl("must end in|must name a .tsv", conditionMessage(condition)))
  lines <- readLines(ptmqtl_path("cases.tsv"), warn = FALSE)
  repeated <- tempfile(fileext = ".tsv")
  writeLines(c(lines[1L], lines[2L], lines[2L]), repeated)
  expect_error(ptmqtl_occupancy(repeated), class = "mzlib_usage_error", contains = "not one per PSM")
  expect_error(ptmqtl_trait_effects(ptmqtl_path("trait_occupancy.tsv"), ptmqtl_path("trait_runs.tsv"), logit_epsilon = 0.5),
    class = "mzlib_usage_error", contains = "above 0 and below 0.5")
})
