# diff_design() and diff_observations(), with factor_spec() and covariate_spec().
#
# Ported from pyMzLib's test_diff.py and test_diff_live.py. The payloads are pyMzLib's recordings,
# verbatim; the sample sheets and peptide tables are mzLib's, shipped in extdata/diff/. Three layers:
#   1. always: the stdin lines each call writes (and nothing for what was not given), what is
#      parsed, refusals before the bridge, and the recorded matrix against model.matrix() run here;
#   2. replayed end to end through the stand-in bridge;
#   3. LIVE: the bridge today returns what was recorded; the design file is the model stats_fit()
#      fits; the responses file is what stats_fit() reads.

diff_path <- function(name) system.file("extdata", "diff", name, package = "mzLibR")
diff_recorded <- function(name) mz$diff_parse_design(recorded_data(name))
observations_recorded <- function(name) mz$diff_parse_observations(recorded_data(name))

# ---------------------------------------------------------------- design: what is sent

test_that("design writes one stdin line per term and leaves out what was not given", {
  request <- mz$diff_build_design_request("sheet.tsv",
    factors = list(factor_spec("condition", "ctrl"), "sex"),
    covariates = list(covariate_spec("age", "years", scale = 10, scaled_unit = "decade", centre = 50),
      covariate_spec("bmi", "kg/m2"), covariate_spec("dose", "mg", centre = 0)),
    interactions = list(c("condition", "sex")), strata = "tissue", batch = "plex", out = NULL)
  expect_identical(request$args, c("diff", "design", "--samples", "sheet.tsv"))
  expect_identical(strsplit(request$stdin, "\n", fixed = TRUE)[[1L]], c(
    "factor\tcondition\tctrl",
    "factor\tsex",
    "covariate\tage\tyears\t10\tdecade\t50",
    "covariate\tbmi\tkg/m2",
    "covariate\tdose\tmg\t\t\t0",
    "interaction\tcondition\tsex",
    "stratum\ttissue",
    "batch\tplex"
  ))
})

test_that("design sends out only when given, and no model sends nothing", {
  with_out <- mz$diff_build_design_request("s.tsv", "malat_load", NULL, NULL, NULL, NULL, "dir")
  expect_identical(utils::tail(with_out$args, 2L), c("--out", "dir"))
  bare <- mz$diff_build_design_request("s.tsv", NULL, NULL, NULL, NULL, NULL, NULL)
  expect_false("--out" %in% bare$args)
  expect_true(is.null(bare$stdin), info = "no model: nothing on stdin, mzLib's empty DesignSpec")
  # A single spec is one term, not a list of its fields.
  one <- mz$diff_build_design_request("s.tsv", factor_spec("malat_load", "500ng"), NULL, NULL, NULL, NULL, NULL)
  expect_identical(one$stdin, "factor\tmalat_load\t500ng")
})

test_that("a malformed model is refused before the bridge", {
  build <- function(...) {
    given <- list(samples = "s.tsv", factors = NULL, covariates = NULL, interactions = NULL, strata = NULL, batch = NULL, out = NULL)
    given[names(list(...))] <- list(...)
    do.call(mz$diff_build_design_request, given)
  }
  expect_error(build(factors = list(3)), class = "mzlib_usage_error", contains = "factor_spec() or a column name")
  expect_error(factor_spec("a\tb"), class = "mzlib_usage_error", contains = "tab or a line break")
  expect_error(factor_spec(""), class = "mzlib_usage_error", contains = "non-empty string")
  expect_error(build(covariates = list(c("age", "years"))), class = "mzlib_usage_error", contains = "covariate_spec()")
  expect_error(covariate_spec("age", "years", scale = Inf), class = "mzlib_usage_error", contains = "finite number")
  expect_error(covariate_spec("age", "years", centre = TRUE), class = "mzlib_usage_error", contains = "finite number")
  expect_error(build(interactions = list(c("a", "b", "c"))), class = "mzlib_usage_error", contains = "pair of factor names")
  expect_error(build(interactions = c("a", "b")), class = "mzlib_usage_error", contains = "list of pairs")
  expect_error(build(batch = ""), class = "mzlib_usage_error", contains = "non-empty string")
  expect_error(diff_design(NULL), class = "mzlib_usage_error", contains = "samples must be")
})

test_that("the specs print what mzLib will be told", {
  expect_identical(factor_spec("condition")$reference_level, NA_character_)
  spec <- covariate_spec("bmi", "kg/m2")
  expect_true(is.na(spec$scale) && is.na(spec$scaled_unit) && is.na(spec$centre), info = "not given: mzLib's default")
  expect_true(grepl("every pair of levels", paste(utils::capture.output(print(factor_spec("condition"))), collapse = "")))
})

# ---------------------------------------------------------------- design: what is parsed

test_that("the MALAT1 design recording parses, its stratum keyed by name", {
  design <- diff_recorded("diff_design_malat.json")
  expect_true(inherits(design, "mzlibr_diff_design"))
  expect_identical(design$factors[[1L]]$name, "malat_load")
  expect_identical(design$factors[[1L]]$reference_level, "500ng")
  expect_identical(c(length(design$covariates), length(design$interactions)), c(0L, 0L))
  expect_true(is.na(design$batch_column))
  expect_identical(c(design$stratum_count, design$sample_count), c(1, 27))
  expect_identical(design$estimability_counts, c("estimable" = 2, "not_estimable:confounded" = 0, "not_estimable:rank_deficient" = 0))
  expect_identical(design$written$file_count, 2)
  stratum <- design$strata[["all"]]
  expect_true(inherits(stratum, "mzlibr_design_stratum"))
  expect_identical(stratum$column_names, c("intercept", "malat_load=125ng", "malat_load=250ng"))
  expect_identical(stratum$fitted_column_names, stratum$column_names)
  expect_identical(dim(stratum$matrix), c(27L, 3L))
  expect_identical(stratum$matrix, stratum$fitted_matrix)
  expect_identical(rownames(stratum$matrix), stratum$sample_ids)
  expect_true(endsWith(stratum$design_file, "design_1.tsv") && endsWith(stratum$contrasts_file, "contrasts_1.tsv"))
  expect_identical(design$records$fitted_weights, I(list(c(0, 1, 0), c(0, 0, 1))))
  expect_identical(design$column_names, c("stratum", "contrast_id", "label", "numerator", "denominator", "covariate",
    "covariate_unit", "covariate_scale", "estimability", "weights", "fitted_weights"))
  expect_identical(names(design$records), design$column_names)
  expect_identical(design$records$label, c("malat_load=125ng vs malat_load=500ng", "malat_load=250ng vs malat_load=500ng"))
})

test_that("the recorded design matrix is R's model.matrix(), run here and as recorded", {
  design <- diff_recorded("diff_design_factor_batch_covariate.json")
  stratum <- design$strata[[1L]]
  recorded_r <- as.matrix(utils::read.delim(diff_path("design_factor_batch_covariate_matrix.tsv"), check.names = FALSE))
  expect_true(max(abs(unname(stratum$matrix) - unname(recorded_r))) < 1e-12)
  sheet <- utils::read.delim(diff_path("design_factor_batch_covariate_samples.tsv"), stringsAsFactors = FALSE)
  sheet <- sheet[match(stratum$sample_ids, sheet[[1L]]), ]
  sheet$condition <- stats::relevel(factor(sheet$condition), ref = "ctrl")
  sheet$batch <- factor(sheet$batch)
  here <- stats::model.matrix(~ condition + batch + I((age - 50) / 10), data = sheet)
  expect_true(max(abs(unname(stratum$matrix) - unname(here))) < 1e-12)
  expect_identical(stratum$covariate_centres, c(age = 50))
  expect_identical(design$covariates[[1L]][c("name", "unit", "scale", "scaled_unit", "centre")],
    list(name = "age", unit = "years", scale = 10, scaled_unit = "decade", centre = 50))
  slope <- design$records[2L, ]
  expect_identical(c(slope$covariate, slope$covariate_unit, slope$covariate_scale),
    c("age", "decade", "years/10, centred at 50 years"))
})

test_that("a confounded design says why, and its contrast has no fitted weights", {
  design <- diff_recorded("diff_design_confounded.json")
  stratum <- design$strata[[1L]]
  expect_identical(stratum$rank, 2)
  expect_identical(stratum$aliased_terms, c("condition", "batch"))
  expect_identical(stratum$fitted_column_names, c("intercept", "condition=B"))
  expect_identical(design$records$estimability, "not_estimable:confounded")
  expect_true(is.null(design$records$fitted_weights[[1L]]), info = "NULL, not NA: no weights exist")
  expect_true(any(grepl("not estimable", design$caveats, fixed = TRUE)))
})

test_that("strata are designed on their own, and say what they cannot run", {
  design <- diff_recorded("diff_design_strata.json")
  expect_identical(names(design$strata), c("organism part=kidney", "organism part=liver", "organism part=muscle"))
  kidney <- design$strata[["organism part=kidney"]]
  expect_identical(kidney$not_run$label, "condition (one level: A)")
  expect_identical(kidney$not_run$reason, "not_in_stratum")
  expect_identical(kidney$column_names, "intercept")
  expect_identical(design$records$contrast_id, c("c1", "c1"))
  expect_identical(design$stratum_columns, "organism part")
  expect_identical(nrow(design$strata[["organism part=liver"]]$not_run), 0L)
})

# ---------------------------------------------------------------- observations: what is sent

test_that("observations write one line per run and send only what was given", {
  runs <- data.frame(file = c("a", "b", "c", "d"), sample = c("S1", "S1", "S1", "S2"),
    fraction = c(NA, 1, NA, NA), technical_replicate = c(NA, NA, 1, 2))
  request <- mz$diff_build_observations_request("p.tsv", runs, "msms_only", NULL, NULL, FALSE, NULL, NULL, NULL)
  expect_identical(request$args, c("diff", "observations", "--peptides", "p.tsv", "--basis", "msms_only"))
  expect_identical(request$stdin, "a\tS1\nb\tS1\t1\nc\tS1\t\t1\nd\tS2\t\t2")
  every <- mz$diff_build_observations_request("p.tsv", data.frame(file = "a", sample = "S1"), "mbr_kept", "g.tsv", "s.tsv",
    TRUE, 50, "o.tsv", "r.tsv")
  expect_identical(every$args, c("diff", "observations", "--peptides", "p.tsv", "--basis", "mbr_kept",
    "--protein-groups", "g.tsv", "--samples", "s.tsv", "--normalize", "--min-peptides", "50", "--out", "o.tsv",
    "--responses-out", "r.tsv"))
  # The result's own runs, file_name and sample_id, are accepted as they are.
  table <- observations_recorded("diff_observations_spike_in.json")
  again <- mz$diff_build_observations_request("p.tsv", table$runs, "msms_only", NULL, NULL, FALSE, NULL, NULL, NULL)
  expect_identical(again$stdin, "09-04-18_EcoliSpikeInSingleShot1x\tspike_1x\t0\t0\n09-04-18_EcoliSpikeInSingleShot2x\tspike_2x\t0\t0")
})

test_that("bad observations arguments are refused before the bridge", {
  build <- function(runs = data.frame(file = "a", sample = "S1"), basis = "msms_only", normalize = FALSE, min_peptides = NULL) {
    mz$diff_build_observations_request("p.tsv", runs, basis, NULL, NULL, normalize, min_peptides, NULL, NULL)
  }
  expect_error(build(basis = "mbr"), class = "mzlib_usage_error", contains = "basis must be one of msms_only, mbr_kept")
  expect_error(build(basis = NULL), class = "mzlib_usage_error", contains = "basis must be one of")
  expect_error(diff_observations("p.tsv", data.frame(file = "a", sample = "S1")), class = "mzlib_usage_error",
    contains = "basis is required")
  expect_error(build(runs = data.frame(file = character(0), sample = character(0))), class = "mzlib_usage_error", contains = "needs the runs")
  expect_error(build(runs = data.frame(file = "a")), class = "mzlib_usage_error", contains = "file and a sample column")
  expect_error(build(runs = "a\tS1"), class = "mzlib_usage_error", contains = "must be a data.frame")
  expect_error(build(runs = data.frame(file = "a", sample = "S1", fraction = 1.5)), class = "mzlib_usage_error",
    contains = "fraction must be a whole number")
  expect_error(build(runs = data.frame(file = "a", sample = "S1", technical_replicate = TRUE)), class = "mzlib_usage_error",
    contains = "technical replicate must be a whole number")
  expect_error(build(runs = data.frame(file = "a b\t", sample = "S1")), class = "mzlib_usage_error", contains = "tab or a line break")
  expect_error(build(normalize = "yes"), class = "mzlib_usage_error", contains = "normalize must be TRUE or FALSE")
  expect_error(build(min_peptides = 50.5), class = "mzlib_usage_error", contains = "min_peptides must be a whole number")
})

# ---------------------------------------------------------------- observations: what is parsed

test_that("the MetaMorpheus recording parses; missing is NA, never 0, and says why", {
  table <- observations_recorded("diff_observations_mm_mbr_kept.json")
  expect_true(inherits(table, "mzlibr_observations"))
  expect_identical(table$basis, "mbr_kept")
  expect_false(table$normalize)
  expect_true(is.na(table$min_peptides))
  expect_identical(c(table$peptide_count, table$sample_count, table$run_count, table$row_count), c(3, 18, 18, 54))
  expect_identical(nrow(table$records), 54L)
  expect_identical(names(table$state_counts), mz$DIFF_OBSERVATION_STATES)
  expect_identical(as.list(table$runs[1L, ]),
    list(file_name = "QE-002106_GM1_a", sample_id = "QE-002106_GM1_a", fraction = 0, technical_replicate = 0))
  alb <- table$protein_groups[table$protein_groups$name == "P02769", ]
  expect_identical(list(alb$q_value, alb$is_decoy, alb$is_contaminant, alb$is_entrapment, alb$genes), list(0, FALSE, TRUE, FALSE, "ALB"))
  expect_identical(table$protein_groups$name[table$protein_groups$is_decoy], "DECOY_P62750")
  expect_true(endsWith(table$protein_groups_file, "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"))
  expect_true(is.na(table$samples_file) && is.null(table$written) && is.null(table$responses_written))
  expect_identical(table$normalization$setting, "none")
  expect_true(all(is.na(table$normalization$shifts)))
  r <- table$records
  expect_true(all(is.na(r$log2_intensity) | r$log2_intensity > 0))
  expect_true(all(r$state %in% mz$DIFF_OBSERVATION_STATES))
  expect_true(all(r$state[!is.na(r$log2_intensity)] %in% c("quantified", "mbr_transferred")))
  expect_identical(r$protein_groups[[1L]], "Q16643")
  expect_identical(r$unique_protein_group[[1L]], "Q16643")
  expect_true(is.list(r$protein_groups), info = "a list column, even when every cell holds one group")
})

test_that("the two bases differ only where a transfer is", {
  kept <- observations_recorded("diff_observations_mm_mbr_kept.json")$records
  only <- observations_recorded("diff_observations_mm_msms_only.json")$records
  differ <- which(!(is.na(kept$log2_intensity) & is.na(only$log2_intensity)) &
    (is.na(kept$log2_intensity) != is.na(only$log2_intensity) | kept$log2_intensity != only$log2_intensity))
  expect_identical(length(differ), 1L)
  expect_identical(c(kept$state[differ], only$state[differ]), c("mbr_transferred", "mbr_transferred"))
  expect_true(is.na(only$log2_intensity[differ]))
  expect_identical(kept$sample[differ], "QE-002118_GM6_a")
})

test_that("the spike-in normalization takes away the spike, and wrote the responses", {
  plain <- observations_recorded("diff_observations_spike_in.json")
  shifted <- observations_recorded("diff_observations_spike_in_normalized.json")
  expect_true(shifted$normalize)
  expect_identical(shifted$min_peptides, 50)
  expect_identical(shifted$normalization$setting, "shared_peptide_median")
  expect_identical(shifted$normalization$reference_set_size, 86)
  shifts <- shifted$normalization$shifts
  expect_true(abs(shifts[["spike_1x"]] + shifts[["spike_2x"]]) < 1e-12)
  expect_identical(c(shifted$responses_written$row_count, shifted$responses_written$sample_count), c(243, 2))
  a <- plain$records
  b <- shifted$records
  expect_identical(is.na(a$log2_intensity), is.na(b$log2_intensity))
  known <- !is.na(a$log2_intensity)
  expect_true(max(abs(b$log2_intensity[known] - (a$log2_intensity[known] - shifts[a$sample[known]]))) < 1e-9)
  expect_identical(plain$normalization$setting, "none", info = "86 shared peptides, under mzLib's floor of 100")
})

test_that("the example inputs ship with the package", {
  for (name in c("malat_dilution_samples.tsv", "design_factor_batch_covariate_samples.tsv", "design_confounded_samples.tsv",
                 "stratified_samples.tsv", "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv",
                 "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv", "EcoliSpikeIn_QuantifiedPeptides.tsv")) {
    expect_true(nzchar(diff_path(name)), info = name)
  }
})

# ---------------------------------------------------------------- end to end, replayed

# The runs of mzLib's MetaMorpheus table: each file its own sample.
diff_mm_runs <- function() {
  header <- strsplit(readLines(diff_path("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"), n = 1L), "\t", fixed = TRUE)[[1L]]
  files <- sub("-calib$", "", sub("^Intensity_", "", header[startsWith(header, "Intensity_")]))
  data.frame(file = files, sample = files, stringsAsFactors = FALSE)
}

DIFF_SPIKE_IN <- data.frame(file = c("09-04-18_EcoliSpikeInSingleShot1x", "09-04-18_EcoliSpikeInSingleShot2x"),
  sample = c("spike_1x", "spike_2x"), stringsAsFactors = FALSE)

test_that("design and observations replay end to end, the model and runs choosing the recording", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  out <- file.path(tempdir(), "malat_design")
  malat <- diff_design(diff_path("malat_dilution_samples.tsv"), factors = list(factor_spec("malat_load", "500ng")), out = out)
  pairs <- diff_design(diff_path("malat_dilution_samples.tsv"), factors = list(factor_spec("malat_load")), out = out)
  expect_identical(c(malat$row_count, pairs$row_count), c(2, 3))
  covariate <- diff_design(diff_path("design_factor_batch_covariate_samples.tsv"),
    factors = list(factor_spec("condition", "ctrl")),
    covariates = list(covariate_spec("age", "years", scale = 10, scaled_unit = "decade", centre = 50)), batch = "batch")
  expect_identical(covariate$batch_column, "batch")
  expect_error(diff_design(diff_path("design_factor_batch_covariate_samples.tsv"), factors = list(factor_spec("condition", "ctrl")),
    batch = "batch"), class = "mzlib_usage_error", contains = "replay bridge")
  stratum <- malat$strata[["all"]]
  fit <- stats_fit("malat_dilution_log2_vs_fluc.tsv", stratum$design_file, c("malat_load=250ng", "malat_load=125ng"))
  both <- stats_fit("malat_dilution_log2_vs_fluc.tsv", stratum$design_file, c("malat_load=250ng", "malat_load=125ng"),
    contrasts = stratum$contrasts_file)
  expect_true(is.null(fit$contrasts) && !is.null(both$contrasts))
  kept <- diff_observations(diff_path("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"), diff_mm_runs(), basis = "mbr_kept",
    protein_groups = diff_path("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"))
  expect_identical(kept$basis, "mbr_kept")
  plain <- diff_observations(diff_path("EcoliSpikeIn_QuantifiedPeptides.tsv"), DIFF_SPIKE_IN, basis = "msms_only")
  shifted <- diff_observations(diff_path("EcoliSpikeIn_QuantifiedPeptides.tsv"), DIFF_SPIKE_IN, basis = "msms_only",
    min_peptides = 50, normalize = TRUE, responses_out = file.path(tempdir(), "spike_in_responses.tsv"))
  expect_identical(c(plain$normalize, shifted$normalize), c(FALSE, TRUE))
  expect_error(diff_observations(diff_path("EcoliSpikeIn_QuantifiedPeptides.tsv"), DIFF_SPIKE_IN[1L, ], basis = "msms_only"),
    class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

diff_live <- function() {
  skip_unless_bridge_has("diff design")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
}

# Enough of a design to compare two answers field by field, paths reduced to file names.
diff_comparable <- function(design) {
  base <- function(path) if (is.na(path)) NA_character_ else sub("^.*[/\\\\]", "", path)
  strata <- lapply(design$strata, function(s) {
    s$design_file <- base(s$design_file)
    s$contrasts_file <- base(s$contrasts_file)
    unclass(s)
  })
  list(strata = strata, records = design$records, caveats = design$caveats,
    estimability_counts = design$estimability_counts, unused_columns = design$unused_columns)
}

test_that("LIVE: the design recordings are still what the bridge says", {
  diff_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  out <- tempfile("malat_design")
  cases <- list(
    diff_design_malat.json = list(diff_path("malat_dilution_samples.tsv"), factors = list(factor_spec("malat_load", "500ng")), out = out),
    diff_design_factor_batch_covariate.json = list(diff_path("design_factor_batch_covariate_samples.tsv"),
      factors = list(factor_spec("condition", "ctrl")),
      covariates = list(covariate_spec("age", "years", scale = 10, scaled_unit = "decade", centre = 50)), batch = "batch"),
    diff_design_confounded.json = list(diff_path("design_confounded_samples.tsv"), factors = list(factor_spec("condition", "A")),
      batch = "batch"),
    diff_design_strata.json = list(diff_path("stratified_samples.tsv"), factors = list(factor_spec("condition", "A")),
      strata = "organism part")
  )
  for (name in names(cases)) {
    live <- do.call(diff_design, cases[[name]])
    expect_identical(diff_comparable(live), diff_comparable(diff_recorded(name)), info = name)
  }
})

test_that("LIVE: the observation recordings are still what the bridge says", {
  diff_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  groups <- diff_path("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
  mm <- diff_path("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv")
  spike <- diff_path("EcoliSpikeIn_QuantifiedPeptides.tsv")
  cases <- list(
    diff_observations_mm_mbr_kept.json = list(mm, diff_mm_runs(), basis = "mbr_kept", protein_groups = groups),
    diff_observations_mm_msms_only.json = list(mm, diff_mm_runs(), basis = "msms_only", protein_groups = groups),
    diff_observations_spike_in.json = list(spike, DIFF_SPIKE_IN, basis = "msms_only"),
    diff_observations_spike_in_normalized.json = list(spike, DIFF_SPIKE_IN, basis = "msms_only", min_peptides = 50,
      normalize = TRUE, responses_out = tempfile(fileext = ".tsv"))
  )
  for (name in names(cases)) {
    live <- do.call(diff_observations, cases[[name]])
    recorded <- observations_recorded(name)
    expect_identical(records_disagree(live$records, recorded$records), "", info = name)
    expect_identical(live$normalization$setting, recorded$normalization$setting, info = name)
    expect_true(isTRUE(all.equal(live$normalization$shifts, recorded$normalization$shifts, tolerance = 1e-9)), info = name)
    expect_identical(live$state_counts, recorded$state_counts, info = name)
    expect_identical(live$caveats, recorded$caveats, info = name)
  }
})

test_that("LIVE: the design file diff_design() writes is the committed design's model", {
  diff_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  stats <- function(name) system.file("extdata", "stats", name, package = "mzLibR")
  design <- diff_design(diff_path("malat_dilution_samples.tsv"), factors = list(factor_spec("malat_load", "500ng")),
    out = tempfile("design"))
  ours <- stats_fit(stats("malat_dilution_log2_vs_fluc.tsv"), design$strata[[1L]]$design_file,
    c("malat_load=250ng", "malat_load=125ng"))
  committed <- stats_fit(stats("malat_dilution_log2_vs_fluc.tsv"), stats("malat_dilution_design.tsv"),
    c("malat_250ng", "malat_125ng"))
  expect_identical(ours$status_counts, committed$status_counts)
  a <- ours$records$estimate
  b <- committed$records$estimate
  expect_identical(is.na(a), is.na(b))
  expect_true(all(abs(a - b)[!is.na(b)] <= 1e-9 * abs(b[!is.na(b)]) + 1e-12))
})

test_that("LIVE: the responses table diff_observations() writes is what stats_fit() reads", {
  diff_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  responses <- tempfile(fileext = ".tsv")
  table <- diff_observations(diff_path("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"), diff_mm_runs(), basis = "mbr_kept",
    responses_out = responses)
  design <- tempfile(fileext = ".tsv")
  writeLines(c("sample\tintercept", paste0(table$samples, "\t1")), design)
  fit <- stats_fit(responses, design, "intercept")
  expect_identical(fit$feature_count, table$peptide_count)
  expect_identical(fit$sample_names, table$samples)
})

test_that("LIVE: a refusal crosses as a usage error", {
  diff_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  sheet <- tempfile(fileext = ".tsv")
  writeLines(c("sample\tcondition", "s1\tA", "s2\t"), sheet)
  expect_error(diff_design(sheet, factors = list(factor_spec("condition", "A"))), class = "mzlib_usage_error",
    contains = "has no value for factor 'condition'")
  expect_error(diff_observations(diff_path("EcoliSpikeIn_QuantifiedPeptides.tsv"), DIFF_SPIKE_IN[1L, ], basis = "msms_only"),
    class = "mzlib_usage_error", contains = "not in the experimental design")
})
