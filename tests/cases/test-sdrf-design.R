# A label-free design from an SDRF (sdrf_design() and its two hand-offs), and isobaric kits.
#
# Ported from pyMzLib's test_sdrf_design.py and test_isobaric.py. The payloads were recorded from
# the real bridge by pyMzLib over three real PRIDE deposits' SDRFs and are shared verbatim. The
# design and every refusal are mzLib's; this layer shapes arguments and projects the table.

design_fixture <- function(name) {
  parsed <- mz$json_parse(paste(readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}
design_both <- c("factor value[genotype]", "factor value[treatment]")
design_valid <- function() mz$sdrf_parse_design(design_fixture("sdrf_design_PXD067622.json"))
design_refused <- function() mz$sdrf_parse_design(design_fixture("sdrf_design_PXD049018.json"))

# ---------------------------------------------------------------- sdrf_design: the payload

test_that("a valid design projects every run, its coordinates 0-based", {
  d <- design_valid()
  expect_true(inherits(d, "mzlibr_sdrf_design"))
  expect_true(d$is_valid)
  expect_identical(d$file_count, 24)
  expect_identical(nrow(d$records), 24L)
  expect_identical(c(length(d$refusals), length(d$notes)), c(0L, 0L))
  expect_identical(d$condition_columns, design_both)
  expect_true(d$condition_columns_declared)
  expect_false(d$searched_files_given)
  expect_identical(d$file_key_column, "comment[data file]")
  expect_true(is.null(d$written))
  expect_identical(length(unique(d$records$condition)), 8L)
  expect_identical(sort(unique(d$records$biological_replicate)), c(0, 1, 2))
  expect_identical(d$records$full_path, paste0(d$records$file_name, ".raw"))
  expect_true(grepl("24 file(s), 8 condition(s)", d$report, fixed = TRUE))
  expect_identical(length(d$caveats), 6L)
})

test_that("the design feeds FlashLFQ unchanged", {
  d <- design_valid()
  spectra <- sdrf_design_spectra(d)
  expect_identical(names(spectra), c("path", "condition", "biological_replicate", "technical_replicate", "fraction"))
  expect_identical(spectra$path[[1L]], "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4.raw")
  expect_identical(spectra$condition[[1L]], "SPRTN-TurboID CA_DMSO (vehicle)")
  expect_identical(unlist(spectra[1L, 3:5], use.names = FALSE), c(0, 0, 0))

  runs <- sdrf_design_run_design(d)
  expect_identical(runs$file_name[[1L]], "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4")
  # FlashLFQ's own stdin renderer accepts every row: 0-based, non-negative whole numbers.
  lines <- mz$flashlfq_design_stdin(runs)
  expect_identical(length(lines), 24L)
  expect_identical(strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]][-1L],
                   c("SPRTN-TurboID CA_DMSO (vehicle)", "0", "0", "0"))
})

test_that("a refusal is a result, and blocks only the hand-off", {
  d <- design_refused()
  expect_false(d$is_valid)
  expect_identical(length(d$refusals), 20L)
  expect_identical(c(d$file_count, nrow(d$records)), c(0, 0))
  expect_true(grepl("MSB67868ABand_01.raw", d$refusals[[1L]], fixed = TRUE))
  expect_true(grepl("'not available'", d$refusals[[1L]], fixed = TRUE))
  # The reasons travel with the error.
  expect_error(sdrf_design_spectra(d), class = "mzlib_usage_error", contains = c("refused", d$refusals[[1L]]))
  expect_error(sdrf_design_run_design(d), class = "mzlib_usage_error", contains = "refused")
  expect_error(sdrf_design_spectra(list()), class = "mzlib_usage_error", contains = "sdrf_design() result")
})

test_that("the renumbering is kept in notes, and the runs match the reference", {
  ranked <- mz$sdrf_parse_design(design_fixture("sdrf_design_studywide.json"))
  expect_true(ranked$is_valid)
  expect_identical(length(ranked$notes), 7L)
  expect_true(any(grepl("22 -> 1, 23 -> 2, 24 -> 3", ranked$notes, fixed = TRUE)))
  expect_identical(ranked$records, design_valid()$records)
})

# ---------------------------------------------------------------- sdrf_design: arguments

test_that("the condition columns travel tab-joined, and none sends no option", {
  request <- mz$sdrf_build_design_request("PXD067622.sdrf.tsv", design_both, NULL, NULL)
  expect_identical(request$args, c("sdrf", "design", "--path", "PXD067622.sdrf.tsv",
                                   "--condition-columns", "factor value[genotype]\tfactor value[treatment]"))
  expect_true(is.null(request$stdin))
  expect_identical(mz$sdrf_build_design_request("PXD067622.sdrf.tsv", NULL, NULL, NULL)$args,
                   c("sdrf", "design", "--path", "PXD067622.sdrf.tsv"))
})

test_that("searched files go on stdin with the flag, and out is passed through", {
  request <- mz$sdrf_build_design_request("a.sdrf.tsv", NULL, c("/data/a.raw", "b.raw"), "ExperimentalDesign.tsv")
  expect_true("--searched-files-stdin" %in% request$args)
  expect_identical(request$stdin, "/data/a.raw\nb.raw\n")
  expect_identical(utils::tail(request$args, 2L), c("--out", "ExperimentalDesign.tsv"))
})

test_that("bad arguments are refused before the bridge", {
  build <- function(...) {
    given <- list(path = "a.sdrf.tsv", condition_columns = NULL, searched_files = NULL, out = NULL)
    given[names(list(...))] <- list(...)
    do.call(mz$sdrf_build_design_request, given)
  }
  expect_error(build(condition_columns = character(0)), class = "mzlib_usage_error", contains = "non-empty")
  expect_error(build(condition_columns = "  "), class = "mzlib_usage_error", contains = "non-blank")
  expect_error(build(condition_columns = "a\tb"), class = "mzlib_usage_error", contains = "tab or newline")
  expect_error(build(condition_columns = 3), class = "mzlib_usage_error", contains = "character vector")
  expect_error(build(searched_files = ""), class = "mzlib_usage_error", contains = "non-blank")
  expect_error(build(searched_files = "a\nb"), class = "mzlib_usage_error", contains = "newline")
  expect_error(build(out = ""), class = "mzlib_usage_error", contains = "out must be")
  expect_error(build(path = ""), class = "mzlib_usage_error")
})

test_that("the design replays end to end, verb check included", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(sdrf_design("PXD067622.sdrf.tsv", condition_columns = design_both)$file_count, 24)
  expect_false(sdrf_design("PXD049018.sdrf.tsv", condition_columns = design_both)$is_valid)
})

# ---------------------------------------------------------------- isobaric_kits

kits_every <- function() mz$isobaric_parse_kits(design_fixture("isobaric_kits.json"))
kits_tmt18 <- function() mz$isobaric_parse_kits(design_fixture("isobaric_kits_TMT18.json"))

test_that("every kit is listed, in mzLib's order, and groups the table", {
  result <- kits_every()
  expect_true(inherits(result, "mzlibr_isobaric_kits"))
  expect_true(is.na(result$kit))
  expect_identical(result$kits$kit, c("TMT6", "TMT10", "TMT11", "TMT16", "TMT18", "iTRAQ4", "iTRAQ8", "diLeu4", "diLeu12"))
  expect_identical(result$kit_count, 9)
  expect_identical(result$record_count, sum(result$kits$channel_count))
  expect_identical(nrow(result$records), as.integer(result$record_count))
  for (i in seq_len(nrow(result$kits))) {
    rows <- result$records[result$records$kit == result$kits$kit[[i]], ]
    expect_identical(nrow(rows), as.integer(result$kits$channel_count[[i]]))
    # 1-based in R, ascending with the m/z.
    expect_identical(rows$channel_index, as.numeric(seq_len(nrow(rows))))
    expect_identical(rows$reporter_ion_mz, sort(rows$reporter_ion_mz))
  }
  expect_identical(result$absolute_tolerance, 0.003)
  expect_identical(length(result$caveats), 5L)
})

test_that("TMTpro 16 is the first sixteen of 18, and iTRAQ8 has no 120", {
  records <- kits_every()$records
  tmt16 <- records$reporter_ion_mz[records$kit == "TMT16"]
  tmt18 <- records$reporter_ion_mz[records$kit == "TMT18"]
  expect_identical(tmt16, tmt18[1:16])
  expect_identical(records$channel_label[records$kit == "iTRAQ8"],
                   c("113", "114", "115", "116", "117", "118", "119", "121"))
})

test_that("one kit, and its matching window", {
  result <- kits_tmt18()
  expect_identical(result$kit, "TMT18")
  expect_identical(result$kits$kit, "TMT18")
  labels <- result$records$channel_label
  expect_identical(labels[1:3], c("126", "127N", "127C"))
  expect_identical(labels[[length(labels)]], "135N")
  expect_true(all(abs(result$records$mz_min - (result$records$reporter_ion_mz - result$absolute_tolerance)) < 1e-9))
  expect_true(all(abs(result$records$mz_max - (result$records$reporter_ion_mz + result$absolute_tolerance)) < 1e-9))
})

test_that("a blank or non-string kit is refused before the bridge", {
  for (kit in list("", "   ", 10, NA_character_, c("TMT6", "TMT10"))) {
    expect_error(isobaric_kits(kit), class = "mzlib_usage_error", contains = "kit must be")
  }
})

test_that("kits replay end to end, and one kit never answers for all", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(isobaric_kits()$kit_count, 9)
  expect_identical(isobaric_kits("TMT18")$record_count, 18)
  expect_error(isobaric_kits("TMT10"), class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

design_live_dir <- function() {
  dir <- Sys.getenv("MZLIBR_SDRF_FIXTURES", "")
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "PXD067622.sdrf.tsv")),
          "no SDRF documents to read (set MZLIBR_SDRF_FIXTURES to pyMzLib's tests/fixtures)")
  normalizePath(dir)
}

test_that("LIVE: the design matches the recording and writes a 1-based file", {
  skip_unless_bridge_has("sdrf design")
  dir <- design_live_dir()
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  out <- file.path(tempdir(), "ExperimentalDesign.tsv")
  on.exit(unlink(out), add = TRUE)
  d <- sdrf_design(file.path(dir, "PXD067622.sdrf.tsv"), condition_columns = design_both, out = out)
  expect_identical(d$records, design_valid()$records)
  expect_identical(d$written$file_count, 24)
  lines <- readLines(out, warn = FALSE)
  expect_identical(lines[[1L]], "FileName\tCondition\tBiorep\tFraction\tTechrep")
  expect_identical(as.numeric(strsplit(lines[[2L]], "\t", fixed = TRUE)[[1L]][3L]),
                   d$records$biological_replicate[[1L]] + 1)
  expect_error(sdrf_design(file.path(dir, "PXD067622.sdrf.tsv"), out = file.path(tempdir(), "design.txt")),
               class = "mzlib_usage_error", contains = ".tsv")
})

test_that("LIVE: kits match the recording and refuse a partial name", {
  skip_unless_bridge_has("isobaric kits")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  expect_identical(isobaric_kits()$records, kits_every()$records)
  expect_identical(isobaric_kits("iTRAQ-4plex on K")$kits$kit, "iTRAQ4")
  expect_error(isobaric_kits("TMT10plex"), class = "mzlib_usage_error", contains = "TMT10plex")
})
