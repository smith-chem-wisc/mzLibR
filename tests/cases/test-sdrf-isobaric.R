# sdrf_isobaric_design(), the projection of mzLib's SdrfIsobaricDesign, and sdrf_restrict(), of
# SdrfSearchScope.Restrict.
#
# Ported from pyMzLib's test_sdrf_isobaric_design.py and test_sdrf_restrict.py. The payloads are
# pyMzLib's recordings from the real bridge over mzLib's own isobaric design fixtures, which the
# live tests read from pyMzLib's tests/fixtures (MZLIBR_SDRF_FIXTURES): PXD008841.subset.sdrf.tsv
# (TMT10, two plexes named only in the file names, 40 rows, 4 files), PXD061609.subset.sdrf.tsv
# (TMTpro 18, three SAX arms of one plex in comment[sample preparation batch]) and
# PXD008841_searched.sdrf.tsv (PXD008841 restricted to a calibrated search by sdrf restrict out =).
#
# What this layer must get right: every option is sent only when given, the condition columns cross
# as ONE tab-joined option and searched files on stdin, a refusal stays a result, the numbers stay
# 1-based (nothing converts them), and a restricted document is an mzlibr_sdrf.

ISOBARIC_POOL <- "TMT_?pool(?<plex>\\d+)"
ISOBARIC_BATCH <- "comment[sample preparation batch]"
ISOBARIC_COLUMNS <- c(
  "full_path", "plex", "plex_id", "sample_name", "channel", "condition",
  "biological_replicate", "fraction", "technical_replicate"
)
RESTRICT_POOL1 <- "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
RESTRICT_POOL5 <- "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
# A calibrated run named by a Windows path, the next one bare, the same run's raw file as well, a
# renamed file whose acquisition is named, and a file the SDRF lacks.
RESTRICT_MESSY <- stats::setNames(
  c(
    paste0("D:\\search\\", RESTRICT_POOL1, "1-calib.mzML"), paste0(RESTRICT_POOL1, "2-calib.mzML"),
    paste0(RESTRICT_POOL1, "2.raw"), "sample_A.mzML", "not_in_sdrf.mzML"
  ),
  c("", "", "", paste0(RESTRICT_POOL5, "1.raw"), "")
)
RESTRICT_CLEAN <- paste0("D:\\search\\", rep(c(RESTRICT_POOL1, RESTRICT_POOL5), each = 2), c(1, 2), "-calib.mzML")
# The same four under this platform's own separator. restrict reduces either separator on any
# platform, but the design readers split only this platform's (Path.GetFileName), so a design given
# D:\search\... on Linux keys on the whole string (SDRF-B5 asks mzLib for one rule). file.path()
# writes "/", which .NET reads as a separator on Windows too.
RESTRICT_NATIVE <- file.path("", "search", paste0(rep(c(RESTRICT_POOL1, RESTRICT_POOL5), each = 2), c(1, 2), "-calib.mzML"))

isobaric_design_of <- function(name) mz$sdrf_parse_isobaric_design(recorded_data(name))

# ---------------------------------------------------------------- isobaric design: what is sent

test_that("a valid design projects every annotated channel", {
  seen <- recorded_bridge_call("sdrf_isobaric_design_PXD008841.json", function() {
    sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10", plex_file_name_pattern = ISOBARIC_POOL)
  })
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c(
    "sdrf", "isobaric-design", "--path", "PXD008841.subset.sdrf.tsv",
    "--kit", "TMT10", "--plex-file-name-pattern", ISOBARIC_POOL
  ))
  expect_true(is.null(seen$stdin))
  expect_identical(seen$timeout, 60)
  d <- seen$result
  expect_true(inherits(d, "mzlibr_sdrf_isobaric_design"))
  expect_true(d$is_valid)
  expect_identical(c(length(d$refusals), length(d$notes)), c(0L, 0L))
  expect_identical(c(d$kit, d$kit_type), c("TMT10", "TMT10"))
  expect_identical(length(d$kit_channels), 10L)
  expect_identical(d$plex_file_name_pattern, ISOBARIC_POOL)
  expect_true(is.na(d$plex_column) && is.na(d$single_plex))
  expect_identical(d$plex_ids, c("1" = 1, "5" = 2))
  expect_identical(d$plex_count, 2)
  expect_identical(c(d$file_count, d$record_count), c(4, 40))
  expect_identical(nrow(d$records), 40L)
  expect_identical(d$file_key_column, "comment[data file]")
  expect_identical(d$condition_columns, "factor value[disease]")
  expect_false(d$condition_columns_declared)
  expect_false(d$searched_files_given)
  expect_true(is.null(d$written))
  expect_identical(d$column_names, ISOBARIC_COLUMNS)
  expect_identical(names(d$records), ISOBARIC_COLUMNS)
  expect_identical(length(d$caveats), 8L)
  expect_true(startsWith(d$caveats[[1L]], "1-BASED"))
  expect_true(grepl("4 file(s), 2 plex(es), kit TMT10", d$report, fixed = TRUE))
})

test_that("the numbers stay 1-based", {
  records <- isobaric_design_of("sdrf_isobaric_design_PXD008841.json")$records
  expect_identical(sort(unique(records$fraction)), c(1, 2))
  expect_identical(unique(records$technical_replicate), 1)
  expect_identical(unique(records$biological_replicate), 1)
  expect_identical(sort(unique(records$plex_id)), c(1, 2))
  first_file <- records[records$full_path == records$full_path[[1L]], ]
  # Reporter m/z order.
  expect_identical(first_file$channel, isobaric_design_of("sdrf_isobaric_design_PXD008841.json")$kit_channels)
  expect_identical(first_file$sample_name[[nrow(first_file)]], "pool")
  row <- records[10L, ]
  expect_identical(list(row$channel, row$sample_name, row$plex_id, row$fraction, row$biological_replicate),
                   list("131N", "pool", 1, 1, 1))
})

test_that("a refusal is a result", {
  seen <- recorded_bridge_call("sdrf_isobaric_design_PXD008841_noplex.json", function() {
    sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10")
  })
  expect_identical(seen$args, c("sdrf", "isobaric-design", "--path", "PXD008841.subset.sdrf.tsv", "--kit", "TMT10"))
  d <- seen$result
  expect_false(d$is_valid)
  expect_true(startsWith(d$refusals[[1L]], "No plex source was declared"))
  expect_identical(d$plex_source, "(not declared)")
  expect_identical(c(d$file_count, d$record_count), c(0, 0))
  # A zero-row frame that keeps its columns, never a row of NA.
  expect_identical(nrow(d$records), 0L)
  expect_identical(names(d$records), ISOBARIC_COLUMNS)
  expect_identical(length(d$plex_ids), 0L)
})

test_that("a declared batch column gives three plexes", {
  seen <- recorded_bridge_call("sdrf_isobaric_design_PXD061609.json", function() {
    sdrf_isobaric_design("PXD061609.subset.sdrf.tsv", kit = "TMT18", plex_column = ISOBARIC_BATCH)
  })
  expect_identical(utils::tail(seen$args, 4L), c("--kit", "TMT18", "--plex-column", ISOBARIC_BATCH))
  d <- seen$result
  expect_true(d$is_valid)
  expect_identical(d$plex_ids, c(HighSalt = 1, LowSalt = 2, NoSAX = 3))
  expect_identical(d$record_count, 108)
  expect_identical(length(d$kit_channels), 18L)
})

test_that("one plex over three arms is refused, every reason at once", {
  d <- isobaric_design_of("sdrf_isobaric_design_PXD061609_oneplex.json")
  expect_false(d$is_valid)
  expect_identical(d$single_plex, "tube1")
  expect_identical(length(d$refusals), 36L)
  expect_true(all(grepl("more than once", d$refusals, fixed = TRUE)))
})

test_that("condition columns and searched files cross as sdrf_design sends them", {
  request <- mz$sdrf_build_isobaric_request("PXD008841_searched.sdrf.tsv", "TMT10", NULL, ISOBARIC_POOL, NULL,
    "factor value[disease]", c("/data/a-calib.mzML", "b-calib.mzML"), "design")
  expect_identical(utils::tail(request$args, 5L), c(
    "--condition-columns", "factor value[disease]", "--searched-files-stdin", "--out", "design"
  ))
  expect_identical(request$stdin, "/data/a-calib.mzML\nb-calib.mzML")
  two <- mz$sdrf_build_isobaric_request("a.sdrf.tsv", NULL, NULL, NULL, NULL,
    c("factor value[genotype]", "factor value[treatment]"), NULL, NULL)
  expect_identical(utils::tail(two$args, 1L), "factor value[genotype]\tfactor value[treatment]")
})

test_that("a restricted document keys on the searched name, and out writes TmtDesign.txt", {
  d <- isobaric_design_of("sdrf_isobaric_design_searched.json")
  expect_identical(d$file_key_column, "comment[searched data file]")
  expect_true(d$searched_files_given)
  # Recorded on Windows from searched files named D:\search\...: compared as strings, never split.
  expect_identical(sort(unique(d$records$full_path)), sort(RESTRICT_CLEAN))
  expect_identical(c(d$written$record_count, d$written$file_count), c(40, 4))
  expect_true(endsWith(d$written$path, "TmtDesign.txt"))
})

test_that("no option is sent unless given", {
  request <- mz$sdrf_build_isobaric_request("PXD008841.subset.sdrf.tsv", NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  expect_identical(request$args, c("sdrf", "isobaric-design", "--path", "PXD008841.subset.sdrf.tsv"))
  expect_true(is.null(request$stdin))
})

test_that("bad isobaric arguments are refused before the bridge", {
  cases <- list(
    list(list(kit = ""), "kit must be a non-blank string"),
    list(list(kit = 10), "kit must be a non-blank string"),
    list(list(kit = c("TMT10", "TMT11")), "kit must be a non-blank string"),
    list(list(plex_column = "  "), "plex_column must be a non-blank string"),
    list(list(plex_file_name_pattern = "a\nb"), "newline"),
    list(list(single_plex = ""), "single_plex must be a non-blank string"),
    list(list(condition_columns = character(0)), "non-empty"),
    list(list(condition_columns = "a\tb"), "tab or newline"),
    list(list(searched_files = ""), "non-blank"),
    list(list(out = ""), "out must name a folder"),
    list(list(path = ""), "path must be")
  )
  for (case in cases) {
    arguments <- utils::modifyList(list(path = "PXD008841.subset.sdrf.tsv"), case[[1L]])
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() do.call(sdrf_isobaric_design, arguments)),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

test_that("an older bridge is named, not spawned, for either SDRF verb", {
  for (call in list(
    function() sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10"),
    function() sdrf_restrict("PXD008841.subset.sdrf.tsv", "x.mzML")
  )) {
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the verb"), call, version = "readers_formats.json"),
      class = "mzlib_usage_error", contains = "pyMzLib 0.6.0"
    )
  }
})

test_that("isobaric designs replay end to end, each recording only for its own call", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_true(sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10", plex_file_name_pattern = ISOBARIC_POOL)$is_valid)
  # No plex source: the refusal, never the design read with the pattern.
  expect_false(sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10")$is_valid)
  expect_identical(sdrf_isobaric_design("PXD061609.subset.sdrf.tsv", kit = "TMT18", plex_column = ISOBARIC_BATCH)$plex_count, 3)
})

# ---------------------------------------------------------------- restrict: what is sent

test_that("a named element sends its acquired file, and only that one", {
  seen <- recorded_bridge_call("sdrf_restrict_PXD008841.json", function() {
    sdrf_restrict("PXD008841.subset.sdrf.tsv", RESTRICT_MESSY, limit = 2)
  })
  expect_identical(seen$args, c("sdrf", "restrict", "--path", "PXD008841.subset.sdrf.tsv", "--limit", "2"))
  expect_identical(seen$stdin, paste(c(
    paste0("D:\\search\\", RESTRICT_POOL1, "1-calib.mzML"),
    paste0(RESTRICT_POOL1, "2-calib.mzML"),
    paste0(RESTRICT_POOL1, "2.raw"),
    paste0("sample_A.mzML\t", RESTRICT_POOL5, "1.raw"),
    "not_in_sdrf.mzML"
  ), collapse = "\n"))
  expect_identical(seen$timeout, 60)
  expect_identical(seen$result$searched_files, unname(RESTRICT_MESSY))
  expect_identical(seen$result$acquired_name_count, 1)
})

test_that("an unnamed vector sends one name per line, and out is passed through", {
  request <- mz$sdrf_build_restrict_request("PXD008841.subset.sdrf.tsv", RESTRICT_CLEAN, "searched.sdrf.tsv", 2, 0)
  expect_identical(utils::tail(request$args, 4L), c("--limit", "2", "--out", "searched.sdrf.tsv"))
  expect_identical(request$stdin, paste(RESTRICT_CLEAN, collapse = "\n"))
  # NA names are unknown acquisitions, as "" is.
  partly <- stats::setNames(c("a.mzML", "b.mzML"), c(NA, "B.raw"))
  expect_identical(mz$sdrf_restrict_lines(partly), c("a.mzML", "b.mzML\tB.raw"))
})

test_that("bad searched files and options are refused before the bridge", {
  cases <- list(
    list(list(searched_files = character(0)), "at least one searched file"),
    list(list(searched_files = ""), "non-blank"),
    list(list(searched_files = NA_character_), "non-blank"),
    list(list(searched_files = "a\tb"), "tab or newline"),
    list(list(searched_files = c("b\n.raw" = "a.mzML")), "tab or newline"),
    list(list(searched_files = list("a.mzML")), "character vector"),
    list(list(searched_files = 42), "character vector"),
    list(list(out = ""), "ending in .tsv"),
    list(list(limit = -1), "limit"),
    list(list(offset = -1), "offset"),
    list(list(path = ""), "path must be")
  )
  for (case in cases) {
    arguments <- utils::modifyList(list(path = "PXD008841.subset.sdrf.tsv", searched_files = RESTRICT_CLEAN), case[[1L]])
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() do.call(sdrf_restrict, arguments)),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

# ---------------------------------------------------------------- restrict: what is projected

restricted <- function() mz$sdrf_parse_restricted(recorded_data("sdrf_restrict_PXD008841.json"))

test_that("the reports are mzLib's", {
  scoped <- restricted()
  expect_identical(c(scoped$input_row_count, scoped$row_count, scoped$returned_count), c(40, 30, 2))
  expect_true(scoped$truncated)
  expect_identical(scoped$offset, 0)
  expect_identical(scoped$searched_without_row, "not_in_sdrf.mzML")
  expect_identical(scoped$dropped_data_files, paste0(RESTRICT_POOL5, "2.raw"))
  expect_identical(scoped$ambiguous, paste0(RESTRICT_POOL1, "2: ", RESTRICT_POOL1, "2-calib.mzML, ", RESTRICT_POOL1, "2.raw"))
  expect_identical(scoped$shared_stem, character(0))
  expect_true(is.null(scoped$written))
  expect_identical(length(scoped$caveats), 6L)
})

test_that("the restricted document is an SDRF document", {
  scoped <- restricted()
  expect_true(inherits(scoped, "mzlibr_restricted_sdrf"))
  expect_true(inherits(scoped, "mzlibr_sdrf"))
  data_file <- sdrf_index_of(scoped, "comment[data file]")
  expect_identical(scoped$columns[[data_file + 1L]], "comment[searched data file]")
  expect_identical(sdrf_value(scoped, "comment[data file]"), rep(paste0(RESTRICT_POOL1, "1.raw"), 2L))
  # The name, no folder.
  expect_identical(sdrf_searched_data_files(scoped), rep(paste0(RESTRICT_POOL1, "1-calib.mzML"), 2L))
  # Padded, never ragged.
  expect_true(all(lengths(scoped$rows) == length(scoped$columns)))
  expect_identical(sdrf_ragged_row_count(scoped), 0L)
})

test_that("out reports the file mzLib wrote and the whole document's rows", {
  written <- mz$sdrf_parse_restricted(recorded_data("sdrf_restrict_PXD008841_out.json"))
  expect_identical(written$written, list(path = "PXD008841_searched.sdrf.tsv", row_count = 40))
  expect_identical(c(written$row_count, written$returned_count), c(40, 2))
})

test_that("a restriction replays end to end, as the help page runs it", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  scoped <- sdrf_restrict("PXD008841.subset.sdrf.tsv", RESTRICT_MESSY, limit = 2)
  expect_identical(scoped$row_count, 30)
  expect_identical(sdrf_restrict("PXD008841.subset.sdrf.tsv", RESTRICT_CLEAN, out = "PXD008841_searched.sdrf.tsv",
                                 limit = 2)$written$row_count, 40)
  # The acquired name is part of the call: without it, the recording does not answer.
  expect_error(sdrf_restrict("PXD008841.subset.sdrf.tsv", unname(RESTRICT_MESSY), limit = 2),
               class = "mzlib_usage_error", contains = "acquired names")
})

# ---------------------------------------------------------------- the pool's written path

test_that("a pool's written path is the file mzLib wrote, passed through as sent", {
  payload <- recorded_data("sdrf_pool_two.json")
  payload$written <- list(path = "C:/out/x.tsv.sdrf.tsv", row_count = 24)
  expect_identical(mz$sdrf_parse_pooled(payload)$written$path, "C:/out/x.tsv.sdrf.tsv")
})

# ---------------------------------------------------------------- against a real mzLib

sdrf_live_fixtures <- function(verb) {
  skip_unless_bridge_has(verb)
  dir <- Sys.getenv("MZLIBR_SDRF_FIXTURES", "")
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "PXD008841.subset.sdrf.tsv")),
          "no isobaric SDRF documents to read (set MZLIBR_SDRF_FIXTURES to pyMzLib's tests/fixtures)")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  normalizePath(dir)
}

test_that("LIVE: the design matches the recording and writes a 1-based TmtDesign.txt", {
  dir <- sdrf_live_fixtures("sdrf isobaric-design")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  folder <- file.path(tempfile("isobaric-"), "new", "design")
  on.exit(unlink(dirname(dirname(folder)), recursive = TRUE), add = TRUE)
  d <- sdrf_isobaric_design(file.path(dir, "PXD008841.subset.sdrf.tsv"), kit = "TMT10",
    plex_file_name_pattern = ISOBARIC_POOL, out = folder)
  recorded <- isobaric_design_of("sdrf_isobaric_design_PXD008841.json")
  expect_identical(d$records, recorded$records)
  expect_identical(d$plex_ids, recorded$plex_ids)
  expect_identical(normalizePath(d$written$path, winslash = "/"), normalizePath(file.path(folder, "TmtDesign.txt"), winslash = "/"))
  expect_identical(c(d$written$file_count, d$written$record_count), c(4, 40))
  lines <- readLines(file.path(folder, "TmtDesign.txt"), warn = FALSE, encoding = "UTF-8")
  expect_identical(strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]], c(
    "File", "Plex", "Sample Name", "TMT Channel", "Condition", "Biological Replicate", "Fraction", "Technical Replicate"
  ))
  first <- d$records[1L, ]
  # The file and the wire carry the same 1-based numbers.
  expect_identical(strsplit(lines[[2L]], "\t", fixed = TRUE)[[1L]], c(
    first$full_path, first$plex, first$sample_name, first$channel, first$condition,
    format(first$biological_replicate), format(first$fraction), format(first$technical_replicate)
  ))
})

test_that("LIVE: an unknown kit is isobaric_kits()'s usage error", {
  dir <- sdrf_live_fixtures("sdrf isobaric-design")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  from_design <- expect_error(sdrf_isobaric_design(file.path(dir, "PXD008841.subset.sdrf.tsv"), kit = "TMT10plex",
    plex_file_name_pattern = ISOBARIC_POOL), class = "mzlib_usage_error")
  from_kits <- expect_error(isobaric_kits("TMT10plex"), class = "mzlib_usage_error")
  expect_identical(conditionMessage(from_design), conditionMessage(from_kits))
})

test_that("LIVE: the restriction matches the recording", {
  dir <- sdrf_live_fixtures("sdrf restrict")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  scoped <- sdrf_restrict(file.path(dir, "PXD008841.subset.sdrf.tsv"), RESTRICT_MESSY, limit = 2)
  recorded <- restricted()
  expect_identical(scoped$rows, recorded$rows)
  expect_identical(scoped$columns, recorded$columns)
  for (name in c("searched_without_row", "dropped_data_files", "ambiguous", "shared_stem")) {
    expect_identical(scoped[[name]], recorded[[name]], info = name)
  }
})

test_that("LIVE: out writes the whole document where mzLib says", {
  dir <- sdrf_live_fixtures("sdrf restrict")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  folder <- file.path(tempfile("restrict-"), "new")
  on.exit(unlink(dirname(folder), recursive = TRUE), add = TRUE)
  scoped <- sdrf_restrict(file.path(dir, "PXD008841.subset.sdrf.tsv"), RESTRICT_CLEAN,
    out = file.path(folder, "searched.tsv"), limit = 1)
  # mzLib adds .sdrf.tsv to a name not ending in it, and says which file it wrote.
  expect_identical(basename(scoped$written$path), "searched.tsv.sdrf.tsv")
  expect_true(file.exists(scoped$written$path))
  expect_false(file.exists(file.path(folder, "searched.tsv")))
  back <- sdrf_read(scoped$written$path)
  expect_identical(c(back$row_count, scoped$returned_count), c(40, 1))
  expect_identical(back$rows, sdrf_read(file.path(dir, "PXD008841_searched.sdrf.tsv"))$rows)
  expect_error(sdrf_restrict(file.path(dir, "PXD008841.subset.sdrf.tsv"), RESTRICT_CLEAN,
    out = file.path(folder, "searched.txt")), class = "mzlib_usage_error", contains = ".tsv")
})

test_that("LIVE: a pool's out names the file mzLib wrote", {
  dir <- sdrf_live_fixtures("sdrf restrict")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  folder <- tempfile("pool-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  pooled <- sdrf_pool(c(one = file.path(dir, "PXD008841.subset.sdrf.tsv")), out = file.path(folder, "x.tsv"), limit = 0)
  expect_identical(basename(pooled$written$path), "x.tsv.sdrf.tsv")
  expect_true(file.exists(pooled$written$path))
})

test_that("LIVE: restrict, then read the design", {
  dir <- sdrf_live_fixtures("sdrf isobaric-design")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  d <- sdrf_isobaric_design(file.path(dir, "PXD008841_searched.sdrf.tsv"), kit = "TMT10",
    plex_file_name_pattern = ISOBARIC_POOL, searched_files = RESTRICT_NATIVE)
  expect_true(d$is_valid, info = d$report)
  expect_identical(d$file_key_column, "comment[searched data file]")
  expect_identical(sort(unique(d$records$full_path)), sort(RESTRICT_NATIVE))
  unrestricted <- sdrf_isobaric_design(file.path(dir, "PXD008841.subset.sdrf.tsv"), kit = "TMT10",
    plex_file_name_pattern = ISOBARIC_POOL, searched_files = RESTRICT_NATIVE)
  expect_false(unrestricted$is_valid)
  expect_true(any(grepl("has no SDRF row", unrestricted$refusals, fixed = TRUE)))
})

test_that("LIVE: a document without a data-file column is a bridge error", {
  sdrf_live_fixtures("sdrf restrict")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  path <- tempfile(fileext = ".sdrf.tsv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("source name\tassay name", "S1\trun 1"), path)
  condition <- expect_error(sdrf_restrict(path, "X.raw"), class = "mzlib_bridge_error")
  expect_identical(condition$error_type, "ArgumentException")
  expect_true(grepl("comment[data file]", conditionMessage(condition), fixed = TRUE))
})
