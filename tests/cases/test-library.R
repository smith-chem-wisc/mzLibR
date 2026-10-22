# library_read() and library_read_many(): mzLib's `.msl` spectral-library reader.
#
# Ported from pyMzLib's test_library.py and test_library_live.py. Every payload is pyMzLib's
# recording from the real bridge, over the two libraries in pyMzLib's tests/fixtures/library/
# (written by the bridge tests' LibraryFixtures), which the live tests read from there, as the SDRF
# ones read its SDRF files (MZLIBR_SDRF_FIXTURES). They are not copied here: an .msl file begins
# with the bytes R CMD check takes for a Windows executable.
#
#   myPrositLib.msl        mzLib's own myPrositLib.msp test library, five Prosit-predicted
#                          precursors, converted by mzLib. library_read_prosit.json
#   synthetic_library.msl  six precursors covering every field: target, decoy, entrapment target
#                          with calibrated minutes, entrapment decoy, a predicted precursor with
#                          every "not available" sentinel, a proteoform. library_read_synthetic.json,
#                          library_read_fragments.json (one precursor's fragment ions),
#                          library_read_index_only.json
#   library_read_many.json both, through --paths-stdin
#
# What this layer must get right: the call it sends, and that it projects mzLib's answer without
# repairing it - an NA stays NA, never a sentinel, and nothing is converted but precursor_index and
# source_index, which are 1-based in R.

library_call <- function(fixture, call) recorded_bridge_call(fixture, call)
library_synthetic <- function() mz$library_parse_records(recorded_data("library_read_synthetic.json"))

# ---------------------------------------------------------------- the call

test_that("a plain read sends only the path", {
  seen <- library_call("library_read_prosit.json", function() library_read("myPrositLib.msl"))
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("library", "read", "--path", "myPrositLib.msl"))
  expect_true(is.null(seen$stdin))
  expect_true(is.null(seen$timeout))
})

test_that("every option is sent under its wire name", {
  args <- mz$library_build_read_args("synthetic_library.msl", TRUE, 0, 1, 2, "ions.tsv")
  expect_identical(args, c(
    "library", "read", "--path", "synthetic_library.msl", "--limit", "1", "--offset", "2",
    "--out", "ions.tsv", "--fragments", "--index-only-threshold", "0"
  ))
  # Past an R integer, as a 4 GiB threshold is, the number still crosses whole.
  expect_identical(utils::tail(mz$library_build_read_args("a.msl", FALSE, 4294967296, NULL, 0, NULL), 1L), "4294967296")
  seen <- library_call("library_read_fragments.json", function() {
    library_read("synthetic_library.msl", fragments = TRUE, limit = 1, timeout = 30)
  })
  expect_identical(seen$timeout, 30)
})

test_that("the threshold is left to mzLib when not given", {
  seen <- library_call("library_read_synthetic.json", function() library_read("synthetic_library.msl"))
  expect_false("--index-only-threshold" %in% seen$args)
  # The echo is mzLib's own default, MslFileTypeHandler.DefaultIndexOnlyThresholdBytes.
  expect_identical(seen$result$index_only_threshold, 1073741824)
})

test_that("bad arguments are refused before the bridge is spawned, the verb check included", {
  cases <- list(
    list(list(fragments = "yes"), "fragments must be TRUE or FALSE"),
    list(list(fragments = NA), "fragments must be TRUE or FALSE"),
    list(list(index_only_threshold = -1), "index_only_threshold"),
    list(list(index_only_threshold = TRUE), "index_only_threshold"),
    list(list(index_only_threshold = 1.5e9 + 0.5), "index_only_threshold"),
    list(list(index_only_threshold = "1GB"), "index_only_threshold"),
    list(list(limit = -1), "limit"),
    list(list(offset = TRUE), "offset"),
    list(list(path = ""), "path must be"),
    list(list(path = c("a.msl", "b.msl")), "library_read_many()")
  )
  for (case in cases) {
    arguments <- utils::modifyList(list(path = "synthetic_library.msl"), case[[1L]])
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() do.call(library_read, arguments),
                           version = "no-such-fixture.json"),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

test_that("read_many sends the list on stdin", {
  seen <- library_call("library_read_many.json", function() {
    library_read_many(c("synthetic_library.msl", "myPrositLib.msl"), threads = 2, on_error = "skip", fragments = TRUE)
  })
  expect_identical(seen$args, c(
    "library", "read", "--paths-stdin", "--threads", "2", "--on-error", "skip", "--fragments"
  ))
  expect_identical(seen$stdin, "synthetic_library.msl\nmyPrositLib.msl\n")
})

test_that("read_many refuses bad arguments before the bridge", {
  cases <- list(
    list(list(paths = character(0)), "at least one path"),
    list(list(paths = c("a.msl", "")), "Path 2"),
    list(list(paths = list("a.msl")), "character vector of paths"),
    list(list(threads = 0), "threads"),
    list(list(on_error = "ignore"), "on_error"),
    list(list(index_only_threshold = "1GB"), "index_only_threshold")
  )
  for (case in cases) {
    arguments <- utils::modifyList(list(paths = "a.msl"), case[[1L]])
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() do.call(library_read_many, arguments),
                           version = "no-such-fixture.json"),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

test_that("an older bridge is named, not spawned, for the library verb", {
  expect_error(
    recorded_bridge_call(function(args, stdin) stop("reached the verb"), function() library_read("a.msl"),
                         version = "readers_formats.json"),
    class = "mzlib_usage_error", contains = "pyMzLib 0.6.0"
  )
})

# ---------------------------------------------------------------- the projection

test_that("the envelope is projected field for field", {
  result <- library_synthetic()
  wire <- recorded_data("library_read_synthetic.json")
  expect_true(inherits(result, "mzlibr_library_records"))
  expect_true(is.na(result$file_type))
  expect_identical(result$reader, "MslLibrary")
  expect_identical(list(result$format_version, result$is_compressed, result$is_index_only), list(3, FALSE, FALSE))
  expect_identical(c(result$target_count, result$decoy_count, result$proteoform_count), c(4, 2, 1))
  expect_identical(c(result$record_count, result$returned_count, result$row_count), c(6, 6, 6))
  expect_identical(result$min_precursor_mz, wire$min_precursor_mz)
  expect_identical(result$max_precursor_mz, wire$max_precursor_mz)
  expect_true(is.na(result$rows_not_read))
  expect_true(is.na(result$retention_time_unit))
  expect_false(result$fragments_included)
  expect_false(result$truncated)
  expect_identical(result$offset, 0)
  expect_true(is.null(result$output))
  expect_identical(result$caveats, mz$wire_strings(wire$caveats))
  expect_identical(result$column_names, mz$wire_strings(wire$column_names))
  expect_identical(names(result$records), result$column_names)
  expect_identical(nrow(result$excluded_fields), 0L)
})

test_that("rows are the stored order, not m/z order, and precursor_index is 1-based", {
  records <- library_synthetic()$records
  expect_identical(records$precursor_index, c(1, 2, 3, 4, 5, 6))
  expect_false(identical(records$precursor_mz, sort(records$precursor_mz)))
})

test_that("labels and retention-time units are per row", {
  records <- library_synthetic()$records
  expect_identical(records$decoy_contaminant_target, c("T", "D", "ET", "ED", "T", "T"))
  expect_identical(records$is_decoy, c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(records$is_entrapment, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
  # Precursor 3 is calibrated minutes; every other is iRT. Nothing is converted.
  expect_identical(records$retention_time_unit, c("iRT", "iRT", "minutes", "iRT", "iRT", "iRT"))
  expect_true(abs(records$retention_time[[3L]] - 21.4) < 1e-4)
})

test_that("mzLib's not-available sentinels arrive as NA", {
  records <- library_synthetic()$records
  predicted <- 5L
  for (name in c("ion_mobility", "q_value", "nce", "protein_accession", "protein_name", "gene_name")) {
    expect_true(is.na(records[[name]][[predicted]]), info = name)
  }
  expect_identical(records$dissociation_type[[predicted]], "Unknown")
  expect_true(abs(records$ion_mobility[[1L]] - 0.95) < 1e-6)
  expect_true(abs(records$q_value[[1L]] - 0.001) < 1e-12)
  expect_identical(records$nce[[1L]], 28)
})

test_that("a library with no proteins names those fields absent", {
  result <- mz$library_parse_records(recorded_data("library_read_prosit.json"))
  expect_identical(result$absent_fields, c("ion_mobility", "protein_accession", "protein_name", "gene_name"))
  expect_true(all(is.na(result$records$protein_accession)))
  # Converted from a text library: mzLib's not-available defaults for what .msp cannot carry.
  expect_true(all(is.na(result$records$nce)))
  expect_identical(result$records$source_type, rep("Predicted", 5L))
  expect_identical(result$records$full_sequence[[4L]], "QSQHM[Common Variable:Oxidation on M]TEVVR")
  expect_identical(result$records$retention_time_unit[1:2], c("iRT", "iRT"))
})

test_that("fragment rows repeat the precursor and keep mzLib's values", {
  result <- mz$library_parse_records(recorded_data("library_read_fragments.json"))
  rows <- result$records
  expect_true(result$fragments_included)
  expect_identical(list(result$returned_count, result$row_count, result$truncated), list(1, 6, TRUE))
  expect_identical(unique(rows$full_sequence), "PEPTIDEK")
  expect_identical(rows$annotation, c("y1", "b2", "y2", "bIy[2-5]", "y3-18.0106", "y4-43.0058"))
  custom <- rows[6L, ]
  expect_true(abs(custom$neutral_loss - (-43.005814)) < 1e-6)
  # mzLib reads 0 back for a custom loss's position.
  expect_true(is.na(custom$residue_position))
  expect_identical(list(rows$secondary_product_type[[4L]], rows$secondary_fragment_number[[4L]]), list("y", 5))
  expect_true(all(is.na(rows$secondary_fragment_number[!rows$is_internal])))
  # No loss is 0, not NA.
  expect_identical(rows$neutral_loss[[1L]], 0)
  expect_identical(max(rows$fragment_intensity), 1)
  expect_identical(unique(rows$fragment_count), 6)
})

test_that("an index-only read says so and gives the same rows", {
  index_only <- mz$library_parse_records(recorded_data("library_read_index_only.json"))
  full <- library_synthetic()$records
  expect_true(index_only$is_index_only)
  expect_identical(index_only$index_only_threshold, 0)
  rownames(full) <- NULL
  expect_identical(index_only$records, full[1:2, , drop = FALSE])
})

test_that("a batch projects one row of files per library", {
  batch <- library_call("library_read_many.json", function() {
    library_read_many(c("synthetic_library.msl", "myPrositLib.msl"))
  })$result
  expect_true(inherits(batch, "mzlibr_library_batch"))
  expect_true(inherits(batch, "mzlibr_batch"))
  expect_identical(c(batch$file_count, batch$read_count, batch$failed_count), c(2, 2, 0))
  expect_identical(c(batch$record_count, batch$returned_count, batch$row_count), c(11, 11, 11))
  expect_identical(batch$column_names[1:2], c("source_index", "source_path"))
  # 1-based, so it indexes files directly.
  expect_identical(batch$records$source_index, c(rep(1, 6), rep(2, 5)))
  expect_identical(batch$records$precursor_index, c(1:6, 1:5) + 0)
  expect_identical(batch$files$record_count, c(6, 5))
  expect_identical(batch$files$reader, c("MslLibrary", "MslLibrary"))
  expect_identical(batch$files$proteoform_count[[1L]], 1)
  expect_identical(batch$files$absent_fields[[2L]], c("ion_mobility", "protein_accession", "protein_name", "gene_name"))
  expect_identical(nrow(mz$bulk_failed(batch)), 0L)
  expect_identical(batch$records$full_sequence[[7L]], "ALAVDGAGKPGAEE")
})

test_that("a skipped library keeps its unknown facts as NA", {
  files <- mz$bulk_parse_files(list(list(
    path = "C:/missing.msl", file_type = NA, reader = NA, format_version = NA, is_compressed = NA,
    is_index_only = NA, target_count = NA, decoy_count = NA, proteoform_count = NA,
    min_precursor_mz = NA, max_precursor_mz = NA, record_count = NA, rows_not_read = NA,
    retention_time_unit = NA, caveats = list(), column_names = list(), absent_fields = list(),
    failed_fields = list(), excluded_fields = list(),
    error = list(kind = "usage", type = "usage", message = "File not found: 'C:/missing.msl'.")
  )))
  expect_identical(files$error_kind, "usage")
  expect_true(is.na(files$record_count) && is.na(files$target_count) && is.na(files$is_index_only))
})

test_that("a table written to out returns no records", {
  wire <- recorded_data("library_read_synthetic.json")
  wire$returned_count <- 0
  wire$row_count <- 0
  wire$columns <- NA
  wire$output <- list(path = "C:/out/library.tsv", format = "tsv", row_count = 6)
  written <- mz$library_parse_records(wire)
  expect_true(is.null(written$records))
  expect_identical(list(written$output$path, written$output$row_count), list("C:/out/library.tsv", 6))
})

test_that("libraries replay end to end, as the help pages run them", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(library_read("synthetic_library.msl")$record_count, 6)
  expect_identical(library_read("synthetic_library.msl", fragments = TRUE, limit = 1)$row_count, 6)
  expect_true(library_read("synthetic_library.msl", index_only_threshold = 0, limit = 2)$is_index_only)
  expect_identical(library_read_many(c("synthetic_library.msl", "myPrositLib.msl"))$record_count, 11)
  # The fragments recording does not answer a call that did not ask for fragments.
  expect_error(library_read("synthetic_library.msl", limit = 1), class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

library_live <- function() {
  skip_unless_bridge_has("library read")
  root <- file.path(Sys.getenv("MZLIBR_SDRF_FIXTURES", ""), "library")
  skip_if(!file.exists(file.path(root, "synthetic_library.msl")),
    "no spectral libraries to read (set MZLIBR_SDRF_FIXTURES to pyMzLib's tests/fixtures)")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  list(synthetic = file.path(root, "synthetic_library.msl"), prosit = file.path(root, "myPrositLib.msl"))
}

# A parsed read without its machine-specific absolute paths, to compare with a recording.
library_without_paths <- function(x) {
  x$path <- NULL
  x$files$path <- NULL
  x$records$source_path <- NULL
  x
}

test_that("LIVE: the recordings are what the bridge says now", {
  at <- library_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  cases <- list(
    list("library_read_prosit.json", function() library_read(at$prosit)),
    list("library_read_synthetic.json", function() library_read(at$synthetic)),
    list("library_read_fragments.json", function() library_read(at$synthetic, fragments = TRUE, limit = 1)),
    list("library_read_index_only.json", function() library_read(at$synthetic, index_only_threshold = 0, limit = 2))
  )
  for (case in cases) {
    recorded <- mz$library_parse_records(recorded_data(case[[1L]]))
    expect_identical(library_without_paths(case[[2L]]()), library_without_paths(recorded), info = case[[1L]])
  }
  live <- library_read_many(c(at$synthetic, at$prosit))
  recorded <- mz$library_parse_batch(recorded_data("library_read_many.json"))
  expect_identical(library_without_paths(live), library_without_paths(recorded))
})

test_that("LIVE: a batch is identical at one and four threads", {
  at <- library_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  one <- library_read_many(c(at$synthetic, at$prosit), fragments = TRUE, threads = 1)
  four <- library_read_many(c(at$synthetic, at$prosit), fragments = TRUE, threads = 4)
  expect_identical(one, four)
  expect_true(one$row_count > one$returned_count)
  expect_identical(one$returned_count, 11)
})

test_that("LIVE: index-only and full loads give the same fragment rows", {
  at <- library_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  full <- library_read(at$synthetic, fragments = TRUE)
  index_only <- library_read(at$synthetic, fragments = TRUE, index_only_threshold = 0)
  expect_identical(c(full$is_index_only, index_only$is_index_only), c(FALSE, TRUE))
  expect_identical(index_only$records, full$records)
})

test_that("LIVE: out writes the table", {
  at <- library_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  target <- file.path(tempfile("library-"), "nested", "ions.tsv")
  on.exit(unlink(dirname(dirname(target)), recursive = TRUE), add = TRUE)
  written <- library_read(at$synthetic, fragments = TRUE, out = target)
  lines <- readLines(target, warn = FALSE, encoding = "UTF-8")
  expect_true(is.null(written$records))
  expect_identical(strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]], written$column_names)
  expect_identical(length(lines) - 1, written$output$row_count)
  expect_identical(written$output$row_count, 14)
})

test_that("LIVE: a text library is a usage error, and a damaged one mzLib's correctness error", {
  at <- library_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  folder <- tempfile("library-bad-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  text <- file.path(folder, "library.msp")
  writeLines("Name: PEPTIDE/2", text)
  expect_error(library_read(text), class = "mzlib_usage_error", contains = ".msl")
  damaged <- file.path(folder, "damaged.msl")
  bytes <- readBin(at$synthetic, "raw", file.size(at$synthetic))
  # Inside the data the CRC-32 covers.
  bytes[[201L]] <- xor(bytes[[201L]], as.raw(0xFF))
  writeBin(bytes, damaged)
  condition <- expect_error(library_read(damaged), class = "mzlib_bridge_error")
  expect_false(inherits(condition, "mzlib_usage_error"))
  expect_true(condition$error_type %in% c("InvalidDataException", "FormatException"))
})

test_that("LIVE: a skipped library is reported and the rest read", {
  at <- library_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  batch <- library_read_many(c(at$synthetic, file.path(tempdir(), "missing.msl")), on_error = "skip")
  expect_identical(c(batch$read_count, batch$failed_count), c(1, 1))
  failed <- mz$bulk_failed(batch)
  expect_identical(failed$error_kind, "usage")
  expect_true(is.na(failed$record_count) && is.na(failed$target_count))
})
