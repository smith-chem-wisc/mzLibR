# Many files in one call, the per-file block, and the mzLib 1.0.592 reader fields.
#
# The `_many` functions are a projection of the bridge's bulk verbs (its BULK.md): the whole list
# goes to ONE bridge process on stdin, with the thread count stated, and comes back as one long
# table plus one row of per-file facts per input. These tests pin that shape against pyMzLib's
# recordings, and pin the argument checks that happen before anything is spawned.

recorded_data <- function(name) {
  parsed <- mz$json_parse(paste(
    readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  ))
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}

# ---------------------------------------------------------------- the request

test_that("a list of paths goes on stdin, one per line, with threads and on_error explicit", {
  request <- mz$bulk_request(c("readers", "read-spectra"), c("a.mzML", " b.mzML "), 2, "skip", NULL, "f")
  expect_identical(request$stdin, "a.mzML\nb.mzML\n")
  expect_identical(
    request$args,
    c("readers", "read-spectra", "--paths-stdin", "--threads", "2", "--on-error", "skip")
  )
  expect_false("--path" %in% request$args)
})

test_that("threads defaults to 1 and -1 means one per core", {
  expect_identical(mz$bulk_check_threads(1), "1")
  expect_identical(mz$bulk_check_threads(-1), "-1")
  for (bad in list(0, -2, 1.5, "2", NA_real_, c(1, 2))) {
    expect_error(mz$bulk_check_threads(bad), class = "mzlib_usage_error", contains = "threads")
  }
})

test_that("on_error is fail or skip and nothing else", {
  expect_identical(mz$bulk_check_on_error("fail"), "fail")
  for (bad in list("ignore", NA_character_, TRUE, c("fail", "skip"))) {
    expect_error(mz$bulk_check_on_error(bad), class = "mzlib_usage_error", contains = "on_error")
  }
})

test_that("an unusable list is refused before anything is spawned", {
  expect_error(mz$bulk_paths_stdin(character(0), "f"), contains = "at least one path")
  expect_error(mz$bulk_paths_stdin(list("a"), "f"), contains = "character vector")
  expect_error(mz$bulk_paths_stdin(c("a", NA), "f"), contains = "Path 2")
  expect_error(mz$bulk_paths_stdin(c("a", "  "), "f"), contains = "Path 2")
  # A line break inside a path would split it in two on the wire.
  expect_error(mz$bulk_paths_stdin("a\nb", "f"), contains = "line break")
})

test_that("every _many function takes the options its single form does, plus the list", {
  expect_true(all(c("paths", "threads", "on_error", "out") %in% names(formals(readers_read_spectra_many))))
  expect_true(all(c("ms_order", "peaks") %in% names(formals(readers_read_spectra_many))))
  expect_true("scores" %in% names(formals(readers_read_matches_many)))
  expect_identical(formals(readers_read_records_many)$threads, 1)
  expect_identical(formals(readers_read_records_many)$on_error, "fail")
  # Windowing belongs to one file (BULK.md section 1), so no _many form takes limit or offset.
  for (fn in list(readers_read_spectra_many, readers_read_records_many, readers_read_results_many)) {
    expect_false(any(c("limit", "offset") %in% names(formals(fn))))
  }
})

# ---------------------------------------------------------------- the batch

spectra_batch <- function() {
  mz$bulk_parse_batch(recorded_data("readers_many_spectra.json"), "readers read-spectra", "mzlibr_read_batch")
}

test_that("a batch has one long table, and one row of facts per input in input order", {
  batch <- spectra_batch()
  expect_true(inherits(batch, "mzlibr_read_batch"))
  expect_true(inherits(batch, "mzlibr_batch"))
  expect_identical(batch$file_count, 3)
  expect_identical(batch$read_count, 2)
  expect_identical(batch$failed_count, 1)
  expect_identical(nrow(batch$files), 3L)
  expect_identical(basename(gsub("\\\\", "/", batch$files$path)), c("sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"))
  expect_identical(names(batch$records)[1:2], c("source_index", "source_path"))
  expect_identical(nrow(batch$records), as.integer(batch$row_count))
})

test_that("source_index is 1-based, so it indexes files directly", {
  batch <- spectra_batch()
  expect_true(all(batch$records$source_index %in% c(1, 3)))
  expect_identical(
    unname(batch$files$path[batch$records$source_index]),
    unname(batch$records$source_path)
  )
  # Rows are grouped by input, in input order.
  expect_false(is.unsorted(batch$records$source_index))
})

test_that("a failed input is a row with its error, and the rest were read", {
  batch <- spectra_batch()
  failed <- mz$bulk_failed(batch)
  expect_identical(nrow(failed), 1L)
  expect_identical(failed$error_kind, "usage")
  expect_true(grepl("not found", failed$error_message))
  # A fact the failed input could not establish is NA, not zero.
  expect_true(is.na(batch$files$record_count[2]))
  expect_identical(batch$files$record_count[c(1, 3)], c(6, 2))
})

test_that("each spectra file's source is flattened into columns of files", {
  batch <- spectra_batch()
  expect_true(all(c("instrument_model", "acquisition_start_time_is_utc") %in% names(batch$files)))
  # MGF records none of it; the failed input has none either.
  expect_true(is.na(batch$files$instrument_model[3]))
  expect_true(is.list(batch$files$caveats))
  expect_true(is.list(batch$files$absent_fields))
})

test_that("a batch prints its inputs, rows and failures", {
  output <- paste(utils::capture.output(print(spectra_batch())), collapse = "\n")
  expect_true(grepl("3 inputs: 2 read, 1 failed", output, fixed = TRUE))
  expect_true(grepl("no-such-run.mzML", output, fixed = TRUE))
})

test_that("a many-file native read keeps each file's rows under its source", {
  batch <- mz$bulk_parse_batch(recorded_data("readers_many_records_mzid.json"), "readers read-records", "mzlibr_read_batch")
  expect_identical(batch$record_count, 24)
  expect_identical(as.numeric(table(batch$records$source_index)), c(12, 12))
  expect_true("q_value" %in% names(batch$records))
})

test_that("identify_many answers per path, with is_quantifiable", {
  batch <- mz$bulk_parse_batch(recorded_data("readers_many_identify.json"), "readers identify", "mzlibr_identify_batch")
  expect_identical(nrow(batch$files), 3L)
  expect_identical(batch$files$file_type[1], "MzIdentML")
  expect_true(is.list(batch$files$views))
  expect_identical(batch$files$error_kind[2], "usage")
})

test_that("a batch converts retention times by each file's own unit", {
  batch <- structure(list(
    verb = "readers read-results",
    files = data.frame(path = c("a", "b"), retention_time_unit = c("minutes", "seconds"), stringsAsFactors = FALSE),
    records = data.frame(source_index = c(1, 1, 2), source_path = c("a", "a", "b"), retention_time = c(1, 2, 120))
  ), class = c("mzlibr_read_batch", "mzlibr_batch"))
  expect_equal(readers_retention_time_in_minutes(batch), c(1, 2, 2))

  batch$files$retention_time_unit[2] <- "unknown"
  expect_error(readers_retention_time_in_minutes(batch), class = "mzlib_usage_error", contains = "'b'")
})

# ---------------------------------------------------------------- the per-file block, one file

test_that("a spectra read carries the run's source", {
  scans <- mz$readers_parse_scan_records(recorded_data("readers_spectra_mzml.json"))
  expect_true(is.list(scans$source))
  expect_true(all(c(
    "instrument_model", "instrument_model_accession", "instrument_serial_number",
    "acquisition_start_time", "acquisition_start_time_is_utc"
  ) %in% names(scans$source)))
  expect_true(is.logical(scans$source$acquisition_start_time_is_utc))
  expect_identical(scans$absent_fields, character(0))
})

test_that("mzIdentML matches carry a q-value, rank and threshold, and say what was skipped", {
  matches <- mz$readers_parse_match_records(recorded_data("readers_matches_mzid.json"))
  expect_true(all(c("q_value", "rank", "pass_threshold") %in% names(matches$records)))
  expect_identical(matches$skipped_count, 0)
  expect_true(is.data.frame(matches$skipped))
  expect_identical(nrow(matches$skipped), 0L)
  # mzIdentML's isDecoy defaults to false when omitted, so the column has no source here.
  expect_true("is_decoy" %in% matches$absent_fields)
  expect_true(all(is.na(matches$records$is_decoy)))
})

test_that("scores make the table long by score, while limit still counts matches", {
  scored <- mz$readers_parse_match_records(recorded_data("readers_matches_mzid_scores.json"))
  expect_true(scored$scores_included)
  expect_identical(scored$returned_count, 1)
  expect_true(scored$row_count > scored$returned_count)
  expect_identical(nrow(scored$records), as.integer(scored$row_count))
  expect_true(all(c("match_index", "score_name", "score_value") %in% names(scored$records)))
})

test_that("a format with no skip list says so with NA, not zero", {
  denovo <- mz$readers_parse_match_records(recorded_data("readers_matches_casanovo.json"))
  expect_true(is.na(denovo$skipped_count))
  expect_true(is.null(denovo$skipped))
  expect_true(all(c("is_decoy", "q_value") %in% denovo$absent_fields))
})

test_that("scores must be TRUE or FALSE", {
  expect_identical(mz$readers_scores_args(TRUE), "--scores")
  expect_identical(mz$readers_scores_args(FALSE), character(0))
  expect_error(mz$readers_scores_args("yes"), class = "mzlib_usage_error")
})

test_that("absent fields are printed, so an NA column never reads as missing data", {
  denovo <- mz$readers_parse_match_records(recorded_data("readers_matches_casanovo.json"))
  output <- paste(utils::capture.output(print(denovo)), collapse = "\n")
  expect_true(grepl("absent from this file", output, fixed = TRUE))
})

# ---------------------------------------------------------------- quantification tables

test_that("protein groups arrive long: one row per group per sample group", {
  groups <- mz$readers_parse_quant_records(recorded_data("readers_protein_groups.json"), "mzlibr_protein_group_records")
  expect_identical(groups$returned_count, 1)
  expect_identical(length(groups$sample_labels), 18L)
  expect_identical(nrow(groups$records), 18L)
  expect_identical(groups$row_count, 18)
  expect_true(all(c("q_value", "sample_label", "spectral_count", "intensity") %in% names(groups$records)))
  # The occupancy cells have no column shape here; the verb that carries them is named.
  expect_true("readers read-occupancy" %in% groups$excluded_fields$verb)
})

test_that("a peptide table with no retention times says so in absent_fields", {
  peptides <- mz$readers_parse_quant_records(
    recorded_data("readers_quantified_peptides.json"), "mzlibr_quantified_peptide_records", "minutes"
  )
  expect_true(all(c("peak_order", "retention_time") %in% peptides$absent_fields))
  expect_true(all(is.na(peptides$records$retention_time)))
  expect_identical(peptides$retention_time_unit, "minutes")
})

test_that("occupancy is one row per site, with a count of truncated cells", {
  sites <- mz$readers_parse_quant_records(recorded_data("readers_occupancy.json"), "mzlibr_occupancy_records")
  expect_identical(sites$truncated_cell_count, 0)
  expect_identical(nrow(sites$records), as.integer(sites$row_count))
  expect_true(all(sites$records$basis %in% c("count", "intensity")))
  output <- paste(utils::capture.output(print(sites)), collapse = "\n")
  expect_true(grepl("mzlibr_occupancy_records", output, fixed = TRUE))
})

# ---------------------------------------------------------------- which verbs the bridge has

version_payload <- function(verbs = NULL) {
  paste0(
    '{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"10.0.8"',
    if (!is.null(verbs)) paste0(',"verbs":[', paste0('"', verbs, '"', collapse = ","), "]") else "",
    '},"error":null}'
  )
}

test_that("the version reports every verb the bridge dispatches, and NA from an older bridge", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  with_bridge_config(option = path, {
    listed <- mzlibr_bridge_version(runner = stub_runner(stdout = version_payload(c("version", "sdrf read")))$run)
    expect_identical(listed$verbs, c("version", "sdrf read"))
    older <- mzlibr_bridge_version(runner = stub_runner(stdout = version_payload())$run)
    expect_identical(older$verbs, NA_character_)
  })
})

test_that("a verb the bridge lacks is refused naming the bridge it needs, and asked once", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  with_bridge_config(option = path, {
    runner <- stub_runner(stdout = version_payload(c("version", "readers read-occupancy")))
    expect_no_warning(mz$bridge_require_verb("readers read-occupancy", "0.2.0", runner = runner$run))
    expect_error(
      mz$bridge_require_verb("readers read-protein-groups", "0.2.0", runner = runner$run),
      class = "mzlib_usage_error", contains = c("pyMzLib 0.2.0", "readers read-protein-groups")
    )
    # The answer is cached per bridge: a second question asks nothing.
    runner$seen$called <- FALSE
    try(mz$bridge_require_verb("readers read-protein-groups", "0.2.0", runner = runner$run), silent = TRUE)
    expect_false(runner$seen$called)
  })
})

test_that("a bridge that lists no verbs at all predates every checked verb", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  with_bridge_config(option = path, {
    expect_error(
      mz$bridge_require_verb("readers read-occupancy", "0.2.0", runner = stub_runner(stdout = version_payload())$run),
      class = "mzlib_usage_error", contains = "0.2.0"
    )
  })
})

test_that("a quantification read replays end to end, verb check included", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  groups <- readers_read_protein_groups("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv", limit = 1)
  expect_true(inherits(groups, "mzlibr_protein_group_records"))
  expect_identical(nrow(groups$records), 18L)
  runs <- readers_read_spectra_many(c("sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"), on_error = "skip")
  expect_identical(runs$failed_count, 1)
  mz$replay_bridge_stop(old)
})

# ---------------------------------------------------------------- against a real bridge

test_that("LIVE: many paths are identified in one call, a missing one skipped", {
  # Needs the bridge pyMzLib 0.2.0 publishes; an older one has no --paths-stdin.
  skip_unless_bridge_has("readers read-occupancy")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  real <- tempfile("mzlibr-bulk-", fileext = ".txt")
  writeLines("not a proteomics file", real)
  on.exit(unlink(real), add = TRUE)
  missing <- file.path(tempdir(), "mzlibr-no-such-file.mzML")

  found <- readers_identify_many(c(real, missing), on_error = "skip")
  expect_identical(found$file_count, 2)
  expect_identical(found$files$file_type[1], "CruxResult")
  expect_identical(found$files$error_kind[2], "usage")
})
