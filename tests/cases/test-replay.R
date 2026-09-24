# The stand-in bridge the help-page examples replay (R/replay.R).
#
# Its whole value is that it refuses a call its recording does not fit - otherwise an example
# could print output recorded for other arguments and nobody would notice. So most of these
# tests are refusals.

test_that("argv splits into the verb words and --name value options", {
  call <- mz$replay_parse_argv(c("readers", "read-spectra", "--path", "a.mzML", "--peaks", "--limit", "3"))
  expect_identical(call$verb, "readers read-spectra")
  expect_identical(call$options$path, "a.mzML")
  expect_true(isTRUE(call$options$peaks))
  expect_identical(call$options$limit, "3")
})

test_that("a path is compared by file name, whichever machine wrote it", {
  expect_identical(mz$replay_base("E:\\data\\run.mzML"), "run.mzML")
  expect_identical(mz$replay_base("/abs/run.mzML"), "run.mzML")
  data <- list(path = "E:\\data\\run.mzML", record_count = 6, returned_count = 6, offset = 0)
  expect_identical(mz$replay_mismatch(data, list(path = "run.mzML")), "")
  expect_true(grepl("not 'other.mzML'", mz$replay_mismatch(data, list(path = "other.mzML")), fixed = TRUE))
})

test_that("limit and offset must reproduce the recorded window", {
  data <- list(path = "run.mzML", record_count = 6, returned_count = 3, offset = 0)
  expect_identical(mz$replay_mismatch(data, list(path = "run.mzML", limit = "3")), "")
  expect_true(nzchar(mz$replay_mismatch(data, list(path = "run.mzML"))))
  expect_true(nzchar(mz$replay_mismatch(data, list(path = "run.mzML", limit = "4"))))
  expect_true(nzchar(mz$replay_mismatch(data, list(path = "run.mzML", limit = "3", offset = "1"))))
})

test_that("a filter or a flag the recording applied must be asked for", {
  filtered <- list(path = "run.mzML", ms_order = 2, record_count = 3, returned_count = 3, offset = 0)
  expect_true(grepl("ms-order", mz$replay_mismatch(filtered, list(path = "run.mzML")), fixed = TRUE))
  expect_identical(mz$replay_mismatch(filtered, list(path = "run.mzML", `ms-order` = "2")), "")

  peaks <- list(path = "run.mzML", peaks_included = TRUE, record_count = 1, returned_count = 1, offset = 0)
  expect_true(grepl("peaks", mz$replay_mismatch(peaks, list(path = "run.mzML")), fixed = TRUE))
  headers <- list(path = "run.mzML", peaks_included = FALSE, record_count = 1, returned_count = 1, offset = 0)
  expect_true(grepl("--peaks given", mz$replay_mismatch(headers, list(path = "run.mzML", peaks = TRUE)), fixed = TRUE))
})

test_that("a many-files call fits only a many-files recording, and the reverse", {
  bulk <- list(files = list(list(path = "/a.mzML"), list(path = "/b.mzML")), read_count = 2)
  one <- list(path = "a.mzML")
  expect_identical(mz$replay_mismatch(bulk, list(`paths-stdin` = TRUE)), "")
  expect_true(grepl("bulk", mz$replay_mismatch(bulk, list())))
  expect_true(grepl("one-document", mz$replay_mismatch(one, list(`paths-stdin` = TRUE))))
  # PRIDE's manifest also has a `files` list, of PRIDE files that carry no path: not bulk.
  pride <- list(accession = "PXD000001", files = list(list(file_name = "a.raw")))
  expect_identical(mz$replay_mismatch(pride, list(accession = "PXD000001")), "")
})

test_that("an echoed scalar option must equal the recording's", {
  data <- list(accession = "PXD000001")
  expect_identical(mz$replay_mismatch(data, list(accession = "PXD000001")), "")
  expect_true(grepl("accession=PXD000001", mz$replay_mismatch(data, list(accession = "PXD000002")), fixed = TRUE))
})

test_that("the staged recordings answer, and a call no recording fits is a usage error", {
  dir <- mz$replay_dir()
  skip_if(!nzchar(dir), "the package's replay recordings are not installed")
  answer <- mz$replay_answer(c("readers", "read-spectra", "--path", "sliced_ethcd.mzML", "--limit", "3"), dir)
  expect_true(startsWith(answer, "{\"ok\":true"))

  refused <- mz$json_parse(mz$replay_answer(c("readers", "read-spectra", "--path", "other.mzML"), dir))
  expect_false(refused$ok)
  expect_identical(refused$error$type, "usage")
  expect_true(grepl("readers_spectra_mzml.json", refused$error$message, fixed = TRUE))

  unknown <- mz$json_parse(mz$replay_answer(c("no", "such-verb"), dir))
  expect_true(grepl("no recording is staged for 'no such-verb'", unknown$error$message, fixed = TRUE))
})

test_that("an example's stand-in bridge runs end to end, and is removed afterwards", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  launcher <- getOption("mzlibr.bridge")
  on.exit(mz$replay_bridge_stop(old), add = TRUE)

  scans <- readers_read_spectra("sliced_ethcd.mzML", limit = 3)
  expect_identical(scans$returned_count, 3)
  expect_error(readers_read_spectra("sliced_ethcd.mzML"), class = "mzlib_usage_error", contains = "replay bridge")

  mz$replay_bridge_stop(old)
  expect_false(file.exists(launcher))
  expect_identical(getOption("mzlibr.bridge"), old$mzlibr.bridge)
})
