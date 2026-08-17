# Reading result files.
#
# The point of this module is to stop someone believing mzLib reads 31 formats into one uniform
# shape. Most of these tests assert that the exceptions are visible.

recorded_formats <- function() {
  payload <- mz$json_parse(paste(
    readLines(fixture_path("readers_formats.json"), warn = FALSE),
    collapse = "\n"
  ))
  mz$readers_parse_formats(payload$data)
}

# ---------------------------------------------------------------- formats

test_that("every format mzLib recognises becomes a row", {
  formats <- recorded_formats()
  expect_identical(nrow(formats), 31L)
  expect_true(is.data.frame(formats))
  expect_false(is.factor(formats$file_type))
})

test_that("exactly four file types are quantifiable", {
  # The precondition for flashlfq_quantify(), and it is a much smaller set than "31 formats"
  # suggests. The number is quoted in ?readers_formats, so it is pinned here.
  #
  # It moved from three to four when mzLib 1.0.585 added DiaNnReport (mzLib #1120), which is
  # how DIA data reaches flashlfq_quantify() at all. This test failing on a bridge repin is
  # the mechanism working: the claim changed, so the name and the number both moved.
  formats <- recorded_formats()
  quantifiable <- formats[formats$is_quantifiable, ]
  expect_identical(nrow(quantifiable), 4L)
  expect_identical(sort(quantifiable$file_type),
                   sort(c("psmtsv", "osmtsv", "MsFraggerPsm", "DiaNnReport")))
})

test_that("most formats have no views at all, and that is a real answer", {
  # 14 of 31. An empty views list means mzLib can parse the file but offers no cross-format
  # projection of it - not that anything failed.
  formats <- recorded_formats()
  viewless <- vapply(formats$views, function(v) length(v) == 0L, logical(1L))
  expect_identical(sum(viewless), 14L)
})

test_that("the view vocabulary is the four documented families", {
  formats <- recorded_formats()
  expect_identical(
    sort(unique(unlist(formats$views))),
    sort(c("ms1_features", "quantifiable", "spectra", "spectral_match"))
  )
})

test_that("an extension does not identify a format on its own", {
  # BrukerD and BrukerTimsTof are both ".d", told apart by what the directory holds, and several
  # formats share ".tsv". Anyone dispatching on extension alone is wrong.
  formats <- recorded_formats()
  expect_true(anyDuplicated(formats$extension) > 0L)
})

test_that("views survive subsetting the formats frame", {
  # A list column for the same reason PRIDE's locations is one.
  formats <- recorded_formats()
  some <- formats[formats$is_quantifiable, ]
  expect_true("views" %in% names(some))
  expect_true(all(vapply(some$views, function(v) "quantifiable" %in% v, logical(1L))))
})

# ---------------------------------------------------------------- identify

test_that("a file with no views says so rather than looking broken", {
  info <- mz$readers_parse_file_info(list(
    path = "/abs/spectra.mzML", file_type = "MzML", extension = ".mzML",
    reader = "MsDataFile", views = list("spectra")
  ))
  expect_true(inherits(info, "mzlibr_file_info"))
  expect_false(info$is_quantifiable)
  expect_identical(info$views, "spectra")

  bare <- mz$readers_parse_file_info(list(
    path = "/abs/x.tsv", file_type = "CruxResult", extension = ".tsv",
    reader = "CruxResultFile", views = list()
  ))
  expect_identical(bare$views, character(0))
  expect_false(bare$is_quantifiable)
  output <- paste(capture.output(print(bare)), collapse = "\n")
  expect_true(grepl("no views", output, fixed = TRUE), info = output)
})

test_that("printing a quantifiable file still warns about trusting it", {
  # is_quantifiable reports what mzLib's interface offers, not that the numbers are comparable.
  # MSFragger is quantifiable by interface and should not be quantified.
  info <- mz$readers_parse_file_info(list(
    path = "/abs/psm.tsv", file_type = "MsFraggerPsm", extension = "psm.tsv",
    reader = "MsFraggerPsmFile", views = list("quantifiable")
  ))
  expect_true(info$is_quantifiable)
  output <- paste(capture.output(print(info)), collapse = "\n")
  expect_true(grepl("before trusting the numbers", output, fixed = TRUE), info = output)
})

test_that("a blank path is refused before anything is spawned", {
  for (path in list("", "   ", NA_character_, NULL, 42)) {
    expect_error(mz$readers_normalise_path(path), class = "mzlib_usage_error")
  }
})

# ---------------------------------------------------------------- read_results

# Deliberately not modifyList(): it recurses into list-valued fields, and for an *unnamed* list
# it iterates over names(value), which is NULL — so the override is silently ignored and the
# base value survives. Four tests here passed against data they had not actually set before that
# was noticed. Replacing wholesale is the only behaviour that means what it looks like.
sample_records <- function(...) {
  base <- list(
    path = "/abs/AllPSMs.psmtsv", file_type = "psmtsv",
    record_count = 4, returned_count = 2, offset = 0, truncated = TRUE,
    retention_time_unit = "minutes", rows_not_read = 1,
    caveats = list("MetaMorpheus retention times are minutes."),
    column_names = list("sequence", "retention_time", "charge"),
    columns = list(
      sequence = list("PEPTIDEK", "ACDEFR"),
      retention_time = list(30.1, NA),
      charge = list(2, 3)
    ),
    output = NULL
  )
  overrides <- list(...)
  for (name in names(overrides)) {
    base[[name]] <- overrides[[name]]
  }
  base
}

test_that("the record columns become a data.frame directly", {
  # pyMzLib returns field -> values because that is what pandas wants; in R the equivalent IS a
  # data.frame, so ResultRecords.records - Python's second, row-wise accessor - has no
  # counterpart here. One object is both.
  records <- mz$readers_parse_records(sample_records())
  expect_true(is.data.frame(records$records))
  expect_identical(nrow(records$records), 2L)
  expect_identical(names(records$records), c("sequence", "retention_time", "charge"))
  expect_identical(records$records$sequence, c("PEPTIDEK", "ACDEFR"))
})

test_that("a null inside a column becomes NA and does not shorten it", {
  # The hazard that runs through this whole package. A shortened column would silently
  # misalign every other column in the table.
  records <- mz$readers_parse_records(sample_records())
  expect_identical(length(records$records$retention_time), 2L)
  expect_true(is.na(records$records$retention_time[2]))
})

test_that("columns of unequal length are refused rather than recycled", {
  # R would recycle silently and produce a table that is wrong rather than absent.
  broken <- sample_records(columns = list(
    sequence = list("A", "B"),
    retention_time = list(1)
  ), column_names = list("sequence", "retention_time"))
  expect_error(mz$readers_parse_records(broken),
    class = "mzlib_protocol_error", contains = "different lengths"
  )
})

test_that("column order follows column_names, not the map's own order", {
  records <- mz$readers_parse_records(sample_records(
    column_names = list("charge", "sequence", "retention_time")
  ))
  expect_identical(names(records$records), c("charge", "sequence", "retention_time"))
})

test_that("truncation and unreadable rows are both surfaced", {
  # A short answer and a complete one must not look alike, and mzLib drops a malformed row
  # silently - so rows_not_read is the only sign that a table is incomplete.
  records <- mz$readers_parse_records(sample_records())
  expect_true(records$truncated)
  expect_identical(records$rows_not_read, 1)
  expect_identical(records$record_count, 4)
  expect_identical(records$returned_count, 2)

  output <- paste(capture.output(print(records)), collapse = "\n")
  expect_true(grepl("truncated", output, fixed = TRUE), info = output)
  expect_true(grepl("partly unreadable", output, fixed = TRUE), info = output)
})

test_that("caveats are printed, because that is where the units live", {
  records <- mz$readers_parse_records(sample_records())
  output <- paste(capture.output(print(records)), collapse = "\n")
  expect_true(grepl("retention times are minutes", output, fixed = TRUE), info = output)
})

test_that("a written table reports where it went and in what format", {
  records <- mz$readers_parse_records(sample_records(
    columns = NULL, column_names = list(), returned_count = 0,
    output = list(path = "/abs/out.tsv", format = "tsv", row_count = 4)
  ))
  expect_true(is.null(records$records))
  expect_identical(records$output$format, "tsv")
  expect_identical(records$output$row_count, 4)
  output <- paste(capture.output(print(records)), collapse = "\n")
  expect_true(grepl("written to /abs/out.tsv", output, fixed = TRUE), info = output)
})

# ---------------------------------------------------------------- retention time units

test_that("minutes pass through and seconds are converted", {
  minutes <- mz$readers_parse_records(sample_records())
  expect_identical(mz$readers_retention_time_in_minutes(minutes)[1], 30.1)

  seconds <- mz$readers_parse_records(sample_records(
    file_type = "MsFraggerPsm", retention_time_unit = "seconds",
    columns = list(
      sequence = list("A", "B"), retention_time = list(1800, 3600), charge = list(2, 2)
    )
  ))
  expect_equal(mz$readers_retention_time_in_minutes(seconds), c(30, 60))
})

test_that("an unknown unit raises rather than guessing", {
  # A silently unconverted time axis is the specific mistake this module exists to prevent: a
  # 60x error in a retention-time comparison reads as a chromatography problem, not a units one.
  unknown <- mz$readers_parse_records(sample_records(retention_time_unit = "unknown"))
  expect_error(mz$readers_retention_time_in_minutes(unknown),
    class = "mzlib_usage_error", contains = "no basis to say what unit"
  )
})

test_that("the converter refuses anything that is not a record view", {
  expect_error(mz$readers_retention_time_in_minutes(list()), class = "mzlib_usage_error")
})

# ---------------------------------------------------------------- argument assembly

test_that("read-results defaults send neither limit nor offset", {
  args <- mz$readers_build_read_args("psm.tsv", NULL, 0, NULL)
  expect_identical(args, c("readers", "read-results", "--path", "psm.tsv"))
})

test_that("limit, offset and out reach the wire in the bridge's spelling", {
  args <- mz$readers_build_read_args("psm.tsv", 100, 50, "out.tsv")
  expect_identical(args[which(args == "--limit") + 1L], "100")
  expect_identical(args[which(args == "--offset") + 1L], "50")
  expect_identical(args[which(args == "--out") + 1L], "out.tsv")
})

test_that("an impossible limit or offset is refused", {
  for (bad in list(0, -1, 1.5, "10")) {
    expect_error(mz$readers_build_read_args("psm.tsv", bad, 0, NULL), class = "mzlib_usage_error")
  }
  for (bad in list(-1, 1.5, "10", NA_real_)) {
    expect_error(mz$readers_build_read_args("psm.tsv", NULL, bad, NULL), class = "mzlib_usage_error")
  }
})

test_that("a large limit is not written in scientific notation", {
  args <- mz$readers_build_read_args("psm.tsv", 1000000, 0, NULL)
  expect_identical(args[which(args == "--limit") + 1L], "1000000")
})

# ---------------------------------------------------------------- against a real mzLib

test_that("LIVE: mzLib still recognises 31 formats, four of them quantifiable", {
  # Enumerated from mzLib itself, so this is the test that notices when the installed version
  # changes what it supports - which is exactly when the numbers in ?readers_formats go stale.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  formats <- readers_formats()
  expect_identical(nrow(formats), 31L)
  expect_identical(sum(formats$is_quantifiable), 4L)
  expect_identical(sum(vapply(formats$views, length, integer(1L)) == 0L), 14L)
})

test_that("LIVE: identify dispatches on extension and does not validate contents", {
  # Worth knowing and easy to get wrong. A plain text file containing nothing proteomic is
  # identified as `CruxResult`, because mzLib dispatches `.txt` there and identify() is
  # deliberately cheap - it resolves the type and stops. It is not a validity check, and a
  # confident-looking file_type is not evidence the file is what it claims.
  #
  # The honest signal is `views`, which is empty here. That is what to branch on.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  scratch <- tempfile("mzlibr-not-a-result-", fileext = ".txt")
  writeLines("this is not a proteomics result file", scratch)
  on.exit(unlink(scratch), add = TRUE)

  info <- readers_identify(scratch)
  expect_identical(info$file_type, "CruxResult")
  expect_identical(info$views, character(0))
  expect_false(info$is_quantifiable)
})

test_that("LIVE: reading a file with no quantifiable view names the views it does have", {
  # The error is the documentation a stuck user reads.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  scratch <- tempfile("mzlibr-spectra-", fileext = ".mzML")
  writeLines("<mzML></mzML>", scratch)
  on.exit(unlink(scratch), add = TRUE)

  condition <- tryCatch(
    {
      readers_read_results(scratch, timeout = 120)
      NULL
    },
    mzlib_error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_error"))
})

# ---------------------------------------------------------------- exhaustive coverage
#
# `readers_read_results()` reaches four of the 31 file types. These four verbs reach the rest. The
# payloads are recorded from the real bridge against real mzLib fixtures, so a wire-shape change
# shows up here as a parse failure rather than as a fixture that agrees with an R file and with
# nothing else.

recorded_payload <- function(name) {
  mz$json_parse(paste(
    readLines(fixture_path(name), warn = FALSE),
    collapse = "\n"
  ))$data
}

recorded_native <- function() {
  mz$readers_parse_native_records(recorded_payload("readers_records_toppic.json"))
}

recorded_features <- function() {
  mz$readers_parse_feature_records(recorded_payload("readers_features_topfd.json"))
}

recorded_matches <- function() {
  mz$readers_parse_match_records(recorded_payload("readers_matches_casanovo.json"))
}

recorded_scans <- function() {
  mz$readers_parse_scan_records(recorded_payload("readers_spectra_mzml.json"))
}

# ---------------------------------------------------------------- read_records

test_that("a format with no cross-format view still becomes a table", {
  # The whole point of the verb. TopPIC has no view at all and was unreadable before it.
  records <- recorded_native()
  expect_identical(records$file_type, "ToppicPrsm")
  expect_identical(records$record_type, "ToppicPrsm")
  expect_identical(length(records$views), 0L)
  expect_true(is.data.frame(records$records))
  expect_true("e_value" %in% names(records$records))
})

test_that("the native table has one row per returned record", {
  records <- recorded_native()
  expect_identical(nrow(records$records), as.integer(records$returned_count))
})

test_that("column names are mzLib's own, so they cross-reference the source", {
  # EValue -> e_value, MIScore -> mi_score, FixedPTMs -> fixed_ptms. A pluralising 's' belongs to
  # the acronym before it; getting that wrong would give fixed_pt_ms.
  records <- recorded_native()
  expect_true("e_value" %in% records$column_names)
  expect_true("mi_score" %in% records$column_names)
  expect_true("fixed_ptms" %in% records$column_names)
})

test_that("a field that could not become a column is named with its reason", {
  # A column that simply vanished is indistinguishable from a field the format does not have.
  records <- recorded_native()
  expect_true(is.data.frame(records$excluded_fields))
  expect_true("alternative_identifications" %in% records$excluded_fields$field)
  expect_true(all(nzchar(records$excluded_fields$reason)))
})

test_that("excluded_fields is a frame with the right columns even when empty", {
  # So `nrow(x$excluded_fields)` works without a NULL check.
  empty <- mz$readers_parse_native_records(list())
  expect_true(is.data.frame(empty$excluded_fields))
  expect_identical(nrow(empty$excluded_fields), 0L)
  expect_identical(sort(names(empty$excluded_fields)), sort(c("field", "type", "reason")))
})

test_that("read_records sends its own verb", {
  args <- mz$readers_build_read_args("a.tsv", NULL, 0, NULL, "read-records")
  expect_identical(args[1:2], c("readers", "read-records"))
})

test_that("a native table cannot be retention-time converted", {
  # Its columns are the format's own and carry no declared unit, so a conversion would be a guess
  # dressed as a conversion.
  records <- recorded_native()
  expect_error(mz$readers_retention_time_in_minutes(records), contains = "no declared")
})

# ---------------------------------------------------------------- read_features

test_that("the feature view has the six interface columns", {
  features <- recorded_features()
  expect_identical(
    features$column_names,
    c(
      "mz", "charge", "retention_time_start", "retention_time_end",
      "intensity", "number_of_isotopes"
    )
  )
})

test_that("the _ms1.feature retention-time unit is unknown, not guessed", {
  # TopFD wrote seconds through v1.6.2 and minutes from v1.7.0 without changing the file type.
  features <- recorded_features()
  expect_identical(features$retention_time_unit, "unknown")
})

test_that("converting an unknown retention-time unit raises rather than guessing", {
  features <- recorded_features()
  expect_error(mz$readers_retention_time_in_minutes(features), contains = "no basis to say")
})

test_that("the feature view says its rows are charge states, not file lines", {
  # record_count exceeds the file's line count for this format, and a caller comparing the two
  # must be told why.
  features <- recorded_features()
  expect_true(any(grepl("CHARGE STATE", features$caveats, fixed = TRUE)))
})

test_that("number_of_isotopes is NA for _ms1.feature rather than zero", {
  # mzLib's single-charge expansion never sets it. Zero would read as "no isotopes found".
  features <- recorded_features()
  expect_true(all(is.na(features$records$number_of_isotopes)))
})

# ---------------------------------------------------------------- read_matches

test_that("the spectral-match view carries identity fields and modifications", {
  matches <- recorded_matches()
  expect_true("accession" %in% matches$column_names)
  expect_true("modifications" %in% matches$column_names)
  expect_true("modification_count" %in% matches$column_names)
})

test_that("Casanovo is_decoy is NA, not FALSE", {
  # De novo sequencing has no target/decoy label at all, so FALSE would let someone filter on a
  # fabricated column - the trap readers_read_results() already refuses for MSFragger.
  matches <- recorded_matches()
  expect_true(all(is.na(matches$records$is_decoy)))
})

test_that("the spectral-match view says nothing in it is FDR-filtered", {
  matches <- recorded_matches()
  expect_true(any(grepl("FDR", matches$caveats, fixed = TRUE)))
})

test_that("the Casanovo scan-number caveat is present", {
  # It is an mzTab index, not an instrument scan number, and joining on it is wrong.
  # Matched case-insensitively: the caveat capitalises INDEX for emphasis, and a fixed = TRUE
  # match on the lowercase form silently never fires.
  matches <- recorded_matches()
  expect_true(any(grepl("index", matches$caveats, ignore.case = TRUE)))
})

# ---------------------------------------------------------------- read_spectra

test_that("scan headers parse into a table with minutes", {
  scans <- recorded_scans()
  expect_true("one_based_scan_number" %in% names(scans$records))
  expect_identical(scans$retention_time_unit, "minutes")
})

test_that("the file's total scan count is reported alongside the filtered count", {
  # So an ms_order filter that matched nothing can never look like an empty file.
  scans <- recorded_scans()
  expect_true(scans$scan_count >= scans$record_count)
})

test_that("peaks are absent unless asked for", {
  scans <- recorded_scans()
  expect_false(scans$peaks_included)
  expect_false("mz" %in% names(scans$records))
  # ...but you are still told how many peaks each scan has.
  expect_true("peak_count" %in% names(scans$records))
})

test_that("read_spectra builds its own verb and window options", {
  built <- mz$readers_build_read_args("run.mzML", 10, 0, NULL, "read-spectra")
  expect_identical(built[1:2], c("readers", "read-spectra"))
  expect_true("--limit" %in% built)
})

test_that("a bad ms_order is refused before the bridge is spawned", {
  for (bad in list(0, -1, 1.5, "2", NA_real_)) {
    expect_error(
      mz$readers_read_spectra("run.mzML", ms_order = bad),
      contains = "ms_order"
    )
  }
})

test_that("a non-logical peaks argument is refused", {
  expect_error(mz$readers_read_spectra("run.mzML", peaks = "yes"), contains = "peaks must be")
  expect_error(mz$readers_read_spectra("run.mzML", peaks = NA), contains = "peaks must be")
})

# ---------------------------------------------------------------- shared argument rules

test_that("every read verb refuses a blank path before spawning anything", {
  for (verb in c("read-records", "read-features", "read-matches", "read-spectra")) {
    expect_error(mz$readers_build_read_args("", NULL, 0, NULL, verb), contains = "non-empty file path")
    expect_error(mz$readers_build_read_args(NA_character_, NULL, 0, NULL, verb), contains = "non-empty")
  }
})

test_that("every read verb refuses a bad limit and a bad offset", {
  for (verb in c("read-records", "read-features", "read-matches", "read-spectra")) {
    expect_error(mz$readers_build_read_args("a.tsv", 0, 0, NULL, verb), contains = "positive whole number")
    expect_error(mz$readers_build_read_args("a.tsv", 1.5, 0, NULL, verb), contains = "positive whole number")
    expect_error(mz$readers_build_read_args("a.tsv", NULL, -1, NULL, verb), contains = "non-negative")
  }
})

test_that("a zero offset is not sent", {
  # A default that is sent explicitly is a default the bridge can later disagree with.
  args <- mz$readers_build_read_args("a.tsv", NULL, 0, NULL, "read-records")
  expect_false("--offset" %in% args)
})

# ---------------------------------------------------------------- printing

test_that("every result type prints without erroring", {
  for (record in list(recorded_native(), recorded_features(), recorded_matches(), recorded_scans())) {
    output <- utils::capture.output(print(record))
    expect_true(length(output) > 0L)
  }
})

test_that("the native print names what could not be projected", {
  output <- paste(utils::capture.output(print(recorded_native())), collapse = "\n")
  expect_true(grepl("could not become columns", output, fixed = TRUE))
  expect_true(grepl("ToppicPrsm", output, fixed = TRUE))
})

test_that("the scan print says peaks were omitted", {
  output <- paste(utils::capture.output(print(recorded_scans())), collapse = "\n")
  expect_true(grepl("peaks not included", output, fixed = TRUE))
})

test_that("the feature print shows the unknown unit rather than hiding it", {
  output <- paste(utils::capture.output(print(recorded_features())), collapse = "\n")
  expect_true(grepl("retention_time_unit: unknown", output, fixed = TRUE))
})

test_that("per-scan peak arrays become a list column, not a length error", {
  # peaks = TRUE returns one mz array and one intensity array PER SCAN. Unlisting those would
  # splice every scan's peaks into one vector and the length check would then reject the table
  # with a message about column lengths - a shape complaint about the correct shape.
  scans <- mz$readers_parse_scan_records(recorded_payload("readers_spectra_peaks.json"))

  expect_true(scans$peaks_included)
  expect_true("mz" %in% names(scans$records))
  expect_identical(nrow(scans$records), as.integer(scans$returned_count))
  # One array per scan, each as long as that scan's reported peak count.
  expect_identical(
    vapply(scans$records$mz, length, integer(1L)),
    as.integer(scans$records$peak_count)
  )
})

test_that("a fabricated zero intensity is disclosed as fabricated", {
  # A within-type schema variant: Apex_intensity is optional and the FLASHDeconv/OpenMS
  # _ms1.feature layout omits it, so mzLib substitutes zero for every feature. A whole column of
  # zeros is indistinguishable from real measurements of nothing.
  #
  # The bridge passes mzLib's value through and SAYS the zero is fabricated. Crossing it as NA
  # makes the wire disagree with mzLib about a number, which ships separately (pyMzLib #28); when
  # that lands, this flips to asserting NA and this package gains its parity port.
  features <- mz$readers_parse_feature_records(recorded_payload("readers_features_flashdeconv.json"))

  expect_true(all(features$records$intensity == 0))
  expect_true(any(grepl("FABRICATED", features$caveats, fixed = TRUE)))
})

test_that("a TopFD feature file carries no fabrication caveat", {
  # The counterpart, and the fixture that proves the caveat above is conditional: TopFD writes
  # Apex_intensity, so its intensities are real and nothing is claimed about them.
  features <- recorded_features()
  expect_true(all(features$records$intensity > 0))
  expect_false(any(grepl("FABRICATED", features$caveats, fixed = TRUE)))
})

test_that("the Casanovo modification caveat describes what mzLib actually does", {
  # An earlier version said modifications were not loaded "because mzLib's file factory does not
  # enable it". The factory does - the parameterless constructor chains to this(TRUE) - and the
  # recorded payload carries populated full sequences alongside the caveat that denied them.
  matches <- recorded_matches()
  expect_false(any(grepl("not loaded", matches$caveats, fixed = TRUE)))
  expect_true(any(grepl("mass shifts", matches$caveats, fixed = TRUE)))
})
