# Sequence conversion: peptidoform_convert().
#
# Ported from pyMzLib's test_peptidoform_convert.py. The payloads were recorded from the real
# bridge by pyMzLib (pyMzLib #75) and are shared verbatim:
#
#   * peptidoform_convert_unimod.json: four probe sequences to Unimod under ReturnNull - two UniProt
#     modifications, one MetaMorpheus one, and one mzLib cannot map;
#   * peptidoform_convert_psmtsv.json / ..._proforma.json: the full sequences of mzLib's own
#     BottomUpExample.psmtsv (readers_results_psmtsv.json) to Unimod and to ProForma.
#
# The accessions are mzLib's. What is pinned here is that mzLibR sends the call mzLib needs and
# projects mzLib's answer row for row, without repairing it.

convert_fixture <- function(name) {
  parsed <- mz$json_parse(paste(readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}
convert_recorded <- function(name) mz$peptidoform_parse_conversions(convert_fixture(name))
psmtsv_full_sequences <- function() {
  mz$wire_strings(convert_fixture("readers_results_psmtsv.json")$columns$full_sequence)
}

CONVERT_PROBE <- c(
  "[UniProt:N-acetylserine on S]SEQK",
  "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
  "PEPM[Common Variable:Oxidation on M]K",
  "PEPK[Made Up:Not a modification on K]R"
)

# ---------------------------------------------------------------- the request

test_that("the request sends mzLib's defaults and the sequences on stdin, unchanged", {
  request <- mz$peptidoform_convert_request(CONVERT_PROBE, "mzLib", "Unimod", "ReturnNull", 1)
  expect_identical(request$args, c(
    "peptidoform", "convert", "--from", "mzLib", "--to", "Unimod",
    "--mode", "ReturnNull", "--threads", "1"
  ))
  expect_identical(request$stdin, paste(CONVERT_PROBE, collapse = "\n"))
  expect_identical(mz$PEPTIDOFORM_CONVERT_SINCE, "0.4.0")
  expect_identical(formals(peptidoform_convert)$source, "mzLib")
  expect_identical(formals(peptidoform_convert)$target, "Unimod")
  expect_identical(formals(peptidoform_convert)$mode, "ReturnNull")
})

test_that("mode is sent in mzLib's spelling, and threads -1 passes through", {
  request <- mz$peptidoform_convert_request("PEPTIDE", "mzLib", "ProForma", "removeincompatibleelements", -1)
  expect_identical(utils::tail(request$args, 4L), c("--mode", "RemoveIncompatibleElements", "--threads", "-1"))
  expect_identical(request$args[[6L]], "ProForma")
  expect_identical(mz$PEPTIDOFORM_CONVERSION_MODES,
                   c("ThrowException", "ReturnNull", "RemoveIncompatibleElements", "UsePrimarySequence"))
})

test_that("bad arguments are refused before the bridge", {
  build <- function(sequences = "PEPTIDE", ...) {
    given <- list(source = "mzLib", target = "Unimod", mode = "ReturnNull", threads = 1)
    given[names(list(...))] <- list(...)
    do.call(mz$peptidoform_convert_request, c(list(sequences), given))
  }
  expect_error(build(character(0)), class = "mzlib_usage_error", contains = "At least one sequence")
  expect_error(build(list("PEPTIDE")), class = "mzlib_usage_error", contains = "character vector")
  expect_error(build(factor("PEPTIDE")), class = "mzlib_usage_error", contains = "character vector")
  expect_error(build(c("PEPTIDE", "  ")), class = "mzlib_usage_error", contains = "blank")
  expect_error(build(c("PEPTIDE", NA)), class = "mzlib_usage_error", contains = "NA")
  expect_error(build("PEP\nTIDE"), class = "mzlib_usage_error", contains = "line break")
  expect_error(build(mode = "Strict"), class = "mzlib_usage_error", contains = "mode must be one of")
  expect_error(build(mode = NULL), class = "mzlib_usage_error", contains = "mode must be one of")
  expect_error(build(source = " "), class = "mzlib_usage_error", contains = "source must be")
  expect_error(build(target = NA_character_), class = "mzlib_usage_error", contains = "target must be")
  expect_error(build(threads = 0), class = "mzlib_usage_error", contains = "threads")
  expect_error(build(threads = -2), class = "mzlib_usage_error", contains = "threads")
  expect_error(build(threads = 1.5), class = "mzlib_usage_error", contains = "threads")
  expect_error(build(threads = TRUE), class = "mzlib_usage_error", contains = "threads")
})

# ---------------------------------------------------------------- the payload

test_that("the probe cases map to their Unimod accessions", {
  result <- convert_recorded("peptidoform_convert_unimod.json")
  expect_true(inherits(result, "mzlibr_sequence_conversions"))
  r <- result$records
  expect_identical(r$input, CONVERT_PROBE)
  expect_identical(r$output[[1L]], "[UNIMOD:1]SEQK")   # N-acetylserine
  expect_identical(r$output[[2L]], "PEPK[UNIMOD:36]R") # N6,N6-dimethyllysine, not Ethyl (UNIMOD:280)
  expect_identical(r$output[[3L]], "PEPM[UNIMOD:35]K") # Oxidation on M
  expect_identical(r$status[1:3], rep("converted", 3L))
})

test_that("an unmappable modification fails its row and names itself", {
  result <- convert_recorded("peptidoform_convert_unimod.json")
  bad <- result$records[result$records$status != "converted", ]
  expect_identical(nrow(bad), 1L)
  expect_identical(bad$input, "PEPK[Made Up:Not a modification on K]R")
  expect_true(is.na(bad$output))
  expect_identical(bad$status, "failed")
  expect_identical(bad$incompatible_items[[1L]], "Made Up:Not a modification on K @3(K)")
  # mzLib's Unimod serializer under ReturnNull records no reason code; projected as NA, not filled in.
  expect_true(is.na(bad$failure_reason))
  expect_identical(result$records$output, c("[UNIMOD:1]SEQK", "PEPK[UNIMOD:36]R", "PEPM[UNIMOD:35]K", NA))
  expect_identical(c(result$record_count, result$converted_count, result$warned_count, result$failed_count),
                   c(4, 3, 0, 1))
})

test_that("the envelope names the formats mzLib registered", {
  result <- convert_recorded("peptidoform_convert_unimod.json")
  expect_identical(c(result$source_format, result$target_format, result$mode), c("mzLib", "Unimod", "ReturnNull"))
  expect_identical(result$source_formats, c("MassShift", "Modomics", "ProForma", "mzLib"))
  expect_identical(result$target_formats, c("Chronologer", "Essential", "MassShift", "ProForma", "Unimod", "mzLib"))
  expect_identical(result$column_names, mz$PEPTIDOFORM_CONVERT_COLUMNS)
  expect_identical(names(result$records), result$column_names)
  expect_identical(length(result$caveats), 6L)
  expect_true(any(grepl("ProForma", result$caveats) & grepl("UniProt", result$caveats)))
})

test_that("the array columns are list columns of character vectors, empty ones kept", {
  r <- convert_recorded("peptidoform_convert_unimod.json")$records
  for (column in c("incompatible_items", "warnings", "errors")) {
    expect_true(is.list(r[[column]]), info = column)
    expect_identical(length(r[[column]]), 4L, info = column)
    expect_true(all(vapply(r[[column]], is.character, logical(1L))), info = column)
  }
  expect_identical(r$incompatible_items[[1L]], character(0))
  expect_true(is.character(r$failure_reason))
})

test_that("real psmtsv sequences convert to Unimod", {
  result <- convert_recorded("peptidoform_convert_psmtsv.json")
  expect_identical(result$records$input, psmtsv_full_sequences())
  expect_identical(result$records$output[[1L]], "YPIEH[UNIMOD:34]GIVTNWDDMEK")
  expect_identical(result$records$output[[5L]], "YPIEH[UNIMOD:34]GIVTNWDDM[UNIMOD:35]EK")
  expect_identical(c(result$converted_count, result$record_count), c(8, 8))
  expect_false(any(grepl("UniProt:", result$records$output, fixed = TRUE)))
})

test_that("ProForma leaves UniProt modifications unresolved", {
  # The mzLib gap the caveats describe (mzLib#1401), pinned so a fix upstream is noticed when the
  # fixture is re-recorded.
  result <- convert_recorded("peptidoform_convert_psmtsv_proforma.json")
  expect_identical(result$target_format, "ProForma")
  expect_identical(result$records$output[[1L]], "YPIEH[UniProt:Tele-methylhistidine on H]GIVTNWDDMEK")
  expect_identical(result$records$output[[3L]], "AYHEQLSVAEITNAC[UNIMOD:4]FEPANQMVK")
  expect_identical(result$records$status[[1L]], "converted")
})

test_that("an older bridge without the new fields still projects", {
  result <- mz$peptidoform_parse_conversions(list(columns = list(input = list("A"), output = list(NA))))
  expect_true(is.na(result$record_count))
  expect_identical(result$caveats, character(0))
  expect_identical(result$source_formats, character(0))
  expect_identical(result$records$output, NA_character_)
  empty <- mz$peptidoform_parse_conversions(list())
  expect_identical(names(empty$records), mz$PEPTIDOFORM_CONVERT_COLUMNS)
  expect_identical(nrow(empty$records), 0L)
})

test_that("columns of different lengths are a protocol error", {
  expect_error(mz$peptidoform_parse_conversions(list(
    column_names = list("input", "output"), columns = list(input = list("A", "B"), output = list("A"))
  )), class = "mzlib_protocol_error")
})

test_that("printing says how many converted, warned and failed", {
  output <- paste(capture.output(print(convert_recorded("peptidoform_convert_unimod.json"))), collapse = "\n")
  expect_true(grepl("4 sequences, mzLib to Unimod (ReturnNull)", output, fixed = TRUE))
  expect_true(grepl("3 converted, 0 with warnings, 1 failed", output, fixed = TRUE))
})

test_that("all three recordings replay end to end, and each answers only its own sequences", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(peptidoform_convert(CONVERT_PROBE)$failed_count, 1)
  psms <- readers_read_results("BottomUpExample.psmtsv")
  expect_identical(psms$records$full_sequence, psmtsv_full_sequences())
  expect_identical(peptidoform_convert(psms$records$full_sequence)$converted_count, 8)
  proforma <- peptidoform_convert(psms$records$full_sequence, target = "ProForma")
  expect_identical(proforma$records$output[[1L]], "YPIEH[UniProt:Tele-methylhistidine on H]GIVTNWDDMEK")
  expect_error(peptidoform_convert(rev(CONVERT_PROBE)), class = "mzlib_usage_error", contains = "replay bridge")
  expect_error(peptidoform_convert(CONVERT_PROBE, mode = "UsePrimarySequence"), class = "mzlib_usage_error",
               contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

test_that("LIVE: the probe cases convert to the recorded Unimod accessions", {
  skip_unless_bridge_has("peptidoform convert")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  result <- peptidoform_convert(CONVERT_PROBE)
  recorded <- convert_recorded("peptidoform_convert_unimod.json")
  expect_identical(result$records, recorded$records)
  expect_identical(result$source_formats, recorded$source_formats)
  expect_identical(result$target_formats, recorded$target_formats)
  expect_identical(result$caveats, recorded$caveats)
})

test_that("LIVE: the psmtsv sequences agree with both recordings, at any thread count", {
  skip_unless_bridge_has("peptidoform convert")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  sequences <- psmtsv_full_sequences()
  expect_identical(peptidoform_convert(sequences, threads = -1)$records,
                   convert_recorded("peptidoform_convert_psmtsv.json")$records)
  expect_identical(peptidoform_convert(sequences, target = "proforma")$records,
                   convert_recorded("peptidoform_convert_psmtsv_proforma.json")$records)
})

test_that("LIVE: the modes drop, keep or refuse what the target cannot write", {
  skip_unless_bridge_has("peptidoform convert")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  bad <- "PEPK[Made Up:Not a modification on K]R"
  for (mode in c("RemoveIncompatibleElements", "UsePrimarySequence")) {
    row <- peptidoform_convert(bad, mode = mode)$records
    expect_identical(row$status, "converted_with_warnings", info = mode)
    expect_false(is.na(row$output), info = mode)
    expect_false(grepl("Made Up", row$output, fixed = TRUE), info = mode)
    expect_true(length(row$incompatible_items[[1L]]) > 0L, info = mode)
  }
  expect_error(peptidoform_convert(c("PEPTIDE", bad), mode = "ThrowException"),
               class = "mzlib_usage_error", contains = "Made Up")
  expect_error(peptidoform_convert("PEPTIDE", source = "NotAFormat"), class = "mzlib_usage_error",
               contains = "mzLib")
})
