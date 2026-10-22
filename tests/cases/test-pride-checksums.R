# PRIDE's checksum list: pride_file_checksums().
#
# Ported from pyMzLib's test_pride.py and test_pride_live.py (the checksum sections). The payload
# was recorded live through the bridge by pyMzLib (mzLib #1427) and is shared verbatim, with the
# PXD000001 manifest and FTP listing the other PRIDE tests use. What this layer must get right: the
# accession is normalised before anything runs, the rows are read strictly, an empty list is an
# answer, and only mzLib's unknown-project wording becomes mzlib_project_not_found.

recorded_checksums <- function() {
  mz$pride_parse_checksums(recorded_data("pride_PXD000001_checksums.json"), "PXD000001")
}
checksum_manifest <- function() mz$pride_parse_manifest(recorded_data("pride_PXD000001_files.json"), "PXD000001")
CHECKSUM_RAW <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw"

# ---------------------------------------------------------------- the payload

test_that("every checksum row is keyed by its file name, in the wire's order", {
  checksums <- recorded_checksums()
  wire <- recorded_data("pride_PXD000001_checksums.json")
  expect_identical(nrow(checksums), 12L)
  expect_identical(wire$file_count, 12)
  expect_identical(names(checksums), c("file_name", "md5", "size_bytes", "project_accession"))
  expect_identical(rownames(checksums), checksums$file_name)
  # The wire sorts by name (ordinal), so the frame is the same on every run.
  expect_identical(checksums$file_name, vapply(wire$files, function(f) f$file_name, character(1L)))
  expect_identical(checksums$file_name, sort(checksums$file_name, method = "radix"))
  expect_false(is.factor(checksums$md5))
})

test_that("a checksum row is typed, and carries its project", {
  raw <- recorded_checksums()[CHECKSUM_RAW, ]
  expect_identical(raw$md5, "c567cecbb9866796e690a4958746f7b6")
  # A double, so a raw file past 2 GiB cannot overflow an integer.
  expect_true(is.double(raw$size_bytes))
  expect_identical(raw$size_bytes, 220475548)
  expect_identical(raw$project_accession, "PXD000001")
  readme <- recorded_checksums()["README.txt", ]
  expect_identical(readme$md5, "492753137013825025e6729c8798664f")
  expect_identical(readme$size_bytes, 1645)
})

test_that("the exact size is not the manifest's, and the manifest's checksum is not the MD5", {
  checksums <- recorded_checksums()
  files <- checksum_manifest()
  gz <- "PRIDE_Exp_Complete_Ac_22134.xml.gz"
  expect_identical(files$file_size_bytes[files$file_name == gz], 10677205)
  expect_identical(checksums[gz, "size_bytes"], 10668000)
  expect_identical(unique(files$checksum), "")
  # The two PRIDE-generated files are the manifest's files without a row.
  expect_identical(
    files$file_name[!files$file_name %in% checksums$file_name],
    c("PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz", "PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz")
  )
  # It names a file the REST manifest omits.
  mzml <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"
  expect_true(mzml %in% checksums$file_name)
  expect_false(mzml %in% files$file_name)
  # Exact sizes, so an exact sum of what the deposited files weigh.
  expect_identical(round(sum(checksums$size_bytes) / 1e9, 2), 1.43)
})

test_that("the generated/ files of the FTP tree have no checksum row", {
  listing <- mz$pride_parse_ftp_files(recorded_data("pride_PXD000001_ftp_files.json"), "PXD000001")
  unchecked <- listing$relative_path[!listing$relative_path %in% recorded_checksums()$file_name]
  expect_identical(unchecked, c(
    "generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz",
    "generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz"
  ))
})

test_that("an existing project with no checksums is zero rows, not an error", {
  empty <- mz$pride_parse_checksums(list(accession = "PXD000001", file_count = 0, files = list()), "PXD000001")
  expect_identical(nrow(empty), 0L)
  expect_identical(names(empty), c("file_name", "md5", "size_bytes", "project_accession"))
})

test_that("a checksum row without its MD5 or size fails rather than defaulting", {
  # Defaulting it to "" or 0 would make every later check of that file blame the download.
  rows <- list(list(file_name = "run1.raw", size_bytes = 10), list(file_name = "run1.raw", md5 = "a"))
  for (row in rows) {
    payload <- list(accession = "PXD000001", file_count = 1, files = list(row))
    expect_error(mz$pride_parse_checksums(payload, "PXD000001"), class = "mzlib_protocol_error")
  }
})

# ---------------------------------------------------------------- the call

test_that("pride_file_checksums sends the verb and the normalised accession, and nothing else", {
  seen <- recorded_bridge_call("pride_PXD000001_checksums.json", function() {
    pride_file_checksums("  pxd000001 ", timeout = 42)
  })
  # The verb check asks the bridge first, then the one call. PRIDE's route is case-sensitive.
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("pride", "checksums", "--accession", "PXD000001"))
  expect_identical(seen$timeout, 42)
  expect_identical(nrow(seen$result), 12L)
})

test_that("pride_file_checksums adds no timeout of its own", {
  # mzLib bounds each PRIDE request itself (#1412), as for the other PRIDE metadata calls.
  expect_true(is.null(formals(pride_file_checksums)$timeout))
  seen <- recorded_bridge_call("pride_PXD000001_checksums.json", function() pride_file_checksums("PXD000001"))
  expect_true(is.null(seen$timeout))
})

test_that("a bad accession is refused before any work, the verb check included", {
  for (accession in list("", "   ", "not-an-accession", 12345, NA_character_, c("PXD000001", "PXD000002"))) {
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() {
        pride_file_checksums(accession)
      }),
      class = "mzlib_usage_error"
    )
  }
})

test_that("an unknown project is mzlib_project_not_found, and a broken list is not", {
  unknown <- bridge_error_answer("MzLibException", "PRIDE Archive has no project with accession 'PXD999999999'.")
  expect_error(
    recorded_bridge_call(unknown, function() pride_file_checksums("PXD999999999")),
    class = "mzlib_project_not_found", contains = "no public project 'PXD999999999'"
  )
  broken <- c(
    "PRIDE Archive's checksum list for 'PXD000001' starts with 'File-Name,MD5', not the expected 'File-Name,File-MD5Checksum,File-Size'.",
    "PRIDE Archive's checksum list for 'PXD000001' lists 'run1.raw' twice."
  )
  for (message in broken) {
    # The project exists and its list is broken: "check for a typo" would point at the wrong mistake.
    condition <- expect_error(
      recorded_bridge_call(bridge_error_answer("MzLibException", message), function() pride_file_checksums("PXD000001")),
      class = "mzlib_bridge_error"
    )
    expect_false(inherits(condition, "mzlib_project_not_found"))
    expect_identical(condition$error_type, "MzLibException")
  }
  expect_error(
    recorded_bridge_call(
      bridge_error_answer("ServiceUnavailable", "PRIDE Archive request failed with status 503"),
      function() pride_file_checksums("PXD000001")
    ),
    class = "mzlib_service_unavailable"
  )
})

test_that("a bridge older than the checksum verb is named, not spawned for it", {
  # readers_formats.json has no `verbs`, like a bridge from before pyMzLib 0.2.0.
  expect_error(
    recorded_bridge_call(function(args, stdin) stop("reached the verb"), function() pride_file_checksums("PXD000001"),
                         version = "readers_formats.json"),
    class = "mzlib_usage_error", contains = "pyMzLib 0.6.0"
  )
})

test_that("checksums replay end to end, as the help page runs them", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(pride_file_checksums("PXD000001")["README.txt", "size_bytes"], 1645)
})

# ---------------------------------------------------------------- against the real archive

test_that("LIVE: the checksum list covers the deposited files, and none of generated/", {
  skip_unless_bridge_has("pride checksums")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  both <- tryCatch(
    list(checksums = pride_file_checksums("PXD000001"), ftp = pride_list_ftp_files("PXD000001")),
    mzlib_service_unavailable = function(e) skip(paste("PRIDE unavailable:", conditionMessage(e)))
  )
  # Relations to the FTP tree rather than counts, so this survives PRIDE re-curating the project.
  expect_true(nrow(both$checksums) > 0L)
  root <- both$ftp$relative_path[!grepl("/", both$ftp$relative_path, fixed = TRUE)]
  expect_true(all(both$checksums$file_name %in% root))
  expect_true(any(startsWith(both$ftp$relative_path, "generated/")))
  expect_true(all(nchar(both$checksums$md5) == 32L & both$checksums$size_bytes > 0))
})

test_that("LIVE: a checksum row matches the bytes a download delivers", {
  skip_unless_bridge_has("pride checksums")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  destination <- tempfile("pride-checksum-")
  on.exit(unlink(destination, recursive = TRUE), add = TRUE)
  fetched <- tryCatch(
    {
      checksums <- pride_file_checksums("PXD000001")
      files <- pride_list_files("PXD000001")
      candidates <- files[files$downloadable & files$file_name %in% checksums$file_name, ]
      chosen <- candidates[which.min(checksums[candidates$file_name, "size_bytes"]), ]
      list(checksums = checksums, chosen = chosen, written = pride_download_files(chosen, destination))
    },
    mzlib_service_unavailable = function(e) skip(paste("PRIDE unavailable:", conditionMessage(e)))
  )
  expected <- fetched$checksums[fetched$chosen$file_name, ]
  expect_identical(file.size(fetched$written[[1L]]), expected$size_bytes)
  expect_identical(unname(tools::md5sum(fetched$written[[1L]])), tolower(expected$md5))
})

test_that("LIVE: an unknown accession raises rather than reporting no checksums", {
  skip_unless_bridge_has("pride checksums")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  condition <- tryCatch(
    {
      pride_file_checksums("PXD999999999")
      NULL
    },
    mzlib_service_unavailable = function(e) skip(paste("PRIDE unavailable:", conditionMessage(e))),
    mzlib_error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_project_not_found"))
})
