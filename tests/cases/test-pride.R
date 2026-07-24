# PRIDE Archive access.
#
# Every number asserted here comes from the recorded PXD000001 manifest, which is the same
# fixture pyMzLib and mzLibRust use. That is deliberate: the documentation quotes these figures,
# and a number in a doc string that no test reproduces is a number that will rot. One already
# did — "37 peptides out of about 200" was stale by 5x and had been copied into a second
# binding before anyone measured it.

recorded_manifest <- function() {
  mz$json_parse(paste(readLines(fixture_path("pride_PXD000001_files.json"), warn = FALSE),
    collapse = "\n"
  ))
}

recorded_files <- function() {
  mz$pride_parse_manifest(recorded_manifest(), "PXD000001")
}

download_args <- function(accession = "PXD000001", dest = "out", category = NULL,
                          extensions = NULL, overwrite = TRUE) {
  mz$pride_build_download_args(accession, dest, category, extensions, overwrite)
}

# ---------------------------------------------------------------- parsing the manifest

test_that("every file in the manifest becomes a row", {
  expect_identical(nrow(recorded_files()), length(recorded_manifest()$files))
})

test_that("the manifest is a plain data.frame with no factors", {
  # `stringsAsFactors` defaulted to TRUE before R 4.0. A package that relies on the default
  # gives factors on old R and characters on new, and the difference only shows up when a user
  # compares a column to a string. Every data.frame in mzLibR passes it explicitly.
  files <- recorded_files()
  expect_true(is.data.frame(files))
  for (column in c("file_name", "category", "checksum", "extension", "project_accession")) {
    expect_true(is.character(files[[column]]), info = column)
    expect_false(is.factor(files[[column]]), info = column)
  }
})

test_that("file sizes are doubles, never integers", {
  # A PRIDE project runs past 2^31 — this one totals 1.44 GB on disk — so an integer column
  # would overflow to NA on real data.
  files <- recorded_files()
  expect_true(is.double(files$file_size_bytes))
  expect_true(all(files$file_size_bytes > 0))
  expect_equal(files$size_mb, files$file_size_bytes / 1e6)
})

test_that("the recorded manifest holds the eight files PRIDE's API publishes", {
  # Eight, while the FTP tree holds thirteen. The number is quoted in ?pride_list_files, so it
  # is pinned here.
  expect_identical(nrow(recorded_files()), 8L)
})

test_that("total_size_bytes agrees with the manifest's own total", {
  files <- recorded_files()
  expect_identical(mz$pride_total_size_bytes(files), sum(files$file_size_bytes))
  expect_identical(mz$pride_total_size_bytes(files), recorded_manifest()$total_size_bytes)
  # 0.51 GB reported against 1.44 GB on disk. Both figures are in the docs; this pins the one
  # mzLibR can actually compute.
  expect_true(mz$pride_total_size_bytes(files) < 0.52e9)
})

test_that("a compressed file's extension is .gz, not what it is compressed from", {
  # The single most common way to get an empty result, and the trap a biologist walked into in
  # the mzLibRust bake-off. Asserted on real data so the doc string cannot drift from it.
  files <- recorded_files()
  mgf <- files[files$file_name == "PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz", ]
  expect_identical(nrow(mgf), 1L)
  expect_identical(mgf$extension, ".gz")
  expect_identical(sum(files$extension == ".mgf"), 0L)
  # And the over-match in the other direction: ".gz" is not selective either.
  expect_true(sum(files$extension == ".gz") > 1L)
})

test_that("extensions are lowercased and keep their dot", {
  files <- recorded_files()
  present <- files$extension[nzchar(files$extension)]
  expect_true(all(startsWith(present, ".")))
  expect_identical(present, tolower(present))
})

test_that("a file with no extension gets an empty string, not a dot", {
  expect_identical(mz$pride_extension("README"), "")
  expect_identical(mz$pride_extension("x.MGF.GZ"), ".gz")
  expect_identical(mz$pride_extension("dir.with.dots/file"), "")
})

test_that("downloadable follows the HTTPS URL", {
  files <- recorded_files()
  expect_identical(files$downloadable, !is.na(files$https_url))
  # And it is never NA, because it is used as a subsetting predicate — `files[files$downloadable, ]`
  # on an NA yields a row of NAs rather than dropping the row.
  expect_false(any(is.na(files$downloadable)))
})

test_that("every row knows which project it came from", {
  # `pride_download_files()` needs this to tell which project to fetch from, so the caller does
  # not have to carry the accession alongside the frame.
  expect_true(all(recorded_files()$project_accession == "PXD000001"))
})

test_that("timestamps parse to POSIXct and unreadable ones become NA", {
  # A caller asked for a file list, not a date. One unreadable timestamp must not fail a whole
  # manifest.
  files <- recorded_files()
  expect_true(inherits(files$submission_date, "POSIXct"))
  expect_false(all(is.na(files$submission_date)))

  expect_true(is.na(mz$pride_parse_timestamp("not-a-date")))
  expect_true(is.na(mz$pride_parse_timestamp("")))
  expect_false(is.na(mz$pride_parse_timestamp("2012-02-07T00:00:00+00:00")))
  expect_false(is.na(mz$pride_parse_timestamp("2012-02-07T00:00:00Z")))
  expect_false(is.na(mz$pride_parse_timestamp("2012-02-07T00:00:00")))
  expect_false(is.na(mz$pride_parse_timestamp("2012-02-07")))
})

test_that("an offset timestamp lands on the right instant", {
  # "+00:00" and "Z" are rewritten before strptime sees them, because %z accepts neither on
  # every platform and the ones it silently fails on differ between Windows and Linux.
  utc <- mz$pride_parse_timestamp("2012-03-13T00:00:00+00:00")
  plus_one <- mz$pride_parse_timestamp("2012-03-13T01:00:00+01:00")
  expect_identical(as.numeric(utc), as.numeric(plus_one))
})

test_that("a valid but unknown accession raises rather than returning zero rows", {
  # PRIDE answers an unknown accession with HTTP 200 and an empty list, not a 404. A zero-row
  # frame would let a typo produce "0 files, done" and a green exit.
  expect_error(
    mz$pride_parse_manifest(list(accession = "PXD999999", files = list()), "PXD999999"),
    class = "mzlib_project_not_found", contains = "typo"
  )
})

# ---------------------------------------------------------------- locations

test_that("locations survive subsetting the data.frame", {
  # This is why locations is a list column and not an attribute: `[.data.frame` drops
  # attributes, so a biologist writing `files[files$size_mb < 5, ]` would silently lose them.
  files <- recorded_files()
  small <- files[files$size_mb < 5, ]
  expect_true(nrow(small) > 0L)
  expect_true("locations" %in% names(small))
  expect_true(is.data.frame(small$locations[[1]]))
})

test_that("the published location order is not stable, which is why https_url exists", {
  # In this one project the mztab lists FTP first and the MGF lists Aspera first. Code that
  # took locations[[1]] would get an unfetchable prd_ascp@... address for some files and a
  # working URL for others. mzLib's TryGetHttpsDownloadUrl searches; this test pins the fact
  # that searching is necessary.
  locations <- mz$pride_locations(recorded_files())
  first_of_each <- locations[!duplicated(locations$file_name), ]
  expect_true(length(unique(first_of_each$name)) > 1L,
    info = paste(unique(first_of_each$name), collapse = " / ")
  )

  aspera <- locations[locations$name == "Aspera Protocol", ]
  expect_true(nrow(aspera) > 0L)
  expect_true(all(grepl("^prd_ascp@", aspera$value)))
  # Whatever the order, the resolved URL is always HTTPS.
  files <- recorded_files()
  expect_true(all(startsWith(files$https_url[!is.na(files$https_url)], "https://")))
})

test_that("pride_locations returns one long frame", {
  locations <- mz$pride_locations(recorded_files())
  expect_true(is.data.frame(locations))
  expect_identical(names(locations), c("file_name", "accession", "name", "value"))
  expect_true(nrow(locations) >= nrow(recorded_files()))
})

# ---------------------------------------------------------------- accessions and paging

test_that("a blank or malformed accession is rejected before any work", {
  for (accession in list("", "   ", "banana", "PXD", "12345", "PXD00", "PXD000001x", "-PXD1")) {
    expect_error(mz$pride_normalise_accession(accession), class = "mzlib_usage_error")
  }
  expect_error(mz$pride_normalise_accession(NULL), class = "mzlib_usage_error")
  expect_error(mz$pride_normalise_accession(c("PXD000001", "PXD000002")), class = "mzlib_usage_error")
})

test_that("accession case and whitespace are normalised", {
  # PRIDE's API is case-sensitive on the accession while category matching is not. Two rules
  # pointing opposite ways is a trap, and this is the one fixable without surprise.
  for (accession in c("pxd000001", "  PXD000001  ", "Pxd000001")) {
    expect_identical(mz$pride_normalise_accession(accession), "PXD000001")
  }
})

test_that("the list args carry the canonical accession and page size", {
  expect_identical(
    mz$pride_build_list_args("  pxd000001 ", 100),
    c("pride", "files", "--accession", "PXD000001", "--page-size", "100")
  )
})

test_that("an impossible page size is refused", {
  expect_error(mz$pride_build_list_args("PXD000001", 0),
    class = "mzlib_usage_error", contains = "must be positive"
  )
  expect_error(mz$pride_build_list_args("PXD000001", -1), class = "mzlib_usage_error")
  expect_error(mz$pride_build_list_args("PXD000001", 2147483648),
    class = "mzlib_usage_error", contains = "larger than the API allows"
  )
  expect_error(mz$pride_build_list_args("PXD000001", 1.5), class = "mzlib_usage_error")
  expect_error(mz$pride_build_list_args("PXD000001", "100"), class = "mzlib_usage_error")
})

test_that("a large page size is not written in scientific notation", {
  # `format(1e6)` is "1e+06", which the bridge's integer parser rejects. Easy to miss because R
  # only switches notation above a threshold.
  args <- mz$pride_build_list_args("PXD000001", 1000000)
  expect_identical(args[6], "1000000")
})

# ---------------------------------------------------------------- download arguments

test_that("download passes the accession and destination", {
  args <- download_args()
  expect_identical(args[1:2], c("pride", "download"))
  expect_true("PXD000001" %in% args)
  expect_true("out" %in% args)
})

test_that("filters are omitted when not asked for", {
  args <- download_args()
  expect_false("--category" %in% args)
  expect_false("--ext" %in% args)
})

test_that("category and extensions reach the wire in the bridge's spelling", {
  args <- download_args(category = "RAW")
  expect_identical(args[which(args == "--category") + 1L], "RAW")

  args <- download_args(extensions = c(".raw", ".mzML"))
  expect_identical(args[which(args == "--ext") + 1L], ".raw,.mzML")
})

test_that("the overwrite flag is not inverted", {
  # The bridge takes --no-overwrite while the R argument is `overwrite`. Getting the polarity
  # wrong re-downloads a project someone was resuming, or skips files they wanted replaced.
  expect_false("--no-overwrite" %in% download_args())
  expect_true("--no-overwrite" %in% download_args(overwrite = FALSE))
})

test_that("a blank destination is refused instead of writing to the working directory", {
  for (dest in list("", "   ", NA_character_, NULL)) {
    expect_error(download_args(dest = dest),
      class = "mzlib_usage_error", contains = "destination directory is required"
    )
  }
})

test_that("a blank category is refused rather than selecting everything", {
  expect_error(download_args(category = "  "),
    class = "mzlib_usage_error", contains = "Omit it to download every category"
  )
})

test_that("a flag-like filter value is refused", {
  # The bridge's parser reads `--a --b` as two flags, so a value beginning with '-' silently
  # discards the option it belonged to and can smuggle in a flag nobody intended.
  for (value in c("--no-overwrite", "-x")) {
    expect_error(download_args(category = value),
      class = "mzlib_usage_error", contains = "may not begin with"
    )
  }
})

test_that("an extension list that names nothing is refused", {
  expect_error(mz$pride_normalise_extensions(c("  ", "")),
    class = "mzlib_usage_error", contains = "names no extensions"
  )
})

test_that("an extension containing a comma is refused", {
  # Commas are the wire's separator, so one inside a value would silently become two filters.
  expect_error(mz$pride_normalise_extensions(".raw,.mzML"),
    class = "mzlib_usage_error", contains = "may not contain a comma"
  )
})

test_that("NULL extensions means no filter", {
  expect_identical(mz$pride_normalise_extensions(NULL), character(0))
  expect_identical(mz$pride_normalise_extensions(character(0)), character(0))
})

test_that("a filter that matched nothing raises, and names the .gz trap", {
  # The error message is the only documentation a stuck user reliably reads, so the explanation
  # goes in it rather than only in ?pride_download.
  expect_error(
    mz$pride_check_filter_matched(character(0), "PXD000001", "RAW", character(0)),
    class = "mzlib_usage_error", contains = c("No file in PXD000001 matched", ".gz")
  )
  expect_error(
    mz$pride_check_filter_matched(character(0), "PXD000001", NULL, ".mgf"),
    class = "mzlib_usage_error", contains = ".gz"
  )
})

test_that("a download with no filter may legitimately write nothing", {
  # With overwrite = FALSE and everything already present, zero written files is the right
  # answer, not an error.
  expect_true(is.null(mz$pride_check_filter_matched(character(0), "PXD000001", NULL, character(0))))
})

# ---------------------------------------------------------------- download_files

test_that("the selection travels on stdin, not argv", {
  # argv has a ceiling of roughly 32 KB and a few thousand file names go straight past it.
  files <- recorded_files()
  selection <- files[files$downloadable, ][1:2, ]
  prepared <- mz$pride_build_download_files_args(selection, "out", TRUE)

  expect_true("--names-from-stdin" %in% prepared$args)
  expect_identical(prepared$stdin, selection$file_name)
  expect_identical(length(prepared$stdin), 2L)
})

test_that("an empty selection is refused rather than reported as success", {
  files <- recorded_files()
  expect_error(
    mz$pride_build_download_files_args(files[0, ], "out", TRUE),
    class = "mzlib_usage_error", contains = "No files selected"
  )
})

test_that("files with no HTTPS location are refused up front", {
  # Better than failing halfway through a multi-gigabyte transfer.
  files <- recorded_files()
  files$downloadable[1] <- FALSE
  files$https_url[1] <- NA_character_
  expect_error(
    mz$pride_build_download_files_args(files, "out", TRUE),
    class = "mzlib_usage_error", contains = c("no HTTPS location", "files$downloadable")
  )
})

test_that("a selection spanning two projects is refused", {
  files <- recorded_files()[1:2, ]
  files$project_accession[2] <- "PXD000002"
  expect_error(
    mz$pride_build_download_files_args(files, "out", TRUE),
    class = "mzlib_usage_error", contains = "must come from one project"
  )
})

test_that("rows carrying no project accession are refused", {
  files <- recorded_files()[1:2, ]
  files$project_accession <- ""
  expect_error(
    mz$pride_build_download_files_args(files, "out", TRUE),
    class = "mzlib_usage_error", contains = "pride_list_files()"
  )
})

test_that("a file name containing a line break is refused, not mis-parsed", {
  # Newline-delimited framing cannot represent it, and mis-parsing would silently select the
  # wrong files. PRIDE has never published such a name, but "never seen it" is not a contract.
  files <- recorded_files()[1:2, ]
  files$file_name[1] <- "two\nlines.raw"
  expect_error(
    mz$pride_build_download_files_args(files, "out", TRUE),
    class = "mzlib_usage_error", contains = "line break"
  )
})

test_that("download_files honours no-overwrite and refuses a blank destination", {
  files <- recorded_files()[1, ]
  expect_true("--no-overwrite" %in% mz$pride_build_download_files_args(files, "out", FALSE)$args)
  expect_error(mz$pride_build_download_files_args(files, "", TRUE), class = "mzlib_usage_error")
})

test_that("something that is not a manifest is refused with a useful message", {
  expect_error(
    mz$pride_build_download_files_args(list(file_name = "a"), "out", TRUE),
    class = "mzlib_usage_error", contains = "data.frame"
  )
  expect_error(
    mz$pride_build_download_files_args(data.frame(file_name = "a"), "out", TRUE),
    class = "mzlib_usage_error", contains = "missing the column"
  )
})

test_that("written paths are read back out of the payload", {
  expect_identical(
    mz$pride_parse_paths(list(paths = list("out/a.raw", "out/b.raw"))),
    c("out/a.raw", "out/b.raw")
  )
  expect_identical(mz$pride_parse_paths(list()), character(0))
  expect_identical(mz$pride_parse_paths(list(paths = list())), character(0))
})

# ---------------------------------------------------------------- against the real archive

test_that("LIVE: PXD000001 lists the files PRIDE's API publishes", {
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  files <- tryCatch(
    pride_list_files("PXD000001"),
    mzlib_service_unavailable = function(e) skip(paste("PRIDE unavailable:", conditionMessage(e)))
  )

  expect_true(nrow(files) > 0L)
  expect_true(all(c("file_name", "file_size_bytes", "https_url", "downloadable") %in% names(files)))
  # The trap, live: whatever else changed, a .gz is still a .gz.
  expect_identical(sum(files$extension == ".mgf"), 0L)
})

test_that("LIVE: an unknown accession raises rather than returning nothing", {
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  condition <- tryCatch(
    {
      pride_list_files("PXD999999")
      NULL
    },
    mzlib_service_unavailable = function(e) skip(paste("PRIDE unavailable:", conditionMessage(e))),
    mzlib_error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_project_not_found"))
})
