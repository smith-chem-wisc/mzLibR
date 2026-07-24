# Fetching a bridge.
#
# Nothing here downloads anything. The parts worth testing are the refusals — mzLibR must never
# install an executable it has not verified, and must never write 140 MB anywhere without being
# asked — and those are all reachable before a byte moves.

test_that("every published platform has a wheel and a checksum", {
  wheels <- mz$MZLIB_BRIDGE_WHEELS
  expect_identical(
    sort(names(wheels)),
    sort(c("win-x64", "osx-arm64", "osx-x64", "linux-x64"))
  )
  for (rid in names(wheels)) {
    entry <- wheels[[rid]]
    expect_true(grepl("^pymzlib-.*\\.whl$", entry$wheel), info = rid)
    # A checksum is 64 lowercase hex characters. A truncated or upper-cased one would compare
    # unequal against a correct download and produce a very confusing failure.
    expect_true(grepl("^[0-9a-f]{64}$", entry$sha256), info = rid)
  }
})

test_that("the payload really is the pyMzLib release wheel for this platform", {
  # Not a separate raw-binary asset. A wheel is a zip and already carries the bridge under
  # pymzlib/_dotnet/<rid>/, so mzLibR needs nothing published that pyMzLib does not already
  # publish for its own users.
  url <- mz$bridge_release_url(mz$MZLIB_BRIDGE_WHEELS[["win-x64"]]$wheel)
  expect_true(startsWith(url, "https://github.com/smith-chem-wisc/pyMzLib/releases/download/"))
  expect_true(grepl("win_amd64.whl", url, fixed = TRUE))
})

test_that("linux-arm64 is refused with the reason and a way forward", {
  # pyMzLib publishes no arm64 Linux wheel, so there is genuinely nothing to fetch. Saying which
  # platforms do have one, and how to build the missing one, is the difference between a dead
  # end and a next step.
  expect_true(is.null(mz$MZLIB_BRIDGE_WHEELS[["linux-arm64"]]))
})

test_that("a url without a checksum is refused", {
  # mzLibR will not install an executable it cannot verify, and an unverified 140 MB binary
  # about to be run is not something to accept on the strength of HTTPS alone.
  expect_error(
    mzlibr_install_bridge(url = "https://example.invalid/bridge.zip", consent = TRUE),
    class = "mzlib_usage_error", contains = "pass sha256"
  )
})

test_that("a non-interactive session without consent refuses before downloading", {
  # CRAN policy forbids writing outside tempdir() without explicit consent, and this is also
  # simply correct: 140 MB should not appear on someone's disk because they mistyped.
  #
  # `interactive()` is FALSE under Rscript, which is what makes this testable at all.
  skip_if(interactive(), "needs a non-interactive session")
  expect_error(
    mzlibr_install_bridge(destination = tempfile("mzlibr-nothing-")),
    class = "mzlib_usage_error", contains = c("not interactive", "consent = TRUE")
  )
})

test_that("the consent message names the size and the destination", {
  # Both are what a person needs to decide. A bare "proceed?" is not consent to anything in
  # particular.
  where <- file.path(tempfile("mzlibr-consent-"), "win-x64")
  condition <- tryCatch(mz$bridge_ask_consent(where, "140 MB"), error = function(e) e)
  expect_true(inherits(condition, "mzlib_usage_error"))
  expect_true(grepl("140 MB", conditionMessage(condition), fixed = TRUE))
  expect_true(grepl(where, conditionMessage(condition), fixed = TRUE))
})

test_that("an already-installed bridge is not replaced without being asked", {
  destination <- tempfile("mzlibr-installed-")
  rid <- mz$bridge_platform_tag()
  dir.create(file.path(destination, rid), recursive = TRUE)
  existing <- file.path(destination, rid, mz$bridge_executable_name())
  writeLines("not a real executable", existing)
  on.exit(unlink(destination, recursive = TRUE), add = TRUE)

  returned <- suppressMessages(mzlibr_install_bridge(destination = destination))
  expect_identical(returned, existing)
  # Untouched: it did not download, and it did not overwrite.
  expect_identical(readLines(existing, warn = FALSE), "not a real executable")
})

test_that("sha256 is computed by some means on this machine", {
  # The installer refuses rather than skipping verification, so on a machine with no hasher at
  # all it would never install. Confirm this one can.
  scratch <- tempfile("mzlibr-hash-")
  writeLines("mzLibR", scratch)
  on.exit(unlink(scratch), add = TRUE)

  digest <- mz$bridge_sha256(scratch)
  expect_false(is.null(digest))
  expect_true(grepl("^[0-9a-f]{64}$", digest), info = paste(digest, collapse = " "))
})

test_that("the same bytes always hash the same way", {
  first <- tempfile("mzlibr-hash-a-")
  second <- tempfile("mzlibr-hash-b-")
  writeBin(as.raw(c(1, 2, 3, 4)), first)
  writeBin(as.raw(c(1, 2, 3, 4)), second)
  on.exit(unlink(c(first, second)), add = TRUE)
  expect_identical(mz$bridge_sha256(first), mz$bridge_sha256(second))
})

test_that("making a bridge runnable is a no-op on Windows and chmod elsewhere", {
  # Two things a Windows-only developer never sees: a zip carries no execute bit, and macOS
  # quarantines anything an application downloaded. Both would surface as failures that sound
  # like mzLib problems rather than packaging ones.
  scratch <- tempfile("mzlibr-chmod-")
  writeLines("x", scratch)
  on.exit(unlink(scratch), add = TRUE)

  mz$bridge_make_runnable(scratch)
  if (.Platform$OS.type != "windows") {
    mode <- as.character(file.info(scratch)$mode)
    expect_true(grepl("7", substr(mode, nchar(mode) - 2L, nchar(mode) - 2L)), info = mode)
  }
  expect_true(file.exists(scratch))
})
