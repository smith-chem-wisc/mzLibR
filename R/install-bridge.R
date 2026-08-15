# Fetching a bridge executable, on the user's explicit instruction and never otherwise.
#
# Python solves the ~140 MB payload with a wheel that carries it; Rust solves it with a build
# script that downloads one. R has neither escape hatch. CRAN's package size limit is about
# 5 MB, CRAN policy forbids writing outside tempdir() without explicit user consent, and it
# forbids downloading during install or .onLoad. So the design is forced and it is different
# from both siblings: an exported function the user calls, which asks first.
#
# The payload source is the raw bridge tarball pyMzLib publishes on every release,
# `mzlib-bridge-<rid>.tar.gz`, alongside a `SHA256SUMS` manifest.
#
# It used to be pyMzLib's release WHEEL. That worked - a wheel is a zip and already carries the
# bridge under `pymzlib/_dotnet/<rid>/` - and it was defensible while nothing neutral was published:
# it meant mzLibR needed no asset pyMzLib was not already producing for its own users. But it left
# an R package reaching through a Python packaging format for no reason except that no other format
# existed, and it was one of the two things that made distribution shared surface rather than
# per-binding packaging (BRIDGE-DISTRIBUTION, 2026-08-08). pyMzLib #31 published the neutral
# artefact; this is mzLibR taking it.
#
# Two concrete gains beyond tidiness. The tarball restores the execute bit, where `utils::unzip()`
# drops it - the reason `bridge_make_runnable()` below exists at all. And it is the payload alone,
# not a whole wheel: about 64-174 MB per platform rather than the same tree wrapped in Python
# packaging metadata this package has no use for.

# The tarball that carries each platform's bridge, and its published SHA-256.
#
# Recorded here rather than fetched, so the checksum being verified does not come down the same
# connection as the thing it is verifying.
#
# The block below is REWRITTEN MECHANICALLY by .github/workflows/bridge-watch.yml, which reads the
# SHA256SUMS asset pyMzLib publishes on every release and regenerates the pins from it. Editing it
# by hand works and is fine for a one-off; the next bump will overwrite the edit. The comment above
# used to say "bumping this means re-recording all four" — that was the state of things before
# pyMzLib #31 published a checksum manifest, and re-recording four digests by hand is exactly the
# chore that manifest exists to end. Keep the markers: the workflow locates the region by them, and
# fails loudly rather than guessing if either is missing.
#
# BEGIN generated bridge pins
MZLIB_BRIDGE_VERSION <- "0.1.0.dev5"

MZLIB_BRIDGE_ASSETS <- list(
  "win-x64" = list(
    asset = "mzlib-bridge-win-x64.tar.gz",
    sha256 = "6ff5ae776c3daa2b228ad78555e026335725396e3b72cfb965534a9476b6ad0e"
  ),
  "osx-arm64" = list(
    asset = "mzlib-bridge-osx-arm64.tar.gz",
    sha256 = "397e1520df30f6cb06fcc27c114ea72fddf7044929826186c4ae9fec0519f4eb"
  ),
  "osx-x64" = list(
    asset = "mzlib-bridge-osx-x64.tar.gz",
    sha256 = "90fc1d23574cdfdccebaadb93a649a0b74a9eda05f298ac03924850d7513a58e"
  ),
  "linux-x64" = list(
    asset = "mzlib-bridge-linux-x64.tar.gz",
    sha256 = "fb6a15eaa1e9ea8cc719a245b2169ec88a9a33e2ad6b9302685ac56ee933741b"
  )
)
# END generated bridge pins

bridge_release_url <- function(asset, version = MZLIB_BRIDGE_VERSION) {
  paste0(
    "https://github.com/smith-chem-wisc/pyMzLib/releases/download/v", version, "/", asset
  )
}

# The SHA-256 of a file, or NULL if this R cannot compute one.
#
# `tools::sha256sum()` is recent; on an older R the system tools are tried, and if none is
# present the caller must decide. What is never done is to skip verification quietly: a 140 MB
# executable that is about to be run is not something to accept on the strength of HTTPS alone.
bridge_sha256 <- function(path) {
  if (exists("sha256sum", where = asNamespace("tools"), inherits = FALSE)) {
    value <- getExportedValue("tools", "sha256sum")(path)
    if (!is.na(value)) {
      return(unname(value))
    }
  }

  candidates <- list(
    list(command = "sha256sum", args = c(shQuote(path))),
    list(command = "shasum", args = c("-a", "256", shQuote(path))),
    list(command = "certutil", args = c("-hashfile", shQuote(path), "SHA256"))
  )
  for (candidate in candidates) {
    if (!nzchar(Sys.which(candidate$command))) {
      next
    }
    output <- tryCatch(
      suppressWarnings(system2(candidate$command, candidate$args, stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0)
    )
    hit <- grep("^[0-9a-fA-F]{64}", trimws(output), value = TRUE)
    if (length(hit) > 0L) {
      return(tolower(substr(trimws(hit[1L]), 1L, 64L)))
    }
    # certutil prints the digest on its own line, sometimes spaced.
    joined <- tolower(gsub("[^0-9a-f]", "", paste(tolower(output), collapse = "")))
    if (nchar(joined) >= 64L) {
      return(substr(joined, 1L, 64L))
    }
  }
  NULL
}

# Ask, when there is somebody to ask.
bridge_ask_consent <- function(destination, size_note) {
  if (!interactive()) {
    stop(mzlib_usage_error(paste0(
      "mzlibr_install_bridge() downloads about ", size_note, " and writes it to\n  ",
      destination, "\n",
      "This session is not interactive, so it cannot ask. Pass consent = TRUE to proceed, ",
      "or set the MZLIB_BRIDGE environment variable to a bridge you already have."
    )))
  }

  cat("mzlibr_install_bridge() will download about ", size_note, " and write it to\n  ",
    destination, "\n", sep = ""
  )
  answer <- readline("Proceed? [y/N] ")
  if (!tolower(trimws(answer)) %in% c("y", "yes")) {
    stop(mzlib_usage_error("Cancelled; nothing was downloaded or written."))
  }
  invisible(TRUE)
}

#' Download a bridge executable into a local cache
#'
#' mzLibR needs a bridge executable and does not ship one: it is about 140 MB, which no CRAN
#' package may carry. This fetches the one built for your platform, verifies it against a
#' recorded SHA-256, and unpacks it where [mzlibr_bridge_path()] will find it.
#'
#' **It is never called for you.** CRAN policy forbids a package downloading anything at install
#' time or writing outside the session's temporary directory without consent, and both are right
#' — so this asks in an interactive session and requires `consent = TRUE` otherwise.
#'
#' If you already have a bridge, you do not need this at all: set `MZLIB_BRIDGE` or
#' `options(mzlibr.bridge=)` to point at it. That is also how you would relink a modified mzLib,
#' which LGPL section 4 requires this package to permit.
#'
#' @param version The pyMzLib release to take the payload from. Every pyMzLib release publishes
#'   the raw bridge for each platform as `mzlib-bridge-<rid>.tar.gz`, together with a `SHA256SUMS`
#'   manifest; this fetches the one matching your platform.
#' @param destination Directory to unpack into. Defaults to R's per-user cache directory for
#'   this package, in the same location `tools::R_user_dir()` would choose.
#' @param consent Set `TRUE` to confirm the download in a non-interactive session. Left as `NA`
#'   an interactive session asks, and a non-interactive one refuses.
#' @param overwrite Whether to replace a bridge already present at the destination.
#' @param url An explicit URL, overriding `version`. When you pass this, pass `sha256` too.
#' @param sha256 The expected SHA-256 of the download, lowercase hex. Required with `url`.
#' @param quiet Suppress progress output.
#' @param timeout Seconds to allow for the download. The default is generous because the payload
#'   is large and R's own default of 60 seconds is not enough for it on most connections.
#'
#' @return The path of the installed bridge executable, invisibly.
#' @seealso [mzlibr_bridge_path()], [mzlibr_bridge_version()]
#' @export
mzlibr_install_bridge <- function(version = MZLIB_BRIDGE_VERSION, destination = NULL,
                                  consent = NA, overwrite = FALSE, url = NULL, sha256 = NULL,
                                  quiet = FALSE, timeout = 1800) {
  rid <- bridge_platform_tag()

  if (is.null(url)) {
    known <- MZLIB_BRIDGE_ASSETS[[rid]]
    if (is.null(known)) {
      stop(mzlib_bridge_not_found(paste0(
        "No published bridge for ", rid, ".\n",
        "pyMzLib publishes one for ", paste(names(MZLIB_BRIDGE_ASSETS), collapse = ", "),
        ".\nBuild one with pyMzLib's pkg/build/publish-bridge.ps1 and point MZLIB_BRIDGE at it, ",
        "or pass url= and sha256= explicitly."
      )))
    }
    url <- bridge_release_url(known$asset, version)
    if (is.null(sha256)) {
      sha256 <- known$sha256
    }
  } else if (is.null(sha256)) {
    stop(mzlib_usage_error(
      "When you pass url=, pass sha256= as well. mzLibR will not install an executable it cannot verify."
    ))
  }

  if (is.null(destination)) {
    destination <- mzlibr_cache_dir()
  }
  target_dir <- file.path(destination, rid)
  target <- file.path(target_dir, bridge_executable_name())

  if (file.exists(target) && !overwrite) {
    if (!quiet) {
      message("A bridge is already installed at ", target, ". Pass overwrite = TRUE to replace it.")
    }
    return(invisible(target))
  }

  if (!isTRUE(consent)) {
    bridge_ask_consent(target_dir, "140 MB")
  }

  archive <- tempfile("mzlib-bridge-", fileext = ".tar.gz")
  on.exit(unlink(archive, recursive = TRUE), add = TRUE)

  # R's download timeout defaults to 60 seconds, which a 140 MB payload will not meet.
  previous <- options(timeout = max(timeout, getOption("timeout", 60)))
  on.exit(options(previous), add = TRUE)

  if (!quiet) {
    message("Downloading ", url)
  }
  status <- tryCatch(
    utils::download.file(url, archive, mode = "wb", quiet = quiet),
    error = function(e) {
      stop(mzlib_service_unavailable(paste0(
        "Could not download the bridge from ", url, ": ", conditionMessage(e)
      )))
    }
  )
  if (!identical(as.integer(status), 0L) || !file.exists(archive)) {
    stop(mzlib_service_unavailable(paste0("Could not download the bridge from ", url, ".")))
  }

  observed <- bridge_sha256(archive)
  if (is.null(observed)) {
    stop(mzlib_usage_error(paste0(
      "This R cannot compute a SHA-256, and no sha256sum, shasum or certutil was found on the ",
      "PATH, so the download cannot be verified. mzLibR will not install an unverified ",
      "executable. Upgrade R, install one of those tools, or download the bridge yourself and ",
      "set MZLIB_BRIDGE."
    )))
  }
  if (!identical(tolower(observed), tolower(sha256))) {
    stop(mzlib_protocol_error(paste0(
      "The download does not match its expected checksum.\n",
      "  expected ", tolower(sha256), "\n",
      "  observed ", tolower(observed), "\n",
      "Nothing was installed. This means the file was corrupted in transit, or is not the file ",
      "mzLibR expected."
    )))
  }

  # The tarball's root IS the payload tree - the executable sits at the top alongside Resources/,
  # the Bruker and timsTOF native libraries, and libmmd. So it unpacks straight into place, with no
  # intermediate directory and no copy step. The wheel needed both, because the payload was buried
  # under `pymzlib/_dotnet/<rid>/` and everything above that was Python packaging.
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  extracted <- tryCatch(
    utils::untar(archive, exdir = target_dir),
    error = function(e) {
      stop(mzlib_protocol_error(paste0(
        "Could not unpack the downloaded bridge: ", conditionMessage(e)
      )))
    }
  )
  if (!identical(as.integer(extracted), 0L)) {
    stop(mzlib_protocol_error(paste0(
      "Unpacking the downloaded bridge failed (untar returned ", extracted, ")."
    )))
  }

  # Checked rather than assumed: the checksum proves the bytes are the ones published, not that
  # they contain what this platform needs. A tarball built for the wrong runtime identifier would
  # pass verification and then fail at the first call with a missing-file error.
  if (!file.exists(target)) {
    stop(mzlib_bridge_not_found(paste0(
      "Unpacked the payload but found no ", bridge_executable_name(), " at ", target,
      ". The archive for ", rid, " does not carry the executable at its root; this is a packaging ",
      "problem upstream."
    )))
  }

  bridge_make_runnable(target)

  if (!quiet) {
    message("Installed ", target)
  }
  invisible(target)
}

# Make a freshly unpacked bridge actually runnable on this platform.
#
# The chmod is now belt and braces rather than load-bearing. It existed because `utils::unzip()`
# drops the execute bit, so the binary extracted from a wheel arrived non-executable and failed
# with a permission error that said nothing about the cause. tar records the mode and R restores
# it - the published payload lists `mzlib-bridge` as `-rwxr-xr-x`, verified on the real asset. It
# is kept anyway: it costs nothing, and it covers a user who points `url=` at an archive built by
# something less careful.
#
# The macOS half is still load-bearing and has nothing to do with archive format. A file
# downloaded by an application inherits a quarantine attribute, so Gatekeeper refuses to run it
# and reports that the developer cannot be verified - which sounds like a security problem with
# mzLib rather than a side effect of downloading.
bridge_make_runnable <- function(path) {
  if (.Platform$OS.type == "windows") {
    return(invisible(NULL))
  }

  Sys.chmod(path, mode = "0755", use_umask = FALSE)

  if (identical(Sys.info()[["sysname"]], "Darwin") && nzchar(Sys.which("xattr"))) {
    # Best effort: if the attribute is absent xattr fails, and that is the good case.
    try(
      suppressWarnings(system2(
        "xattr", c("-d", "com.apple.quarantine", shQuote(path)),
        stdout = FALSE, stderr = FALSE
      )),
      silent = TRUE
    )
  }
  invisible(NULL)
}
