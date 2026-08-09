#!/usr/bin/env Rscript

# Regenerates the pinned bridge version and digests in R/install-bridge.R from a pyMzLib
# SHA256SUMS manifest.
#
# Before pyMzLib #31 there was no manifest, so the four digests in install-bridge.R were recorded
# by hand and the comment above them said as much: "bumping this means re-recording all four."
# That is the chore this replaces. #31 publishes SHA256SUMS alongside the wheels and the raw bridge
# tarballs on every release, so the digests can be read from the same authority that produced the
# files rather than transcribed from a browser.
#
# Run by .github/workflows/bridge-watch.yml, and runnable by hand:
#
#   Rscript scripts/regen-bridge-pins.R SHA256SUMS R/install-bridge.R 0.1.0.dev4
#
# It rewrites only the region between the BEGIN/END generated bridge pins markers, and refuses to
# do anything at all if it cannot find both — a script that silently appends when its anchor moves
# is worse than one that stops.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: regen-bridge-pins.R <SHA256SUMS> <install-bridge.R> <version>", call. = FALSE)
}
sums_path <- args[[1L]]
target_path <- args[[2L]]
version <- args[[3L]]

# The platform wheel that carries each runtime identifier's bridge. These tags are set by
# pyMzLib's wheels.yml build matrix; if that matrix changes its plat_tag values, this table is what
# needs to change with it, and the "no digest for" error below is what will say so.
TAGS <- c(
  "win-x64"   = "win_amd64",
  "osx-arm64" = "macosx_12_0_arm64",
  "osx-x64"   = "macosx_12_0_x86_64",
  "linux-x64" = "manylinux_2_28_x86_64"
)

# --- read the manifest -----------------------------------------------------------------------

lines <- trimws(readLines(sums_path, warn = FALSE))
lines <- lines[nzchar(lines)]

# `sha256sum` writes "<digest>  <name>", and "<digest> *<name>" when it was run in binary mode.
# Accept both rather than assuming which side produced the file.
matched <- regmatches(lines, regexec("^([0-9a-fA-F]{64})[[:space:]]+[*]?(.+)$", lines))
matched <- Filter(function(x) length(x) == 3L, matched)
if (length(matched) == 0L) {
  stop("No 'sha256  filename' lines found in ", sums_path, call. = FALSE)
}

digests <- stats::setNames(
  tolower(vapply(matched, `[`, character(1L), 2L)),
  vapply(matched, `[`, character(1L), 3L)
)

# --- build the replacement block -------------------------------------------------------------

entry <- function(rid) {
  wheel <- sprintf("pymzlib-%s-py3-none-%s.whl", version, TAGS[[rid]])
  if (!wheel %in% names(digests)) {
    stop(
      "No digest for ", wheel, " in ", sums_path, ".\n",
      "The manifest lists: ", paste(names(digests), collapse = ", "),
      call. = FALSE
    )
  }
  sprintf(
    '  "%s" = list(\n    wheel = "%s",\n    sha256 = "%s"\n  )',
    rid, wheel, digests[[wheel]]
  )
}

# One element per LINE, not per entry. `entry()` returns a four-line chunk, and pasting those into
# a single string would give a character vector whose elements span newlines: writeLines() flattens
# it to exactly the right file, so the output looks correct, but the `identical()` no-op check
# below then compares a 4-element vector against the 16-element one readLines() returns and reports
# a change on every run. Caught only by running it — the generated file was byte-identical and the
# script still said it had repinned.
entries <- paste(vapply(names(TAGS), entry, character(1L)), collapse = ",\n")
block <- c(
  "# BEGIN generated bridge pins",
  sprintf('MZLIB_BRIDGE_VERSION <- "%s"', version),
  "",
  "MZLIB_BRIDGE_WHEELS <- list(",
  strsplit(entries, "\n", fixed = TRUE)[[1L]],
  ")",
  "# END generated bridge pins"
)

# --- splice it in ----------------------------------------------------------------------------

source_lines <- readLines(target_path, warn = FALSE)
begin <- which(source_lines == "# BEGIN generated bridge pins")
end <- which(source_lines == "# END generated bridge pins")

if (length(begin) != 1L || length(end) != 1L || end <= begin) {
  stop(
    "Expected exactly one BEGIN and one END 'generated bridge pins' marker in ", target_path,
    " (found ", length(begin), " and ", length(end), "). Nothing was written.",
    call. = FALSE
  )
}

updated <- c(
  source_lines[seq_len(begin - 1L)],
  block,
  source_lines[seq.int(end + 1L, length(source_lines))]
)

if (identical(updated, source_lines)) {
  cat("Already pinned to", version, "- nothing to do.\n")
  quit(status = 0L)
}

writeLines(updated, target_path)
cat("Repinned", target_path, "to pyMzLib", version, "\n")
