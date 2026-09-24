# Copy the per-verb specs from a bridge checkout into docs/specs/, and record where they came from.
#
#   Rscript scripts/sync-specs.R --from ../bridge/design/verbs
#   Rscript scripts/sync-specs.R --check        # CI: the .json copies still match the .yaml
#
# The specs (one YAML file per wire verb: params, units, result fields, what NA means, errors,
# caveats, citations, each binding's spelling) are written once, in the bridge repository, so that
# pyMzLib, mzLibRust and mzLibR cannot drift apart on facts. That repository is private and this
# one is public, so CI here cannot read it. The copy in docs/specs/ is what scripts/build-man.R
# renders into man/, what scripts/name-parity.R lints the help text against, and what
# scripts/stage-replay.R reads to pick the recordings the examples replay.
#
# **Never edit docs/specs/ by hand.** A fact that is wrong is wrong in the bridge: fix it there, then
# run this. Every binding then picks up the correction, and every binding's lint flags the pages
# that repeated the error.
#
# Why a .json beside every .yaml. R has no YAML reader in base, and this package has no
# dependencies - not for users, and not for developers running build-man.R. So this script, which
# only a maintainer syncing specs ever runs, is the one place that needs the `yaml` package: it
# writes each spec a second time as JSON, which the package's own reader (R/json.R) already parses.
# `--check` re-derives every .json from its .yaml and fails on a difference, so the two cannot
# drift; CI installs `yaml` for that job only.
#
# Only specs the bridge has COMMITTED are copied (git ls-files), and SOURCE records the commit, so
# a reviewer can tie every vendored spec to a bridge commit rather than to someone's working tree.

args <- commandArgs(trailingOnly = TRUE)
dest <- file.path("docs", "specs")
if (!dir.exists("R") || !file.exists("DESCRIPTION")) {
  stop("run this from the package root")
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "scripts/sync-specs.R needs the 'yaml' package: install.packages(\"yaml\"). ",
    "Nothing else in mzLibR does."
  )
}

# ---------------------------------------------------------------- YAML in, canonical JSON out

# Every YAML sequence becomes an R list, never an atomic vector, so `[a]` stays an array of one
# rather than collapsing into the scalar `a` - the distinction the specs' `views: [spectra]` needs.
read_spec <- function(path) {
  yaml::read_yaml(
    path,
    fileEncoding = "UTF-8",
    handlers = list(
      seq = function(x) as.list(x),
      map = function(x) if (length(x) == 0L) structure(list(), names = character(0)) else x
    )
  )
}

json_string <- function(x) {
  x <- enc2utf8(x)
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

json_scalar <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) {
    return("null")
  }
  if (is.logical(x)) {
    return(if (x) "true" else "false")
  }
  if (is.integer(x)) {
    return(as.character(x))
  }
  if (is.numeric(x)) {
    if (x == round(x) && abs(x) < 1e15) {
      return(formatC(x, format = "d", big.mark = ""))
    }
    return(formatC(x, digits = 15, format = "g"))
  }
  json_string(as.character(x))
}

# Two-space indentation and keys in the spec's own order, so a diff of a .json reads like a diff
# of its .yaml.
to_json <- function(x, indent = "") {
  inner <- paste0(indent, "  ")
  if (is.list(x)) {
    if (!is.null(names(x))) {
      if (length(x) == 0L) {
        return("{}")
      }
      items <- vapply(seq_along(x), function(i) {
        paste0(inner, json_string(names(x)[i]), ": ", to_json(x[[i]], inner))
      }, character(1L))
      return(paste0("{\n", paste(items, collapse = ",\n"), "\n", indent, "}"))
    }
    if (length(x) == 0L) {
      return("[]")
    }
    items <- vapply(x, function(item) paste0(inner, to_json(item, inner)), character(1L))
    return(paste0("[\n", paste(items, collapse = ",\n"), "\n", indent, "]"))
  }
  if (length(x) > 1L) {
    return(to_json(as.list(x), indent))
  }
  json_scalar(x)
}

write_utf8 <- function(text, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(paste0(text, "\n"))), con)
}

read_bytes <- function(path) {
  if (!file.exists(path)) {
    return(raw(0))
  }
  readBin(path, what = "raw", n = file.info(path)$size)
}

json_for <- function(yaml_path) {
  to_json(read_spec(yaml_path))
}

# ---------------------------------------------------------------- --check

if (identical(args, "--check")) {
  specs <- sort(list.files(dest, pattern = "[.]yaml$", full.names = TRUE))
  problems <- character(0)
  for (spec in specs) {
    target <- sub("[.]yaml$", ".json", spec)
    want <- charToRaw(enc2utf8(paste0(json_for(spec), "\n")))
    if (!identical(read_bytes(target), want)) {
      problems <- c(problems, paste0(basename(target), ": differs from ", basename(spec)))
    }
  }
  strays <- setdiff(
    sub("[.]json$", "", list.files(dest, pattern = "[.]json$")),
    sub("[.]yaml$", "", basename(specs))
  )
  if (length(strays) > 0L) {
    problems <- c(problems, paste0(strays, ".json: no .yaml beside it"))
  }
  if (length(problems) > 0L) {
    cat(paste0("  ! ", problems, "\n"), sep = "")
    cat("docs/specs/*.json is stale. Run: Rscript scripts/sync-specs.R --from <bridge>/design/verbs\n")
    quit(status = 1L)
  }
  cat(length(specs), "specs: every .json matches its .yaml\n")
  quit(status = 0L)
}

# ---------------------------------------------------------------- --from

from <- which(args == "--from")
if (length(from) != 1L || from == length(args)) {
  stop("usage: Rscript scripts/sync-specs.R --from <bridge>/design/verbs  |  --check")
}
source_dir <- normalizePath(args[from + 1L], winslash = "/", mustWork = TRUE)

git <- function(...) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(source_dir), ...), stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
  if (!is.null(attr(out, "status"))) character(0) else out
}

tracked <- git("ls-files", "--", "*.yaml")
names_wanted <- if (length(tracked) > 0L) {
  sort(basename(tracked))
} else {
  sort(list.files(source_dir, pattern = "[.]yaml$"))
}
if (length(names_wanted) == 0L) {
  stop("no *.yaml specs in ", source_dir, "; is that a bridge's design/verbs?")
}
untracked <- setdiff(list.files(source_dir, pattern = "[.]yaml$"), names_wanted)
for (name in untracked) {
  cat("skipped ", name, " (not committed in the bridge)\n", sep = "")
}

dir.create(dest, recursive = TRUE, showWarnings = FALSE)
for (stale in setdiff(list.files(dest, pattern = "[.](yaml|json)$"),
                      c(names_wanted, sub("[.]yaml$", ".json", names_wanted)))) {
  unlink(file.path(dest, stale))
  cat("removed ", stale, " (no longer in the bridge)\n", sep = "")
}

for (name in names_wanted) {
  target <- file.path(dest, name)
  before <- read_bytes(target)
  file.copy(file.path(source_dir, name), target, overwrite = TRUE)
  if (!identical(before, read_bytes(target))) {
    cat(if (length(before) > 0L) "updated " else "added   ", name, "\n", sep = "")
  }
  write_utf8(json_for(target), sub("[.]yaml$", ".json", target))
}

commit <- git("rev-parse", "HEAD")
commit <- if (length(commit) == 1L) commit else "unknown (not a git checkout)"
dirty <- git("status", "--porcelain", "--untracked-files=no", "--", ".")
write_utf8(paste(c(
  "# Written by scripts/sync-specs.R. Do not edit; re-run the script.",
  "repository: trishorts/bridge (private)",
  "path: design/verbs",
  paste0("commit: ", commit),
  paste0("uncommitted_changes: ", if (length(dirty) > 0L) "yes" else "no"),
  paste0("synced: ", format(Sys.Date())),
  paste0("specs: ", length(names_wanted))
), collapse = "\n"), file.path(dest, "SOURCE"))

cat(length(names_wanted), " specs from ", substr(commit, 1L, 12L),
  if (length(dirty) > 0L) " (with uncommitted changes)" else "", "\n", sep = "")
cat("Now run: R CMD INSTALL . && Rscript scripts/build-man.R\n")
