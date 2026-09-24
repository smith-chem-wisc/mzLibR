# Stage the recordings the help-page examples replay into inst/replay/.
#
#   Rscript scripts/stage-replay.R            # (re)write inst/replay/
#   Rscript scripts/stage-replay.R --check    # CI: fail if that would change anything
#
# The examples run under `R CMD check` against a stand-in bridge (R/replay.R) that answers from
# recordings of the real one. Examples run from the INSTALLED package, where tests/ does not exist,
# so the recordings they use are copied into inst/replay/, byte for byte from tests/fixtures/ -
# which holds pyMzLib's fixtures, byte for byte. Only the ones an example may use are shipped.
#
# Which recording may answer which verb comes from the verb specs' `examples` (docs/specs/), the
# same list pyMzLib's doctests replay, plus REPLAY_EXTRA below for recordings no spec names.
# inst/replay/TABLE records the result, one `verb<TAB>fixture` line each.
#
# Base R only, apart from reading the specs, which needs mzLibR installed (scripts/spec-facts.R).

suppressMessages(library(mzLibR))
source(file.path("scripts", "spec-facts.R"))

check <- identical(commandArgs(trailingOnly = TRUE), "--check")
fixtures <- file.path("tests", "fixtures")
dest <- file.path("inst", "replay")

# Recordings an example needs that no spec lists as an example come from REPLAY_EXTRA, kept per
# module in scripts/deviations/<module>.R, each with the reason.

rows <- list()
for (spec in load_specs()) {
  examples <- spec$examples
  if (!is.list(examples)) next
  for (example in examples) {
    rows[[length(rows) + 1L]] <- c(verb = spec$verb, fixture = spec_text(example$fixture))
  }
}
rows <- c(rows, REPLAY_EXTRA)

table <- unique(do.call(rbind, rows))
present <- file.exists(file.path(fixtures, table[, "fixture"]))
for (missing in unique(table[!present, "fixture"])) {
  cat("not staged: ", missing, " (not in tests/fixtures yet)\n", sep = "")
}
table <- table[present, , drop = FALSE]
# method = "radix" sorts in the C locale, so the table is the same on every machine.
table <- table[order(table[, "verb"], table[, "fixture"], method = "radix"), , drop = FALSE]

table_text <- c(
  "# Written by scripts/stage-replay.R from the verb specs' examples. Do not edit; re-run it.",
  "# verb<TAB>recording: the recordings the stand-in bridge (R/replay.R) may answer each verb with.",
  paste(table[, "verb"], table[, "fixture"], sep = "\t")
)
wanted <- unique(table[, "fixture"])

read_bytes <- function(path) {
  if (!file.exists(path)) raw(0) else readBin(path, what = "raw", n = file.info(path)$size)
}

problems <- character(0)
if (check) {
  have <- setdiff(list.files(dest), "TABLE")
  for (extra in setdiff(have, wanted)) {
    problems <- c(problems, paste0(extra, ": staged but no example uses it"))
  }
  for (name in wanted) {
    if (!identical(read_bytes(file.path(dest, name)), read_bytes(file.path(fixtures, name)))) {
      problems <- c(problems, paste0(name, ": differs from tests/fixtures (or is missing)"))
    }
  }
  current <- if (file.exists(file.path(dest, "TABLE"))) readLines(file.path(dest, "TABLE"), warn = FALSE) else character(0)
  if (!identical(current, table_text)) {
    problems <- c(problems, "TABLE: differs from what the specs ask for")
  }
  if (length(problems) > 0L) {
    cat(paste0("  ! ", problems, "\n"), sep = "")
    cat("inst/replay/ is stale. Run: Rscript scripts/stage-replay.R\n")
    quit(status = 1L)
  }
  cat(length(wanted), "recordings staged for", length(unique(table[, "verb"])), "verbs, all current\n")
  quit(status = 0L)
}

dir.create(dest, recursive = TRUE, showWarnings = FALSE)
for (extra in setdiff(setdiff(list.files(dest), "TABLE"), wanted)) {
  unlink(file.path(dest, extra))
  cat("removed ", extra, "\n", sep = "")
}
for (name in wanted) {
  file.copy(file.path(fixtures, name), file.path(dest, name), overwrite = TRUE)
}
con <- file(file.path(dest, "TABLE"), open = "wb")
writeBin(charToRaw(paste0(paste(table_text, collapse = "\n"), "\n")), con)
close(con)
cat(length(wanted), "recordings staged for", length(unique(table[, "verb"])), "verbs\n")
