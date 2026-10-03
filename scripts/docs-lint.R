# Hold the articles and help pages to the promises the docs make about themselves.
#
#   Rscript scripts/docs-lint.R    # CI: exits non-zero on any finding
#
# The R counterpart of pyMzLib's tests/test_docs_lint.py:
#
#   * Every R chunk in an article runs, or says why not. The articles are vignettes, so R CMD check
#     runs every chunk against the replay bridge. A chunk with `eval = FALSE` must open with a
#     "# Not run: <reason>" comment, and an R block that knitr would not run (```r) is not allowed.
#   * Every article opens with a question -> function -> mzLib table, so a reader can find the call
#     for their question before reading the prose.
#   * Every article ends with what to cite, rendered from the specs by scripts/render-cite.R.
#   * No counts of formats or verbs, and no mzLib versions, in prose: articles, the README and the
#     help-page text in R/. "36 formats" was once written in seven places across the bindings and
#     went stale in all of them at one mzLib release. A count belongs in executed output; a version
#     belongs in NEWS.md.
#
# Base R only.

if (!dir.exists("R") || !file.exists("DESCRIPTION")) {
  stop("run this from the package root")
}

COUNT <- paste0(
  "\\b([0-9]+|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|nineteen)",
  "\\s+(file types|file formats|formats|wire verbs|verbs|readers)\\b"
)
VERSION <- "\\bmzLib\\s+v?1\\.0\\.[0-9]+"
CITE_BEGIN <- "<!-- BEGIN cite:"

problems <- character(0)
note <- function(path, line, message) {
  problems <<- c(problems, paste0(path, ":", line, ": ", message))
}

# Lines outside fenced blocks and HTML comments, with their numbers: what a reader reads as claims.
prose_lines <- function(lines) {
  inside <- FALSE
  keep <- logical(length(lines))
  for (i in seq_along(lines)) {
    if (grepl("^\\s*```", lines[i])) {
      inside <- !inside
      next
    }
    keep[i] <- !inside && !grepl("^\\s*<!--", lines[i])
  }
  which(keep)
}

lint_prose <- function(path, lines, numbers) {
  for (i in numbers) {
    for (pattern in c(COUNT, VERSION)) {
      found <- regmatches(lines[i], gregexpr(pattern, lines[i], ignore.case = identical(pattern, COUNT), perl = TRUE))[[1L]]
      for (hit in found) {
        note(path, i, paste0(if (identical(pattern, COUNT)) "a count of formats or verbs" else "an mzLib version",
          " in prose: '", hit, "'"))
      }
    }
  }
}

# ---------------------------------------------------------------- the articles

for (path in sort(list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE))) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  # The YAML header ends at the second `---`; the table must come within 40 lines of it.
  dashes <- which(lines == "---")
  body_start <- if (length(dashes) >= 2L) dashes[2L] + 1L else 1L
  head <- lines[body_start:min(length(lines), body_start + 39L)]
  if (!any(grepl("^\\| question \\|", head))) {
    note(path, body_start, "open with a '| question | function | mzLib |' table within 40 lines")
  }
  if (!any(startsWith(lines, CITE_BEGIN))) {
    note(path, length(lines), "end with a 'What to cite' section (scripts/render-cite.R)")
  }

  for (i in grep("^\\s*```", lines)) {
    fence <- trimws(lines[i])
    if (grepl("^```\\s*[rR]\\s*$", fence)) {
      note(path, i, "an R block knitr will not run: write ```{r} so it runs, or ```text for output")
    }
    if (grepl("^```\\{r.*eval\\s*=\\s*FALSE", fence)) {
      first <- if (i < length(lines)) trimws(lines[i + 1L]) else ""
      if (!startsWith(first, "# Not run:")) {
        note(path, i, "an eval = FALSE chunk must open with '# Not run: <reason>'")
      }
    }
  }
  lint_prose(path, lines, prose_lines(lines))
}

# ---------------------------------------------------------------- README and help-page text

readme <- readLines("README.md", warn = FALSE, encoding = "UTF-8")
lint_prose("README.md", readme, prose_lines(readme))

for (path in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  roxygen <- which(startsWith(lines, "#'") & !grepl("^#' *(\\\\dontshow|@examples)", lines))
  lint_prose(path, lines, roxygen)
}

if (length(problems) > 0L) {
  cat(paste0("  ! ", problems, "\n"), sep = "")
  cat(length(problems), "docs finding(s). See scripts/docs-lint.R for why each rule exists.\n")
  quit(status = 1L)
}
cat("docs lint: every article runs, opens with its question table and ends with what to cite;",
  "no counts or versions in prose\n")
