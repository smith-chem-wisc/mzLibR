# Generate man/*.Rd from the `#'` blocks in R/*.R.
#
# roxygen2 is an excellent tool and it is also a dependency, which this package does not have —
# not for users, not for developers. But hand-written .Rd files sitting apart from the code are
# how documentation rots, and PLAN.md section 8.3 is emphatic that rotting documentation is the
# expensive kind. So the `#'` blocks stay next to the functions they describe, and this converts
# them.
#
# It implements only the subset mzLibR actually uses: @param, @return, @section, @seealso,
# @export, @examples, and the title/description/details paragraph convention. It is not roxygen2
# and does not try to be; anything it does not understand is an error rather than a silent
# omission.
#
# Run, from the package root, with mzLibR installed:
#   Rscript scripts/build-man.R

suppressMessages(library(mzLibR))

sources <- list.files("R", pattern = "[.]R$", full.names = TRUE)
if (length(sources) == 0L) {
  stop("run this from the package root")
}
dir.create("man", showWarnings = FALSE)

KNOWN_TAGS <- c("param", "return", "section", "seealso", "export", "examples", "noRd")

# ---------------------------------------------------------------- inline markup

# Turn the Markdown used in the `#'` blocks into Rd markup.
#
# The escaping has to happen first and only on the characters Rd treats specially, because doing
# it after would mangle the markup this function itself emits.
rd_inline <- function(text) {
  text <- gsub("\\", "\\\\", text, fixed = TRUE)
  # Percent begins a comment in Rd, and this documentation is full of them — "a 63% under-count",
  # "loses about 30%". Unescaped, the rest of the line silently disappears.
  text <- gsub("%", "\\%", text, fixed = TRUE)
  text <- gsub("{", "\\{", text, fixed = TRUE)
  text <- gsub("}", "\\}", text, fixed = TRUE)

  # `[fn()]` and `[fn]` become links; do this before code spans so the brackets are still there.
  text <- gsub("\\[([a-zA-Z_][a-zA-Z0-9_.]*)\\(\\)\\]", "\\\\code\\{\\\\link\\{\\1\\}\\}", text)
  text <- gsub("\\[([a-zA-Z_][a-zA-Z0-9_.]*)\\]", "\\\\code\\{\\\\link\\{\\1\\}\\}", text)

  text <- gsub("\\*\\*([^*]+)\\*\\*", "\\\\strong\\{\\1\\}", text)
  text <- gsub("`([^`]+)`", "\\\\code\\{\\1\\}", text)
  text
}

# ---------------------------------------------------------------- block parsing

# Split a block's lines into the leading prose and the tagged sections.
parse_block <- function(lines, target) {
  prose <- character(0)
  tags <- list()
  current <- NULL

  for (line in lines) {
    tag_match <- regmatches(line, regexec("^@([a-zA-Z]+)\\s*(.*)$", line))[[1]]
    if (length(tag_match) == 3L) {
      tag <- tag_match[2]
      if (!tag %in% KNOWN_TAGS) {
        stop("unknown tag @", tag, " in the block for ", target)
      }
      current <- list(tag = tag, lines = tag_match[3])
      tags[[length(tags) + 1L]] <- current
      next
    }
    if (is.null(current)) {
      prose <- c(prose, line)
    } else {
      tags[[length(tags)]]$lines <- c(tags[[length(tags)]]$lines, line)
    }
  }

  list(prose = prose, tags = tags)
}

# Group lines into paragraphs on blank lines.
paragraphs <- function(lines) {
  out <- list()
  buffer <- character(0)
  for (line in c(lines, "")) {
    if (!nzchar(trimws(line))) {
      if (length(buffer) > 0L) {
        out[[length(out) + 1L]] <- paste(buffer, collapse = "\n")
        buffer <- character(0)
      }
    } else {
      buffer <- c(buffer, line)
    }
  }
  out
}

# ---------------------------------------------------------------- usage

# The usage line, taken from the installed function so it cannot disagree with the code.
# The usage block, taken from the installed function so it cannot disagree with the code.
#
# `deparse()`'s own line breaks are kept rather than collapsed. A long signature on one line is
# an "Rd line widths" NOTE from R CMD check, and worse, it is truncated in the PDF manual —
# which silently hides the tail of exactly the signatures long enough to need reading.
# Split a parameter list on the commas that separate parameters, ignoring commas nested inside a
# default value such as `c("a", "b")`.
split_parameters <- function(text) {
  characters <- strsplit(text, "", fixed = TRUE)[[1L]]
  depth <- 0L
  quoted <- FALSE
  parts <- character(0)
  buffer <- character(0)

  for (character in characters) {
    if (quoted) {
      buffer <- c(buffer, character)
      if (character == "\"") quoted <- FALSE
      next
    }
    if (character == "\"") {
      quoted <- TRUE
    } else if (character %in% c("(", "[", "{")) {
      depth <- depth + 1L
    } else if (character %in% c(")", "]", "}")) {
      depth <- depth - 1L
    } else if (character == "," && depth == 0L) {
      parts <- c(parts, paste(buffer, collapse = ""))
      buffer <- character(0)
      next
    }
    buffer <- c(buffer, character)
  }
  trimws(c(parts, paste(buffer, collapse = "")))
}

# The usage block, taken from the installed function so it cannot disagree with the code.
#
# Wrapped by hand rather than left to `deparse()`, whose `width.cutoff` is a soft target that a
# single long parameter name goes straight past. Over 90 characters is an "Rd line widths" NOTE
# from R CMD check, and worse, the line is truncated in the PDF manual — silently hiding the
# tail of exactly the signatures long enough to need reading.
usage_for <- function(name) {
  fn <- get(name, envir = asNamespace("mzLibR"))
  signature <- paste(deparse(args(fn), width.cutoff = 500L), collapse = " ")
  signature <- sub("^function\\s*\\(", "", signature)
  signature <- trimws(sub("\\)\\s*NULL\\s*$", "", signature))

  method <- regmatches(name, regexec("^([a-zA-Z0-9.]+?)\\.(mzlibr_[a-z_]+)$", name))[[1]]
  head <- if (length(method) == 3L) {
    paste0("\\method{", method[2], "}{", method[3], "}(")
  } else {
    paste0(name, "(")
  }

  parameters <- split_parameters(signature)
  parameters <- parameters[nzchar(parameters)]
  if (length(parameters) == 0L) {
    return(paste0(head, ")"))
  }

  lines <- character(0)
  current <- head
  for (index in seq_along(parameters)) {
    piece <- paste0(parameters[index], if (index < length(parameters)) "," else ")")
    candidate <- if (identical(current, head) || endsWith(current, "(")) {
      paste0(current, piece)
    } else {
      paste0(current, " ", piece)
    }
    if (nchar(candidate) > 88L && !endsWith(current, "(")) {
      lines <- c(lines, current)
      current <- paste0("  ", piece)
    } else {
      current <- candidate
    }
  }
  paste(c(lines, current), collapse = "\n")
}

# ---------------------------------------------------------------- emit

written <- character(0)

for (source in sources) {
  lines <- readLines(source, warn = FALSE)
  block <- character(0)

  for (index in seq_along(lines)) {
    line <- lines[index]

    if (grepl("^#'", line)) {
      block <- c(block, sub("^#'\\s?", "", line))
      next
    }

    definition <- regmatches(line, regexec("^([a-zA-Z0-9._]+)\\s*<-\\s*function", line))[[1]]
    if (length(block) > 0L && length(definition) == 2L) {
      name <- definition[2]
      parsed <- parse_block(block, name)
      tag_names <- vapply(parsed$tags, function(entry) entry$tag, character(1L))

      if (!"export" %in% tag_names || "noRd" %in% tag_names) {
        block <- character(0)
        next
      }

      chunks <- paragraphs(parsed$prose)
      if (length(chunks) == 0L) {
        stop("no title in the block for ", name)
      }
      title <- rd_inline(gsub("\n", " ", chunks[[1]]))
      description <- if (length(chunks) >= 2L) chunks[[2]] else chunks[[1]]
      details <- if (length(chunks) > 2L) chunks[-(1:2)] else list()

      out <- c(
        "% Generated by scripts/build-man.R from the roxygen-style block in the source.",
        "% Do not edit by hand; edit R/ and re-run.",
        paste0("\\name{", name, "}"),
        paste0("\\alias{", name, "}"),
        paste0("\\title{", title, "}"),
        "\\usage{",
        usage_for(name),
        "}"
      )

      params <- Filter(function(entry) entry$tag == "param", parsed$tags)
      if (length(params) > 0L) {
        out <- c(out, "\\arguments{")
        for (entry in params) {
          # Split on the first run of whitespace of the FIRST line rather than of the whole
          # block: `[\s\S]` is a Perl idiom and R's default engine does not read it that way,
          # so a multi-line @param quietly failed to match at all.
          first <- entry$lines[1L]
          if (!grepl("^\\S+\\s", first)) {
            stop("malformed @param in ", name, ": ", substr(first, 1L, 60L))
          }
          parameter <- sub("^(\\S+)\\s.*$", "\\1", first)
          body <- paste(c(sub("^\\S+\\s+", "", first), entry$lines[-1L]), collapse = "\n")
          out <- c(out, paste0("\\item{", parameter, "}{", rd_inline(body), "}"), "")
        }
        out <- c(out, "}")
      }

      out <- c(out, "\\description{", rd_inline(description), "}")

      if (length(details) > 0L) {
        out <- c(out, "\\details{", rd_inline(paste(unlist(details), collapse = "\n\n")), "}")
      }

      for (entry in parsed$tags) {
        if (entry$tag == "section") {
          first <- entry$lines[1L]
          if (!grepl(":", first, fixed = TRUE)) {
            stop("malformed @section in ", name, ": the heading must end with a colon")
          }
          heading <- sub(":.*$", "", first)
          body <- paste(c(sub("^[^:]*:\\s*", "", first), entry$lines[-1L]), collapse = "\n")
          out <- c(
            out,
            paste0("\\section{", rd_inline(heading), "}{"),
            rd_inline(body),
            "}"
          )
        }
        if (entry$tag == "return") {
          out <- c(out, "\\value{", rd_inline(paste(entry$lines, collapse = "\n")), "}")
        }
        if (entry$tag == "seealso") {
          out <- c(out, "\\seealso{", rd_inline(paste(entry$lines, collapse = "\n")), "}")
        }
        if (entry$tag == "examples") {
          out <- c(out, "\\examples{", paste(entry$lines, collapse = "\n"), "}")
        }
      }

      writeLines(out, file.path("man", paste0(name, ".Rd")))
      written <- c(written, name)
      block <- character(0)
      next
    }

    if (!grepl("^\\s*$", line)) {
      block <- character(0)
    }
  }
}

# Every exported object must be documented or `R CMD check` complains, and an undocumented
# export is a real gap rather than a formality.
exported <- getNamespaceExports("mzLibR")
undocumented <- setdiff(exported, written)
if (length(undocumented) > 0L) {
  stop("exported but undocumented: ", paste(sort(undocumented), collapse = ", "))
}

cat("wrote", length(written), "Rd files\n")
