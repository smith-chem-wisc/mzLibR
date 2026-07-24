# Generate docs/name-parity.md by comparing mzLibR's names against pyMzLib's and the wire's.
#
# pyMzLib is the parent binding; mzLibRust derives from it and mzLibR from both. Names therefore
# flow downward, and any place mzLibR spells something differently is either a defect or a
# decision that has to be written down. This script finds them mechanically so the answer does
# not depend on anyone remembering.
#
# Two halves, because there are two ways to get a name wrong:
#
#   1. **Functions and their parameters**, checked against pyMzLib's module functions *and* its
#      class methods — an R package has no per-module namespace and no methods, so
#      `Digest.truncated` becomes `digest_truncated()`, and that flattening is the only renaming
#      mzLibR is entitled to do.
#   2. **data.frame columns**, checked against the field names on the wire, which are themselves
#      the snake_case of mzLib's own names. This is the half that matters most day to day: a
#      column is what a user types.
#
# Run:
#   Rscript scripts/name-parity.R > docs/name-parity.md

args <- commandArgs(trailingOnly = TRUE)
python_src <- if (length(args) >= 1L) {
  args[1L]
} else {
  "E:/CodeReview/pyMzLib/code/pyMzLib/pkg/python/src/pymzlib"
}
fixtures <- if (length(args) >= 2L) args[2L] else "tests/fixtures"

modules <- c("pride", "peptidoform", "flashlfq")

# ---------------------------------------------------------------- the Python side

# Parameter names out of a `def ...(...)` signature beginning at `start`.
python_parameters <- function(lines, start) {
  depth <- 0L
  body <- character(0)
  for (index in seq(start, length(lines))) {
    body <- c(body, lines[index])
    opens <- lengths(regmatches(lines[index], gregexpr("(", lines[index], fixed = TRUE)))
    closes <- lengths(regmatches(lines[index], gregexpr(")", lines[index], fixed = TRUE)))
    depth <- depth + opens - closes
    if (depth <= 0L) break
  }

  text <- paste(body, collapse = "\n")
  inner <- sub("^\\s*def [a-z_][a-z0-9_]*\\(", "", text)
  inner <- sub("\\)[^)]*$", "", inner)
  parameters <- trimws(strsplit(inner, ",\n|,(?![^\\[]*\\])", perl = TRUE)[[1]])
  parameters <- sub("[:=].*$", "", parameters)
  parameters <- trimws(gsub("[*]", "", parameters))
  parameters[nzchar(parameters) & !parameters %in% c("self", "cls")]
}

# Every callable in a module, keyed `module.function` or `Class.method`.
python_signatures <- function(path, module) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  current_class <- NA_character_

  for (index in seq_along(lines)) {
    line <- lines[index]
    if (grepl("^class [A-Za-z_]", line)) {
      current_class <- sub("^class ([A-Za-z_][A-Za-z0-9_]*).*$", "\\1", line)
      next
    }
    if (grepl("^def [a-z_][a-z0-9_]*\\(", line)) {
      name <- sub("^def ([a-z_][a-z0-9_]*)\\(.*$", "\\1", line)
      out[[paste0(module, ".", name)]] <- python_parameters(lines, index)
      next
    }
    if (grepl("^    def [a-z_][a-z0-9_]*\\(", line) && !is.na(current_class)) {
      name <- sub("^    def ([a-z_][a-z0-9_]*)\\(.*$", "\\1", line)
      if (startsWith(name, "_")) next
      out[[paste0(current_class, ".", name)]] <- python_parameters(lines, index)
    }
  }
  out
}

parent <- list()
for (module in modules) {
  path <- file.path(python_src, paste0(module, ".py"))
  if (file.exists(path)) {
    parent <- c(parent, python_signatures(path, module))
  }
}

# ---------------------------------------------------------------- the mapping

# R flattens pyMzLib's classes away, so the correspondence has to be stated. Anything not named
# here is reported as unmapped rather than quietly passing.
mapping <- c(
  pride_list_files = "pride.list_files",
  pride_download = "pride.download",
  pride_download_files = "pride.download_files",
  pride_total_size_bytes = "pride.total_size_bytes",
  peptidoform_fragments = "peptidoform.fragments",
  digest_truncated = "Digest.truncated",
  digest_modified_peptides = "Digest.modified_peptides",
  census_explain = "ModificationCensus.explain",
  census_excluded = "ModificationCensus.excluded",
  peptide_mz = "Peptide.mz",
  flashlfq_quantify = "flashlfq.quantify",
  flashlfq_peptide_count = "FlashLfqResults.peptide_count",
  flashlfq_protein_count = "FlashLfqResults.protein_count",
  flashlfq_mbr_peak_count = "FlashLfqResults.mbr_peak_count",
  flashlfq_mbr_peaks = "FlashLfqResults.mbr_peaks",
  flashlfq_mbr_rescued_peptide_count = "FlashLfqResults.mbr_rescued_peptide_count"
)

# Deliberately absent from pyMzLib, with the reason. Anything here is a considered addition, not
# an oversight.
additions <- c(
  digest_distinct_base_sequences = "from mzLibRust (Digest::distinct_base_sequences); pyMzLib has no equivalent",
  digest_fragments_by_series = "from mzLibRust (Digest::fragments_by_series); pyMzLib has only fragment_count",
  pride_locations = "R-only: unnests the `locations` list column, which is pyMzLib's PrideFile.locations field",
  mzlibr_bridge_path = "transport; pyMzLib's equivalent is private (_bridge)",
  mzlibr_bridge_version = "transport; pyMzLib's equivalent is private (_bridge)",
  mzlibr_install_bridge = paste(
    "R-only: Python ships the payload inside the wheel and Rust downloads it from build.rs.",
    "CRAN allows neither, so the download has to be a function the user calls."
  )
)

# pyMzLib callables with no mzLibR counterpart, with the reason.
omissions <- c(
  "PrideFile.size_mb" = "a column, `size_mb`, not a function",
  "PrideFile.extension" = "a column, `extension`, not a function",
  "PrideFile.downloadable" = "a column, `downloadable`, not a function",
  "PrideFile.as_dict" = "meaningless in R: the data.frame is already the record",
  "Peptide.is_modified" = "a column comparison, `modification_count > 0`",
  "Peptide.intensity" = "a row of the long `peptides` frame",
  "Peptide.detection_type" = "a column of the long `peptides` frame",
  "ProteinGroup.intensity" = "a row of the long `proteins` frame",
  "Peak.is_mbr" = "a column comparison, `detection_type == \"MBR\"`",
  "Digest.fragment_count" = "`nrow(digest$fragments)`; not promoted to a function because a bare total folds in the spurious ETD y series — see ?digest_fragments_by_series"
)

# ---------------------------------------------------------------- the R side

suppressMessages(library(mzLibR))
exported <- sort(getNamespaceExports("mzLibR"))
exported <- exported[!grepl("^print\\.", exported)]

r_parameters <- function(name) {
  formal_names <- names(formals(get(name, envir = asNamespace("mzLibR"))))
  formal_names[formal_names != "..."]
}

# ---------------------------------------------------------------- report: functions

cat("# Name parity with pyMzLib and the wire\n\n")
cat("Generated by `scripts/name-parity.R`. Do not edit by hand; re-run it.\n\n")
cat("pyMzLib is the parent binding. Names flow pyMzLib -> mzLibRust -> mzLibR.\n")
cat("R has no per-module namespace and no methods, so `pride.list_files` becomes\n")
cat("`pride_list_files()` and `Digest.truncated` becomes `digest_truncated()`. That flattening\n")
cat("is the only renaming mzLibR is entitled to do.\n\n")

cat("## Functions and parameters\n\n")
cat("| mzLibR | pyMzLib | parameters |\n|---|---|---|\n")

problems <- character(0)

for (name in exported) {
  ours <- r_parameters(name)

  if (name %in% names(additions)) {
    cat("| `", name, "()` | *none* | ", additions[[name]], " |\n", sep = "")
    next
  }

  mapped <- mapping[[name]]
  if (is.null(mapped)) {
    problems <- c(problems, paste0(name, ": exported by mzLibR but not mapped to a parent"))
    cat("| `", name, "()` | **UNMAPPED** | ",
      paste0("`", ours, "`", collapse = ", "), " |\n",
      sep = ""
    )
    next
  }

  theirs <- parent[[mapped]]
  if (is.null(theirs)) {
    problems <- c(problems, paste0(name, ": mapped to ", mapped, ", which was not found"))
    cat("| `", name, "()` | `", mapped, "` **NOT FOUND** | |\n", sep = "")
    next
  }

  # The receiver differs by necessity: a Python method takes `self`, an R function takes the
  # object as its first argument. That first argument is excluded from the comparison.
  ours_compared <- if (grepl("^[A-Z]", mapped)) ours[-1L] else ours
  only_ours <- setdiff(ours_compared, theirs)
  only_theirs <- setdiff(theirs, ours_compared)

  if (length(only_ours) == 0L && length(only_theirs) == 0L) {
    cat("| `", name, "()` | `", mapped, "` | ", length(theirs), " identical |\n", sep = "")
  } else {
    detail <- character(0)
    if (length(only_ours) > 0L) {
      detail <- c(detail, paste0("only mzLibR: ", paste(only_ours, collapse = ", ")))
    }
    if (length(only_theirs) > 0L) {
      detail <- c(detail, paste0("only pyMzLib: ", paste(only_theirs, collapse = ", ")))
    }
    problems <- c(problems, paste0(name, ": ", paste(detail, collapse = "; ")))
    cat("| `", name, "()` | `", mapped, "` | **", paste(detail, collapse = "; "), "** |\n", sep = "")
  }
}

cat("\n### pyMzLib callables with no mzLibR function\n\n")
cat("| pyMzLib | why not |\n|---|---|\n")
for (name in names(omissions)) {
  cat("| `", name, "` | ", omissions[[name]], " |\n", sep = "")
}

# ---------------------------------------------------------------- report: columns

cat("\n## data.frame columns against the wire\n\n")
cat("The wire's field names are the snake_case of mzLib's own. A column marked *derived* is not\n")
cat("on the wire; every other column must match a wire field exactly.\n\n")

read_fixture <- function(name) {
  path <- file.path(fixtures, name)
  if (!file.exists(path)) {
    return(NULL)
  }
  asNamespace("mzLibR")$json_parse(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

report_columns <- function(label, columns, wire_fields, derived) {
  cat("\n### `", label, "`\n\n", sep = "")
  cat("| column | on the wire |\n|---|---|\n")
  for (column in columns) {
    verdict <- if (column %in% wire_fields) {
      "yes"
    } else if (column %in% names(derived)) {
      paste0("*derived* — ", derived[[column]])
    } else {
      problems <<- c(problems, paste0(label, "$", column, ": not a wire field and not declared derived"))
      "**UNDECLARED**"
    }
    cat("| `", column, "` | ", verdict, " |\n", sep = "")
  }
}

pride <- read_fixture("pride_PXD000001_files.json")
if (!is.null(pride)) {
  files <- asNamespace("mzLibR")$pride_parse_manifest(pride, "PXD000001")
  report_columns("pride_list_files()", names(files), names(pride$files[[1]]), c(
    size_mb = "file_size_bytes / 1e6 (pyMzLib PrideFile.size_mb)",
    extension = "from file_name (pyMzLib PrideFile.extension)",
    downloadable = "https_url is present (pyMzLib PrideFile.downloadable)",
    project_accession = "stamped on by pride_list_files so a selection knows its project"
  ))
}

peptidoform <- read_fixture("peptidoform_P02768_small.json")
if (!is.null(peptidoform)) {
  digest <- asNamespace("mzLibR")$peptidoform_parse(peptidoform)
  report_columns("peptidoform_fragments()$peptides", names(digest$peptides),
    names(peptidoform$peptides[[1]]),
    c(peptide_index = "row key joining peptides, fragments and modifications")
  )
  report_columns("peptidoform_fragments()$fragments", names(digest$fragments),
    names(peptidoform$peptides[[1]]$fragments[[1]]),
    c(peptide_index = "row key")
  )
  report_columns("peptidoform_fragments()$modifications", names(digest$modifications),
    names(peptidoform$peptides[[2]]$modifications[[1]]),
    c(peptide_index = "row key")
  )
}

quant <- read_fixture("flashlfq_small.json")
if (!is.null(quant)) {
  results <- asNamespace("mzLibR")$flashlfq_parse(quant)
  report_columns("flashlfq_quantify()$spectra_files", names(results$spectra_files),
    names(quant$spectra_files[[1]]), character(0)
  )
  report_columns("flashlfq_quantify()$peptides", names(results$peptides),
    names(quant$peptides[[1]]),
    c(
      file_name = "key of the `intensities` map, unnested to long form",
      intensity = "value of the `intensities` map",
      detection_type = "value of the `detection_types` map"
    )
  )
  report_columns("flashlfq_quantify()$proteins", names(results$proteins),
    names(quant$proteins[[1]]),
    c(
      file_name = "key of the `intensities` map, unnested to long form",
      intensity = "value of the `intensities` map"
    )
  )
  report_columns("flashlfq_quantify()$peaks", names(results$peaks),
    names(quant$peaks[[1]]), character(0)
  )
}

# ---------------------------------------------------------------- verdict

cat("\n## Anything needing a reason\n\n")
if (length(problems) == 0L) {
  cat("None. Every parameter matches its parent, and every column is either a wire field or a\n")
  cat("declared derivation.\n")
} else {
  for (line in problems) {
    cat("- ", line, "\n", sep = "")
  }
  quit(status = 1L)
}
