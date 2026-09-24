# Generate docs/name-parity.md: mzLibR's names against pyMzLib's, its columns against the wire, and
# its help pages against the bridge's verb specs.
#
# pyMzLib is the parent binding; mzLibRust derives from it and mzLibR from both. Names therefore
# flow downward, and any place mzLibR spells something differently is either a defect or a
# decision that has to be written down. This script finds them mechanically so the answer does
# not depend on anyone remembering.
#
# Three halves, because there are three ways to get it wrong:
#
#   1. **Functions and their parameters**, checked against pyMzLib's module functions *and* its
#      class methods - an R package has no per-module namespace and no methods, so
#      `Digest.truncated` becomes `digest_truncated()`, and that flattening is the only renaming
#      mzLibR is entitled to do. Every pyMzLib callable must be mapped, or omitted with a reason.
#   2. **data.frame columns**, checked against the field names on the wire, which are themselves
#      the snake_case of mzLib's own names. A column is what a user types.
#   3. **The verb specs** (docs/specs/, vendored from the bridge). For every verb mzLibR projects:
#      each spec parameter is an argument (or a declared deviation), each `@param` for a parameter
#      with a unit names that unit, each returned field with a unit is named in `@return` with its
#      unit, and each field the spec says the R object carries really is in it - checked by
#      parsing the verb's recorded fixture. A deviation that names nothing in its spec fails too.
#
# Run from the package root, with mzLibR installed and man/ current:
#   Rscript scripts/name-parity.R <pyMzLib>/pkg/python/src/pymzlib > docs/name-parity.md
#
# with pyMzLib checked out at the ref named in scripts/pymzlib-ref - the release whose bridge this
# package pins, or, while a port runs ahead of that release, the pyMzLib commit it ports.
#
# CI checks pyMzLib out at the ref in scripts/pymzlib-ref and runs this with `git diff --exit-code`,
# so the page cannot go stale and a parity break fails the pull request that causes it.

args <- commandArgs(trailingOnly = TRUE)
python_src <- if (length(args) >= 1L) {
  args[1L]
} else {
  "E:/CodeReview/pyMzLib/code/pyMzLib/pkg/python/src/pymzlib"
}
# What the page says it compared against: the pyMzLib ref CI checks out (scripts/pymzlib-ref), so a
# local run and CI write the same page.
python_label <- if (length(args) >= 2L) {
  args[2L]
} else if (file.exists(file.path("scripts", "pymzlib-ref"))) {
  paste0("`", trimws(readLines(file.path("scripts", "pymzlib-ref"), warn = FALSE)[1L]), "`")
} else {
  "a local checkout"
}
fixtures <- file.path("tests", "fixtures")
if (!dir.exists(python_src)) {
  stop("no pyMzLib source at ", python_src, " - pass <pyMzLib>/pkg/python/src/pymzlib")
}

modules <- c("pride", "peptidoform", "flashlfq", "readers", "sdrf", "proteins")

suppressMessages(library(mzLibR))
mz <- asNamespace("mzLibR")
source(file.path("scripts", "spec-facts.R"))

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
  parameters[nzchar(parameters) & !parameters %in% c("self", "cls", "/")]
}

# Every public callable in a module, keyed `module.function` or `Class.method`. Private classes
# and methods (a leading underscore) are the module's business, not its API.
python_signatures <- function(path, module) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  out <- list()
  current_class <- NA_character_

  for (index in seq_along(lines)) {
    line <- lines[index]
    if (grepl("^class [A-Za-z_]", line)) {
      current_class <- sub("^class ([A-Za-z_][A-Za-z0-9_]*).*$", "\\1", line)
      next
    }
    if (grepl("^def [a-z_][a-z0-9_]*\\(", line)) {
      current_class <- NA_character_
      name <- sub("^def ([a-z_][a-z0-9_]*)\\(.*$", "\\1", line)
      if (startsWith(name, "_")) next
      out[[paste0(module, ".", name)]] <- python_parameters(lines, index)
      next
    }
    if (grepl("^    def [a-z_][a-z0-9_]*\\(", line) && !is.na(current_class)) {
      name <- sub("^    def ([a-z_][a-z0-9_]*)\\(.*$", "\\1", line)
      if (startsWith(name, "_") || startsWith(current_class, "_")) next
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
  pride_list_ftp_files = "pride.list_ftp_files",
  pride_download = "pride.download",
  pride_download_files = "pride.download_files",
  pride_total_size_bytes = "pride.total_size_bytes",
  pride_approximate_total_size_bytes = "pride.approximate_total_size_bytes",
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
  flashlfq_mbr_rescued_peptide_count = "FlashLfqResults.mbr_rescued_peptide_count",
  readers_formats = "readers.formats",
  readers_identify = "readers.identify",
  readers_read_results = "readers.read_results",
  readers_read_records = "readers.read_records",
  readers_read_features = "readers.read_features",
  readers_read_matches = "readers.read_matches",
  readers_read_spectra = "readers.read_spectra",
  sdrf_read = "sdrf.read",
  sdrf_pool = "sdrf.pool",
  sdrf_index_of = "SdrfDocument.index_of",
  sdrf_indexes_of = "SdrfDocument.indexes_of",
  sdrf_value = "SdrfDocument.value",
  sdrf_all = "SdrfDocument.all",
  sdrf_records = "SdrfDocument.records",
  sdrf_has_repeated_columns = "SdrfDocument.has_repeated_columns",
  sdrf_ragged_row_count = "SdrfDocument.ragged_row_count",
  sdrf_source_documents = "PooledSdrf.source_documents"
)

# Parameters that differ by necessity rather than by accident, with the reason. Keyed by the R
# function; each entry names the parameters on each side that the comparison should not count.
parameter_reasons <- list(
  sdrf_pool = list(
    r = character(0), py = character(0),
    why = "R names the documents vector's elements with their labels; pyMzLib takes a {path: label} map. Same parameter, `documents`, both sides"
  )
)

# Deliberately absent from pyMzLib, with the reason. Anything here is a considered addition, not
# an oversight.
additions <- c(
  digest_distinct_base_sequences = "from mzLibRust (Digest::distinct_base_sequences); pyMzLib has no equivalent",
  digest_fragments_by_series = "from mzLibRust (Digest::fragments_by_series); pyMzLib has only fragment_count",
  pride_locations = "R-only: unnests the `locations` list column, which is pyMzLib's PrideFile.locations field",
  readers_retention_time_in_minutes = paste(
    "one function with a `column` argument where pyMzLib has one property per column",
    "(ResultRecords.retention_time_in_minutes, FeatureRecords.retention_time_start_in_minutes and",
    "_end_in_minutes), because R dispatches on the object's class and Python does not"
  ),
  mzlibr_bridge_path = "transport; pyMzLib's equivalent is private (_bridge)",
  mzlibr_bridge_version = "transport; pyMzLib's equivalent is pymzlib.bridge_version, outside the modules compared here",
  # Says why THIS package needs the function, not what the other bindings do instead. The earlier
  # wording ("Rust downloads it from build.rs") was a claim about another repository, which nothing
  # here can test and which goes false when that repository changes with nothing changing here -
  # the mistake mzLibRust #16 was opened to stop making.
  mzlibr_install_bridge = paste(
    "no pyMzLib counterpart: a wheel carries the payload, so Python never has to fetch one.",
    "CRAN forbids downloading at install time and writing outside tempdir() without consent,",
    "so here the download has to be a function the user calls."
  )
)

# pyMzLib callables with no mzLibR counterpart, with the reason. Every pyMzLib callable must be
# mapped above or listed here; one that is neither fails the check, which is how a new pyMzLib
# function gets noticed on the day it lands.
omissions <- c(
  "PrideFile.size_mb" = "a column, `size_mb`, not a function",
  "PrideFile.extension" = "a column, `extension`, not a function",
  "PrideFile.downloadable" = "a column, `downloadable`, not a function",
  "PrideFile.as_dict" = "meaningless in R: the data.frame is already the record",
  "PrideFtpFile.approximate_size_mb" = "a column, `approximate_size_mb`, not a function",
  "PrideFtpFile.extension" = "a column, `extension`, not a function",
  "PrideFtpFile.as_dict" = "meaningless in R: the data.frame is already the record",
  "PrideProjectSearchResult.matched_fields" = "parity debt: `pride search` arrives with the mzLib 1.0.592 port",
  "PrideProjectSearchResult.as_dict" = "parity debt: `pride search` arrives with the mzLib 1.0.592 port",
  "pride.search" = "parity debt: `pride search` arrives with the mzLib 1.0.592 port",
  "flashlfq.median_polish" = "parity debt: `quant median-polish` arrives with the mzLib 1.0.592 port",
  "Peptide.is_modified" = "a column comparison, `modification_count > 0`",
  "Peptide.intensity" = "a row of the long `peptides` frame",
  "Peptide.detection_type" = "a column of the long `peptides` frame",
  "ProteinGroup.intensity" = "a row of the long `proteins` frame",
  "Peak.is_mbr" = "a column comparison, `detection_type == \"MBR\"`",
  "Digest.fragment_count" = "`nrow(digest$fragments)`; not promoted to a function because a bare total folds in the spurious ETD y series — see ?digest_fragments_by_series",
  "Format.is_quantifiable" = "a column of readers_formats(), `is_quantifiable`",
  "FileInfo.is_quantifiable" = "an element of readers_identify(), `is_quantifiable`",
  "ResultRecords.retention_time_in_minutes" = "readers_retention_time_in_minutes()",
  "FeatureRecords.retention_time_start_in_minutes" = "readers_retention_time_in_minutes(column = \"retention_time_start\")",
  "FeatureRecords.retention_time_end_in_minutes" = "readers_retention_time_in_minutes(column = \"retention_time_end\")",
  "ScanRecords.total_ion_current" = "a column of `records`, `total_ion_current`"
)

# ---------------------------------------------------------------- the R side

exported <- sort(getNamespaceExports("mzLibR"))
exported <- exported[!grepl("^print\\.", exported)]
functions <- exported[vapply(exported, function(n) is.function(get(n, envir = mz)), logical(1L))]

r_parameters <- function(name) {
  formal_names <- names(formals(get(name, envir = mz)))
  formal_names[formal_names != "..."]
}

# ---------------------------------------------------------------- report: functions

cat("# Name parity with pyMzLib, the wire and the verb specs\n\n")
cat("Generated by `scripts/name-parity.R`. Do not edit by hand; re-run it.\n\n")
cat("Compared against pyMzLib at ", python_label, ", and the verb specs vendored in `docs/specs/` ",
  "from bridge commit `", substr(spec_source_commit(), 1L, 12L), "`.\n\n", sep = "")
cat("pyMzLib is the parent binding. Names flow pyMzLib -> mzLibRust -> mzLibR.\n")
cat("R has no per-module namespace and no methods, so `pride.list_files` becomes\n")
cat("`pride_list_files()` and `Digest.truncated` becomes `digest_truncated()`. That flattening\n")
cat("is the only renaming mzLibR is entitled to do.\n\n")

cat("## Functions and parameters\n\n")
cat("| mzLibR | pyMzLib | parameters |\n|---|---|---|\n")

problems <- character(0)

for (name in functions) {
  ours <- r_parameters(name)

  if (name %in% names(additions)) {
    cat("| `", name, "()` | *none* | ", additions[[name]], " |\n", sep = "")
    next
  }

  mapped <- if (name %in% names(mapping)) mapping[[name]] else NULL
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
  if (!name %in% names(parent)) {
    problems <- c(problems, paste0(name, ": listed as omitted, but pyMzLib has no such callable"))
  }
  cat("| `", name, "` | ", omissions[[name]], " |\n", sep = "")
}
unaccounted <- setdiff(names(parent), c(unname(mapping), names(omissions)))
for (name in unaccounted) {
  problems <- c(problems, paste0(name, ": pyMzLib callable with no mzLibR function and no reason given"))
  cat("| `", name, "` | **UNACCOUNTED FOR** |\n", sep = "")
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
  text <- rawToChar(readBin(path, what = "raw", n = file.info(path)$size))
  Encoding(text) <- "UTF-8"
  data <- mz$json_parse(text)
  if (is.list(data) && !is.null(data[["ok"]]) && "data" %in% names(data)) data[["data"]] else data
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
  files <- mz$pride_parse_manifest(pride, "PXD000001")
  report_columns("pride_list_files()", names(files), names(pride$files[[1]]), c(
    size_mb = "file_size_bytes / 1e6 (pyMzLib PrideFile.size_mb)",
    extension = "from file_name (pyMzLib PrideFile.extension)",
    downloadable = "https_url is present (pyMzLib PrideFile.downloadable)",
    project_accession = "stamped on by pride_list_files so a selection knows its project"
  ))
}

ftp <- read_fixture("pride_ftp_PXD000001.json")
if (!is.null(ftp)) {
  listing <- mz$pride_parse_ftp_files(ftp, "PXD000001")
  report_columns("pride_list_ftp_files()", names(listing), names(ftp$files[[1]]), c(
    approximate_size_mb = "approximate_size_bytes / 1e6 (pyMzLib PrideFtpFile.approximate_size_mb)",
    extension = "from file_name (pyMzLib PrideFtpFile.extension)"
  ))
}

peptidoform <- read_fixture("peptidoform_P02768_small.json")
if (!is.null(peptidoform)) {
  digest <- mz$peptidoform_parse(peptidoform)
  report_columns("peptidoform_fragments()$peptides", names(digest$peptides),
    names(peptidoform$peptides[[1]]),
    c(peptide_index = "row key joining peptides, fragments and modifications")
  )
  report_columns("peptidoform_fragments()$fragments", names(digest$fragments),
    names(peptidoform$peptides[[1]]$fragments[[1]]),
    c(peptide_index = "row key")
  )
  with_modifications <- Filter(function(p) length(p$modifications) > 0L, peptidoform$peptides)
  if (length(with_modifications) > 0L) {
    report_columns("peptidoform_fragments()$modifications", names(digest$modifications),
      names(with_modifications[[1]]$modifications[[1]]),
      c(peptide_index = "row key")
    )
  }
}

quant <- read_fixture("flashlfq_small.json")
if (!is.null(quant)) {
  results <- mz$flashlfq_parse(quant)
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

formats_payload <- read_fixture("readers_formats.json")
if (!is.null(formats_payload)) {
  formats <- mz$readers_parse_formats(formats_payload)
  report_columns("readers_formats()", names(formats),
    names(formats_payload$formats[[1]]),
    c(is_quantifiable = "\"quantifiable\" is among views (pyMzLib Format.is_quantifiable)")
  )
}

# ---------------------------------------------------------------- report: the verb specs

cat("\n## Help pages against the verb specs\n\n")
cat("For each wire verb mzLibR projects: every spec parameter is an argument or a declared\n")
cat("deviation, and its `@param` names its unit; every returned field that has a unit is named in\n")
cat("`@return` with its unit; and every field the spec says the R object carries is in the object\n")
cat("parsed from the verb's recorded fixture. Deviations are declared, with reasons, in\n")
cat("`scripts/spec-facts.R`.\n\n")

rd_cache <- new.env()
rd_for <- function(name) {
  if (is.null(rd_cache[[name]])) {
    rd_cache[[name]] <- tools::parse_Rd(file.path("man", paste0(name, ".Rd")))
  }
  rd_cache[[name]]
}

rd_flat <- function(x) {
  gsub("\\s+", " ", paste(unlist(x), collapse = ""))
}

rd_argument_text <- function(name, argument) {
  rd <- rd_for(name)
  for (section in rd) {
    if (identical(attr(section, "Rd_tag"), "\\arguments")) {
      for (item in section) {
        if (identical(attr(item, "Rd_tag"), "\\item") && identical(rd_flat(item[[1L]]), argument)) {
          return(rd_flat(item[[2L]]))
        }
      }
    }
  }
  NA_character_
}

rd_value_text <- function(name) {
  rd <- rd_for(name)
  for (section in rd) {
    if (identical(attr(section, "Rd_tag"), "\\value")) {
      return(rd_flat(section))
    }
  }
  ""
}

# Is `field` mentioned in `text` with its unit close by, before or after?
mentioned_with_unit <- function(text, field, unit) {
  hits <- gregexpr(field, text, fixed = TRUE)[[1L]]
  if (hits[1L] == -1L) {
    return(NA)
  }
  any(vapply(hits, function(at) {
    window <- substr(text, max(1L, at - 120L), at + nchar(field) + 160L)
    mentions_unit(window, unit)
  }, logical(1L)))
}

# How to parse each verb's recorded fixture into the R object, for the fields-present check.
PARSE_FIXTURE <- list(
  "version" = list("bridge_version.json", function(d) {
    text <- rawToChar(readBin(file.path(fixtures, "bridge_version.json"), "raw", 1e6))
    old <- options(mzlibr.bridge = file.path(fixtures, "bridge_version.json"))
    on.exit(options(old))
    mzlibr_bridge_version(runner = function(exe, args, stdin = NULL, timeout = NULL) {
      list(stdout = paste0("{\"ok\":true,\"data\":", text, "}"), stderr = "", status = 0L, timed_out = FALSE)
    })
  }),
  "readers formats" = list("readers_formats.json", function(d) mz$readers_parse_formats(d)),
  "readers identify" = list("readers_identify_mzid.json", function(d) mz$readers_parse_file_info(d)),
  "readers read-results" = list("readers_results_fragger.json", function(d) mz$readers_parse_records(d)),
  "readers read-records" = list("readers_records_toppic.json", function(d) mz$readers_parse_native_records(d)),
  "readers read-features" = list("readers_features_topfd.json", function(d) mz$readers_parse_feature_records(d)),
  "readers read-matches" = list("readers_matches_casanovo.json", function(d) mz$readers_parse_match_records(d)),
  "readers read-spectra" = list("readers_spectra_mzml.json", function(d) mz$readers_parse_scan_records(d)),
  "sdrf read" = list("sdrf_read_PXD000070.json", function(d) mz$sdrf_parse_document(d)),
  "sdrf pool" = list("sdrf_pool_two.json", function(d) mz$sdrf_parse_pooled(d)),
  "pride files" = list("pride_PXD000001_files.json", function(d) mz$pride_parse_manifest(d, "PXD000001")),
  "pride ftp-files" = list("pride_ftp_PXD000001.json", function(d) mz$pride_parse_ftp_files(d, "PXD000001")),
  "peptidoform fragments" = list("peptidoform_P02768_small.json", function(d) mz$peptidoform_parse(d)),
  "quant flashlfq" = list("flashlfq_small.json", function(d) mz$flashlfq_parse(d))
)

# Every name reachable in the R object: its elements, the columns of any data.frame among them,
# and `a$b` paths one level into a list element (the peptidoform census).
object_names <- function(object) {
  out <- names(object)
  if (is.null(out)) out <- character(0)
  for (n in names(object)) {
    element <- object[[n]]
    if (is.data.frame(element)) {
      out <- c(out, names(element))
    }
    if (is.list(element) && !is.null(names(element))) {
      out <- c(out, paste0(n, "$", names(element)))
    }
  }
  unique(out)
}

# The functions projecting each spec, with the page mode each is rendered in.
spec_pages <- function(spec) {
  pages <- list()
  main <- spec_r_name(spec)
  if (!is.na(main) && main %in% functions) pages[[main]] <- ""
  bulk <- spec_r_bulk_name(spec)
  if (!is.na(bulk) && bulk %in% functions) pages[[bulk]] <- "bulk"
  selection <- spec_text(spec$bindings$r$selection, NA_character_)
  if (!is.na(selection) && selection %in% functions) pages[[selection]] <- "selection"
  pages
}

all_spec_keys <- function(spec) {
  params <- c(spec$params, spec$result$bulk$params)
  fields <- c(spec$result$envelope_fields, spec$result$bulk$envelope_fields,
    if (is.list(spec$result$columns)) spec$result$columns else list())
  keys <- c(
    paste0("param.", vapply(params, function(p) p$name, character(1L))),
    paste0("field.", vapply(fields, function(f) f$wire, character(1L)))
  )
  for (table in names(spec_tables(spec))) {
    keys <- c(keys, paste0("table.", table),
      paste0("field.", vapply(spec_tables(spec)[[table]], field_key, character(1L))))
  }
  unique(keys)
}

specs <- load_specs()
spec_rows <- character(0)
spelling <- character(0)

for (key in names(R_DEVIATIONS)) {
  verb <- sub(" \\((bulk|selection)\\)$", "", key)
  if (!verb %in% names(specs)) {
    problems <- c(problems, paste0("R_DEVIATIONS names '", key, "', which has no spec"))
    next
  }
  for (entry in setdiff(names(R_DEVIATIONS[[key]]), all_spec_keys(specs[[verb]]))) {
    problems <- c(problems, paste0("R_DEVIATIONS[['", key, "']] names ", entry, ", which its spec does not have"))
  }
}
for (verb in c(names(R_NAMES), names(R_TABLE))) {
  if (!verb %in% names(specs)) {
    problems <- c(problems, paste0("scripts/spec-facts.R names verb '", verb, "', which has no spec"))
  }
}

for (spec in specs) {
  pages <- spec_pages(spec)
  since_r <- spec_text(spec$since$mzlibr, "")
  declared <- spec_text(spec$bindings$r$name, "")
  if (spec$verb %in% names(R_NAMES)) {
    spelling <- c(spelling, paste0(
      "- `", spec$verb, "`: the spec says `", declared, "`, mzLibR has `", R_NAMES[[spec$verb]], "`"
    ))
  }
  if (length(pages) == 0L) {
    spec_rows <- c(spec_rows, paste0(
      "| `", spec$verb, "` | `", declared, "` | ", if (nzchar(since_r)) since_r else "—",
      " | *not projected yet* | | |"
    ))
    next
  }

  checked_params <- 0L
  checked_fields <- 0L
  verb_problems <- character(0)

  for (fn in names(pages)) {
    page_spec <- spec
    page_spec$.page <- pages[[fn]]
    bulk <- identical(pages[[fn]], "bulk")
    formals_here <- r_parameters(fn)

    for (param in spec_params(page_spec, bulk)) {
      r_name <- spec_r_param(page_spec, param$name)
      if (is.na(r_name)) next
      checked_params <- checked_params + 1L
      if (!r_name %in% formals_here) {
        verb_problems <- c(verb_problems, paste0(fn, "(): spec parameter --", param$name, " is not an argument (expected `", r_name, "`)"))
        next
      }
      if (!spec_null(param$unit)) {
        text <- rd_argument_text(fn, r_name)
        if (is.na(text) || !mentions_unit(text, spec_text(param$unit))) {
          verb_problems <- c(verb_problems, paste0(
            fn, "(): @param ", r_name, " does not state its unit (", spec_text(param$unit), ")"
          ))
        }
      }
    }

    value <- rd_value_text(fn)
    fields <- c(spec_fields(page_spec, "envelope_fields", bulk), spec_fields(page_spec, "columns", bulk),
      unlist(unname(spec_tables(page_spec)), recursive = FALSE))
    for (field in fields) {
      deviation <- spec_deviation(page_spec, "field", field_key(field))
      if (!is.null(deviation) && (is.na(deviation$r) || isTRUE(deviation$pending))) next
      if (spec_null(field$unit)) next
      r_name <- sub("^.*\\$", "", spec_r_field(page_spec, field_key(field)))
      checked_fields <- checked_fields + 1L
      found <- mentioned_with_unit(value, paste0(r_name), spec_text(field$unit))
      if (is.na(found)) {
        verb_problems <- c(verb_problems, paste0(fn, "(): @return does not mention `", r_name, "` (", spec_text(field$unit), ")"))
      } else if (!found) {
        verb_problems <- c(verb_problems, paste0(fn, "(): @return mentions `", r_name, "` without its unit (", spec_text(field$unit), ")"))
      }
    }
  }

  # The fields the spec says the R object carries, against the object itself.
  check <- PARSE_FIXTURE[[spec$verb]]
  fixture_note <- "no fixture"
  if (!is.null(check) && file.exists(file.path(fixtures, check[[1L]]))) {
    object <- check[[2L]](read_fixture(check[[1L]]))
    present <- object_names(object)
    fields <- c(spec_fields(spec, "envelope_fields"), spec_fields(spec, "columns"),
      unlist(unname(spec_tables(spec)), recursive = FALSE))
    missing <- character(0)
    for (field in fields) {
      deviation <- spec_deviation(spec, "field", field_key(field))
      if (!is.null(deviation) && (is.na(deviation$r) || isTRUE(deviation$pending))) next
      r_name <- spec_r_field(spec, field_key(field))
      if (identical(r_name, "columns") || !spec_null(field$present_when)) next
      if (!r_name %in% present) missing <- c(missing, r_name)
    }
    for (m in missing) {
      verb_problems <- c(verb_problems, paste0(
        spec$verb, ": the spec's field `", m, "` is not in the R object parsed from ", check[[1L]],
        " - project it, or declare it in R_DEVIATIONS"
      ))
    }
    fixture_note <- paste0("`", check[[1L]], "`")
  }

  problems <- c(problems, verb_problems)
  spec_rows <- c(spec_rows, paste0(
    "| `", spec$verb, "` | ", paste0("`", names(pages), "()`", collapse = ", "), " | ",
    if (nzchar(since_r)) since_r else "—", " | ", checked_params, " params, ", checked_fields,
    " fields with units | ", fixture_note, " | ",
    if (length(verb_problems) == 0L) "ok" else paste0("**", length(verb_problems), " problem(s)**"), " |"
  ))
}

cat("| verb | mzLibR | spec `since.mzlibr` | checked | fields checked against | verdict |\n")
cat("|---|---|---|---|---|---|\n")
cat(paste0(spec_rows, "\n"), sep = "")

cat("\n### Where mzLibR's name differs from the spec's `bindings.r`\n\n")
if (length(spelling) == 0L) {
  cat("None.\n")
} else {
  cat(paste0(spelling, "\n"), sep = "")
}

# ---------------------------------------------------------------- verdict

cat("\n## Anything needing a reason\n\n")
if (length(problems) == 0L) {
  cat("None. Every parameter matches its parent, every column is either a wire field or a\n")
  cat("declared derivation, and every help page states the units its spec gives.\n")
} else {
  for (line in problems) {
    cat("- ", line, "\n", sep = "")
  }
  quit(status = 1L)
}
