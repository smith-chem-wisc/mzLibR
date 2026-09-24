# Shared reading of the vendored verb specs in docs/specs/.
#
# Sourced by scripts/build-man.R (which renders each spec's fact tables into the man page of the
# function that projects the verb), scripts/name-parity.R (which lints the hand-written @param and
# @return text against them) and scripts/stage-replay.R (which picks the recordings the examples
# replay), so the three agree on what a spec says and how its wire names become R names.
#
# Base R only. The specs are read from the .json copy scripts/sync-specs.R writes beside each
# .yaml, with this package's own JSON reader, so none of these scripts needs the `yaml` package.
# They need mzLibR INSTALLED, like build-man.R always has.

SPEC_DIR <- file.path("docs", "specs")

# How the spec's language-neutral error kinds surface in R: a condition class, dispatched on by
# tryCatch(). The bridge reports `usage` and `ServiceUnavailable` by name; anything else is a
# correctness failure carrying mzLib's own exception type in `error_type`. See ?mzlib_error.
R_CONDITION <- c(
  usage = "mzlib_usage_error",
  service_unavailable = "mzlib_service_unavailable",
  correctness = "mzlib_bridge_error"
)

# Where mzLibR's function name differs from the spec's `bindings.r.name`. Every entry is reported
# back to the bridge so the spec can be corrected; until then the lint uses this one.
R_NAMES <- c(
  # Transport functions carry the package prefix, like mzlibr_bridge_path() and
  # mzlibr_install_bridge(), so a user can find all three together.
  "version" = "mzlibr_bridge_version"
)

# Where the R object keeps the verb's table, when it has one. The wire's `columns` map becomes a
# data.frame under this name.
R_TABLE <- c(
  "readers formats" = "(the returned data.frame itself)",
  "readers read-results" = "records",
  "readers read-records" = "records",
  "readers read-features" = "records",
  "readers read-matches" = "records",
  "readers read-spectra" = "records",
  "pride files" = "(the returned data.frame itself)",
  "pride ftp-files" = "(the returned data.frame itself)",
  "peptidoform fragments" = "peptides"
)

# A field of the verb that is still on the wire at the spec's version but that this package does
# not return yet. Listed as such on the man page, so a page never claims a field it lacks.
pending <- function(why = "arrives with the mzLib 1.0.592 port") {
  list(r = NA, pending = TRUE, why = why)
}

not_here <- function(why) {
  list(r = NA, why = why)
}

as_r <- function(r, why) {
  list(r = r, why = why)
}

ONE_PATH_ERROR <- not_here("always null for one path: a file that cannot be read raises instead")
RECORDS <- as_r("records", "the table is a data.frame, so it is `records`")

# Every place mzLibR deliberately projects a spec param or result field under another name, or not
# as a named element at all. **Nothing else may be skipped by the lint**, and the lint fails on an
# entry here that names no param or field of the spec, so this list cannot go stale silently.
#
# Key: the spec's `verb`, or `"<verb> (bulk)"` / `"<verb> (selection)"` for entries that apply to
# one page only (they are consulted first). Value: a list of `"param.<wire name>"`,
# `"field.<wire name>"`, `"field.<table>.<wire name>"` (a field of `result.tables`) and
# `"table.<name>"` entries, each `list(r = <the R name, or NA when it is not a named element>,
# why = <the reason>)`. An R name may be a path such as `census$sites`.
R_DEVIATIONS <- list(
  "version" = list(
    "field.verbs" = pending()
  ),
  "readers formats" = list(
    "field.format_count" = not_here("the result is the formats data.frame; this is nrow() of it"),
    "field.formats" = not_here("the result is this list itself, as a data.frame, one row per format")
  ),
  "readers identify" = list(
    "field.error" = ONE_PATH_ERROR
  ),
  "readers read-results" = list(
    "field.columns" = RECORDS,
    "field.reader" = pending(),
    "field.absent_fields" = pending(),
    "field.failed_fields" = pending(),
    "field.excluded_fields" = pending(),
    "field.error" = ONE_PATH_ERROR
  ),
  "readers read-records" = list(
    "field.columns" = RECORDS,
    "field.retention_time_unit" = pending(),
    "field.caveats" = pending(),
    "field.absent_fields" = pending(),
    "field.skipped_count" = pending(),
    "field.skipped" = pending(),
    "field.rows_not_read" = pending(),
    "field.error" = ONE_PATH_ERROR
  ),
  "readers read-features" = list(
    "field.columns" = RECORDS,
    "field.reader" = pending(),
    "field.rows_not_read" = pending(),
    "field.absent_fields" = pending(),
    "field.failed_fields" = pending(),
    "field.excluded_fields" = pending(),
    "field.error" = ONE_PATH_ERROR
  ),
  "readers read-matches" = list(
    "param.scores" = pending(),
    "field.columns" = RECORDS,
    "field.reader" = pending(),
    "field.skipped_count" = pending(),
    "field.skipped" = pending(),
    "field.rows_not_read" = pending(),
    "field.retention_time_unit" = pending(),
    "field.absent_fields" = pending(),
    "field.failed_fields" = pending(),
    "field.excluded_fields" = pending(),
    "field.scores_included" = pending(),
    "field.row_count" = pending(),
    "field.q_value" = pending(),
    "field.rank" = pending(),
    "field.pass_threshold" = pending(),
    "field.match_index" = pending(),
    "field.score_name" = pending(),
    "field.score_value" = pending(),
    "field.error" = ONE_PATH_ERROR
  ),
  "readers read-spectra" = list(
    "field.columns" = RECORDS,
    "field.source" = pending(),
    "field.rows_not_read" = pending(),
    "field.absent_fields" = pending(),
    "field.failed_fields" = pending(),
    "field.excluded_fields" = pending(),
    "field.error" = ONE_PATH_ERROR
  ),
  # SDRF's column names repeat, so a document keeps its header as `columns` - a character vector
  # that may hold duplicates - and its cells as `rows`, rather than as a data.frame.
  "sdrf read" = list(
    "field.column_names" = as_r("columns", "the header, duplicates kept, is `columns`")
  ),
  "sdrf pool" = list(
    "param.stdin" = as_r("documents", "one stdin line per element of `documents`; its name is the label"),
    "field.column_names" = as_r("columns", "the header, duplicates kept, is `columns`")
  ),
  "pride files" = list(
    "field.accession" = as_r("project_accession", "stamped on every row, so a selection knows its project"),
    "field.file_count" = not_here("the result is the files data.frame; this is nrow() of it"),
    "field.total_size_bytes" = not_here("pride_total_size_bytes() of the result"),
    "field.files" = not_here("the result is this list itself, as a data.frame, one row per file")
  ),
  "pride ftp-files" = list(
    "field.accession" = not_here("an FTP listing is not a download selection, so no row carries it"),
    "field.file_count" = not_here("the result is the files data.frame; this is nrow() of it"),
    "field.approximate_total_size_bytes" = not_here("pride_approximate_total_size_bytes() of the result"),
    "field.files" = not_here("the result is this list itself, as a data.frame, one row per file")
  ),
  "pride download" = list(
    "param.dest" = as_r("destination", "spelled out"),
    "param.ext" = as_r("extensions", "spelled out; a character vector"),
    "param.no-overwrite" = as_r("overwrite", "stated positively: overwrite = FALSE sends --no-overwrite"),
    "param.names-from-stdin" = not_here("pride_download_files() sends a selection, and sets this itself"),
    "param.stdin" = not_here("pride_download_files() sends the selected file names"),
    "field.accession" = not_here("the result is the character vector of paths"),
    "field.destination_directory" = not_here("the result is the character vector of paths"),
    "field.downloaded_count" = not_here("length() of the result"),
    "field.paths" = not_here("the result is this character vector itself")
  ),
  "pride download (selection)" = list(
    "param.accession" = not_here("taken from the rows' project_accession"),
    "param.category" = not_here("a selection is already filtered, by `[`"),
    "param.ext" = not_here("a selection is already filtered, by `[`"),
    "param.names-from-stdin" = not_here("always set: the selection travels on stdin"),
    "param.stdin" = as_r("files", "the file_name of each selected row, one per stdin line")
  ),
  "peptidoform fragments" = list(
    "param.no-modifications" = as_r("modifications", "stated positively: modifications = FALSE sends --no-modifications"),
    "param.max-mods" = as_r("max_modifications", "spelled out"),
    "field.annotated_modification_sites" = as_r("census$sites", "gathered into the census; see census_explain()"),
    "field.annotated_modifications_loaded" = as_r("census$applied", "gathered into the census; see census_explain()"),
    "field.uniprot_annotated_features" = as_r("census$annotated", "gathered into the census; see census_explain()"),
    "field.unresolved_modifications" = as_r("census$unresolved", "gathered into the census; see census_excluded()"),
    "field.uniprot_features_by_type" = as_r("census$by_type", "gathered into the census, as a data.frame"),
    "field.peptide_count" = not_here("nrow() of `peptides`"),
    "field.modifications" = as_r("modifications", "its own long data.frame, joined to `peptides` on peptide_index"),
    "field.fragments" = as_r("fragments", "its own long data.frame, joined to `peptides` on peptide_index")
  ),
  "quant flashlfq" = list(
    "param.stdin" = as_r("spectra", "one stdin line per run, with its design columns"),
    "param.ppm" = as_r("ppm_tolerance", "mzLib's parameter name"),
    "param.isotope-ppm" = as_r("isotope_ppm_tolerance", "mzLib's parameter name"),
    "param.mbr" = as_r("match_between_runs", "mzLib's parameter name"),
    "param.mbr-ppm" = as_r("mbr_ppm_tolerance", "mzLib's parameter name"),
    "param.mbr-q" = as_r("mbr_q_value_threshold", "mzLib's parameter name"),
    "param.shared-peptides" = as_r("use_shared_peptides_for_protein_quant", "mzLib's parameter name"),
    "param.bayesian" = as_r("bayesian_protein_quant", "mzLib's parameter name"),
    "param.use-pep-q" = as_r("use_pep_q_value", "mzLib's parameter name"),
    "param.threads" = as_r("max_threads", "mzLib's parameter name; the default differs, see its argument"),
    "param.out" = as_r("output_directory", "mzLib's parameter name"),
    "field.peptide_count" = not_here("flashlfq_peptide_count() of the result"),
    "field.protein_count" = not_here("flashlfq_protein_count() of the result"),
    "field.peptides.intensities" = as_r("intensity", "unnested: one row per peptide per run, the run in file_name"),
    "field.peptides.detection_types" = as_r("detection_type", "unnested beside intensity"),
    "field.proteins.intensities" = as_r("intensity", "unnested: one row per protein group per sample, the sample in file_name")
  )
)

# ---------------------------------------------------------------- reading

spec_null <- function(x) {
  is.null(x) || (length(x) == 1L && !is.list(x) && is.na(x))
}

spec_text <- function(x, default = "") {
  if (spec_null(x)) default else paste(as.character(unlist(x)), collapse = " ")
}

spec_read_json <- function(path) {
  size <- file.info(path)$size
  text <- rawToChar(readBin(path, what = "raw", n = size))
  Encoding(text) <- "UTF-8"
  asNamespace("mzLibR")$json_parse(text)
}

# Every vendored spec, keyed by its wire verb, each with `.file` set to its file stem.
load_specs <- function(dir = SPEC_DIR) {
  files <- sort(list.files(dir, pattern = "[.]json$", full.names = TRUE))
  specs <- lapply(files, function(path) {
    spec <- spec_read_json(path)
    spec$.file <- sub("[.]json$", "", basename(path))
    spec
  })
  names(specs) <- vapply(specs, function(s) s$verb, character(1L))
  specs
}

spec_by_file <- function(specs, stem) {
  hit <- Filter(function(s) identical(s$.file, stem), specs)
  if (length(hit) != 1L) {
    stop("no vendored spec docs/specs/", stem, ".json - run scripts/sync-specs.R")
  }
  hit[[1L]]
}

spec_source_commit <- function(dir = SPEC_DIR) {
  path <- file.path(dir, "SOURCE")
  if (!file.exists(path)) {
    return("unknown")
  }
  line <- grep("^commit:", readLines(path, warn = FALSE), value = TRUE)
  if (length(line) == 0L) "unknown" else trimws(sub("^commit:", "", line[1L]))
}

# The R function that projects the verb, per the spec (or R_NAMES where mzLibR differs).
spec_r_name <- function(spec) {
  if (spec$verb %in% names(R_NAMES)) {
    return(R_NAMES[[spec$verb]])
  }
  spec_text(spec$bindings$r$name, NA_character_)
}

spec_r_bulk_name <- function(spec) {
  spec_text(spec$bindings$r$bulk, NA_character_)
}

# `ms-order` becomes `ms_order`, unless the deviations say otherwise. NA: not a formal.
spec_r_param <- function(spec, wire) {
  deviation <- spec_deviation(spec, "param", wire)
  if (!is.null(deviation)) deviation$r else gsub("-", "_", wire, fixed = TRUE)
}

spec_r_field <- function(spec, wire) {
  deviation <- spec_deviation(spec, "field", wire)
  # A table field is keyed `<table>.<wire>`; wire names themselves never contain a dot.
  if (!is.null(deviation)) deviation$r else sub("^.*[.]", "", wire)
}

# The deviation for one param, field or table, looking first at the page's own list
# ("<verb> (bulk)" or "<verb> (selection)") when the spec is rendered for such a page.
spec_deviation <- function(spec, kind, wire) {
  key <- paste0(kind, ".", wire)
  page <- spec$.page
  if (!is.null(page) && nzchar(page)) {
    own <- R_DEVIATIONS[[paste0(spec$verb, " (", page, ")")]][[key]]
    if (!is.null(own)) {
      return(own)
    }
  }
  R_DEVIATIONS[[spec$verb]][[key]]
}

spec_params <- function(spec, bulk = FALSE) {
  params <- spec$params
  if (!is.list(params)) {
    params <- list()
  }
  if (!bulk) {
    return(params)
  }
  # The bulk page: the _many form takes a list of paths instead of `path`, drops the window (BULK.md
  # section 1 refuses offset and limit with --paths-stdin), and adds the bulk options.
  kept <- Filter(function(p) !p$name %in% c("path", "paths-stdin", "offset", "limit", "out", "threads", "on-error"), params)
  extra <- spec$result$bulk$params
  if (!is.list(extra)) {
    extra <- Filter(function(p) p$name %in% c("paths-stdin", "threads", "on-error"), params)
  }
  c(kept, extra)
}

spec_fields <- function(spec, section = c("envelope_fields", "columns"), bulk = FALSE) {
  section <- match.arg(section)
  result <- spec$result
  fields <- if (bulk && section == "envelope_fields" && is.list(result$bulk$envelope_fields)) {
    result$bulk$envelope_fields
  } else {
    result[[section]]
  }
  if (!is.list(fields)) {
    return(list())
  }
  # A column present only in the --paths-stdin shape (source_index, source_path) belongs to the
  # bulk page and not the one-path page.
  Filter(function(f) {
    when <- spec_text(f$present_when)
    !(identical(when, "--paths-stdin") && !bulk)
  }, fields)
}

# `result.tables`, each field carrying `.key` = "<table>.<wire>" for R_DEVIATIONS lookups.
spec_tables <- function(spec) {
  tables <- spec$result$tables
  if (!is.list(tables) || is.null(names(tables))) {
    return(list())
  }
  out <- lapply(names(tables), function(name) {
    fields <- tables[[name]]
    if (!is.list(fields)) {
      return(list())
    }
    lapply(fields, function(f) {
      f$.key <- paste0(name, ".", f$wire)
      f
    })
  })
  names(out) <- names(tables)
  out
}

# ---------------------------------------------------------------- units, for the lint

# What counts as mentioning a spec unit in prose. "Rows to skip" mentions `rows`; "skip this many"
# does not, and fails, exactly as pyMzLib's lint does.
unit_patterns <- function(unit) {
  unit <- trimws(unit)
  known <- list(
    "min" = c("minute", "min"),
    "s" = c("second"),
    "ms" = c("millisecond", "ms"),
    "Da" = c("dalton", "Da"),
    "m/z" = c("m/z"),
    "V" = c("volt", "V"),
    "charge" = c("charge"),
    "retention_time_unit" = c("retention_time_unit")
  )
  if (unit %in% names(known)) {
    return(known[[unit]])
  }
  head <- trimws(sub("\\s*\\(.*$", "", unit))
  if (grepl("^intensity", unit)) {
    return(c("intensity", "instrument units"))
  }
  if (grepl("^engine-defined", unit)) {
    return(c("engine", "score_name"))
  }
  if (grepl("^PSMs", unit)) {
    return(c("PSM"))
  }
  # A count: "rows" is mentioned by "row" as much as by "rows".
  stem <- if (nchar(head) > 3L && grepl("s$", head)) sub("s$", "", head) else head
  stem
}

mentions_unit <- function(text, unit) {
  patterns <- unit_patterns(unit)
  any(vapply(patterns, function(p) {
    if (nchar(p) <= 2L) {
      grepl(paste0("(^|[^A-Za-z])", gsub("([/.])", "\\\\\\1", p), "([^A-Za-z]|$)"), text)
    } else {
      grepl(p, text, ignore.case = TRUE, fixed = FALSE)
    }
  }, logical(1L)))
}

# ---------------------------------------------------------------- rendering into Rd

MZLIB_BLOB <- "https://github.com/smith-chem-wisc/mzLib/blob"

# Escape spec text for Rd, turning `code` into \code{} and **bold** into \strong{}. Unlike
# build-man.R's rd_inline(), square brackets are NOT links here: spec prose says `files[i].error`,
# and a link to a topic called `i` would be nonsense.
rd_spec <- function(text) {
  text <- gsub("\\s+", " ", trimws(text))
  text <- gsub("\\", "\\\\", text, fixed = TRUE)
  text <- gsub("%", "\\%", text, fixed = TRUE)
  text <- gsub("{", "\\{", text, fixed = TRUE)
  text <- gsub("}", "\\}", text, fixed = TRUE)
  text <- gsub("\\*\\*([^*]+)\\*\\*", "\\\\strong\\{\\1\\}", text)
  # A code span holding a quote goes in \samp{}: \code{} is parsed as R, where a lone quote opens a
  # string that swallows the closing brace.
  text <- gsub("`([^`]*['\"][^`]*)`", "\\\\samp\\{\\1\\}", text)
  text <- gsub("`([^`]+)`", "\\\\code\\{\\1\\}", text)
  # A bare wire option such as --paths-stdin: in running Rd text `--` renders as a dash.
  text <- gsub("(?<![{\\w-])(--[a-z][a-z-]*)", "\\\\code\\{\\1\\}", text, perl = TRUE)
  text
}

# Whether mzLibR exports a function of this name, so a page links only to what exists.
r_exports <- function(name) {
  !is.na(name) && name %in% getNamespaceExports("mzLibR")
}

# \code{} is parsed as R, where a lone quote opens a string and swallows the closing brace - so
# text with a quote in it ("mzLib's ...") goes in \samp{}, which is verbatim and renders the same.
rd_code <- function(text) {
  text <- gsub("`", "", text, fixed = TRUE)
  macro <- if (grepl("['\"]", text)) "\\samp{" else "\\code{"
  inner <- gsub("\\", "\\\\", text, fixed = TRUE)
  inner <- gsub("%", "\\%", inner, fixed = TRUE)
  inner <- gsub("{", "\\{", inner, fixed = TRUE)
  inner <- gsub("}", "\\}", inner, fixed = TRUE)
  paste0(macro, gsub("\\s+", " ", trimws(inner)), "}")
}

rd_default <- function(param) {
  if (isTRUE(param$required)) {
    return("required")
  }
  default <- param$default
  if (spec_null(default)) {
    return("absent")
  }
  if (is.logical(default)) {
    return(rd_code(if (default) "TRUE" else "FALSE"))
  }
  if (is.list(default)) {
    return(rd_code(paste(unlist(default), collapse = ", ")))
  }
  rd_code(format(default))
}

rd_param_item <- function(spec, param) {
  wire <- param$name
  r_name <- spec_r_param(spec, wire)
  label <- if (is.na(r_name)) {
    paste0("wire \\code{--", rd_spec(wire), "}")
  } else if (identical(r_name, wire)) {
    rd_code(r_name)
  } else {
    paste0(rd_code(r_name), " (wire \\code{--", rd_spec(wire), "})")
  }
  facts <- c(
    rd_spec(spec_text(param$type)),
    if (!spec_null(param$unit)) paste0("in \\strong{", rd_spec(spec_text(param$unit)), "}"),
    if (isTRUE(param$required)) "required" else paste0("default ", rd_default(param)),
    if (!spec_null(param$range)) paste0("range ", rd_code(spec_text(param$range)))
  )
  body <- paste0(paste(facts, collapse = "; "), ". ", rd_spec(spec_text(param$doc)))
  deviation <- spec_deviation(spec, "param", wire)
  if (!is.null(deviation) && is.na(deviation$r)) {
    body <- paste0(body, " \\emph{Not an argument here: ", rd_spec(deviation$why), ".}")
  }
  paste0("\\item{", label, "}{", body, "}")
}

field_key <- function(field) {
  if (is.null(field$.key)) field$wire else field$.key
}

rd_field_item <- function(spec, field) {
  wire <- field$wire
  r_name <- spec_r_field(spec, field_key(field))
  label <- if (!is.na(r_name) && !identical(r_name, wire)) {
    paste0(rd_code(r_name), " (wire \\code{", rd_spec(wire), "})")
  } else {
    rd_code(wire)
  }
  null_text <- if (isTRUE(field$nullable)) {
    paste0("\\code{NA} when ", rd_spec(spec_text(field$null_means, "the value is missing")))
  } else {
    "never \\code{NA}"
  }
  facts <- c(
    rd_spec(spec_text(field$type)),
    if (!spec_null(field$unit)) paste0("in \\strong{", rd_spec(spec_text(field$unit)), "}"),
    null_text
  )
  body <- paste0(paste(facts, collapse = "; "), ". ", rd_spec(spec_text(field$doc)))
  when <- spec_text(field$present_when)
  if (nzchar(when) && !identical(when, "--paths-stdin")) {
    body <- paste0(body, " \\emph{Present only with \\code{", rd_spec(gsub("-", "_", when)), "}.}")
  }
  paste0("\\item{", label, "}{", body, "}")
}

# Split fields into the ones this function returns and the ones it does not (yet).
split_projected <- function(spec, fields) {
  projected <- Filter(function(f) {
    d <- spec_deviation(spec, "field", field_key(f))
    is.null(d) || !isTRUE(d$pending)
  }, fields)
  pending <- Filter(function(f) isTRUE(spec_deviation(spec, "field", field_key(f))$pending), fields)
  list(projected = projected, pending = pending)
}

rd_section <- function(title, lines) {
  c(paste0("\\section{", title, "}{"), lines, "}")
}

rd_describe <- function(items) {
  c("\\describe{", items, "}")
}

# The fact sections of one man page, in the shared reference-page order (bridge
# design/verbs/README.md). Rd fixes the position of \value, \references, \seealso and \examples
# itself, so those land where R puts them; everything else follows the shared order.
render_spec_sections <- function(spec, bulk = FALSE, commit = spec_source_commit()) {
  out <- character(0)
  source_note <- paste0(
    "Generated from the bridge's verb spec \\code{", spec$.file, ".yaml} (bridge commit \\code{",
    substr(commit, 1L, 12L), "}) by \\code{scripts/build-man.R}; the spec owns these facts, ",
    "and all three bindings render the same ones."
  )

  wraps <- spec$wraps
  if (is.list(wraps) && length(wraps) > 0L) {
    items <- vapply(wraps, function(w) {
      pin <- spec_text(w$pin)
      path <- spec_text(w$path)
      symbol <- rd_code(spec_text(w$symbol))
      if (nzchar(pin) && nzchar(path)) {
        paste0("\\item \\href{", MZLIB_BLOB, "/", pin, "/", path, "}{", symbol, "} in \\code{",
          rd_spec(path), "} at mzLib \\code{", substr(pin, 1L, 8L), "}")
      } else {
        paste0("\\item ", symbol)
      }
    }, character(1L))
    out <- c(out, rd_section("Wraps", c(
      paste0("Wire verb \\code{", rd_spec(spec$verb), "}", if (bulk) " with \\code{--paths-stdin}" else "",
        ". ", source_note),
      "", "\\itemize{", items, "}"
    )))
  }

  params <- spec_params(spec, bulk)
  out <- c(out, rd_section(
    "Parameters: units, ranges and defaults",
    if (length(params) == 0L) {
      "This verb takes no parameters."
    } else {
      rd_describe(vapply(params, function(p) rd_param_item(spec, p), character(1L)))
    }
  ))

  envelope <- split_projected(spec, spec_fields(spec, "envelope_fields", bulk))
  columns_raw <- spec$result$columns
  columns <- split_projected(spec, spec_fields(spec, "columns", bulk))
  if (length(envelope$projected) > 0L) {
    out <- c(out, rd_section("Returned fields", c(
      "Each field with its type, its unit, and what \\code{NA} means when it is \\code{NA}.",
      "", rd_describe(vapply(envelope$projected, function(f) rd_field_item(spec, f), character(1L)))
    )))
  }
  table <- if (spec$verb %in% names(R_TABLE)) R_TABLE[[spec$verb]] else NULL
  if (length(columns$projected) > 0L) {
    heading <- if (is.null(table)) "Row fields" else if (startsWith(table, "(")) {
      "Columns"
    } else {
      paste0("Columns of \\code{", table, "}")
    }
    out <- c(out, rd_section(heading, c(
      if (bulk) "The first two columns say which input each row came from; rows are grouped by input, in input order, whatever \\code{threads} is." else character(0),
      rd_describe(vapply(columns$projected, function(f) rd_field_item(spec, f), character(1L)))
    )))
  } else if (identical(columns_raw, "per-format")) {
    rules <- spec$result$cell_rules
    out <- c(out, rd_section(paste0("Columns of \\code{", if (is.null(table)) "records" else table, "}"), c(
      "The columns are \\strong{per format}: this file's own mzLib record fields, named in \\code{column_names}. Every cell follows these rules:",
      "", "\\itemize{",
      vapply(if (is.list(rules)) rules else list(), function(r) paste0("\\item ", rd_spec(spec_text(r))), character(1L)),
      "}"
    )))
  }
  # Named sub-tables (FlashLFQ's peptides, proteins and peaks; a protein database's GO terms): each
  # field keyed `<table>.<wire>` in R_DEVIATIONS, because the same wire name recurs across tables.
  tables <- spec_tables(spec)
  table_pending <- list()
  for (table_name in names(tables)) {
    parts <- split_projected(spec, tables[[table_name]])
    table_pending <- c(table_pending, parts$pending)
    if (length(parts$projected) == 0L) next
    deviation <- spec_deviation(spec, "table", table_name)
    shown <- if (!is.null(deviation) && !is.na(deviation$r)) deviation$r else table_name
    out <- c(out, rd_section(paste0("Fields of \\code{", rd_spec(shown), "}"), c(
      if (!is.null(deviation)) c(paste0(rd_spec(deviation$why), "."), "") else character(0),
      rd_describe(vapply(parts$projected, function(f) rd_field_item(spec, f), character(1L)))
    )))
  }
  pending <- c(envelope$pending, columns$pending, table_pending)
  if (length(pending) > 0L) {
    out <- c(out, rd_section("On the wire but not projected yet", c(
      "The bridge sends these, and this version of mzLibR does not return them yet:",
      "", rd_describe(vapply(pending, function(f) {
        paste0("\\item{", rd_code(f$wire), "}{", rd_spec(spec_deviation(spec, "field", field_key(f))$why), ".}")
      }, character(1L)))
    )))
  }

  errors <- spec$errors
  if (is.list(errors) && length(errors) > 0L) {
    items <- vapply(errors, function(e) {
      kind <- spec_text(e$kind)
      class <- R_CONDITION[kind]
      when <- spec_text(e$when)
      paste0("\\item{", if (is.na(class)) rd_code(kind) else paste0("\\code{", class, "}"),
        " (", rd_spec(kind), ")}{",
        if (identical(when, "never")) "Never raised by this verb." else rd_spec(when), "}")
    }, character(1L))
    out <- c(out, rd_section("Errors", c(
      paste0(
        "Each is an R condition carrying the class shown and \\code{mzlib_error}; see ",
        "\\code{\\link{mzlib_error}}.",
        if (bulk) {
          ""
        } else if (r_exports(spec_r_bulk_name(spec))) {
          paste0(
            " A condition that mentions \\code{paths-stdin}, \\code{threads} or \\code{on-error} ",
            "comes only from \\code{\\link{", spec_r_bulk_name(spec), "}()}."
          )
        } else if (!is.na(spec_r_bulk_name(spec))) {
          paste0(
            " A condition that mentions \\code{paths-stdin}, \\code{threads} or \\code{on-error} ",
            "belongs to the verb's many-files form."
          )
        } else {
          ""
        }
      ),
      "", rd_describe(items)
    )))
  }

  caveats <- spec$caveats
  if (is.list(caveats) && length(caveats) > 0L) {
    out <- c(out, rd_section("Caveats", c(
      "\\itemize{", vapply(caveats, function(c) paste0("\\item ", rd_spec(spec_text(c))), character(1L)), "}"
    )))
  }

  performance <- spec_text(spec$performance)
  performance <- c(
    if (nzchar(performance)) rd_spec(performance) else character(0),
    paste(
      "Every call starts one bridge process, which costs a .NET start-up before any work.",
      if (bulk) {
        "This bulk form pays it once for the whole list, which is the reason it exists: never loop the one-path function over many files."
      } else if (r_exports(spec_r_bulk_name(spec))) {
        paste0("For many files call \\code{\\link{", spec_r_bulk_name(spec), "}()} once rather than looping this function: one process, one start-up, and the thread count stated on the wire.")
      } else {
        ""
      }
    )
  )
  out <- c(out, rd_section("Performance", performance))

  spelling <- function(lang) {
    entry <- spec$bindings[[lang]]
    name <- spec_text(entry$name)
    if (!nzchar(name)) {
      return("not yet")
    }
    text <- paste0("\\code{", rd_spec(name), "}")
    options <- spec_text(entry$options)
    if (nzchar(options)) {
      text <- paste0(text, " with \\code{", rd_spec(options), "}")
    }
    bulk_name <- spec_text(entry$bulk)
    if (nzchar(bulk_name)) {
      text <- paste0(text, "; many files: \\code{", rd_spec(bulk_name), "}")
    }
    text
  }
  out <- c(out, rd_section("Same verb in other bindings", c(
    "\\itemize{",
    paste0("\\item Python (pyMzLib): ", spelling("py")),
    paste0("\\item Rust (mzLibRust): ", spelling("rust")),
    paste0("\\item R (mzLibR): \\code{", spec_r_name(spec), "}",
      if (r_exports(spec_r_bulk_name(spec))) paste0("; many files: \\code{\\link{", spec_r_bulk_name(spec), "}}") else ""),
    "}"
  )))

  since <- spec$since
  since_text <- function(key) {
    value <- since[[key]]
    if (spec_null(value)) "not yet shipped" else rd_spec(format(value))
  }
  out <- c(out, rd_section("Since", c(
    paste0(
      "Wire protocol ", since_text("protocol"), "; pyMzLib ", since_text("pymzlib"),
      "; mzLibRust ", since_text("mzlibrust"), "; mzLibR ", since_text("mzlibr"), "."
    )
  )))

  questions <- spec$open_questions
  if (is.list(questions) && length(questions) > 0L) {
    out <- c(out, rd_section("Not yet verified", c(
      "The spec records these as open. They are listed rather than hidden:",
      "", "\\itemize{",
      vapply(questions, function(q) paste0("\\item ", rd_spec(spec_text(q))), character(1L)),
      "}"
    )))
  }

  out
}

# `\references{}` from the spec's DOIs, or NULL when it cites none.
render_spec_references <- function(spec) {
  cites <- spec$cite
  if (!is.list(cites) || length(cites) == 0L) {
    return(NULL)
  }
  c("\\references{", "\\itemize{", vapply(cites, function(c) {
    doi <- spec_text(c$doi)
    # \href rather than \doi{}, which needs R 3.6 while DESCRIPTION promises 3.5.
    paste0("\\item \\href{https://doi.org/", doi, "}{doi:", rd_spec(doi), "}: ", rd_spec(spec_text(c$`for`)))
  }, character(1L)), "}", "}")
}
