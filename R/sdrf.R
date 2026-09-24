# SDRF-Proteomics experimental-design files: read one, or pool several into one table.
#
# Every other reader here answers *what did the search find*. SDRF answers *what was searched* -
# which sample, which organism part, which replicate - and that is the half a caller needs to
# group results across experiments.
#
# Ported from pyMzLib's `sdrf.py`, which decided the verbs, the wire fields and the caveats. What
# changes here is only the projection into R's idiom: `NA` rather than `None`, 1-based positions
# rather than 0-based, an S3 class with a `print` method rather than a dataclass, and a pool's
# labels carried as the NAMES of a character vector of paths rather than as a {path: label} map.
#
# **This module is row-major, and every other reader here is columnar.** SDRF column names are
# data, and they repeat, so `columns` is a character vector that may hold duplicates and `rows` is
# a list of character vectors. A data.frame would force the names unique and the rows rectangular,
# and both are losses; `sdrf_records()` offers one anyway, for documents where neither applies.

# The column `sdrf_pool()` adds to record which document each row came from.
SDRF_SOURCE_DOCUMENT_COLUMN <- "comment[source document]"

# ---------------------------------------------------------------- parsing

# One wire row: a list of strings. A cell is never null on the wire, but if one ever were it
# would become NA in place rather than shortening the row.
sdrf_parse_row <- function(row) {
  if (!is.list(row) || length(row) == 0L) {
    return(character(0))
  }
  vapply(row, function(cell) as.character(cell)[1L], character(1L), USE.NAMES = FALSE)
}

sdrf_parse_fields <- function(data) {
  rows <- data[["rows"]]
  list(
    path = as.character(wire_field(data, "path", "character", "")),
    columns = readers_parse_strings(data, "column_names"),
    rows = if (is.list(rows)) lapply(rows, sdrf_parse_row) else list(),
    row_count = as.numeric(wire_field(data, "row_count", "numeric", 0)),
    returned_count = as.numeric(wire_field(data, "returned_count", "numeric", 0)),
    offset = as.numeric(wire_field(data, "offset", "numeric", 0)),
    truncated = isTRUE(data[["truncated"]]),
    caveats = readers_parse_strings(data, "caveats")
  )
}

sdrf_parse_document <- function(data) {
  structure(sdrf_parse_fields(data), class = "mzlibr_sdrf")
}

sdrf_parse_pooled <- function(data) {
  written <- data[["written"]]
  written <- if (is.list(written)) {
    list(
      path = as.character(wire_field(written, "path", "character", NA_character_)),
      row_count = as.numeric(wire_field(written, "row_count", "numeric", NA_real_))
    )
  } else {
    NULL
  }

  fields <- sdrf_parse_fields(data)
  # A pooled table is not a file that was read.
  fields$path <- ""
  structure(
    c(
      fields,
      list(
        document_count = as.numeric(wire_field(data, "document_count", "numeric", 0)),
        paths = readers_parse_strings(data, "paths"),
        labels = readers_parse_strings(data, "labels"),
        written = written
      )
    ),
    class = c("mzlibr_pooled_sdrf", "mzlibr_sdrf")
  )
}

# ---------------------------------------------------------------- argument assembly

sdrf_window_args <- function(limit, offset) {
  args <- character(0)
  if (!is.null(limit)) {
    # Zero is allowed, as it is by the bridge: a header-only answer.
    if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) || limit != round(limit) ||
      limit < 0) {
      stop(mzlib_usage_error(paste0(
        "limit must be a non-negative whole number, or NULL for every row; got ",
        paste(deparse(limit), collapse = " "), "."
      )))
    }
    args <- c(args, "--limit", formatC(limit, format = "d"))
  }
  if (!is.numeric(offset) || length(offset) != 1L || is.na(offset) || offset != round(offset) ||
    offset < 0) {
    stop(mzlib_usage_error(paste0(
      "offset must be a non-negative whole number; got ",
      paste(deparse(offset), collapse = " "), "."
    )))
  }
  if (offset > 0) {
    args <- c(args, "--offset", formatC(offset, format = "d"))
  }
  args
}

sdrf_build_read_args <- function(path, limit, offset) {
  c("sdrf", "read", "--path", readers_normalise_path(path), sdrf_window_args(limit, offset))
}

# The arguments and the stdin for a pool, validated before anything is spawned.
#
# Documents cross on stdin, one per line, as `path[\tlabel]`: the bridge keeps named options in a
# dictionary, so a repeated `--path` would silently keep only the last.
sdrf_build_pool_request <- function(documents, out, limit, offset) {
  if (!is.character(documents)) {
    stop(mzlib_usage_error(paste0(
      "documents must be a character vector of paths, optionally named by label; got ",
      class(documents)[1L], "."
    )))
  }
  if (length(documents) == 0L) {
    stop(mzlib_usage_error("At least one SDRF document is required."))
  }

  labels <- names(documents)
  labelled <- !is.null(labels)
  paths <- trimws(unname(documents))

  lines <- character(length(paths))
  for (i in seq_along(paths)) {
    path <- paths[[i]]
    if (is.na(path) || !nzchar(path)) {
      stop(mzlib_usage_error("A document path may not be blank."))
    }
    label <- if (labelled) trimws(labels[[i]]) else ""
    if (labelled && (is.na(label) || !nzchar(label))) {
      stop(mzlib_usage_error(paste0(
        "The label for '", path, "' is blank. Name every document, or pass an unnamed vector to ",
        "accept mzLib's path-derived default for all of them."
      )))
    }
    if (grepl("[\t\r\n]", path) || grepl("[\t\r\n]", label)) {
      stop(mzlib_usage_error(paste0(
        "A path or label contains a tab or newline, which the bridge uses to separate them: '",
        path, "' / '", label, "'."
      )))
    }
    lines[[i]] <- if (labelled) paste0(path, "\t", label) else path
  }

  args <- c("sdrf", "pool", sdrf_window_args(limit, offset))
  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  list(args = args, stdin = paste(lines, collapse = "\n"))
}

# ---------------------------------------------------------------- the public surface

#' Read one SDRF-Proteomics experimental-design file
#'
#' Every other reader here answers what a search **found**. SDRF answers what was **searched** -
#' which sample, which organism part, which replicate, which instrument settings - which is the
#' half you need to compare results across experiments.
#'
#' **Use this, not [readers_read_records()], for SDRF.** `readers_read_records()` joins each row's
#' cells into one semicolon-separated string, and SDRF's own `NT=...;AC=...` grammar puts
#' semicolons inside cells, so that string cannot be split back apart.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param limit Maximum rows to return. `NULL`, the default, returns all of them.
#' @param offset Rows to skip.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf`: `path`, `columns` (a character vector that may repeat), `rows` (a
#'   list of character vectors, one per row), `row_count` (the whole document), `returned_count`,
#'   `offset`, `truncated` and `caveats`.
#'
#' @section The shape is row-major, by necessity:
#'
#' SDRF column names are data, and they **repeat**: 649 files in the curated corpus carry
#' `comment[modification parameters]` more than once, up to eight times in one file. A data.frame
#' would force those names unique. So `columns` keeps the duplicates, `rows` holds the cells, and
#' position links them. Use [sdrf_value()] for the first cell under a name and [sdrf_all()] for
#' every one.
#'
#' **Rows are ragged.** A row may carry fewer cells than `columns` has entries - PXD059974 in
#' mzLib's own fixtures has a 46-column header with 17 of its 22 rows carrying 42 cells - and mzLib
#' preserves that so the file round-trips. [sdrf_value()] gives `NA` for a position a row does not
#' reach.
#'
#' **Cells are raw strings, never interpreted.** The key=value grammar arrives exactly as written,
#' because it cannot be told apart from a `comment[file uri]` cell whose pre-signed URL contains
#' `Signature=` and `Expires=`.
#'
#' **A reserved word is a real value.** `"not available"` and `"not applicable"` mean the
#' experiment stated an absence. `NA` means the document has no such column. Do not collapse them.
#'
#' @section Reading is not validating:
#'
#' This reads whatever the file holds and makes no claim that it is correct. mzLib's judgements
#' about a document are their own functions: [sdrf_validate()] checks the specification's
#' structural rules, [sdrf_assess()] asks whether the sample columns describe a design at all,
#' [sdrf_samples()] merges each sample's rows with its age parsed, and [sdrf_lint()] finds
#' concepts several documents wrote inconsistently.
#'
#' @seealso [sdrf_pool()], [sdrf_value()], [sdrf_all()], [sdrf_records()], [sdrf_validate()]
#' @spec sdrf.read
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' doc
#' sdrf_value(doc, "characteristics[organism]")
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_read <- function(path, limit = NULL, offset = 0, timeout = 60) {
  args <- sdrf_build_read_args(path, limit, offset)
  sdrf_parse_document(bridge_invoke(args, timeout = timeout))
}

#' Merge several SDRF documents into one analysis table
#'
#' Columns are the union of every document's, ordered by SDRF's own block structure, and a name
#' that repeats is carried at the highest multiplicity any single document used, so nothing is
#' dropped. A cell a document did not have is filled with the reserved word `"not available"`, and
#' a `comment[source document]` column records which document each row came from.
#'
#' @param documents A character vector of paths. **Name it** to choose each document's provenance
#'   label: `c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv")`. Name every element
#'   or none.
#' @param out Write the merged document here as SDRF. The **whole** document is written regardless
#'   of `limit` and `offset`.
#' @param limit Maximum rows to return. `NULL`, the default, returns all of them.
#' @param offset Rows to skip.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_pooled_sdrf`, which is also an `mzlibr_sdrf`: everything [sdrf_read()]
#'   returns, plus `document_count`, `paths`, `labels`, and `written` (a list of `path` and
#'   `row_count`, or `NULL` when `out` was not given). `row_count` counts the rows of the whole
#'   pooled document; `returned_count` the rows returned, starting `offset` rows in.
#'
#' @section Give your documents names:
#'
#' With an unnamed vector, mzLib falls back to `containing-folder/file-stem`, which depends on
#' where the files happen to sit - so the same documents pooled from a different directory produce
#' a different table. The returned `caveats` say so when that fallback was used.
#'
#' @section The result is an analysis table, not something to deposit:
#'
#' Source name, assay name and label together are unique within one document, but two
#' experiments may both have a `"Sample 1"`, so a pooled table will usually violate SDRF's
#' uniqueness rule. Use [sdrf_source_documents()] as part of any key.
#'
#' @seealso [sdrf_read()], [sdrf_source_documents()]
#' @spec sdrf.pool
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' pooled <- sdrf_pool(c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv"),
#'   limit = 4)
#' pooled
#' sdrf_source_documents(pooled)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_pool <- function(documents, out = NULL, limit = NULL, offset = 0, timeout = 60) {
  request <- sdrf_build_pool_request(documents, out, limit, offset)
  sdrf_parse_pooled(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

# ---------------------------------------------------------------- accessors

sdrf_check_document <- function(doc) {
  if (!inherits(doc, "mzlibr_sdrf")) {
    stop(mzlib_usage_error("Expected an SDRF document from sdrf_read() or sdrf_pool()."))
  }
}

sdrf_check_column <- function(column) {
  if (!is.character(column) || length(column) != 1L || is.na(column)) {
    stop(mzlib_usage_error("column must be a single column name."))
  }
}

#' Where a column first sits in an SDRF document
#'
#' Comparison is exact and case-**sensitive**, matching the SDRF specification and mzLib.
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#' @param column A column name.
#'
#' @return The 1-based position of the first column with this name, or `NA`.
#' @seealso [sdrf_indexes_of()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' sdrf_index_of(doc, "source name")
#' sdrf_index_of(doc, "no such column")
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_index_of <- function(doc, column) {
  positions <- sdrf_indexes_of(doc, column)
  if (length(positions) == 0L) NA_integer_ else positions[[1L]]
}

#' Every position a column name occupies in an SDRF document
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#' @param column A column name.
#'
#' @return The 1-based positions carrying this name, in document order; empty when the column is
#'   absent.
#' @seealso [sdrf_index_of()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' sdrf_indexes_of(doc, "comment[modification parameters]")
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_indexes_of <- function(doc, column) {
  sdrf_check_document(doc)
  sdrf_check_column(column)
  which(doc$columns == column)
}

#' The first cell under a column, one per row
#'
#' `NA` means **the document does not have this column, or the row is too short to reach it**. It
#' does not mean "empty": the SDRF reserved words `"not available"` and `"not applicable"` are real
#' values an experiment chose to write, and they come back as themselves.
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#' @param column A column name.
#'
#' @return A character vector with one element per returned row.
#' @seealso [sdrf_all()] for a column that repeats.
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' sdrf_value(doc, "characteristics[organism part]")
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_value <- function(doc, column) {
  i <- sdrf_index_of(doc, column)
  vapply(doc$rows, function(row) {
    if (is.na(i) || i > length(row)) NA_character_ else row[[i]]
  }, character(1L))
}

#' Every cell under a column, one list element per row
#'
#' The accessor for a multi-cardinality column such as `comment[modification parameters]`, which
#' legitimately repeats - up to eight times in one corpus file. Positions a row is too short to
#' reach are skipped, so each element holds only cells that exist.
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#' @param column A column name.
#'
#' @return A list with one character vector per returned row.
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' sdrf_all(doc, "comment[modification parameters]")[[1]]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_all <- function(doc, column) {
  positions <- sdrf_indexes_of(doc, column)
  lapply(doc$rows, function(row) row[positions[positions <= length(row)]])
}

#' An SDRF document as a data.frame, when its shape allows one
#'
#' The convenient shape, and a **lossy** one when a column name repeats: each name becomes one
#' column holding its **last** occurrence a row reaches, which is what pyMzLib's `records` does.
#' [sdrf_has_repeated_columns()] says whether that applies. A position a row is too short to
#' reach is `NA`.
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#'
#' @return A data.frame with one row per returned row and one column per distinct name, in order
#'   of first appearance. Names are kept verbatim (`check.names = FALSE`).
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' sdrf_has_repeated_columns(doc)   # TRUE: the data.frame keeps the last of each name
#' records <- sdrf_records(doc)
#' records[, c("source name", "characteristics[organism]")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_records <- function(doc) {
  sdrf_check_document(doc)
  distinct <- unique(doc$columns)
  columns <- lapply(distinct, function(name) {
    positions <- which(doc$columns == name)
    vapply(doc$rows, function(row) {
      reached <- positions[positions <= length(row)]
      if (length(reached) == 0L) NA_character_ else row[[max(reached)]]
    }, character(1L))
  })
  names(columns) <- distinct
  as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE, optional = TRUE)
}

#' Whether an SDRF document repeats a column name
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#' @return `TRUE` or `FALSE`. When `TRUE`, [sdrf_records()] loses cells.
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' doc <- sdrf_read("PXD000070.sdrf.tsv")
#' sdrf_has_repeated_columns(doc)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_has_repeated_columns <- function(doc) {
  sdrf_check_document(doc)
  anyDuplicated(doc$columns) > 0L
}

#' How many rows are shorter than the header
#'
#' @param doc An [sdrf_read()] or [sdrf_pool()] result.
#' @return The number of returned rows carrying fewer cells than there are columns.
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' ragged <- sdrf_read("PXD059974.sdrf.tsv")
#' sdrf_ragged_row_count(ragged)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_ragged_row_count <- function(doc) {
  sdrf_check_document(doc)
  sum(lengths(doc$rows) < length(doc$columns))
}

#' Which document each pooled row came from
#'
#' @param pooled An [sdrf_pool()] result.
#' @return The `comment[source document]` label of each returned row.
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' pooled <- sdrf_pool(c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv"),
#'   limit = 4)
#' sdrf_source_documents(pooled)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_source_documents <- function(pooled) {
  if (!inherits(pooled, "mzlibr_pooled_sdrf")) {
    stop(mzlib_usage_error("Expected a pooled SDRF table from sdrf_pool()."))
  }
  sdrf_value(pooled, SDRF_SOURCE_DOCUMENT_COLUMN)
}

# ---------------------------------------------------------------- print methods

sdrf_print_body <- function(x) {
  cat("  ", length(x$columns), " columns",
    if (sdrf_has_repeated_columns(x)) " (some names repeat - use sdrf_all())" else "", "\n",
    sep = ""
  )
  cat("  ", format(x$row_count), " rows in the document, ", format(x$returned_count), " returned",
    if (isTRUE(x$offset > 0)) paste0(" from offset ", format(x$offset)) else "", "\n",
    sep = ""
  )
  if (isTRUE(x$truncated)) {
    cat("  ! truncated - rows were left behind\n")
  }
  for (caveat in x$caveats) {
    cat("  ! ", caveat, "\n", sep = "")
  }
  invisible(x)
}

#' Print an SDRF document
#'
#' @param x An [sdrf_read()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf <- function(x, ...) {
  cat("<mzlibr_sdrf> ", basename(x$path), "\n", sep = "")
  sdrf_print_body(x)
}

#' Print a pooled SDRF table
#'
#' @param x An [sdrf_pool()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_pooled_sdrf <- function(x, ...) {
  cat("<mzlibr_pooled_sdrf> ", format(x$document_count), " documents: ",
    paste(x$labels, collapse = ", "), "\n",
    sep = ""
  )
  sdrf_print_body(x)
  if (!is.null(x$written)) {
    cat("  written to ", x$written$path, " (", format(x$written$row_count), " rows)\n", sep = "")
  }
  invisible(x)
}
