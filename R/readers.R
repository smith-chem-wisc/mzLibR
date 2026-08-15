# Reading proteomics result files: what a file is, what can be done with it, and its records.
#
# mzLib recognises 31 result-file types written by a dozen tools - MetaMorpheus, MSFragger,
# TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur, DIA-NN, FlashLFQ - and dispatches
# each to a parser it maintains. This module asks it what a path is.
#
# The temptation is to describe mzLib as reading 31 formats into one uniform shape. It does not,
# and the whole design of this module is about not letting anyone believe it does. The formats
# fall into disjoint families and 14 of the 31 belong to no family at all, so an empty `views` is
# a real and common answer rather than a failure.

# The view that matters most: the cross-format record projection, and the input
# `flashlfq_quantify()` accepts. Exactly four of the 31 file types have it.
READERS_QUANTIFIABLE <- "quantifiable"

# ---------------------------------------------------------------- parsing

readers_parse_views <- function(entry) {
  views <- entry[["views"]]
  if (!is.list(views) || length(views) == 0L) {
    return(character(0))
  }
  vapply(views, function(view) as.character(view)[1L], character(1L), USE.NAMES = FALSE)
}

readers_parse_formats <- function(data) {
  entries <- data[["formats"]]
  if (!is.list(entries)) {
    entries <- list()
  }
  if (length(entries) == 0L) {
    return(data.frame(
      file_type = character(0), extension = character(0), reader = character(0),
      is_quantifiable = logical(0), stringsAsFactors = FALSE
    ))
  }

  views <- lapply(entries, readers_parse_views)
  formats <- data.frame(
    file_type = vapply(entries, wire_field, character(1L), "file_type", "character", NA_character_),
    extension = vapply(entries, wire_field, character(1L), "extension", "character", NA_character_),
    reader = vapply(entries, wire_field, character(1L), "reader", "character", NA_character_),
    is_quantifiable = vapply(views, function(v) READERS_QUANTIFIABLE %in% v, logical(1L)),
    stringsAsFactors = FALSE
  )
  # A list column, because a format may have several views or none, and a filtered frame must
  # keep them - the same reason PRIDE's `locations` is one.
  formats$views <- views
  formats
}

readers_parse_file_info <- function(data) {
  structure(
    list(
      path = as.character(wire_field(data, "path", "character", NA_character_)),
      file_type = as.character(wire_field(data, "file_type", "character", NA_character_)),
      extension = as.character(wire_field(data, "extension", "character", NA_character_)),
      reader = as.character(wire_field(data, "reader", "character", NA_character_)),
      views = readers_parse_views(data),
      is_quantifiable = READERS_QUANTIFIABLE %in% readers_parse_views(data)
    ),
    class = "mzlibr_file_info"
  )
}

# The `columns` map becomes a data.frame directly.
#
# pyMzLib returns field -> list of values because that is what pandas accepts; in R the
# equivalent *is* a data.frame, so `ResultRecords.records` - the row-wise view Python needs as a
# second accessor - has no counterpart here. One object is both.
readers_parse_records_table <- function(data) {
  columns <- data[["columns"]]
  names_in_order <- data[["column_names"]]
  order <- if (is.list(names_in_order) && length(names_in_order) > 0L) {
    vapply(names_in_order, function(n) as.character(n)[1L], character(1L), USE.NAMES = FALSE)
  } else if (is.list(columns)) {
    names(columns)
  } else {
    character(0)
  }

  if (!is.list(columns) || length(order) == 0L) {
    return(NULL)
  }

  built <- lapply(order, function(name) {
    values <- columns[[name]]
    if (!is.list(values)) {
      return(rep(NA, 0L))
    }
    # Every value individually, so a `null` in the middle of a numeric column becomes NA rather
    # than shortening the column - the same hazard as everywhere else in this package.
    unlisted <- lapply(values, function(value) {
      if (is.null(value) || (length(value) == 1L && is.logical(value) && is.na(value))) NA else value
    })
    simplified <- unlist(unlisted, use.names = FALSE)
    if (is.null(simplified)) NA else simplified
  })
  names(built) <- order

  lengths_seen <- vapply(built, length, integer(1L))
  if (length(unique(lengths_seen)) > 1L) {
    stop(mzlib_protocol_error(paste0(
      "The record columns have different lengths (",
      paste(paste0(order, "=", lengths_seen), collapse = ", "),
      "), so they cannot form a table."
    )))
  }

  as.data.frame(built, stringsAsFactors = FALSE, optional = TRUE)
}

readers_parse_records <- function(data) {
  caveats <- data[["caveats"]]
  caveats <- if (is.list(caveats) && length(caveats) > 0L) {
    vapply(caveats, function(c) as.character(c)[1L], character(1L), USE.NAMES = FALSE)
  } else {
    character(0)
  }

  column_names <- data[["column_names"]]
  column_names <- if (is.list(column_names) && length(column_names) > 0L) {
    vapply(column_names, function(n) as.character(n)[1L], character(1L), USE.NAMES = FALSE)
  } else {
    character(0)
  }

  output <- data[["output"]]
  written <- if (is.list(output)) {
    list(
      path = as.character(wire_field(output, "path", "character", NA_character_)),
      format = as.character(wire_field(output, "format", "character", NA_character_)),
      row_count = as.numeric(wire_field(output, "row_count", "numeric", NA_real_))
    )
  } else {
    NULL
  }

  structure(
    list(
      path = as.character(wire_field(data, "path", "character", NA_character_)),
      file_type = as.character(wire_field(data, "file_type", "character", NA_character_)),
      record_count = as.numeric(wire_field(data, "record_count", "numeric", NA_real_)),
      returned_count = as.numeric(wire_field(data, "returned_count", "numeric", NA_real_)),
      offset = as.numeric(wire_field(data, "offset", "numeric", 0)),
      truncated = isTRUE(data[["truncated"]]),
      retention_time_unit = as.character(
        wire_field(data, "retention_time_unit", "character", "unknown")
      ),
      rows_not_read = as.numeric(wire_field(data, "rows_not_read", "numeric", NA_real_)),
      caveats = caveats,
      column_names = column_names,
      records = readers_parse_records_table(data),
      output = written
    ),
    class = "mzlibr_result_records"
  )
}

# ---------------------------------------------------------------- argument assembly

readers_normalise_path <- function(path, name = "path") {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    stop(mzlib_usage_error(paste0(name, " must be a single non-empty file path.")))
  }
  trimws(path)
}

readers_build_read_args <- function(path, limit, offset, out) {
  args <- c("readers", "read-results", "--path", readers_normalise_path(path))

  if (!is.null(limit)) {
    if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) || limit != round(limit) ||
      limit < 1) {
      stop(mzlib_usage_error(paste0(
        "limit must be a positive whole number, or NULL for every record; got ",
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

  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  args
}

# ---------------------------------------------------------------- the public surface

#' Every file type mzLib can recognise
#'
#' Enumerated from mzLib itself rather than from a list maintained here, so it reflects the
#' installed version and cannot go stale.
#'
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return A data.frame with one row per format: `file_type`, `extension`, `reader`,
#'   `is_quantifiable`, and a `views` list column.
#'
#' @section Views, and why most formats have none:
#'
#' It is tempting to read "31 formats" as "31 formats in one uniform shape". They are not. The
#' formats fall into disjoint families, and **14 of the 31 belong to none of them** — an empty
#' `views` is a real and common answer, meaning mzLib can parse the file but offers no
#' cross-format projection of it.
#'
#' - `"quantifiable"` — the cross-format record view, and the input [flashlfq_quantify()]
#'   accepts. **Exactly four file types have it**: MetaMorpheus `psmtsv` and `osmtsv`,
#'   `MsFraggerPsm`, and DIA-NN `DiaNnReport`.
#' - `"ms1_features"` — deconvolved MS1 features (TopFD `_ms1.feature`, Dinosaur).
#' - `"spectra"` — the file is spectra, not results.
#' - `"spectral_match"` — identifications that share no file-level interface.
#'
#' Note also that `extension` is **not unique**: `BrukerD` and `BrukerTimsTof` are both `.d`,
#' told apart by what the directory holds, and several formats share `.tsv`.
#'
#' @seealso [readers_identify()]
#' @export
readers_formats <- function(timeout = 60) {
  readers_parse_formats(bridge_invoke(c("readers", "formats"), timeout = timeout))
}

#' Identify a result file without parsing its contents
#'
#' Cheap by design: mzLib resolves the type and stops, so identifying a million-row file costs no
#' more than identifying an empty one.
#'
#' It is not, however, *pure*. mzLib disambiguates a bare `.tsv` by reading its first line, a
#' `.mztab` by its first five, and a Bruker `.d` by which analysis file the directory holds, so
#' an unreadable file raises rather than returning a guess.
#'
#' @param path Path to a result file.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_file_info`: `path`, `file_type`, `extension`, `reader`, `views` and
#'   `is_quantifiable`.
#'
#' @section Ask this before quantifying:
#'
#' `is_quantifiable` is the precondition for [flashlfq_quantify()] — when `FALSE`, mzLib can
#' still read the file, it simply has no uniform view, and quantification would fail on it.
#'
#' **But `TRUE` is not permission.** It reports what mzLib's *interface* offers, not that the
#' numbers are comparable. `MsFraggerPsm` is quantifiable by interface and should not be
#' quantified: among other things its retention times are in seconds while MetaMorpheus's are in
#' minutes, and mzLib does not normalise them. See [readers_read_results()] and the `caveats` it
#' returns.
#'
#' @seealso [readers_formats()], [readers_read_results()]
#' @export
readers_identify <- function(path, timeout = 60) {
  args <- c("readers", "identify", "--path", readers_normalise_path(path))
  readers_parse_file_info(bridge_invoke(args, timeout = timeout))
}

#' Read a result file into the uniform record view
#'
#' Only the three file types offering the `"quantifiable"` view can be read this way — check
#' [readers_identify()] first, or catch the error, which names the views the file does have.
#'
#' @param path Path to a MetaMorpheus `.psmtsv` / `.osmtsv` or an MSFragger `psm.tsv`.
#' @param limit Maximum records to return. `NULL`, the default, returns all of them.
#'
#'   **There is no default row limit**, deliberately. A result file can carry a million rows, and
#'   truncating by default would mean the ordinary call returns a table that looks complete and
#'   is not. `truncated` reports whether anything was left behind.
#' @param offset Records to skip.
#'
#'   **This is a window, not a cursor.** mzLib materialises the whole file on every call — its
#'   readers look lazy and are not — so paging re-reads and re-parses the file once per page. A
#'   loop over pages is quadratic. For a large file use `out` instead, in one call.
#' @param out Write the records to this path as a **tab-separated** table and return only a
#'   summary, instead of carrying them back in the envelope. The intended path for large files.
#'
#'   Tab-separated, not comma-separated, because these fields contain commas: MSFragger's mapped
#'   proteins are a comma-separated list inside a single field. Read it with
#'   `read.delim(path)`, not `read.csv(path)`.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely. A large file legitimately
#'   takes a while.
#'
#' @return An `mzlibr_result_records`. `records` is a data.frame of the record view, or `NULL`
#'   when `out` was given. `record_count` counts the **whole file** regardless of `limit` and
#'   `offset`; `returned_count` counts what came back.
#'
#' @section Two fields to read before you trust the table:
#'
#' `rows_not_read` — data rows that did not become records. **mzLib drops a malformed row
#' silently**, so a non-zero value here means the file is partly unreadable and your table is
#' incomplete. `NA` when the count could not be established.
#'
#' `caveats` — what the uniform view cannot be trusted to mean for this format, each citing the
#' mzLib source it came from. Empty for some formats and not for others. This is where you learn
#' that, e.g., TopPIC retention times are seconds while MetaMorpheus's and MSFragger's are minutes.
#'
#' Relatedly, `retention_time_unit` is per format and mzLib does **not** normalise it. Convert
#' with [readers_retention_time_in_minutes()] rather than by hand.
#'
#' @seealso [readers_identify()], [readers_retention_time_in_minutes()]
#' @export
readers_read_results <- function(path, limit = NULL, offset = 0, out = NULL, timeout = NULL) {
  args <- readers_build_read_args(path, limit, offset, out)
  readers_parse_records(bridge_invoke(args, timeout = timeout))
}

#' Retention times in minutes, whatever unit the format wrote
#'
#' The conversion you would otherwise write by hand, using `retention_time_unit`.
#'
#' **Raises when the unit is `"unknown"` rather than guessing.** A silently unconverted time axis
#' is the specific mistake this module exists to prevent: TopPIC still writes seconds while
#' MetaMorpheus and MSFragger write minutes (mzLib normalises MSFragger since PR #1116, but not
#' TopPIC), and a 60x error in a retention-time comparison looks like a chromatography problem
#' rather than a units problem.
#'
#' @param records A [readers_read_results()] result.
#' @return A numeric vector of retention times in minutes.
#' @export
readers_retention_time_in_minutes <- function(records) {
  if (!inherits(records, "mzlibr_result_records")) {
    stop(mzlib_usage_error("records must be a readers_read_results() result."))
  }
  if (is.null(records$records) || !"retention_time" %in% names(records$records)) {
    return(numeric(0))
  }

  values <- as.numeric(records$records$retention_time)
  if (identical(records$retention_time_unit, "minutes")) {
    return(values)
  }
  if (identical(records$retention_time_unit, "seconds")) {
    return(values / 60)
  }
  stop(mzlib_usage_error(paste0(
    "Cannot convert retention time for '", records$file_type,
    "': mzLib gives no basis to say what unit it is in. Inspect the values against scan numbers ",
    "before comparing them."
  )))
}

#' Print a file identification
#'
#' @param x A [readers_identify()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_file_info <- function(x, ...) {
  cat("<mzlibr_file_info> ", basename(x$path), "\n", sep = "")
  cat("  ", x$file_type, " (", x$extension, "), read by ", x$reader, "\n", sep = "")
  if (length(x$views) == 0L) {
    cat("  no views - mzLib can read this file but offers no cross-format projection of it\n")
  } else {
    cat("  views: ", paste(x$views, collapse = ", "), "\n", sep = "")
  }
  if (isTRUE(x$is_quantifiable)) {
    cat("  quantifiable: yes - but see ?readers_identify before trusting the numbers\n")
  }
  invisible(x)
}

#' Print a record view
#'
#' @param x A [readers_read_results()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_result_records <- function(x, ...) {
  cat("<mzlibr_result_records> ", basename(x$path), " (", x$file_type, ")\n", sep = "")
  cat("  ", format(x$record_count), " records in the file, ",
    format(x$returned_count), " returned",
    if (x$offset > 0) paste0(" from offset ", format(x$offset)) else "", "\n",
    sep = ""
  )
  if (isTRUE(x$truncated)) {
    cat("  ! truncated - records were left behind\n")
  }
  # A silently dropped row is the failure this warning exists for.
  if (!is.na(x$rows_not_read) && x$rows_not_read > 0) {
    cat("  ! ", format(x$rows_not_read),
      " data row(s) did not become records; this file is partly unreadable\n",
      sep = ""
    )
  }
  cat("  retention_time_unit: ", x$retention_time_unit, "\n", sep = "")
  for (caveat in x$caveats) {
    cat("  ! ", caveat, "\n", sep = "")
  }
  if (!is.null(x$output)) {
    cat("  written to ", x$output$path, " (", x$output$format, ", ",
      format(x$output$row_count), " rows)\n",
      sep = ""
    )
  }
  invisible(x)
}
