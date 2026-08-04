# Reading proteomics result files: what a file is, what can be done with it, and its records.
#
# mzLib recognises 29 result-file types written by a dozen tools - MetaMorpheus, MSFragger,
# TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur, FlashLFQ - and dispatches
# each to a parser it maintains. This module asks it what a path is.
#
# The temptation is to describe mzLib as reading 29 formats into one uniform shape. It does not,
# and the whole design of this module is about not letting anyone believe it does. The formats
# fall into disjoint families and 13 of the 29 belong to no family at all, so an empty `views` is
# a real and common answer rather than a failure.

# The view that matters most: the cross-format record projection, and the input
# `flashlfq_quantify()` accepts. Exactly three of the 29 file types have it.
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
    cells <- lapply(values, function(value) {
      if (is.null(value) || (length(value) == 1L && is.logical(value) && is.na(value))) NA else value
    })

    # A cell that is itself an array stays one, as a list column.
    #
    # `readers_read_spectra(peaks = TRUE)` returns one mz array and one intensity array PER SCAN,
    # so the cell is a vector rather than a scalar. Unlisting those would splice every scan's
    # peaks into one long vector, and the length check below would then reject the whole table
    # with a message about column lengths - reporting a shape problem for what is really the
    # correct shape. A list column is what R has for this.
    if (any(vapply(cells, length, integer(1L)) > 1L)) {
      return(I(cells))
    }

    simplified <- unlist(cells, use.names = FALSE)
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

readers_build_read_args <- function(path, limit, offset, out, verb = "read-results") {
  args <- c("readers", verb, "--path", readers_normalise_path(path))

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
#' It is tempting to read "29 formats" as "29 formats in one uniform shape". They are not. The
#' formats fall into disjoint families, and **13 of the 29 belong to none of them** — an empty
#' `views` is a real and common answer, meaning mzLib can parse the file but offers no
#' cross-format projection of it.
#'
#' - `"quantifiable"` — the cross-format record view, and the input [flashlfq_quantify()]
#'   accepts. **Exactly three file types have it**: MetaMorpheus `psmtsv` and `osmtsv`, and
#'   `MsFraggerPsm`.
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
#' @param records A [readers_read_results()], [readers_read_features()] or
#'   [readers_read_spectra()] result.
#' @param column Which retention-time column to convert. Defaults to `retention_time` for a record
#'   or scan table and to `retention_time_start` for a feature table, which has two. Name the other
#'   explicitly — `"retention_time_end"` — to convert it.
#' @return A numeric vector of retention times in minutes.
#'
#' @section When it raises:
#'
#' Whenever the unit is `"unknown"`, which is not hypothetical: it is the honest answer for TopFD
#' `_ms1.feature`, where seconds became minutes at v1.7.0 without the file type changing. Guessing
#' there is what mzLib's own deconvolution code does, and what this deliberately does not.
#'
#' @export
readers_retention_time_in_minutes <- function(records, column = NULL) {
  known <- c(
    "mzlibr_result_records", "mzlibr_feature_records", "mzlibr_scan_records",
    "mzlibr_native_records"
  )
  if (!any(inherits(records, known, which = TRUE) > 0L)) {
    stop(mzlib_usage_error(paste0(
      "records must be a readers_read_results(), readers_read_features(), ",
      "readers_read_spectra() or readers_read_records() result."
    )))
  }
  # A native record table has no declared unit - its columns are the format's own - so converting
  # one would be a guess dressed as a conversion.
  if (inherits(records, "mzlibr_native_records")) {
    stop(mzlib_usage_error(paste0(
      "readers_read_records() returns this format's own columns, which carry no declared ",
      "retention-time unit, so they cannot be converted. Use readers_read_results(), ",
      "readers_read_features() or readers_read_spectra(), which do declare one."
    )))
  }

  if (is.null(column)) {
    column <- if (inherits(records, "mzlibr_feature_records")) {
      "retention_time_start"
    } else {
      "retention_time"
    }
  }
  if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    stop(mzlib_usage_error("column must be a single column name, or NULL for the default."))
  }

  if (is.null(records$records) || !column %in% names(records$records)) {
    return(numeric(0))
  }

  values <- as.numeric(records$records[[column]])
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

# ---------------------------------------------------------------- exhaustive coverage
#
# `readers_read_results()` projects `IQuantifiableResultFile`, which three of the 29 file types
# implement. The four verbs below reach the rest, in two different ways because the gap has two
# shapes.
#
# `readers_read_records()` reads ANY of the 29 by projecting each format's own record type, so
# its columns are deliberately not uniform - a TopPIC file gives TopPIC's own 36. The other three
# project the remaining cross-format views, and are uniform in the way `read_results` is.
#
# Ported from pyMzLib's `readers.py`, which decided the verbs, the wire fields and the caveats;
# what changes here is only the projection into R's idiom - a data.frame rather than a map of
# arrays, `NA` rather than `None`, an S3 class with a `print` method rather than a dataclass.

# Shared by the four parsers: the character vector behind a wire list of strings.
readers_parse_strings <- function(data, name) {
  values <- data[[name]]
  if (!is.list(values) || length(values) == 0L) {
    return(character(0))
  }
  vapply(values, function(v) as.character(v)[1L], character(1L), USE.NAMES = FALSE)
}

# Shared by the four parsers: the `output` block, or NULL when the table came back inline.
readers_parse_output <- function(data) {
  output <- data[["output"]]
  if (!is.list(output)) {
    return(NULL)
  }
  list(
    path = as.character(wire_field(output, "path", "character", NA_character_)),
    format = as.character(wire_field(output, "format", "character", NA_character_)),
    row_count = as.numeric(wire_field(output, "row_count", "numeric", NA_real_))
  )
}

# The fields every read verb reports, so the four parsers cannot drift on them.
readers_parse_common <- function(data) {
  list(
    path = as.character(wire_field(data, "path", "character", NA_character_)),
    file_type = as.character(wire_field(data, "file_type", "character", NA_character_)),
    record_count = as.numeric(wire_field(data, "record_count", "numeric", NA_real_)),
    returned_count = as.numeric(wire_field(data, "returned_count", "numeric", NA_real_)),
    offset = as.numeric(wire_field(data, "offset", "numeric", 0)),
    truncated = isTRUE(data[["truncated"]]),
    column_names = readers_parse_strings(data, "column_names"),
    records = readers_parse_records_table(data),
    output = readers_parse_output(data)
  )
}

readers_parse_native_records <- function(data) {
  excluded <- data[["excluded_fields"]]
  # A data.frame rather than a list of lists: this is a table of (field, type, reason) and R users
  # will want to filter it. Empty with the right columns when nothing was excluded, so
  # `nrow(x$excluded_fields)` works without a NULL check.
  excluded <- if (is.list(excluded) && length(excluded) > 0L) {
    data.frame(
      field = vapply(excluded, wire_field, character(1L), "field", "character", NA_character_),
      type = vapply(excluded, wire_field, character(1L), "type", "character", NA_character_),
      reason = vapply(excluded, wire_field, character(1L), "reason", "character", NA_character_),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      field = character(0), type = character(0), reason = character(0),
      stringsAsFactors = FALSE
    )
  }

  structure(
    c(
      readers_parse_common(data),
      list(
        reader = as.character(wire_field(data, "reader", "character", NA_character_)),
        record_type = as.character(wire_field(data, "record_type", "character", NA_character_)),
        views = readers_parse_views(data),
        excluded_fields = excluded,
        failed_fields = readers_parse_strings(data, "failed_fields")
      )
    ),
    class = "mzlibr_native_records"
  )
}

readers_parse_feature_records <- function(data) {
  structure(
    c(
      readers_parse_common(data),
      list(
        retention_time_unit = as.character(
          wire_field(data, "retention_time_unit", "character", "unknown")
        ),
        caveats = readers_parse_strings(data, "caveats")
      )
    ),
    class = "mzlibr_feature_records"
  )
}

readers_parse_match_records <- function(data) {
  structure(
    c(
      readers_parse_common(data),
      list(caveats = readers_parse_strings(data, "caveats"))
    ),
    class = "mzlibr_match_records"
  )
}

readers_parse_scan_records <- function(data) {
  structure(
    c(
      readers_parse_common(data),
      list(
        reader = as.character(wire_field(data, "reader", "character", NA_character_)),
        scan_count = as.numeric(wire_field(data, "scan_count", "numeric", NA_real_)),
        ms_order = as.numeric(wire_field(data, "ms_order", "numeric", NA_real_)),
        peaks_included = isTRUE(data[["peaks_included"]]),
        retention_time_unit = as.character(
          wire_field(data, "retention_time_unit", "character", "minutes")
        ),
        caveats = readers_parse_strings(data, "caveats")
      )
    ),
    class = "mzlibr_scan_records"
  )
}

#' Read any file mzLib recognises, into that format's own fields
#'
#' The exhaustive verb: if [readers_identify()] succeeds on a path, this reads it. All 29 file
#' types, including the 13 that belong to no cross-format view at all - TopPIC, Crux, MSFragger's
#' peptide and protein tables, the FlashDeconv formats - which no other function here can touch.
#'
#' @param path Path to any file mzLib recognises. A Bruker `.d` directory is also accepted.
#' @param limit Maximum records to return. `NULL`, the default, returns all of them.
#' @param offset Records to skip. A window, not a cursor - see [readers_read_results()].
#' @param out Write a **tab-separated** table here and return only a summary. Read it with
#'   `read.delim()`, not `read.csv()`.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_native_records`. `records` is a data.frame of this format's own fields, or
#'   `NULL` when `out` was given.
#'
#' @section The columns are not uniform, by design:
#'
#' They are this format's own mzLib record fields, under mzLib's names in `snake_case`: a TopPIC
#' file gives 36 columns, a Crux file 23, an experiment annotation 5. Read `column_names`, and use
#' [readers_read_results()], [readers_read_features()] or [readers_read_matches()] when you need
#' columns that mean the same thing across formats.
#'
#' Because the names are mzLib's own, they are **cross-referenceable against the mzLib source**: a
#' column called `e_value` is `ToppicPrsm$EValue`, and `record_type` names the class to look in.
#'
#' @section Nothing is silently dropped:
#'
#' `excluded_fields` is a data.frame of the fields that could **not** become columns, each with the
#' reason - a nested object or a dictionary has no faithful column shape, and inventing one would
#' mean publishing a schema mzLib does not have. A column that simply vanished would be
#' indistinguishable from a field the format does not have.
#'
#' `failed_fields` names the fields that **raised** while being read, with the exception type.
#' Several mzLib properties are computed and assume a UniProt-style FASTA header - Crux's and
#' MsPathFinderT's `accession` are both `protein_id` split on a pipe - so on other databases they
#' throw. Those cells arrive `NA` rather than taking the whole read down, but a failure must not
#' look like missing data.
#'
#' Note that mzLib's documented `-1` "absent" sentinel is **not** mapped to `NA` here, unlike in
#' [readers_read_results()]. In a format's own columns `-1` is frequently a real measurement - a
#' mass difference, a delta, TopPIC's `feature_score` - and nulling those would destroy data.
#'
#' @seealso [readers_identify()], [readers_read_results()]
#' @export
readers_read_records <- function(path, limit = NULL, offset = 0, out = NULL, timeout = NULL) {
  args <- readers_build_read_args(path, limit, offset, out, "read-records")
  readers_parse_native_records(bridge_invoke(args, timeout = timeout))
}

#' Read deconvolved MS1 features, in the cross-format `ms1_features` view
#'
#' Two file types offer it: TopFD/FLASHDeconv `_ms1.feature` and Dinosaur `.feature.tsv`. A file
#' without the view raises, with a message naming the views it does have.
#'
#' @param path Path to an `_ms1.feature` or Dinosaur `.feature.tsv`.
#' @param limit Maximum features to return. `NULL`, the default, returns all of them.
#' @param offset Features to skip.
#' @param out Write a tab-separated table here and return only a summary.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_feature_records`. `records` is a data.frame with `mz`, `charge`,
#'   `retention_time_start`, `retention_time_end`, `intensity` and `number_of_isotopes`.
#'
#' @section One row is not one line of the file, for `_ms1.feature`:
#'
#' An `_ms1.feature` row is a deconvolved **neutral mass spanning a charge range**, and mzLib
#' expands it into one single-charge feature per charge in that range. A hundred-feature file can
#' read as a thousand rows. **Dinosaur is one-for-one.** Either way [readers_read_records()] gives
#' you the file's own rows.
#'
#' `intensity` is the **apex** intensity, not the sum over the feature. Both formats carry a summed
#' intensity column too, and [readers_read_records()] has it.
#'
#' @section The retention-time unit is genuinely unknown for `_ms1.feature`:
#'
#' TopFD wrote seconds through v1.6.2 and minutes from v1.7.0 - *within the same file type*, with
#' nothing in the file to say which. mzLib normalises neither, and its own deconvolution code
#' resorts to a heuristic (divide by 60 if the largest end time exceeds 500). This package will not
#' launder a guess into a stated fact, so [readers_retention_time_in_minutes()] raises rather than
#' converting. Dinosaur reports minutes and converts without complaint.
#'
#' @seealso [readers_read_records()], [readers_retention_time_in_minutes()]
#' @export
readers_read_features <- function(path, limit = NULL, offset = 0, out = NULL, timeout = NULL) {
  args <- readers_build_read_args(path, limit, offset, out, "read-features")
  readers_parse_feature_records(bridge_invoke(args, timeout = timeout))
}

#' Read identifications, in the cross-format `spectral_match` view
#'
#' Four file types offer it: MsPathFinderT's targets, decoys and combined results, and Casanovo's
#' `.mztab`. These are the identification formats that share no *file*-level interface, so
#' [readers_read_results()] cannot reach them.
#'
#' @param path Path to an MsPathFinderT `_IcTarget.tsv` / `_IcDecoy.tsv` / `_IcTDA.tsv`, or a
#'   Casanovo `.mztab`.
#' @param limit Maximum matches to return. `NULL`, the default, returns all of them.
#' @param offset Matches to skip.
#' @param out Write a tab-separated table here and return only a summary.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_match_records`. `records` is a data.frame with `file_name_without_extension`,
#'   `one_based_scan_number`, `base_sequence`, `full_sequence`, `accession`, `is_decoy`,
#'   `modifications` and `modification_count`.
#'
#' @section Nothing here is FDR-filtered, and there is nothing to filter on:
#'
#' mzLib's `ISpectralMatch` carries identity fields only. Every one of these formats records an
#' E-value or q-value somewhere; [readers_read_records()] will give you those columns. Filter
#' before you report.
#'
#' @section Two is_decoy traps, both reported in caveats:
#'
#' **MsPathFinderT** infers decoys from the protein *name* - mzLib reports a decoy when the name
#' starts with `XXX`. A database whose decoys carry a different prefix reads entirely as targets.
#'
#' **Casanovo** is de novo and writes no target/decoy label at all. mzLib leaves the field at its
#' default `FALSE` and never assigns it, so `FALSE` would mean *unknown*; it arrives as `NA`
#' instead - the rule [readers_read_results()] already applies to MSFragger.
#'
#' @seealso [readers_read_records()]
#' @export
readers_read_matches <- function(path, limit = NULL, offset = 0, out = NULL, timeout = NULL) {
  args <- readers_build_read_args(path, limit, offset, out, "read-matches")
  readers_parse_match_records(bridge_invoke(args, timeout = timeout))
}

#' Read the scans of a spectra file: headers always, peaks on request
#'
#' Seven file types offer the `spectra` view: `.mzML`, `.mgf`, `_ms1.msalign`, `_ms2.msalign`,
#' Thermo `.raw`, Bruker `.d` and timsTOF `.d`.
#'
#' Retention times here **are** in minutes for every format - mzLib's spectra readers convert at
#' the boundary, unlike its result-file readers, which pass the tool's own unit through untouched.
#'
#' @param path Path to a spectra file. A Bruker `.d` directory is also accepted.
#' @param limit Maximum scans to return. `NULL`, the default, returns all of them.
#' @param offset Scans to skip, applied **after** `ms_order`.
#' @param ms_order Keep only scans at this MS level - `1` for survey scans, `2` for fragment scans.
#'   `NULL`, the default, keeps every scan. Applied before `offset` and `limit`, so
#'   `ms_order = 2, limit = 10` means the first ten MS2 scans rather than the MS2 scans among the
#'   first ten.
#' @param peaks Include the `mz` and `intensity` arrays. `FALSE` by default, and worth leaving so
#'   unless you need them: a scan header is tens of bytes and its peak list is thousands, and a
#'   mid-size mzML holds tens of thousands of scans. `peak_count` still reports how many peaks each
#'   scan has.
#' @param out Write a tab-separated table here and return only a summary. With `peaks = TRUE` each
#'   cell holds a `;`-joined list.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely. A large `.raw` legitimately
#'   takes a while.
#'
#' @return An `mzlibr_scan_records`. `scan_count` is the file's total **before** any `ms_order`
#'   filter, reported alongside `record_count` so a filter that matched nothing can never look like
#'   an empty file.
#'
#' @section Two of the seven need Windows:
#'
#' Bruker `.d` and timsTOF `.d` are read through vendor native libraries (`baf2sql`, `timsdata`)
#' and are **Windows-x64 only**. Thermo `.raw` uses managed vendor assemblies and works everywhere.
#' msalign files hold **deconvolved neutral masses**, not raw m/z - do not re-deconvolve them.
#'
#' @seealso [readers_read_records()]
#' @export
readers_read_spectra <- function(path, limit = NULL, offset = 0, ms_order = NULL, peaks = FALSE,
                                 out = NULL, timeout = NULL) {
  args <- readers_build_read_args(path, limit, offset, out, "read-spectra")

  if (!is.null(ms_order)) {
    if (!is.numeric(ms_order) || length(ms_order) != 1L || is.na(ms_order) ||
      ms_order != round(ms_order) || ms_order < 1) {
      stop(mzlib_usage_error(paste0(
        "ms_order must be a whole number of 1 or more, or NULL for every scan; got ",
        paste(deparse(ms_order), collapse = " "), "."
      )))
    }
    args <- c(args, "--ms-order", formatC(ms_order, format = "d"))
  }

  if (!is.logical(peaks) || length(peaks) != 1L || is.na(peaks)) {
    stop(mzlib_usage_error("peaks must be TRUE or FALSE."))
  }
  if (isTRUE(peaks)) {
    args <- c(args, "--peaks")
  }

  readers_parse_scan_records(bridge_invoke(args, timeout = timeout))
}

# ---------------------------------------------------------------- print methods

# Shared by the four print methods: the counts line, the truncation flag, the caveats and the
# written-table line. Written once so a verb added later cannot quietly print less.
readers_print_body <- function(x, unit_label = NULL, caveats = x$caveats) {
  cat("  ", format(x$record_count), " records in the file, ",
    format(x$returned_count), " returned",
    if (isTRUE(x$offset > 0)) paste0(" from offset ", format(x$offset)) else "", "\n",
    sep = ""
  )
  if (isTRUE(x$truncated)) {
    cat("  ! truncated - records were left behind\n")
  }
  if (!is.null(unit_label)) {
    cat("  retention_time_unit: ", unit_label, "\n", sep = "")
  }
  for (caveat in caveats) {
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

#' Print a native record table
#'
#' @param x A [readers_read_records()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_native_records <- function(x, ...) {
  cat("<mzlibr_native_records> ", basename(x$path), " (", x$file_type, ")\n", sep = "")
  cat("  ", length(x$column_names), " columns from ", x$record_type,
    if (length(x$views) == 0L) {
      " - no cross-format view"
    } else {
      paste0(" - views: ", paste(x$views, collapse = ", "))
    },
    "\n",
    sep = ""
  )
  # This verb reports no caveats: its columns are the format's own, so there is no uniform meaning
  # for them to fail to have. What it reports instead is what could not be projected.
  readers_print_body(x, caveats = character(0))
  if (nrow(x$excluded_fields) > 0L) {
    cat("  ", nrow(x$excluded_fields), " field(s) could not become columns: ",
      paste(x$excluded_fields$field, collapse = ", "), "\n",
      sep = ""
    )
  }
  # A field that threw is not missing data, and must not be read as such.
  for (failed in x$failed_fields) {
    cat("  ! field failed to read: ", failed, "\n", sep = "")
  }
  invisible(x)
}

#' Print an MS1 feature table
#'
#' @param x A [readers_read_features()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_feature_records <- function(x, ...) {
  cat("<mzlibr_feature_records> ", basename(x$path), " (", x$file_type, ")\n", sep = "")
  readers_print_body(x, x$retention_time_unit)
}

#' Print a spectral-match table
#'
#' @param x A [readers_read_matches()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_match_records <- function(x, ...) {
  cat("<mzlibr_match_records> ", basename(x$path), " (", x$file_type, ")\n", sep = "")
  readers_print_body(x)
}

#' Print a scan table
#'
#' @param x A [readers_read_spectra()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_scan_records <- function(x, ...) {
  cat("<mzlibr_scan_records> ", basename(x$path), " (", x$file_type, ", ", x$reader, ")\n", sep = "")
  cat("  ", format(x$scan_count), " scans in the file",
    if (!is.na(x$ms_order)) paste0(", filtered to MS", format(x$ms_order)) else "", "\n",
    sep = ""
  )
  readers_print_body(x, x$retention_time_unit)
  if (!isTRUE(x$peaks_included)) {
    cat("  peaks not included - pass peaks = TRUE for the mz and intensity arrays\n")
  }
  invisible(x)
}
