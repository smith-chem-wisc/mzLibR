# Many inputs in one bridge call: the `_many` functions, and the batch object they return.
#
# The bridge's rule (design/PARALLELISM.md, design/BULK.md in the bridge repository) is that
# parallelism is expressed once, on the wire, and never as a loop in a binding. A `_many` function
# therefore hands its WHOLE list to ONE bridge process, on stdin, with the thread count stated
# explicitly; the bridge reads `threads` inputs at once and returns one long table. **No loop over
# files happens in this file, and none may be added**: one process per file would pay .NET
# start-up once per file, and would leave nobody owning the total thread count.
#
# The spelling is the same in all three bindings: the one-path function is unchanged, and the
# many-path form is a separate `<verb>_many(paths, threads = 1, on_error = "fail")`.
#
# What comes back, in R's idiom:
#
#   * `records` - the one long table, as a data.frame. Its first two columns say which input each
#     row came from: `source_index` and `source_path`. **`source_index` is 1-based here** (the
#     wire's is 0-based), so `batch$files[batch$records$source_index, ]` finds each row's file.
#   * `files` - one row per input, in input order: the per-file facts a one-path read reports at
#     its top level, under the same names. List-valued facts (caveats, absent_fields, ...) are list
#     columns; `error` is split into `error_kind`, `error_type` and `error_message`, `NA` for an
#     input that was read; a nested object such as a spectra file's `source` is flattened into its
#     own columns.
#   * every other top-level field of the envelope, under its wire name.
#
# The per-file block is parsed generically - every field the bridge sends becomes a column - so a
# field the bridge adds later arrives without a change here (CONTRACT.md: an added field is
# optional forever).

# ---------------------------------------------------------------- arguments

bulk_check_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
    threads != round(threads) || (threads < 1 && threads != -1)) {
    stop(mzlib_usage_error(paste0(
      "threads must be a whole number of 1 or more, or -1 for one per core; got ",
      paste(deparse(threads), collapse = " "), "."
    )))
  }
  formatC(threads, format = "d")
}

bulk_check_on_error <- function(on_error) {
  if (!is.character(on_error) || length(on_error) != 1L || is.na(on_error) ||
    !on_error %in% c("fail", "skip")) {
    stop(mzlib_usage_error(paste0(
      "on_error must be \"fail\" or \"skip\"; got ", paste(deparse(on_error), collapse = " "), "."
    )))
  }
  on_error
}

# The list of paths as the bridge's stdin: one per line, validated before anything is spawned.
bulk_paths_stdin <- function(paths, fn) {
  if (!is.character(paths)) {
    stop(mzlib_usage_error(paste0(
      fn, "() takes a character vector of paths; got ", class(paths)[1L], "."
    )))
  }
  if (length(paths) == 0L) {
    stop(mzlib_usage_error(paste0(fn, "() needs at least one path.")))
  }
  paths <- trimws(unname(paths))
  for (i in seq_along(paths)) {
    if (is.na(paths[[i]]) || !nzchar(paths[[i]])) {
      stop(mzlib_usage_error(paste0("Path ", i, " of the list is blank or NA.")))
    }
    # The list travels one path per line, so a line break inside a path would split it in two.
    if (grepl("[\r\n]", paths[[i]])) {
      stop(mzlib_usage_error(paste0("Path ", i, " of the list contains a line break: '", paths[[i]], "'.")))
    }
  }
  paste0(paste(paths, collapse = "\n"), "\n")
}

# `words` is the verb, e.g. c("readers", "read-spectra"). Returns list(args, stdin).
bulk_request <- function(words, paths, threads, on_error, out = NULL, fn) {
  stdin <- bulk_paths_stdin(paths, fn)
  args <- c(
    words, "--paths-stdin",
    "--threads", bulk_check_threads(threads),
    "--on-error", bulk_check_on_error(on_error)
  )
  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  list(args = args, stdin = stdin)
}

# ---------------------------------------------------------------- parsing

# A wire value that is null, however it arrived.
wire_null <- function(value) {
  is.null(value) || (length(value) == 1L && !is.list(value) && is.na(value))
}

# A list of strings as a character vector (empty for null).
wire_strings <- function(value) {
  if (!is.list(value) || length(value) == 0L) {
    return(character(0))
  }
  vapply(value, function(v) if (wire_null(v)) NA_character_ else as.character(v)[1L], character(1L),
    USE.NAMES = FALSE
  )
}

# A list of flat objects as a data.frame, one row per object, columns in first-seen order. Every
# cell a row lacks is NA. `empty` names the columns of the zero-row frame.
wire_objects <- function(value, empty = character(0)) {
  if (!is.list(value) || length(value) == 0L) {
    columns <- rep(list(character(0)), length(empty))
    names(columns) <- empty
    return(as.data.frame(columns, stringsAsFactors = FALSE, optional = TRUE))
  }
  keys <- unique(unlist(lapply(value, names)))
  columns <- lapply(keys, function(key) {
    cells <- lapply(value, function(entry) {
      cell <- entry[[key]]
      if (wire_null(cell)) NA else cell
    })
    if (any(vapply(cells, function(c) is.list(c) || length(c) != 1L, logical(1L)))) {
      return(I(cells))
    }
    unlist(cells, use.names = FALSE)
  })
  names(columns) <- keys
  as.data.frame(columns, stringsAsFactors = FALSE, optional = TRUE)
}

# The `excluded_fields` block: field, type, reason, and the verb that carries the field instead -
# always those four columns, so `nrow()` and `$verb` work without a NULL check.
wire_excluded <- function(value) {
  frame <- wire_objects(value, empty = c("field", "type", "reason", "verb"))
  for (column in setdiff(c("field", "type", "reason", "verb"), names(frame))) {
    frame[[column]] <- rep(NA_character_, nrow(frame))
  }
  frame[, c("field", "type", "reason", "verb", setdiff(names(frame), c("field", "type", "reason", "verb"))), drop = FALSE]
}

# A files[] entry's `error`, or NULL when the input was read.
wire_error <- function(value) {
  if (!is.list(value) || length(value) == 0L) {
    return(NULL)
  }
  list(
    kind = wire_field(value, "kind", "character", NA_character_),
    type = wire_field(value, "type", "character", NA_character_),
    message = wire_field(value, "message", "character", NA_character_)
  )
}

# One cell of the files data.frame, from one field of one files[] entry.
bulk_file_cell <- function(value) {
  if (wire_null(value)) {
    return(NA)
  }
  value
}

# The `files` array as a data.frame: one row per input, in input order.
#
# Scalars become columns; lists of strings become list columns of character vectors; a list of
# objects (excluded_fields, skipped) becomes a list column of data.frames; `error` becomes three
# columns; any other object (a spectra file's `source`) is flattened into its own columns.
bulk_parse_files <- function(entries) {
  if (!is.list(entries) || length(entries) == 0L) {
    return(data.frame(path = character(0), stringsAsFactors = FALSE))
  }

  flat <- lapply(entries, function(entry) {
    row <- list()
    for (key in names(entry)) {
      value <- entry[[key]]
      if (identical(key, "error")) {
        error <- wire_error(value)
        row$error_kind <- if (is.null(error)) NA_character_ else error$kind
        row$error_type <- if (is.null(error)) NA_character_ else error$type
        row$error_message <- if (is.null(error)) NA_character_ else error$message
      } else if (identical(key, "source_index")) {
        row$source_index <- if (wire_null(value)) NA_real_ else as.numeric(value) + 1
      } else if (is.list(value) && !is.null(names(value))) {
        for (inner in names(value)) row[[inner]] <- bulk_file_cell(value[[inner]])
      } else if (identical(key, "source") && wire_null(value)) {
        # A read input whose reader built no source description, or a failed one: its flattened
        # columns stay NA, filled in below from the other rows.
        next
      } else if (is.list(value)) {
        is_objects <- length(value) > 0L && all(vapply(value, is.list, logical(1L)))
        row[[key]] <- list(if (is_objects || key %in% c("excluded_fields", "skipped")) {
          if (identical(key, "excluded_fields")) wire_excluded(value) else wire_objects(value)
        } else {
          wire_strings(value)
        })
      } else {
        row[[key]] <- bulk_file_cell(value)
      }
    }
    row
  })

  keys <- unique(unlist(lapply(flat, names)))
  columns <- lapply(keys, function(key) {
    cells <- lapply(flat, function(row) if (key %in% names(row)) row[[key]] else NA)
    listy <- vapply(cells, is.list, logical(1L))
    if (any(listy)) {
      return(I(lapply(cells, function(c) if (is.list(c)) c[[1L]] else if (is.na(c)) NULL else c)))
    }
    unlist(cells, use.names = FALSE)
  })
  names(columns) <- keys
  as.data.frame(columns, stringsAsFactors = FALSE, optional = TRUE)
}

# The long table of a batch, with `source_index` made 1-based.
bulk_parse_table <- function(data) {
  table <- readers_parse_records_table(data)
  if (!is.null(table) && "source_index" %in% names(table)) {
    table$source_index <- as.numeric(table$source_index) + 1
  }
  table
}

# Every top-level field of a batch envelope, in R's shape. `class` is the batch's own class, which
# precedes the shared "mzlibr_batch".
bulk_parse_batch <- function(data, verb, class) {
  out <- list(verb = verb)
  for (key in names(data)) {
    value <- data[[key]]
    if (identical(key, "columns")) {
      out$records <- bulk_parse_table(data)
    } else if (identical(key, "files")) {
      out$files <- bulk_parse_files(value)
    } else if (identical(key, "output")) {
      out$output <- readers_parse_output(data)
    } else if (wire_null(value)) {
      out[key] <- list(NA)
    } else if (is.list(value) && !is.null(names(value))) {
      out[[key]] <- value
    } else if (is.list(value)) {
      out[[key]] <- if (length(value) > 0L && all(vapply(value, is.list, logical(1L)))) {
        wire_objects(value)
      } else {
        wire_strings(value)
      }
    } else {
      out[[key]] <- value
    }
  }
  if (is.null(out$files)) {
    out$files <- bulk_parse_files(list())
  }
  if (!"records" %in% names(out)) {
    out["records"] <- list(NULL)
  }
  structure(out, class = c(class, "mzlibr_batch"))
}

# The inputs of a batch that could not be read, under on_error = "skip".
bulk_failed <- function(batch) {
  files <- batch$files
  if (!"error_kind" %in% names(files)) {
    return(files[0L, , drop = FALSE])
  }
  files[!is.na(files$error_kind), , drop = FALSE]
}

# ---------------------------------------------------------------- printing

#' Print a batch read from many files
#'
#' @param x A `_many` result, such as [readers_read_spectra_many()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_batch <- function(x, ...) {
  cat("<", class(x)[1L], "> ", x$verb, "\n", sep = "")
  count <- function(name) if (is.null(x[[name]]) || is.na(x[[name]])) NA else format(x[[name]])
  cat("  ", count("file_count"), " inputs: ", count("read_count"), " read, ",
    count("failed_count"), " failed", if (!is.null(x$on_error)) paste0(" (on_error = \"", x$on_error, "\")") else "",
    "\n",
    sep = ""
  )
  if (!is.null(x$records)) {
    cat("  ", nrow(x$records), " rows in `records`",
      if (!is.null(x$record_count)) paste0(", ", format(x$record_count), " records") else "", "\n",
      sep = ""
    )
  }
  failed <- bulk_failed(x)
  for (i in seq_len(nrow(failed))) {
    cat("  ! ", basename(failed$path[[i]]), ": ", failed$error_message[[i]], "\n", sep = "")
  }
  if (!is.null(x$output)) {
    cat("  written to ", x$output$path, " (", x$output$format, ", ",
      format(x$output$row_count), " rows)\n",
      sep = ""
    )
  }
  invisible(x)
}
