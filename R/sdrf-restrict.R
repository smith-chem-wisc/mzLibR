# SDRF-Proteomics restricted to the files one search read: which rows did my search read?
#
#   sdrf_restrict()             mzLib's SdrfSearchScope.Restrict: the SDRF cut to a search's files,
#                               joining X-calib.mzML to the row naming X.raw, and what could not
#                               be placed
#   sdrf_searched_data_files()  the searched file each row was joined to
#
# Ported from pyMzLib's sdrf.py (restrict, RestrictedSdrf), which decided the verb, the wire fields
# and the caveats. The restricted document is an `mzlibr_sdrf`, so every accessor of R/sdrf.R works
# on it, and its searched files are a character vector NAMED by the acquired file where you know it
# (the spec's bindings.r note), as sdrf_pool() names its documents by label.

# The column sdrf_restrict() writes the searched file's name into.
SDRF_SEARCHED_DATA_FILE_COLUMN <- "comment[searched data file]"

# The searched files as stdin lines, `searched` or `searched<TAB>acquired`.
sdrf_restrict_lines <- function(searched_files) {
  if (!is.character(searched_files)) {
    stop(mzlib_usage_error(paste0(
      "searched_files must be a character vector of paths, named by the acquired file where you ",
      "know it, e.g. c(\"X-calib.mzML\", \"X.raw\" = \"sample_A.mzML\"); got ", class(searched_files)[1L], "."
    )))
  }
  if (length(searched_files) == 0L) {
    stop(mzlib_usage_error("sdrf_restrict() needs at least one searched file."))
  }
  searched <- trimws(unname(searched_files))
  acquired <- names(searched_files)
  if (is.null(acquired)) acquired <- rep("", length(searched))
  acquired <- ifelse(is.na(acquired), "", trimws(acquired))
  lines <- character(length(searched))
  for (i in seq_along(searched)) {
    if (is.na(searched[[i]]) || !nzchar(searched[[i]])) {
      stop(mzlib_usage_error(paste0("Every searched file must be non-blank; element ", i, " is blank or NA.")))
    }
    for (text in c(searched[[i]], acquired[[i]])) {
      if (grepl("[\t\r\n]", text)) {
        stop(mzlib_usage_error(paste0("A searched or acquired file contains a tab or newline: ", deparse(text), ".")))
      }
    }
    lines[[i]] <- if (nzchar(acquired[[i]])) paste0(searched[[i]], "\t", acquired[[i]]) else searched[[i]]
  }
  lines
}

sdrf_build_restrict_request <- function(path, searched_files, out, limit, offset) {
  lines <- sdrf_restrict_lines(searched_files)
  args <- c("sdrf", "restrict", "--path", readers_normalise_path(path), sdrf_window_args(limit, offset))
  if (!is.null(out)) {
    if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(trimws(out))) {
      stop(mzlib_usage_error("out must be a file path ending in .tsv, e.g. 'searched.sdrf.tsv'."))
    }
    args <- c(args, "--out", trimws(out))
  }
  # No trailing newline: system2(input =) ends the last line itself.
  list(args = args, stdin = paste(lines, collapse = "\n"))
}

sdrf_parse_restricted <- function(data) {
  written <- data[["written"]]
  structure(
    c(
      sdrf_parse_fields(data),
      list(
        searched_files = wire_strings(data[["searched_files"]]),
        acquired_name_count = as.numeric(wire_field(data, "acquired_name_count", "numeric", 0)),
        input_row_count = as.numeric(wire_field(data, "input_row_count", "numeric", NA_real_)),
        searched_without_row = wire_strings(data[["searched_without_row"]]),
        dropped_data_files = wire_strings(data[["dropped_data_files"]]),
        ambiguous = wire_strings(data[["ambiguous"]]),
        shared_stem = wire_strings(data[["shared_stem"]]),
        written = if (is.list(written)) {
          list(
            path = as.character(wire_field(written, "path", "character", NA_character_)),
            row_count = as.numeric(wire_field(written, "row_count", "numeric", NA_real_))
          )
        } else {
          NULL
        }
      )
    ),
    class = c("mzlibr_restricted_sdrf", "mzlibr_sdrf")
  )
}

#' Restrict an SDRF to the files one search read, and report what could not be placed
#'
#' Answers **which rows of this SDRF did my search read?** A search reads `X-calib.mzML` while the
#' SDRF names the acquired `X.raw`. Calls mzLib's `SdrfSearchScope.Restrict`, which joins them **by
#' stem**, ignoring case, after removing the derivative suffixes MetaMorpheus writes (`-calib`,
#' `-averaged`, in any order and number). Each kept row is carried cell for cell, with the searched
#' file's name in a new `comment[searched data file]` column directly after `comment[data file]`;
#' nothing is renumbered.
#'
#' **Restrict first, then read a design.** [sdrf_design()] and [sdrf_isobaric_design()] match
#' `searched_files` by exact name, and key on `comment[searched data file]` when a document has it,
#' so a restricted document joins the search's own file names.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param searched_files The files the search read, as paths or bare names, a character vector.
#'   **Name an element by the acquired file it came from** when you know it -
#'   `c("X-calib.mzML", "X.raw" = "sample_A.mzML")` - and that overrides the suffix rule for it;
#'   leave the others unnamed. Sent on stdin, one line each.
#' @param out Also write the **whole** restricted document here as SDRF, whatever `limit` and
#'   `offset` are. Must end in `.tsv`; name it `.sdrf.tsv`, because mzLib adds that to a name that
#'   does not end in it. `written$path` names the file written.
#' @param limit Maximum rows to return. `NULL`, the default, returns all of them.
#' @param offset Rows to skip.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_restricted_sdrf`, which is also an `mzlibr_sdrf`, so [sdrf_value()],
#'   [sdrf_all()] and [sdrf_records()] work on it. `path` is the SDRF restricted, as given;
#'   `columns` the restricted header, which may repeat names; `rows` one character vector of raw
#'   cells per returned row, in the input's row order - a kept row shorter than the header is padded
#'   with empty cells, so unlike [sdrf_read()] the rows are never ragged. `row_count` counts the rows
#'   kept, in the **whole** restricted document, whatever `limit` and `offset`; `returned_count` the
#'   rows returned, starting `offset` rows in; `truncated` whether rows were left behind;
#'   `input_row_count` the rows of the SDRF before the restriction.
#'
#'   `searched_files` lists the searched files as sent, in order, and `acquired_name_count` counts
#'   the files you named an acquired file for. What the restriction could not place:
#'   `searched_without_row`, the searched files (file names) the SDRF has no row for - reported,
#'   never invented; `dropped_data_files`, the acquired files the SDRF describes that this search did
#'   not read (a row naming no file is dropped too, and not listed); `ambiguous`, acquisitions two or
#'   more searched files claimed, as `"stem: searched, searched"` - the first in name order kept the
#'   rows, and the others have none and are not in `searched_without_row`; and `shared_stem`, searched
#'   files whose rows name two or more **different** acquired files (one run listed as `X.raw` and
#'   `X.mzML`), as `"searched: acquired, acquired"`, every such row kept. `written` is a list (`path`,
#'   the file mzLib wrote, and `row_count`, the rows of the whole document), or `NULL` without `out`.
#'   `caveats` says how the join works and what it cannot tell you.
#'
#' @section Paths, on any platform:
#'
#' The searched files may be full paths; the join uses the file name. A path written with the other
#' platform's separator (`D:\search\x.mzML` on Linux) is reduced here too, but the design readers
#' split only the running platform's separator, so give [sdrf_design()] and [sdrf_isobaric_design()]
#' paths as this platform writes them (`file.path()` does).
#'
#' @seealso [sdrf_searched_data_files()], [sdrf_isobaric_design()], [sdrf_design()]
#' @spec sdrf.restrict
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # Two calibrated runs, the raw file of one of them as well, a renamed file whose acquisition you
#' # know, and one the SDRF does not describe:
#' pool1 <- "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
#' pool5 <- "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
#' searched <- c(paste0("D:\\\\search\\\\", pool1, "1-calib.mzML"), paste0(pool1, "2-calib.mzML"),
#'   paste0(pool1, "2.raw"), "sample_A.mzML", "not_in_sdrf.mzML")
#' names(searched) <- c("", "", "", paste0(pool5, "1.raw"), "")
#' scoped <- sdrf_restrict("PXD008841.subset.sdrf.tsv", searched, limit = 2)
#' scoped
#' c(input = scoped$input_row_count, kept = scoped$row_count, returned = scoped$returned_count)
#' sdrf_searched_data_files(scoped)
#' scoped$searched_without_row
#' scoped$dropped_data_files
#' writeLines(strwrap(scoped$ambiguous, 80))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_restrict <- function(path, searched_files, out = NULL, limit = NULL, offset = 0, timeout = 60) {
  request <- sdrf_build_restrict_request(path, searched_files, out, limit, offset)
  bridge_require_verb("sdrf restrict", SDRF_ISOBARIC_RESTRICT_SINCE)
  sdrf_parse_restricted(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' The searched file each row of a restricted SDRF was joined to
#'
#' The `comment[searched data file]` cell of each row: the searched file's name, with no folder.
#'
#' @param doc An [sdrf_restrict()] result, or an SDRF read back from one it wrote.
#' @return A character vector, one element per row; `NA` for a row the document does not reach, or
#'   every row when the document has no such column.
#' @seealso [sdrf_restrict()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' pool1 <- "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
#' pool5 <- "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
#' searched <- c(paste0("D:\\\\search\\\\", pool1, "1-calib.mzML"), paste0(pool1, "2-calib.mzML"),
#'   paste0(pool1, "2.raw"), "sample_A.mzML", "not_in_sdrf.mzML")
#' names(searched) <- c("", "", "", paste0(pool5, "1.raw"), "")
#' sdrf_searched_data_files(sdrf_restrict("PXD008841.subset.sdrf.tsv", searched, limit = 2))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_searched_data_files <- function(doc) {
  sdrf_value(doc, SDRF_SEARCHED_DATA_FILE_COLUMN)
}

#' Print a restricted SDRF
#'
#' @param x An [sdrf_restrict()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_restricted_sdrf <- function(x, ...) {
  cat("<mzlibr_restricted_sdrf> ", basename(x$path), ": ", format(x$row_count), " of ",
    format(x$input_row_count), " rows kept for ", length(x$searched_files), " searched file(s)\n",
    sep = ""
  )
  sdrf_print_body(x)
  report <- function(label, values) {
    if (length(values) > 0L) cat("  ", label, ": ", paste(values, collapse = "; "), "\n", sep = "")
  }
  report("searched without a row", x$searched_without_row)
  report("acquired but not searched", x$dropped_data_files)
  report("ambiguous", x$ambiguous)
  report("shared stem", x$shared_stem)
  if (!is.null(x$written)) {
    cat("  written to ", x$written$path, " (", format(x$written$row_count), " rows)\n", sep = "")
  }
  invisible(x)
}
