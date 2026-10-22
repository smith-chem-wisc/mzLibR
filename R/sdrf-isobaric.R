# SDRF-Proteomics as an isobaric (TMT, TMTpro, iTRAQ, DiLeu) design: can this file drive a
# quantification, channel by channel?
#
#   sdrf_isobaric_design()   mzLib's SdrfIsobaricDesign.Read: one row per file and annotated
#                            channel, MetaMorpheus's TmtDesign.txt on request, or every reason
#                            there is no design
#
# Ported from pyMzLib's sdrf.py (isobaric_design, SdrfIsobaricDesign, DesignedChannel), which
# decided the verb, the wire fields and the caveats. The rows are a data.frame, `records`.
# **Every number in it stays 1-based in R, as the SDRF writes it** - fraction, technical and
# biological replicate, plex_id - because they are the SDRF's own numbers, not row indexes (the
# spec's bindings.r note). That is the opposite of sdrf_design(), whose label-free numbers are
# 0-based: the split is mzLib's (SdrfIsobaricDesign against SpectraFileInfo), and mzLibR converts
# neither.

# The pyMzLib release whose bridge first dispatches `sdrf isobaric-design` and `sdrf restrict`.
SDRF_ISOBARIC_RESTRICT_SINCE <- "0.6.0"

# `--<name> value` for a string option, nothing for NULL; a blank value is refused.
sdrf_optional_text <- function(value, name, example) {
  if (is.null(value)) {
    return(character(0))
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(mzlib_usage_error(paste0(
      name, " must be a non-blank string or NULL, e.g. '", example, "'; got ",
      paste(deparse(value), collapse = " "), "."
    )))
  }
  if (grepl("[\r\n]", value)) {
    stop(mzlib_usage_error(paste0(name, " may not contain a newline: ", deparse(value), ".")))
  }
  c(paste0("--", gsub("_", "-", name, fixed = TRUE)), value)
}

sdrf_build_isobaric_request <- function(path, kit, plex_column, plex_file_name_pattern, single_plex,
                                        condition_columns, searched_files, out) {
  args <- c(
    "sdrf", "isobaric-design", "--path", readers_normalise_path(path),
    sdrf_optional_text(kit, "kit", "TMT10"),
    sdrf_optional_text(plex_column, "plex_column", "comment[plex]"),
    sdrf_optional_text(plex_file_name_pattern, "plex_file_name_pattern", "TMT_?pool(?<plex>\\d+)"),
    sdrf_optional_text(single_plex, "single_plex", "plex1")
  )
  if (!is.null(condition_columns)) {
    names <- sdrf_design_names(condition_columns, "condition_columns", "[\t\r\n]")
    args <- c(args, "--condition-columns", paste(names, collapse = "\t"))
  }
  stdin <- NULL
  if (!is.null(searched_files)) {
    files <- sdrf_design_names(searched_files, "searched_files", "[\r\n]")
    args <- c(args, "--searched-files-stdin")
    # No trailing newline: system2(input =) ends the last line itself.
    stdin <- paste(files, collapse = "\n")
  }
  if (!is.null(out)) {
    if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(trimws(out))) {
      stop(mzlib_usage_error("out must name a folder; TmtDesign.txt is written inside it, e.g. 'design'."))
    }
    args <- c(args, "--out", trimws(out))
  }
  list(args = args, stdin = stdin)
}

sdrf_parse_isobaric_design <- function(data) {
  text <- function(name) wire_field(data, name, "character", NA_character_)
  channels <- data[["kit_channels"]]
  ids <- data[["plex_ids"]]
  written <- data[["written"]]
  structure(
    list(
      path = text("path"),
      is_valid = isTRUE(data[["is_valid"]]),
      kit = text("kit"),
      kit_type = text("kit_type"),
      kit_channels = if (is.list(channels)) wire_strings(channels) else NULL,
      plex_column = text("plex_column"),
      plex_file_name_pattern = text("plex_file_name_pattern"),
      single_plex = text("single_plex"),
      plex_source = text("plex_source"),
      plex_ids = proteins_parse_counts(if (is.list(ids)) ids else list()),
      plex_count = as.numeric(wire_field(data, "plex_count", "numeric", NA_real_)),
      file_key_column = text("file_key_column"),
      condition_columns = wire_strings(data[["condition_columns"]]),
      condition_columns_declared = isTRUE(data[["condition_columns_declared"]]),
      searched_files_given = isTRUE(data[["searched_files_given"]]),
      file_count = as.numeric(wire_field(data, "file_count", "numeric", NA_real_)),
      record_count = as.numeric(wire_field(data, "record_count", "numeric", NA_real_)),
      refusals = wire_strings(data[["refusals"]]),
      notes = wire_strings(data[["notes"]]),
      report = text("report"),
      column_names = wire_strings(data[["column_names"]]),
      records = sdrf_design_table(data),
      written = if (is.list(written)) lapply(written, function(v) if (wire_null(v)) NA else v) else NULL,
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_sdrf_isobaric_design"
  )
}

#' Read an isobaric (TMT, TMTpro, iTRAQ, DiLeu) design out of an SDRF, or every reason it cannot be read
#'
#' Answers **can this SDRF drive an isobaric quantification, and which sample is in which channel of
#' which file?** Calls mzLib's `SdrfIsobaricDesign.Read`. It runs MetaMorpheus's own checks on a TMT
#' design - one plex, fraction and technical replicate per file; a channel described once per plex;
#' no sample twice at one fraction in a plex - and **refuses rather than repairs**, reporting every
#' failure at once.
#'
#' **It guesses nothing:** not the kit, because one kit's labels can be a subset of another's, and
#' not the plex, which you declare with exactly one of `plex_column`, `plex_file_name_pattern` and
#' `single_plex`. Giving none, or several, is mzLib's refusal.
#'
#' **The numbers are 1-based, as the SDRF writes them.** [sdrf_design()]'s are 0-based. That split
#' is mzLib's, and mzLibR converts neither.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param kit The kit the search used, named as [isobaric_kits()] takes it (`"TMT10"`, `"tmt18"`,
#'   `"TMT6-plex"`, `"TMT11 on K"`). `NULL` sends none, and mzLib then refuses: the kit is never
#'   guessed from the labels.
#' @param plex_column The column each row's plex is read from, by its exact name.
#' @param plex_file_name_pattern A .NET regular expression matched, ignoring case, against each row's
#'   file name; the plex is its group named `plex`, else its first group, e.g.
#'   `"TMT_?pool(?<plex>\\d+)"` (one backslash in the pattern, two in R source).
#' @param single_plex The name of the one plex every file belongs to.
#' @param condition_columns The `factor value[...]` columns that make up the condition, exactly as
#'   for [sdrf_design()]. `NULL` uses the document's only factor column.
#' @param searched_files The files the search will read, exactly as for [sdrf_design()]: each must
#'   be named **exactly** by some row. A search that read `X-calib.mzML` for an SDRF naming `X.raw`
#'   matches nothing; [sdrf_restrict()] the SDRF first.
#' @param out A **folder**, created when missing. mzLib writes MetaMorpheus's `TmtDesign.txt`
#'   (1-based) inside it, replacing one already there, but only when the design is valid.
#'   MetaMorpheus reads it only from the folder of the first spectra file.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_isobaric_design`. `is_valid` is `TRUE` when there are no refusals.
#'   **A refusal is an answer, not an error**: `refusals` then lists every reason at once, in mzLib's
#'   words, and `records` has no rows.
#'
#'   `kit` is the kit name exactly as you gave it, and `kit_type` mzLib's name for the kit it
#'   resolved to (`"TMT10"`, `"TMT18"`, `"iTRAQ4"`, ...), as [isobaric_kits()] names it; both `NA`
#'   when no kit was given. `kit_channels` lists every channel of the kit in reporter m/z order (a
#'   channel here and absent from `records` is one the SDRF does not annotate), `NULL` without a kit.
#'   `plex_column`, `plex_file_name_pattern` and `single_plex` echo the plex source you gave, `NA`
#'   for the others; `plex_source` is where mzLib read the plex from, as its report prints it,
#'   `"(not declared)"` when no single source could be used. `plex_ids` is a named vector, plex name
#'   to its number, 1 to N in case-insensitive name order (the order MetaMorpheus numbers plexes
#'   in), empty when refused; `plex_count` counts the plexes, 0 when refused. `file_key_column` is
#'   the column files were keyed on - `comment[searched data file]` when the SDRF has it, else
#'   `comment[data file]` - or `NA` when it has neither. `condition_columns` lists the columns the
#'   condition was built from; `condition_columns_declared` and `searched_files_given` say whether
#'   you passed them. `file_count` counts the files in the design and `record_count` the rows of
#'   `records`, one per file and annotated channel; both 0 when refused. `notes` says what was
#'   dropped or filled on the way; `report` is mzLib's account of all of it, to print. `written` is a
#'   list (`path` of the `TmtDesign.txt`, `file_count`, `record_count`), or `NULL` when `out` was not
#'   given or the design was refused. `caveats` says what the design does and does not cover.
#'
#'   `records` is a data.frame, files in SDRF row order, then each file's channels in reporter m/z
#'   order: `full_path` (the searched path you passed for the file, else the file name the SDRF
#'   gives); `plex`, as the declared source gives it; `plex_id`, the plex's number in `plex_ids`;
#'   `sample_name`, the SDRF `source name` of the sample in the channel; `channel`, as the kit spells
#'   it (`"127N"`); `condition`, the condition columns' values joined with `_`; and
#'   `biological_replicate` (as written, never ranked), `fraction` (1 when the SDRF has no fraction
#'   column) and `technical_replicate` (1 when it has no technical-replicate column), **all
#'   1-based**.
#'
#' @section Restrict first, then read the design:
#'
#' A search reads `X-calib.mzML` where the SDRF names the acquired `X.raw`, so `searched_files` of
#' the search's own names match nothing. [sdrf_restrict()] joins the two and writes the searched
#' name into `comment[searched data file]`, and this function keys on that column when a document
#' has it.
#'
#' @seealso [sdrf_restrict()], [sdrf_design()], [isobaric_kits()]
#' @spec sdrf.isobaric-design
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # PXD008841 writes no plex column; the plex survives only in the file names, spelled two ways
#' # (TMTpool1, TMT_pool5), so one pattern reads both:
#' d <- sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10",
#'   plex_file_name_pattern = "TMT_?pool(?<plex>\\\\d+)")
#' d
#' d$plex_ids
#' d$records[10, c("channel", "sample_name", "plex_id", "fraction", "biological_replicate")]
#'
#' # Declare no plex source and mzLib refuses rather than guess one:
#' refused <- sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10")
#' refused$is_valid
#' writeLines(strwrap(refused$refusals[[1]], 90))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_isobaric_design <- function(path, kit = NULL, plex_column = NULL, plex_file_name_pattern = NULL,
                                 single_plex = NULL, condition_columns = NULL, searched_files = NULL,
                                 out = NULL, timeout = 60) {
  request <- sdrf_build_isobaric_request(
    path, kit, plex_column, plex_file_name_pattern, single_plex, condition_columns, searched_files, out
  )
  bridge_require_verb("sdrf isobaric-design", SDRF_ISOBARIC_RESTRICT_SINCE)
  sdrf_parse_isobaric_design(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Print an isobaric design
#'
#' @param x An [sdrf_isobaric_design()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_isobaric_design <- function(x, ...) {
  if (isTRUE(x$is_valid)) {
    cat("<mzlibr_sdrf_isobaric_design> ", basename(x$path), ": ", format(x$file_count), " files, ",
      format(x$plex_count), " plex(es), kit ", x$kit_type, ", ", format(x$record_count),
      " channel rows (1-based numbers)\n",
      sep = ""
    )
  } else {
    cat("<mzlibr_sdrf_isobaric_design> ", basename(x$path), ": REFUSED, ", length(x$refusals),
      " reason(s)\n",
      sep = ""
    )
    for (reason in utils::head(x$refusals, 3L)) cat("  ! ", reason, "\n", sep = "")
    if (length(x$refusals) > 3L) cat("  ... see $refusals\n")
  }
  cat("  plex from: ", x$plex_source, "\n", sep = "")
  if (!is.na(x$file_key_column)) cat("  files keyed on: ", x$file_key_column, "\n", sep = "")
  if (length(x$notes) > 0L) cat("  ", length(x$notes), " note(s) - read $notes\n", sep = "")
  if (!is.null(x$written)) cat("  written: ", x$written$path, "\n", sep = "")
  invisible(x)
}
