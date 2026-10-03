# SDRF-Proteomics as a label-free experimental design: can this file drive a quantification?
#
#   sdrf_design()              mzLib's SdrfLabelFreeDesign.Read: the design MetaMorpheus and
#                              FlashLFQ take, or every reason it was refused
#   sdrf_design_spectra()      that design as flashlfq_quantify() takes `spectra`
#   sdrf_design_run_design()   that design as flashlfq_median_polish() takes `design`
#
# Ported from pyMzLib's sdrf.py (design, SdrfDesign.spectra, SdrfDesign.run_design). The runs are
# a data.frame, `records`. Unlike every other position in this module, `biological_replicate`,
# `technical_replicate` and `fraction` stay 0-based in R: they are design coordinates passed
# straight to the quant functions, which take them 0-based, not row indexes (the spec's
# bindings.r note). The two helpers only rename and select columns; the design is mzLib's.

# The pyMzLib release whose bridge first dispatches `sdrf design`.
SDRF_DESIGN_SINCE <- "0.3.0"

# The design coordinates, in the order the quant functions take them.
SDRF_DESIGN_FIELDS <- c("condition", "biological_replicate", "technical_replicate", "fraction")

sdrf_design_names <- function(values, what, separator) {
  if (!is.character(values) || length(values) == 0L) {
    stop(mzlib_usage_error(paste0(
      what, " must be a non-empty character vector, or NULL; got ",
      paste(deparse(values), collapse = " "), "."
    )))
  }
  for (value in values) {
    if (is.na(value) || !nzchar(trimws(value))) {
      stop(mzlib_usage_error(paste0("Every entry of ", what, " must be non-blank; got a blank or NA.")))
    }
    if (grepl(separator, value)) {
      stop(mzlib_usage_error(paste0("An entry of ", what, " contains a tab or newline: '", value, "'.")))
    }
  }
  unname(values)
}

sdrf_build_design_request <- function(path, condition_columns, searched_files, out) {
  args <- c("sdrf", "design", "--path", readers_normalise_path(path))
  if (!is.null(condition_columns)) {
    names <- sdrf_design_names(condition_columns, "condition_columns", "[\t\r\n]")
    args <- c(args, "--condition-columns", paste(names, collapse = "\t"))
  }
  stdin <- NULL
  if (!is.null(searched_files)) {
    files <- sdrf_design_names(searched_files, "searched_files", "[\r\n]")
    args <- c(args, "--searched-files-stdin")
    stdin <- paste0(paste(files, collapse = "\n"), "\n")
  }
  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  list(args = args, stdin = stdin)
}

# The runs table. A refused design sends every column empty, which must stay a zero-row frame
# with its columns, not one row of NA.
sdrf_design_table <- function(data) {
  columns <- data[["columns"]]
  empty <- !is.list(columns) || all(vapply(columns, function(v) !is.list(v) || length(v) == 0L, logical(1L)))
  if (empty) {
    names <- wire_strings(data[["column_names"]])
    return(as.data.frame(structure(rep(list(logical(0)), length(names)), names = names),
      stringsAsFactors = FALSE, optional = TRUE
    ))
  }
  sdrf_analysis_table(data)
}

sdrf_parse_design <- function(data) {
  records <- sdrf_design_table(data)
  file_key <- data[["file_key_column"]]
  written <- data[["written"]]
  structure(
    list(
      path = wire_field(data, "path", "character", NA_character_),
      is_valid = isTRUE(data[["is_valid"]]),
      file_key_column = if (wire_null(file_key)) NA_character_ else as.character(file_key),
      condition_columns = wire_strings(data[["condition_columns"]]),
      condition_columns_declared = isTRUE(data[["condition_columns_declared"]]),
      searched_files_given = isTRUE(data[["searched_files_given"]]),
      file_count = as.numeric(wire_field(data, "file_count", "numeric", NA_real_)),
      refusals = wire_strings(data[["refusals"]]),
      notes = wire_strings(data[["notes"]]),
      report = wire_field(data, "report", "character", NA_character_),
      column_names = wire_strings(data[["column_names"]]),
      records = records,
      written = if (is.list(written)) lapply(written, function(v) if (wire_null(v)) NA else v) else NULL,
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_sdrf_design"
  )
}

#' Read a label-free experimental design out of an SDRF, or every reason it cannot be read
#'
#' Calls mzLib's `SdrfLabelFreeDesign.Read`. It runs every check MetaMorpheus's own design validator
#' runs and **refuses rather than repairs**: a factor written `not available`, two files claiming the
#' same replicate and fraction, a document with several factor columns and no declaration - each is
#' a refusal, and all of them are reported at once. It relabels in only two ways, both recorded in
#' `notes`: study-wide biological replicate numbers are ranked within each condition (`22 -> 1`),
#' and rows for files the search does not read are dropped.
#'
#' **A refusal is an answer, not an error.** When `is_valid` is `FALSE`, `refusals` lists every
#' reason and `records` has no rows, so a design MetaMorpheus would reject cannot be used by
#' mistake. mzLib refuses because MetaMorpheus skips quantification with only a warning when its
#' design file is invalid.
#'
#' @param path Path to a `.sdrf.tsv` file.
#' @param condition_columns The `factor value[...]` columns that make up the condition, by exact
#'   (case-sensitive) name, in the order to join them with `_`. `NULL` uses the document's only
#'   factor column; a document with several is then refused, because joining all of them could
#'   split a condition on a nuisance factor.
#' @param searched_files The files the search will read, as paths or bare names. Each must be named
#'   **exactly** (case and extension) by one SDRF row; rows for other files are dropped and reported
#'   in `notes`, and the paths you give become `full_path`. `NULL` takes the SDRF's own file names.
#' @param out Also write MetaMorpheus's `ExperimentalDesign.tsv` (1-based) here. Must end in `.tsv`.
#'   Written only when the design is valid. MetaMorpheus finds the file only when it is named
#'   `ExperimentalDesign.tsv` and sits beside the spectra.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sdrf_design`: `path`; `is_valid`; `file_key_column`, the column runs were keyed
#'   on (`NA` when the SDRF names no file, itself a refusal); `condition_columns` used, in join
#'   order; `condition_columns_declared` and `searched_files_given`; `file_count` files (runs) in the
#'   design, 0 when refused; `refusals` and `notes` in mzLib's words; `report`, mzLib's account of
#'   all of it; `written` (`path`, `file_count`), or `NULL`; `caveats`; and `records`, one row per
#'   run in SDRF row order: `full_path`, `file_name` (the run key, `Intensity_<file_name>`),
#'   `condition`, and `biological_replicate`, `technical_replicate` and `fraction`, **0-based**, as
#'   the quant functions take them. `ExperimentalDesign.tsv` writes the same numbers plus one.
#'
#' @section From SDRF to FlashLFQ:
#'
#' [sdrf_design_spectra()] gives the design as [flashlfq_quantify()] takes `spectra`, and
#' [sdrf_design_run_design()] as [flashlfq_median_polish()] takes `design`, so the SDRF drives the
#' quantification with no hand-written design. Both refuse a refused design.
#'
#' Label-free only: an isobaric (TMT, iTRAQ) SDRF needs a channel design this does not produce. See
#' [isobaric_kits()] for the channels themselves.
#'
#' @seealso [sdrf_design_spectra()], [sdrf_validate()], [sdrf_samples()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' both <- c("factor value[genotype]", "factor value[treatment]")
#' d <- sdrf_design("PXD067622.sdrf.tsv", condition_columns = both)
#' d
#' head(d$records[, c("file_name", "condition", "biological_replicate")], 4)
#'
#' # A refusal lists every reason at once - here one per run, as the treatment is "not available":
#' refused <- sdrf_design("PXD049018.sdrf.tsv", condition_columns = both)
#' refused$is_valid
#' length(refused$refusals)
#' writeLines(refused$refusals[[1]])
#'
#' # Study-wide replicate numbers are ranked within each condition, and the mapping is kept:
#' ranked <- sdrf_design("PXD067622_studywide.sdrf.tsv", condition_columns = both)
#' writeLines(ranked$notes[[2]])
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec sdrf.design
#' @export
sdrf_design <- function(path, condition_columns = NULL, searched_files = NULL, out = NULL,
                        timeout = 60) {
  request <- sdrf_build_design_request(path, condition_columns, searched_files, out)
  bridge_require_verb("sdrf design", SDRF_DESIGN_SINCE)
  sdrf_parse_design(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

sdrf_design_require_valid <- function(design, what) {
  if (!inherits(design, "mzlibr_sdrf_design")) {
    stop(mzlib_usage_error("design must be an sdrf_design() result."))
  }
  if (!isTRUE(design$is_valid)) {
    stop(mzlib_usage_error(paste0(
      "The design was refused, so there is no ", what, " to give. mzLib's reasons:\n  ",
      paste(design$refusals, collapse = "\n  ")
    )))
  }
}

#' A label-free design as `flashlfq_quantify()` takes its `spectra`
#'
#' One row per run: `path` (the design's `full_path`) and the four design fields, 0-based - so the
#' SDRF drives FlashLFQ with no hand-written design. Pass `searched_files` to [sdrf_design()] when
#' the SDRF names files without the directory you keep them in.
#'
#' @param design An [sdrf_design()] result.
#' @return A data.frame with `path`, `condition`, `biological_replicate`, `technical_replicate` and
#'   `fraction`, in SDRF row order.
#' @seealso [sdrf_design_run_design()], [flashlfq_quantify()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' d <- sdrf_design("PXD067622.sdrf.tsv",
#'   condition_columns = c("factor value[genotype]", "factor value[treatment]"))
#' head(sdrf_design_spectra(d), 3)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_design_spectra <- function(design) {
  sdrf_design_require_valid(design, "spectra table")
  out <- design$records[, c("full_path", SDRF_DESIGN_FIELDS), drop = FALSE]
  names(out)[1L] <- "path"
  rownames(out) <- NULL
  out
}

#' A label-free design as `flashlfq_median_polish()` takes its `design`
#'
#' One row per run, keyed by `file_name` - the `Intensity_<file_name>` column of a
#' `QuantifiedPeptides.tsv` - with the four design fields, 0-based.
#'
#' @param design An [sdrf_design()] result.
#' @return A data.frame with `file_name`, `condition`, `biological_replicate`,
#'   `technical_replicate` and `fraction`, in SDRF row order.
#' @seealso [sdrf_design_spectra()], [flashlfq_median_polish()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' d <- sdrf_design("PXD067622.sdrf.tsv",
#'   condition_columns = c("factor value[genotype]", "factor value[treatment]"))
#' head(sdrf_design_run_design(d), 3)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
sdrf_design_run_design <- function(design) {
  sdrf_design_require_valid(design, "run design")
  out <- design$records[, c("file_name", SDRF_DESIGN_FIELDS), drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Print a label-free design
#'
#' @param x An [sdrf_design()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sdrf_design <- function(x, ...) {
  if (isTRUE(x$is_valid)) {
    cat("<mzlibr_sdrf_design> ", basename(x$path), ": ", format(x$file_count), " runs in ",
      length(unique(x$records$condition)), " condition(s)\n",
      sep = ""
    )
  } else {
    cat("<mzlibr_sdrf_design> ", basename(x$path), ": REFUSED, ", length(x$refusals),
      " reason(s)\n",
      sep = ""
    )
    for (reason in utils::head(x$refusals, 3L)) cat("  ! ", reason, "\n", sep = "")
    if (length(x$refusals) > 3L) cat("  ... see $refusals\n")
  }
  if (!is.na(x$file_key_column)) cat("  files keyed on: ", x$file_key_column, "\n", sep = "")
  if (length(x$condition_columns) > 0L) {
    cat("  condition from: ", paste(x$condition_columns, collapse = " + "), "\n", sep = "")
  }
  if (length(x$notes) > 0L) cat("  ", length(x$notes), " note(s) - read $notes\n", sep = "")
  if (!is.null(x$written)) cat("  written: ", x$written$path, "\n", sep = "")
  invisible(x)
}
