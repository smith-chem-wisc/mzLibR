# Converting full sequences from one notation to another, with mzLib.
#
#   peptidoform_convert()   SequenceConversionService.Default (ProForma registered) .GetConverter
#                           .Convert, once per input: for example a MetaMorpheus full sequence to
#                           Unimod accessions
#
# Ported from pyMzLib's `peptidoform.convert()` (pyMzLib #75). Every output and every status is
# mzLib's: this file checks the arguments, sends the sequences one per stdin line, and projects the
# answer into a data.frame, `records`, one row per input in input order. It maps no modification
# itself, and repairs nothing mzLib returns - including the ProForma target's UniProt gap, which
# the caveats describe (mzLib#1401).

# The pyMzLib release whose bridge first dispatches `peptidoform convert` (its spec's
# since.pymzlib).
PEPTIDOFORM_CONVERT_SINCE <- "0.4.0"

# mzLib's SequenceConversionHandlingMode names, the values `mode` accepts.
PEPTIDOFORM_CONVERSION_MODES <- c(
  "ThrowException", "ReturnNull", "RemoveIncompatibleElements", "UsePrimarySequence"
)

# The table's columns, and which of them hold one character vector per row.
PEPTIDOFORM_CONVERT_COLUMNS <- c(
  "input", "output", "status", "failure_reason", "incompatible_items", "warnings", "errors"
)
PEPTIDOFORM_CONVERT_LIST_COLUMNS <- c("incompatible_items", "warnings", "errors")

# ---------------------------------------------------------------- argument shaping

# The sequences as stdin text, unchanged, refusing what would shift the rows.
#
# The bridge skips blank lines and splits on line breaks, so either would put every later result
# on the wrong row. Refused here rather than silently re-aligned.
peptidoform_convert_stdin <- function(sequences) {
  if (!is.character(sequences)) {
    stop(mzlib_usage_error(paste0(
      "sequences must be a character vector of full sequences; got ", class(sequences)[1L], "."
    )))
  }
  if (length(sequences) == 0L) {
    stop(mzlib_usage_error(
      "At least one sequence is required, e.g. \"[UniProt:N-acetylserine on S]SEQK\"."
    ))
  }
  for (i in seq_along(sequences)) {
    value <- sequences[[i]]
    if (is.na(value)) {
      stop(mzlib_usage_error(paste0("sequences[", i, "] is NA; every entry is one result row.")))
    }
    if (!nzchar(trimws(value))) {
      stop(mzlib_usage_error(paste0("sequences[", i, "] is blank; every entry is one result row.")))
    }
    if (grepl("[\r\n]", value)) {
      stop(mzlib_usage_error(paste0("sequences[", i, "] contains a line break: '", value, "'.")))
    }
  }
  # No trailing newline: system2(input =) ends the last line itself (see stats_adjust_stdin()).
  paste(unname(sequences), collapse = "\n")
}

peptidoform_convert_format <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(mzlib_usage_error(paste0(
      name, " must be a single non-empty format name; got ", paste(deparse(value), collapse = " "), "."
    )))
  }
  value
}

# mzLib's own spelling of `mode`, matched case-insensitively, as the bridge matches it.
peptidoform_convert_mode <- function(mode) {
  hit <- if (is.character(mode) && length(mode) == 1L && !is.na(mode)) {
    PEPTIDOFORM_CONVERSION_MODES[tolower(PEPTIDOFORM_CONVERSION_MODES) == tolower(mode)]
  } else {
    character(0)
  }
  if (length(hit) != 1L) {
    stop(mzlib_usage_error(paste0(
      "mode must be one of ", paste(PEPTIDOFORM_CONVERSION_MODES, collapse = ", "),
      " (mzLib's SequenceConversionHandlingMode); got ", paste(deparse(mode), collapse = " "), "."
    )))
  }
  hit
}

peptidoform_convert_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) || threads != round(threads) ||
    threads == 0 || threads < -1) {
    stop(mzlib_usage_error(paste0(
      "threads must be 1 or more, or -1 for every core; got ", paste(deparse(threads), collapse = " "), "."
    )))
  }
  formatC(threads, format = "d")
}

peptidoform_convert_request <- function(sequences, source, target, mode, threads) {
  stdin <- peptidoform_convert_stdin(sequences)
  args <- c(
    "peptidoform", "convert",
    "--from", peptidoform_convert_format(source, "source"),
    "--to", peptidoform_convert_format(target, "target"),
    "--mode", peptidoform_convert_mode(mode),
    "--threads", peptidoform_convert_threads(threads)
  )
  list(args = args, stdin = stdin)
}

# ---------------------------------------------------------------- parsing

# The wire table as a data.frame. `incompatible_items`, `warnings` and `errors` are ALWAYS list
# columns of character vectors, character(0) where mzLib recorded none, so a row with nothing to
# say keeps its place; the other columns are character, NA for null.
peptidoform_convert_table <- function(data) {
  columns <- data[["columns"]]
  order <- wire_strings(data[["column_names"]])
  if (length(order) == 0L) {
    order <- if (is.list(columns) && length(columns) > 0L) names(columns) else PEPTIDOFORM_CONVERT_COLUMNS
  }
  built <- lapply(order, function(name) {
    values <- if (is.list(columns)) columns[[name]] else NULL
    if (!is.list(values)) values <- list()
    if (name %in% PEPTIDOFORM_CONVERT_LIST_COLUMNS || any(vapply(values, is.list, logical(1L)))) {
      return(I(lapply(values, function(v) if (is.list(v)) wire_strings(v) else character(0))))
    }
    vapply(values, function(v) if (wire_null(v)) NA_character_ else as.character(v)[1L], character(1L),
      USE.NAMES = FALSE
    )
  })
  names(built) <- order
  lengths_seen <- vapply(built, length, integer(1L))
  if (length(unique(lengths_seen)) > 1L) {
    stop(mzlib_protocol_error(paste0(
      "The table's columns have different lengths (",
      paste(paste0(order, "=", lengths_seen), collapse = ", "), "), so they cannot form a table."
    )))
  }
  as.data.frame(built, stringsAsFactors = FALSE, optional = TRUE)
}

peptidoform_parse_conversions <- function(data) {
  count <- function(name) as.numeric(wire_field(data, name, "numeric", NA_real_))
  structure(
    list(
      source_format = wire_field(data, "source_format", "character", NA_character_),
      target_format = wire_field(data, "target_format", "character", NA_character_),
      mode = wire_field(data, "mode", "character", NA_character_),
      source_formats = wire_strings(data[["source_formats"]]),
      target_formats = wire_strings(data[["target_formats"]]),
      record_count = count("record_count"),
      converted_count = count("converted_count"),
      warned_count = count("warned_count"),
      failed_count = count("failed_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = peptidoform_convert_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_sequence_conversions"
  )
}

# ---------------------------------------------------------------- the public surface

#' Convert full sequences to another notation with mzLib, one row per input
#'
#' Hands each sequence to mzLib's `SequenceConversionService` (with ProForma registered) and returns
#' what mzLib made of it. The main use is MetaMorpheus or mzLib full sequences to Unimod
#' accessions: `[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]SEQK`. Nothing is converted
#' in R or in the bridge; every output and every status is mzLib's.
#'
#' **Choose Unimod, not ProForma, for UniProt-sourced modifications.** mzLib's ProForma target does
#' not resolve them and writes them back under their mzLib name, with status `"converted"`
#' (mzLib#1401). In ProForma output, treat any bracket that is not a `UNIMOD:` term as unresolved.
#' The `pro_forma` column [readers_read_records()] gives a `.psmtsv` comes from the same serializer
#' and has the same gap.
#'
#' @param sequences A character vector of full sequences, one result row each, in order, duplicates
#'   included, sent to mzLib exactly as given. An `NA`, a blank entry or one with a line break is
#'   refused, because it would shift every later row. Split an ambiguous MetaMorpheus full sequence
#'   (`|`-joined candidates) first: mzLib joins the candidates into one sequence, with only a
#'   warning.
#' @param source The notation the sequences are in: one of mzLib's registered source formats (the
#'   result's `source_formats`), matched case-insensitively. Sent as the wire's `--from`.
#' @param target The notation to write: one of mzLib's registered target formats (the result's
#'   `target_formats`), matched case-insensitively. Sent as the wire's `--to`.
#' @param mode mzLib's `SequenceConversionHandlingMode`, case-insensitive: `"ReturnNull"` (a
#'   sequence mzLib cannot convert is a `"failed"` row), `"RemoveIncompatibleElements"` or
#'   `"UsePrimarySequence"` (what the target cannot write is dropped and the row is
#'   `"converted_with_warnings"`), or `"ThrowException"` (the first such sequence, in input order,
#'   fails the whole call with an `mzlib_usage_error` naming it, and nothing is returned).
#' @param threads Sequences converted at once, or `-1` for every core. The rows are identical, in
#'   input order, at any value.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_sequence_conversions`:
#'
#'   - `records`, a data.frame with one row per input, in input order: `input`, the sequence as
#'     sent; `output`, mzLib's converted sequence, `NA` when the row failed; `status`,
#'     `"converted"` (an output and nothing recorded against it), `"converted_with_warnings"` (an
#'     output, but mzLib recorded a warning, an error or an incompatible item) or `"failed"` (no
#'     output); `failure_reason`, mzLib's `ConversionFailureReason` or `NA` when it recorded none -
#'     which a failed row can be: the Unimod target under `"ReturnNull"` names only the
#'     incompatible items; and three list columns of character vectors, `character(0)` when
#'     empty: `incompatible_items` (what the target could not write, as mzLib describes it,
#'     `"Made Up:Not a modification on K @3(K)"`), `warnings` and `errors`.
#'   - `source_format` and `target_format`, spelled as mzLib registered them; `mode`, the handling
#'     mode used; `source_formats` and `target_formats`, every notation mzLib has registered.
#'   - `record_count` sequences, of which `converted_count`, `warned_count` and `failed_count` have
#'     each status; `column_names`; and `caveats`, what a status does and does not promise.
#'
#'   The failed and warned rows are `records[records$status != "converted", ]`.
#'
#' @seealso [peptidoform_fragments()], [readers_read_results()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # UniProt and MetaMorpheus modifications to Unimod accessions, and one mzLib cannot map:
#' sequences <- c(
#'   "[UniProt:N-acetylserine on S]SEQK",
#'   "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
#'   "PEPM[Common Variable:Oxidation on M]K",
#'   "PEPK[Made Up:Not a modification on K]R"
#' )
#' result <- peptidoform_convert(sequences)
#' result
#' result$records[, c("output", "status")]
#' failed <- result$records[result$records$status != "converted", ]
#' failed$incompatible_items[[1]]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec peptidoform.convert
#' @export
peptidoform_convert <- function(sequences, source = "mzLib", target = "Unimod", mode = "ReturnNull",
                                threads = 1, timeout = 120) {
  request <- peptidoform_convert_request(sequences, source, target, mode, threads)
  bridge_require_verb("peptidoform convert", PEPTIDOFORM_CONVERT_SINCE)
  peptidoform_parse_conversions(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

# ---------------------------------------------------------------- printing

#' Print sequence conversions
#'
#' @param x A [peptidoform_convert()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_sequence_conversions <- function(x, ...) {
  cat("<mzlibr_sequence_conversions> ", format(x$record_count), " sequences, ", x$source_format,
    " to ", x$target_format, " (", x$mode, ")\n",
    sep = ""
  )
  cat("  ", format(x$converted_count), " converted, ", format(x$warned_count), " with warnings, ",
    format(x$failed_count), " failed\n",
    sep = ""
  )
  invisible(x)
}
