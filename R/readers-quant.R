# Quantification tables: MetaMorpheus protein groups, FlashLFQ peptides, and PTM site occupancy.
#
# mzLib 1.0.592 reads the tables MetaMorpheus and FlashLFQ write after quantifying (mzLib #1347).
# Each is a wide table with one column per sample, and the bridge sends it LONG - one row per record
# per sample, with the sample's values as columns - so nothing arrives as a dictionary and every
# row is a plain data.frame row. `limit` and `offset` still count records (groups, peptides), not
# rows, so a window never splits one group's samples.
#
# Ported from pyMzLib's readers.py (read_protein_groups, read_quantified_peptides, read_occupancy
# and their `_many` forms). These verbs are new in the bridge pyMzLib 0.2.0 publishes, so each
# asks the bridge whether it has the verb before calling it (bridge_require_verb()), and an older
# bridge is reported as old rather than answering "Unknown command".

# The pyMzLib release whose bridge first has these three verbs (their specs' since.pymzlib).
READERS_QUANT_VIEWS_SINCE <- "0.2.0"

readers_read_quant <- function(verb, class, path, limit, offset, out, timeout, retention_default) {
  args <- readers_build_read_args(path, limit, offset, out, verb)
  bridge_require_verb(paste("readers", verb), READERS_QUANT_VIEWS_SINCE)
  readers_parse_quant_records(bridge_invoke(args, timeout = timeout), class, retention_default)
}

readers_read_quant_many <- function(verb, paths, out, threads, on_error, timeout, fn) {
  request <- bulk_request(c("readers", verb), paths, threads, on_error, out, fn)
  bridge_require_verb(paste("readers", verb), READERS_QUANT_VIEWS_SINCE)
  bulk_parse_batch(
    bridge_invoke(request$args, stdin = request$stdin, timeout = timeout),
    paste("readers", verb), "mzlibr_read_batch"
  )
}

#' Read a MetaMorpheus protein-group table, one row per group per sample group
#'
#' MetaMorpheus writes `AllQuantifiedProteinGroups.tsv` wide: one intensity column and one
#' spectral-count column per sample group. This reads it long - one row per protein group per
#' sample group - so the per-sample values are ordinary columns and nothing arrives as a map.
#'
#' **The table is unfiltered**, as MetaMorpheus writes it: decoys, contaminants and groups above 1%
#' FDR are all rows. Filter on `q_value` and `decoy_contaminant_target` before you count anything.
#'
#' @param path Path to a MetaMorpheus `AllQuantifiedProteinGroups.tsv`
#'   ([readers_identify()] reports `MetaMorpheusQuantifiedProteinGroups`).
#' @param limit Maximum protein groups to return - groups, not rows; each gives one row per sample
#'   group. `NULL`, the default, returns all of them.
#' @param offset Protein groups to skip - groups, not rows. A window, not a cursor: mzLib parses the
#'   whole file on every call.
#' @param out Write the table to this path as tab-separated text and return only a summary.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_protein_group_records`. `record_count` counts the protein groups in the whole
#'   file; `returned_count` the groups returned, starting `offset` groups in; `row_count` the rows
#'   in `records`, which is groups times sample groups. `sample_labels` are the file's sample-group
#'   labels in header order, verbatim.
#'
#'   `records` has one row per group per sample group: `q_value` a fraction from 0 to 1;
#'   `spectral_count` in PSMs; and `intensity` in the instrument's intensity units, `NA` where the
#'   cell was blank - not quantified in that sample group, which is not zero. `gene` and
#'   `organism` are `NA` where the file has no such column, and `absent_fields` then names them.
#'   `rows_not_read` is `NA`: every line is a group, or the file fails to read.
#'
#' @section Sample labels are verbatim:
#'
#' A label such as `QE-002106_GM1_a-calib` is what MetaMorpheus wrote in the header. Condition,
#' replicate and channel cannot be recovered from it reliably, so none is parsed out: map the
#' labels to your design yourself.
#'
#' @section The occupancy columns are elsewhere:
#'
#' MetaMorpheus's per-sample PTM occupancy cells have no column shape and are listed in
#' `excluded_fields`, whose `verb` column names `readers read-occupancy`: use
#' [readers_read_occupancy()], which reads them one row per site.
#'
#' @seealso [readers_read_protein_groups_many()], [readers_read_occupancy()],
#'   [readers_read_quantified_peptides()]
#' @spec readers.read-protein-groups
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' groups <- readers_read_protein_groups("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
#'   limit = 1)
#' groups
#' head(groups$records[, c("protein_group_name", "q_value", "sample_label", "intensity")])
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
readers_read_protein_groups <- function(path, limit = NULL, offset = 0, out = NULL, timeout = NULL) {
  readers_read_quant(
    "read-protein-groups", "mzlibr_protein_group_records", path, limit, offset, out, timeout,
    NA_character_
  )
}

#' Read many protein-group tables into one long table, in one bridge call
#'
#' The many-files form of [readers_read_protein_groups()]: the whole list goes to ONE bridge
#' process, which reads `threads` files at once and returns one table, rows grouped by file in
#' input order whatever `threads` is. Never loop the one-file function over files instead.
#'
#' @param paths A character vector of paths to MetaMorpheus `AllQuantifiedProteinGroups.tsv` files.
#' @param out Write the long table here as tab-separated text, one file at a time, and return only
#'   a summary.
#' @param threads Files read at once. `1`, the default, holds one whole file in memory at a time;
#'   `-1` uses one per core. The table does not depend on it.
#' @param on_error `"fail"`, the default, stops at the first file that cannot be read, with its
#'   error. `"skip"` records the failure in that file's row of `files` and reads the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`. `records` is the long table: `source_index`, the 1-based
#'   position of the row's file in `paths`, and `source_path`, then the same columns as
#'   [readers_read_protein_groups()] - `q_value` a fraction, `spectral_count` in PSMs, `intensity`
#'   in the instrument's intensity units. `files` is a data.frame with one row per input: its
#'   `path`, `file_type`, `reader`, `record_count` in groups, `sample_labels`, `caveats`,
#'   `absent_fields`, `failed_fields`, `excluded_fields` and, for an input that could not be read,
#'   `error_kind`, `error_type` and `error_message`. `file_count`, `read_count` and `failed_count`
#'   count files; `record_count` counts groups over the files read; `row_count` counts rows.
#'
#' @seealso [readers_read_protein_groups()]
#' @spec readers.read-protein-groups bulk
#' @examples
#' \dontrun{
#' # Your own files: no recording of a many-file read of this table exists to replay.
#' batch <- readers_read_protein_groups_many(c("search_1/AllQuantifiedProteinGroups.tsv", "search_2/AllQuantifiedProteinGroups.tsv"),
#'   threads = 2, on_error = "skip")
#' batch$files[, c("path", "record_count", "error_message")]
#' }
#' @export
readers_read_protein_groups_many <- function(paths, out = NULL, threads = 1, on_error = "fail",
                                             timeout = NULL) {
  readers_read_quant_many(
    "read-protein-groups", paths, out, threads, on_error, timeout,
    "readers_read_protein_groups_many"
  )
}

#' Read a FlashLFQ peptide table, one row per peptide per sample
#'
#' Reads FlashLFQ's `QuantifiedPeptides.tsv`, or the `AllQuantifiedPeptides.tsv` MetaMorpheus
#' writes with it, long: one row per peptide per sample, with that sample's intensity and detection
#' type as columns.
#'
#' **An intensity of 0 is not a measurement.** FlashLFQ writes a literal 0 for a peptide it did not
#' quantify in a sample, and mzLib reads it as 0; `detection_type` - `"MSMS"`, `"MBR"`,
#' `"NotDetected"` and so on - is what tells the two apart.
#'
#' @param path Path to a FlashLFQ `QuantifiedPeptides.tsv` or MetaMorpheus
#'   `AllQuantifiedPeptides.tsv` ([readers_identify()] reports `FlashLFQQuantifiedPeptide`).
#' @param limit Maximum peptides to return - peptides, not rows; each gives one row per sample.
#'   `NULL`, the default, returns all of them.
#' @param offset Peptides to skip - peptides, not rows.
#' @param out Write the table to this path as tab-separated text and return only a summary.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_quantified_peptide_records`. `record_count` counts the peptides in the whole
#'   file; `returned_count` the peptides returned, starting `offset` peptides in; `row_count` the
#'   rows in `records`, peptides times samples. `sample_labels` are the file's sample labels in
#'   header order.
#'
#'   `records` has one row per peptide per sample: `intensity` in the instrument's intensity units
#'   (0 is not a measurement - see above); `detection_type`; and `retention_time` in minutes,
#'   written only by IsoTracker. For any other file `retention_time` and `peak_order` are `NA` in
#'   every row and `absent_fields` names them. `rows_not_read` is `NA`: every line is a peptide, or
#'   the file fails to read.
#'
#' @seealso [readers_read_quantified_peptides_many()], [readers_read_protein_groups()],
#'   [flashlfq_quantify()]
#' @spec readers.read-quantified-peptides
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' peptides <- readers_read_quantified_peptides("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv",
#'   limit = 1)
#' peptides$absent_fields   # this file has no retention times: they are NA, not missing rows
#' head(peptides$records[, c("sequence", "sample_label", "intensity", "detection_type")])
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
readers_read_quantified_peptides <- function(path, limit = NULL, offset = 0, out = NULL,
                                             timeout = NULL) {
  readers_read_quant(
    "read-quantified-peptides", "mzlibr_quantified_peptide_records", path, limit, offset, out,
    timeout, "minutes"
  )
}

#' Read many peptide tables into one long table, in one bridge call
#'
#' The many-files form of [readers_read_quantified_peptides()]: one bridge process for the whole
#' list, `threads` files at a time, rows grouped by file in input order.
#'
#' @param paths A character vector of paths to FlashLFQ or MetaMorpheus peptide tables.
#' @param out Write the long table here as tab-separated text and return only a summary.
#' @param threads Files read at once. `1`, the default, or `-1` for one per core. The table does
#'   not depend on it.
#' @param on_error `"fail"`, the default, or `"skip"` to record an unreadable file in `files` and
#'   read the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`: `records`, the long table, begins with `source_index` (1-based,
#'   into `paths`) and `source_path`, then the columns of [readers_read_quantified_peptides()] -
#'   `intensity` in the instrument's intensity units, `retention_time` in minutes. `files` has one
#'   row per input with its per-file facts. `file_count`, `read_count` and `failed_count` count
#'   files; `record_count` counts peptides over the files read and `returned_count` the peptides
#'   returned; `row_count` counts rows.
#'
#' @seealso [readers_read_quantified_peptides()]
#' @spec readers.read-quantified-peptides bulk
#' @examples
#' \dontrun{
#' # Your own files: no recording of a many-file read of this table exists to replay.
#' batch <- readers_read_quantified_peptides_many(c("search_1/AllQuantifiedPeptides.tsv", "search_2/AllQuantifiedPeptides.tsv"),
#'   threads = 2, on_error = "skip")
#' batch$files[, c("path", "record_count", "error_message")]
#' }
#' @export
readers_read_quantified_peptides_many <- function(paths, out = NULL, threads = 1,
                                                  on_error = "fail", timeout = NULL) {
  readers_read_quant_many(
    "read-quantified-peptides", paths, out, threads, on_error, timeout,
    "readers_read_quantified_peptides_many"
  )
}

#' Read PTM site occupancy from a MetaMorpheus protein-group table
#'
#' MetaMorpheus writes each protein group's modification-site occupancy into two cells per sample
#' group - by PSM count and by intensity - in a compact text form. This reads them as one row per
#' group, sample group, basis and modified site.
#'
#' **Read each basis for what it is exact in.** On `basis == "count"` rows `numerator` and
#' `denominator` are exact PSM counts and `fraction` is rounded; on `basis == "intensity"` rows
#' `fraction` is exact and the two intensities are rounded.
#'
#' @param path Path to a MetaMorpheus `AllQuantifiedProteinGroups.tsv`.
#' @param limit Maximum protein groups to read - groups, not rows; a group with no modified sites
#'   gives no rows. `NULL`, the default, reads all of them.
#' @param offset Protein groups to skip - groups, not rows.
#' @param out Write the table to this path as tab-separated text and return only a summary.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_occupancy_records`. `record_count` counts the protein groups in the whole
#'   file; `returned_count` the groups read, starting `offset` groups in; `row_count` the rows in
#'   `records`. `truncated_cell_count` counts occupancy cells MetaMorpheus cut short, in cells.
#'
#'   `records` has one row per group, sample group, basis and site: `fraction` from 0 to 1, and
#'   `numerator` and `denominator` in PSMs on count rows or in intensity on intensity rows.
#'   `position` is the 1-based residue; `cell_is_truncated` marks a row from a cut cell.
#'   `rows_not_read` is `NA`: rows are not counted for this table.
#'
#' @section entity_index is not an accession index:
#'
#' An entity with no modified sites is dropped from the cell, so `entity_index` cannot be zipped
#' with the accessions in `protein_group_name`.
#'
#' @seealso [readers_read_occupancy_many()], [readers_read_protein_groups()]
#' @spec readers.read-occupancy
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' sites <- readers_read_occupancy("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
#' sites
#' counted <- sites$records[sites$records$basis == "count", ]
#' head(counted[, c("protein_group_name", "position", "modification", "numerator", "denominator")])
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
readers_read_occupancy <- function(path, limit = NULL, offset = 0, out = NULL, timeout = NULL) {
  readers_read_quant(
    "read-occupancy", "mzlibr_occupancy_records", path, limit, offset, out, timeout, NA_character_
  )
}

#' Read the PTM site occupancy of many protein-group tables, in one bridge call
#'
#' The many-files form of [readers_read_occupancy()]: one bridge process for the whole list,
#' `threads` files at a time, rows grouped by file in input order.
#'
#' @param paths A character vector of paths to MetaMorpheus `AllQuantifiedProteinGroups.tsv` files.
#' @param out Write the long table here as tab-separated text and return only a summary.
#' @param threads Files read at once. `1`, the default, or `-1` for one per core. The table does
#'   not depend on it.
#' @param on_error `"fail"`, the default, or `"skip"` to record an unreadable file in `files` and
#'   read the rest.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_read_batch`: `records` begins with `source_index` (1-based, into `paths`)
#'   and `source_path`, then the columns of [readers_read_occupancy()] - `fraction` from 0 to 1,
#'   `numerator` and `denominator` in PSMs or intensity by basis. `files` has one row per input,
#'   including its `truncated_cell_count`. `file_count`, `read_count` and `failed_count` count
#'   files; `record_count` counts groups over the files read and `returned_count` the groups
#'   returned; `row_count` counts rows.
#'
#' @seealso [readers_read_occupancy()]
#' @spec readers.read-occupancy bulk
#' @examples
#' \dontrun{
#' # Your own files: no recording of a many-file read of this table exists to replay.
#' batch <- readers_read_occupancy_many(c("search_1/AllQuantifiedProteinGroups.tsv", "search_2/AllQuantifiedProteinGroups.tsv"),
#'   threads = 2, on_error = "skip")
#' batch$files[, c("path", "record_count", "error_message")]
#' }
#' @export
readers_read_occupancy_many <- function(paths, out = NULL, threads = 1, on_error = "fail",
                                        timeout = NULL) {
  readers_read_quant_many(
    "read-occupancy", paths, out, threads, on_error, timeout, "readers_read_occupancy_many"
  )
}

# ---------------------------------------------------------------- print methods

readers_print_quant <- function(x, what) {
  cat("<", class(x)[1L], "> ", basename(x$path), " (", x$file_type, ")\n", sep = "")
  cat("  ", length(x$sample_labels), " samples; ", format(x$row_count), " rows for ",
    format(x$returned_count), " ", what, "\n",
    sep = ""
  )
  readers_print_body(x, unit_label = if (!is.na(x$retention_time_unit)) x$retention_time_unit)
}

#' Print a protein-group table
#'
#' @param x A [readers_read_protein_groups()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_protein_group_records <- function(x, ...) {
  readers_print_quant(x, "protein groups")
}

#' Print a quantified-peptide table
#'
#' @param x A [readers_read_quantified_peptides()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_quantified_peptide_records <- function(x, ...) {
  readers_print_quant(x, "peptides")
}

#' Print a site-occupancy table
#'
#' @param x A [readers_read_occupancy()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_occupancy_records <- function(x, ...) {
  readers_print_quant(x, "protein groups")
  if (isTRUE(x$truncated_cell_count > 0)) {
    cat("  ! ", format(x$truncated_cell_count), " occupancy cell(s) were cut short by the writer\n", sep = "")
  }
  invisible(x)
}
