# Gene Ontology for protein groups: what each group does, every member kept.
#
#   proteins_annotate_go()  mzLib's GoGroupAnnotator over a stored MetaMorpheus protein-group table:
#                           one row per (group, GO term) that ANY member holds, directly or through
#                           an ancestor, naming the members that carry it
#   proteins_update_go()    the one function that downloads a go.obo, on purpose
#
# Ported from pyMzLib's `proteins.py` (annotate_go, update_go). What changes is only the projection
# into R: the annotation table is a data.frame, `annotations`; a cell holding several accessions or
# evidence codes is a list column of character vectors, empty ones included; `evidence_by_member`
# is a list column of named lists; and the provenance blocks are plain lists. Nothing is filtered
# or collapsed here: the consensus and direct-only views are filters the caller applies to rows.

# Every `annotation_status` mzLib writes (GoAnnotationTsv.StatusName). A row with a GO term is
# always "annotated"; the other three mark a group's single term-less row and say why.
PROTEINS_ANNOTATION_STATUSES <- c("annotated", "no_go_terms", "no_entry", "contaminant")

# The pyMzLib release whose bridge first dispatches `proteins annotate-go` and `update-go`.
PROTEINS_GO_SINCE <- "0.3.0"

# The annotation table's columns whose cells are arrays of strings, from the verb's spec.
PROTEINS_GO_LIST_COLUMNS <- c(
  "accession_used", "accession_direct", "accession_inherited", "evidence", "entrapment_members"
)

# ---------------------------------------------------------------- argument shaping

# One path, as clean text; NULL when it is optional and absent.
proteins_go_path <- function(value, name, required = FALSE) {
  if (is.null(value)) {
    if (required) {
      stop(mzlib_usage_error(paste0(name, " is required.")))
    }
    return(NULL)
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(mzlib_usage_error(paste0(name, " must be a single non-empty path.")))
  }
  trimws(value)
}

# A row count for --limit / --offset: a whole number, zero or more.
proteins_go_count <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value != round(value) ||
    value < 0) {
    stop(mzlib_usage_error(paste0(
      name, " must be a whole number, zero or greater; got ", paste(deparse(value), collapse = " "), "."
    )))
  }
  formatC(value, format = "d")
}

proteins_build_go_request <- function(groups, database, go_obo, category_map, skip_unknown_go_ids,
                                      out, categories_out, limit, offset) {
  args <- c(
    "proteins", "annotate-go",
    "--groups", proteins_go_path(groups, "groups", required = TRUE),
    "--database", proteins_go_path(database, "database", required = TRUE),
    "--go-obo", proteins_go_path(go_obo, "go_obo", required = TRUE)
  )
  optional <- list(
    "--category-map" = proteins_go_path(category_map, "category_map"),
    "--out" = proteins_go_path(out, "out"),
    "--categories-out" = proteins_go_path(categories_out, "categories_out")
  )
  for (option in names(optional)) {
    if (!is.null(optional[[option]])) args <- c(args, option, optional[[option]])
  }
  if (!is.logical(skip_unknown_go_ids) || length(skip_unknown_go_ids) != 1L ||
    is.na(skip_unknown_go_ids)) {
    stop(mzlib_usage_error("skip_unknown_go_ids must be TRUE or FALSE."))
  }
  if (isTRUE(skip_unknown_go_ids)) {
    args <- c(args, "--skip-unknown-go-ids")
  }
  if (!is.null(limit)) {
    args <- c(args, "--limit", proteins_go_count(limit, "limit"))
  }
  c(args, "--offset", proteins_go_count(offset, "offset"))
}

# ---------------------------------------------------------------- parsing

# The annotation table, typed by the spec rather than guessed from the first rows: an array column
# is always a list column of character vectors, so a term-less row's `accession_used` is
# character(0), and `evidence_by_member` is always a list column of named lists.
proteins_parse_go_table <- function(block) {
  order <- wire_strings(block[["column_names"]])
  columns <- block[["columns"]]
  if (!is.list(columns)) columns <- list()
  built <- lapply(order, function(name) {
    values <- columns[[name]]
    if (!is.list(values)) values <- list()
    if (name %in% PROTEINS_GO_LIST_COLUMNS) {
      return(I(lapply(values, wire_strings)))
    }
    if (identical(name, "evidence_by_member")) {
      return(I(lapply(values, function(members) {
        if (!is.list(members) || length(members) == 0L) {
          return(structure(list(), names = character(0)))
        }
        lapply(members, wire_strings)
      })))
    }
    if (length(values) == 0L) {
      return(logical(0))
    }
    unlist(lapply(values, function(v) if (wire_null(v)) NA else v), use.names = FALSE)
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

# A flat provenance object, null fields as NA and `caveats` as a character vector.
proteins_parse_go_object <- function(value) {
  if (!is.list(value) || is.null(names(value))) {
    return(NULL)
  }
  out <- lapply(value, function(v) {
    if (wire_null(v)) NA else if (is.list(v)) wire_strings(v) else v
  })
  out
}

proteins_parse_go_categories <- function(value) {
  if (!is.list(value) || is.null(names(value))) {
    return(NULL)
  }
  list(
    map_name = wire_field(value, "map_name", "character", NA_character_),
    map_version = wire_field(value, "map_version", "character", NA_character_),
    source_file_name = wire_field(value, "source_file_name", "character", NA_character_),
    sha256 = wire_field(value, "sha256", "character", NA_character_),
    anchor_count = proteins_count(value, "anchor_count"),
    row_count = proteins_count(value, "row_count"),
    column_names = wire_strings(value[["column_names"]]),
    records = proteins_parse_go_table(value)
  )
}

proteins_parse_go_annotations <- function(data) {
  header <- data[["header"]]
  header <- if (is.list(header) && length(header) > 0L) {
    vapply(header, function(v) if (wire_null(v)) NA_character_ else as.character(v), character(1L))
  } else {
    structure(character(0), names = character(0))
  }
  structure(
    list(
      groups_file = wire_field(data, "groups_file", "character", NA_character_),
      groups_file_sha256 = wire_field(data, "groups_file_sha256", "character", NA_character_),
      table_row_count = proteins_count(data, "table_row_count"),
      decoy_group_count = proteins_count(data, "decoy_group_count"),
      group_count = proteins_count(data, "group_count"),
      annotation_database = proteins_parse_go_object(data[["annotation_database"]]),
      go = proteins_parse_go_object(data[["go"]]),
      skip_unknown_go_ids = isTRUE(data[["skip_unknown_go_ids"]]),
      unresolved_go_ids = wire_strings(data[["unresolved_go_ids"]]),
      header = header,
      row_count = proteins_count(data, "row_count"),
      returned_count = proteins_count(data, "returned_count"),
      offset = proteins_count(data, "offset"),
      truncated = isTRUE(data[["truncated"]]),
      caveats = wire_strings(data[["caveats"]]),
      written = proteins_parse_go_object(data[["written"]]),
      categories = proteins_parse_go_categories(data[["categories"]]),
      categories_written = proteins_parse_go_object(data[["categories_written"]]),
      column_names = wire_strings(data[["column_names"]]),
      annotations = proteins_parse_go_table(data)
    ),
    class = "mzlibr_go_annotations"
  )
}

proteins_parse_go_update <- function(data) {
  previous <- data[["previous_sha256"]]
  structure(
    list(
      go_obo_file = wire_field(data, "go_obo_file", "character", NA_character_),
      url = wire_field(data, "url", "character", NA_character_),
      existed_before = isTRUE(data[["existed_before"]]),
      previous_sha256 = if (wire_null(previous)) NA_character_ else as.character(previous),
      changed = isTRUE(data[["changed"]]),
      go = proteins_parse_go_object(data[["go"]]),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_go_update"
  )
}

# ---------------------------------------------------------------- the public surface

#' Annotate MetaMorpheus protein groups with Gene Ontology terms, keeping every member
#'
#' Wraps mzLib's `GoGroupAnnotator` over a stored protein-group table
#' (`ProteinGroupFromTsv.ToGoAnnotationGroups`). For each group it returns one row per GO term that
#' **any** member holds - directly, or by propagation up `is_a` and `part_of` - and each row names
#' the members that carry it. Nothing is collapsed: MetaMorpheus picks no leading protein, so the
#' table is the union, and the views people usually want are filters on its rows.
#'
#' **Every non-decoy group gets at least one row.** A group with no term gets a single row whose
#' `annotation_status` says why: `"no_go_terms"` (its members have none), `"no_entry"` (a member is
#' not in the database - annotate against the database the search used), or `"contaminant"`.
#' Decoy groups get none.
#'
#' **Pin the ontology.** Terms and their ancestors change between GO releases, so a result means
#' something only relative to one go.obo. This reads the file you name and never downloads one:
#' fetch a release on purpose with [proteins_update_go()], keep the file, and record `go$sha256`
#' with your results.
#'
#' @param groups A MetaMorpheus protein-group table: `AllQuantifiedProteinGroups.tsv`,
#'   `AllProteinGroups.tsv` (a search without quantification), or one file's
#'   `<file>_ProteinGroups.tsv`.
#' @param database The UniProt XML (`.xml` or `.xml.gz`) whose GO terms annotate the members: the
#'   proteome the search used. A FASTA is refused, because it carries no GO and every group would
#'   read `"no_go_terms"` whatever the proteins are.
#' @param go_obo A go.obo file. It must exist; nothing is fetched.
#' @param category_map Optional: your own term-to-category map in mzLib's format
#'   (`#!category_map_format 1`, `#!map_name`, `#!map_version`, then `category`, `subcategory` and
#'   `anchor_go_id` columns). mzLib ships no vocabulary. Adds `categories`.
#' @param skip_unknown_go_ids `FALSE` (the default) fails when the database cites a GO id the
#'   ontology release lacks - usually a UniProt release newer than the go.obo - and names every
#'   missing id. `TRUE` drops each such id and lists it in `unresolved_go_ids`.
#' @param out Write the **whole** table here with mzLib's own `GoAnnotationTsv` writer, provenance
#'   header included, whatever `limit` and `offset` are. Must end in `.tsv`. For a large run pair it
#'   with `limit = 0`, so only the summary comes back.
#' @param categories_out Write the category table here with mzLib's `GoCategoryTsv` writer. Must
#'   end in `.tsv`, and needs `category_map`.
#' @param limit At most this many rows (rows, not groups) in `annotations`. `NULL`, the default,
#'   returns them all. Never shortens `out`.
#' @param offset Rows to skip before the window. Default 0.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait: a whole proteome XML takes a
#'   while.
#'
#' @return An `mzlibr_go_annotations`:
#'
#'   - `annotations`, a data.frame with one row per (group, term) in the returned window, under
#'     mzLib's own `GoAnnotationTsv` column names and order (`column_names`). `accession_used`,
#'     `accession_direct`, `accession_inherited`, `evidence` and `entrapment_members` are list
#'     columns of character vectors; `evidence_by_member` is a list column of named lists, member to
#'     its own evidence codes. `go_id`, `go_name`, `aspect`, `inherited` and `propagated` are `NA`
#'     on a term-less row. `n_members` and `n_with` count proteins (members); `n_with` is 0 on a
#'     term-less row. `q_value` is the group's q-value, a fraction (0 to 1), never filtered on.
#'   - `group_count` groups annotated (every non-decoy group, each once); `table_row_count` rows in
#'     the table read, decoys included; `decoy_group_count` decoy groups skipped.
#'   - `row_count` rows in the whole table; `returned_count` rows in `annotations`; `offset` rows
#'     skipped; `truncated` is `TRUE` when `limit` or `offset` left rows out.
#'   - Provenance: `groups_file` and `groups_file_sha256`; `annotation_database` (`path`,
#'     `file_type`, `reader`, `protein_count`, `sha256` of the decompressed database, `caveats`);
#'     `go` (`source_file_name`, `sha256`, `release`, `term_count`); and `header`, mzLib's `#!`
#'     header as a named character vector. **Its counters count groups at q <= 0.01, not rows.**
#'   - `skip_unknown_go_ids`, `unresolved_go_ids`, `caveats`.
#'   - `written` (`path`, `row_count`) and `categories_written` (`path`), or `NULL` when not asked.
#'   - `categories`, or `NULL` without `category_map`: `map_name`, `map_version`,
#'     `source_file_name`, `sha256`, `anchor_count`, `row_count`, and `records`, one row per (term,
#'     category, subcategory). `subcategory` is `NA` where the term reaches only the category's own
#'     anchor; a term under no anchor has no row. Join to `annotations` on `go_id`.
#'
#' @section The views are filters, not options:
#'
#' - **Consensus** - every member carries the term: `n_with == n_members`.
#' - **Direct annotations only**: `!propagated`.
#' - **No borrowing**: `!inherited`. An isoform (`P04406-2`) or sequence variant (`P04406_A20T`)
#'   absent from the database takes its entry's terms and is listed in `accession_inherited`; an
#'   isoform can differ from its entry in cellular component, which is why the row says so.
#'
#' @seealso [proteins_update_go()], [proteins_read()] for a database's own GO terms
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' go <- proteins_annotate_go("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
#'   go_obo = "go-pxd036557.obo", category_map = "organelle_map.tsv")
#' go
#' go$header[c("status_annotated", "status_contaminant")]
#' histones <- go$annotations[go$annotations$protein_group == "P0C0S5|Q71UI9", ]
#' c(rows = nrow(histones), consensus = sum(histones$n_with == histones$n_members))
#' head(go$categories$records)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec proteins.annotate-go
#' @export
proteins_annotate_go <- function(groups, database, go_obo, category_map = NULL,
                                 skip_unknown_go_ids = FALSE, out = NULL, categories_out = NULL,
                                 limit = NULL, offset = 0, timeout = NULL) {
  args <- proteins_build_go_request(
    groups, database, go_obo, category_map, skip_unknown_go_ids, out, categories_out, limit, offset
  )
  bridge_require_verb("proteins annotate-go", PROTEINS_GO_SINCE)
  proteins_parse_go_annotations(bridge_invoke(args, timeout = timeout))
}

#' Download the current Gene Ontology release to a go.obo, on purpose
#'
#' Wraps mzLib's `Loaders.UpdateGeneOntology`. The whole go.obo (tens of megabytes) is streamed from
#' GO's PURL, which always serves the **current** release. When a file is already at the path, it is
#' kept beside the new one as `go.obo.<yyyyMMdd-HHmmss-fff>` if the download differs, and left alone
#' if it is the same, so earlier runs stay reproducible. A failed download leaves any existing file
#' untouched.
#'
#' This is the only function in mzLibR that fetches a go.obo. [proteins_annotate_go()] never does,
#' so the release a result was computed against is always a file you chose to keep.
#'
#' @param go_obo Where to write, e.g. `"go.obo"`. Its folder must exist.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait; mzLib itself gives up after two
#'   minutes without data.
#'
#' @return An `mzlibr_go_update`: `go_obo_file` written; `url` fetched from; `existed_before`;
#'   `previous_sha256` of the file that was there, `NA` when there was none; `changed`, whether the
#'   file on disk is now different; `go`, the release now at the path (`source_file_name`,
#'   `sha256`, `release`, `term_count`); and `caveats`, which name the backup kept.
#'
#' @seealso [proteins_annotate_go()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' update <- proteins_update_go("go.obo")
#' update
#' update$go[c("release", "term_count")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec proteins.update-go
#' @export
proteins_update_go <- function(go_obo, timeout = NULL) {
  path <- proteins_go_path(go_obo, "go_obo", required = TRUE)
  bridge_require_verb("proteins update-go", PROTEINS_GO_SINCE)
  proteins_parse_go_update(bridge_invoke(c("proteins", "update-go", "--go-obo", path), timeout = timeout))
}

# ---------------------------------------------------------------- printing

#' Print a Gene Ontology annotation
#'
#' @param x A [proteins_annotate_go()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_go_annotations <- function(x, ...) {
  cat("<mzlibr_go_annotations> ", format(x$group_count), " groups, ", format(x$row_count),
    " (group, term) rows", if (isTRUE(x$truncated)) paste0(", ", format(x$returned_count), " returned") else "",
    "\n",
    sep = ""
  )
  if (!is.null(x$go)) {
    cat("  go.obo: ", x$go$source_file_name,
      if (!is.na(x$go$release)) paste0(", ", x$go$release) else ", no data-version",
      ", sha256 ", substr(x$go$sha256, 1L, 12L), "\n",
      sep = ""
    )
  }
  statuses <- paste0("status_", PROTEINS_ANNOTATION_STATUSES)
  shown <- x$header[intersect(statuses, names(x$header))]
  if (length(shown) > 0L) {
    cat("  groups at q <= ", x$header[["counter_q_value_max"]], ": ",
      paste0(sub("^status_", "", names(shown)), " ", shown, collapse = ", "), "\n",
      sep = ""
    )
  }
  if (!is.null(x$categories)) {
    cat("  categories: ", x$categories$map_name, " v", x$categories$map_version, ", ",
      format(x$categories$row_count), " rows\n",
      sep = ""
    )
  }
  if (!is.null(x$written)) cat("  written: ", x$written$path, "\n", sep = "")
  if (length(x$unresolved_go_ids) > 0L) {
    cat("  ! dropped GO ids the release lacks: ", paste(x$unresolved_go_ids, collapse = ", "), "\n", sep = "")
  }
  for (caveat in x$caveats) {
    cat("  ! ", caveat, "\n", sep = "")
  }
  invisible(x)
}

#' Print a Gene Ontology update
#'
#' @param x A [proteins_update_go()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_go_update <- function(x, ...) {
  cat("<mzlibr_go_update> ", x$go_obo_file, if (isTRUE(x$changed)) " changed" else " unchanged",
    "\n",
    sep = ""
  )
  if (!is.null(x$go)) {
    cat("  release: ", if (!is.na(x$go$release)) x$go$release else "(no data-version)", ", ",
      format(x$go$term_count), " terms\n",
      sep = ""
    )
  }
  for (caveat in x$caveats) {
    cat("  ! ", caveat, "\n", sep = "")
  }
  invisible(x)
}
