# A residue across species, through Ensembl Compara's gene-tree alignment (mzLib #1465).
#
#   proteins_map_residues()   which residue of a protein in one species sits in the same column of
#                             Compara's gene-tree alignment as residue p of a protein in another,
#                             or why none does
#
# Ported from pyMzLib's `proteins.py` (map_residues, ResidueMappings, ComparaAlignment,
# AlignedProteins), which decided the verb, the wire fields and the caveats. What changes here is
# only the projection into R's idiom, as for proteins_build_orthology() beside it: the queries are a
# data.frame or a list of `list(source, position, target)`; the gene sets a data.frame with a
# `species` column; the alignment and gene-tree files lists; the tables data.frames, with
# `alignment_index` 1-based as every position in R is; and `NA` where pyMzLib has `None`.
#
# Local files only, and no repair: positions reach mzLib unchecked (mzLib owns the rule that one
# outside its protein fails the call), every outcome is mzLib's ComparaColumnOutcome name, and a gap
# stays NA - never the nearest residue.

# The pyMzLib release whose bridge first dispatches `genes map-residues`.
PROTEINS_MAP_RESIDUES_SINCE <- "0.7.0"

# Every outcome proteins_map_residues() can report: mzLib's ComparaColumnOutcome names, verbatim and
# in mzLib's order. Every query gets exactly one.
PROTEINS_COLUMN_OUTCOMES <- c(
  "Aligned", "GapInTarget", "DifferentAlignments", "SourceNotInAlignment", "TargetNotInAlignment"
)

# ---------------------------------------------------------------- argument shaping

# One query as list(source, position, target), from a data.frame row or a list element.
proteins_query_parts <- function(queries) {
  if (is.data.frame(queries)) {
    if (ncol(queries) != 3L) {
      stop(mzlib_usage_error(paste0(
        "queries as a data.frame must have three columns, source_protein_id, position and ",
        "target_protein_id; got ", ncol(queries), "."
      )))
    }
    return(lapply(seq_len(nrow(queries)), function(i) {
      list(queries[[1L]][[i]], queries[[2L]][[i]], queries[[3L]][[i]])
    }))
  }
  if (!is.list(queries)) {
    stop(mzlib_usage_error(paste0(
      "queries must be a data.frame (source_protein_id, position, target_protein_id) or a list of ",
      "list(source_protein_id, position, target_protein_id); got ", class(queries)[1L], "."
    )))
  }
  queries
}

# The queries as stdin lines, `source<TAB>position<TAB>target`, refusing what cannot travel. The
# position is sent as given: mzLib, not mzLibR, decides whether it lies inside its protein.
proteins_query_lines <- function(queries) {
  parts <- proteins_query_parts(queries)
  lines <- character(length(parts))
  for (i in seq_along(parts)) {
    query <- parts[[i]]
    if (!is.list(query) || length(query) != 3L) {
      stop(mzlib_usage_error(paste0(
        "Query ", i, " must be list(source_protein_id, position, target_protein_id); got ",
        paste(deparse(query), collapse = " "), "."
      )))
    }
    position <- query[[2L]]
    if (!is.numeric(position) || length(position) != 1L || is.na(position) || !is.finite(position) ||
      position != round(position)) {
      stop(mzlib_usage_error(paste0(
        "Query ", i, "'s position must be a whole number (a 1-based residue number); got ",
        paste(deparse(position), collapse = " "), "."
      )))
    }
    ids <- list(source = query[[1L]], target = query[[3L]])
    for (name in names(ids)) {
      value <- ids[[name]]
      if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
        stop(mzlib_usage_error(paste0(
          "Query ", i, "'s ", name, " must be a protein id; got ", paste(deparse(value), collapse = " "), "."
        )))
      }
      if (grepl("[\t\r\n]", value)) {
        stop(mzlib_usage_error(paste0("Query ", i, "'s ", name, " contains a tab or line break: '", value, "'.")))
      }
    }
    # sprintf, not formatC(format = "d"): any whole number is sent, mzLib refuses what it must.
    lines[[i]] <- paste(trimws(ids$source), sprintf("%.0f", position), trimws(ids$target), sep = "\t")
  }
  if (length(lines) == 0L) {
    stop(mzlib_usage_error(
      "queries is empty. Give at least one list(source_protein_id, position, target_protein_id)."
    ))
  }
  lines
}

proteins_build_map_residues_request <- function(queries, alignment, gene_trees, gtfs, gene_sets, proteins) {
  query_lines <- proteins_query_lines(queries)
  species_lines <- proteins_species_lines(gtfs, gene_sets, required = FALSE)
  for (required in c("alignment", "gene_trees")) {
    if (is.null(if (required == "alignment") alignment else gene_trees)) {
      stop(mzlib_usage_error(paste0(
        required, " is required: ",
        if (required == "alignment") {
          "Compara's peptide alignments, Compara.<release>.protein_default.aa.fasta.gz."
        } else {
          "Compara's <collection>.GeneTree_content.default.e<release>.txt.gz."
        }
      )))
    }
  }
  if (!is.logical(proteins) || length(proteins) != 1L || is.na(proteins)) {
    stop(mzlib_usage_error(paste0(
      "proteins must be TRUE or FALSE; got ", paste(deparse(proteins), collapse = " "), "."
    )))
  }
  args <- c(
    "genes", "map-residues",
    "--alignment", proteins_single_path(alignment, "alignment"),
    "--gene-trees", proteins_single_path(gene_trees, "gene_trees")
  )
  if (proteins) {
    args <- c(args, "--proteins")
  }
  # The species lines (none: mzLib keeps every gene), then `--`, then one query per line.
  stdin <- paste0(paste(c(species_lines, PROTEINS_SECTION, query_lines), collapse = "\n"), "\n")
  list(args = args, stdin = stdin)
}

# ---------------------------------------------------------------- parsing

# A `{row_count, column_names, columns}` block as a list with `records`, or NULL for a null block.
proteins_parse_aligned_proteins <- function(block) {
  if (!is.list(block) || is.null(names(block))) {
    return(NULL)
  }
  records <- wire_array_table(block)
  if ("alignment_index" %in% names(records)) {
    records$alignment_index <- as.numeric(records$alignment_index) + 1
  }
  list(
    row_count = proteins_count(block, "row_count"),
    column_names = wire_strings(block[["column_names"]]),
    records = records
  )
}

proteins_parse_residue_mappings <- function(data) {
  alignment <- if (is.list(data[["alignment"]])) data[["alignment"]] else list()
  trees <- if (is.list(data[["gene_trees"]])) data[["gene_trees"]] else list()
  restricted <- trees[["restricted_to_gene_sets"]]
  outcomes <- proteins_parse_counts(data[["outcome_counts"]])
  structure(
    list(
      alignment = list(
        source_file_name = wire_field(alignment, "source_file_name", "character", NA_character_),
        sha256 = wire_field(alignment, "sha256", "character", NA_character_),
        release = wire_field(alignment, "release", "character", NA_character_),
        collection = wire_field(alignment, "collection", "character", NA_character_),
        alignment_count = proteins_count(alignment, "alignment_count"),
        protein_count = proteins_count(alignment, "protein_count"),
        restricted_to_gene_tree_content = wire_field(alignment, "restricted_to_gene_tree_content", "character",
          NA_character_)
      ),
      gene_trees = list(
        source_file_name = wire_field(trees, "source_file_name", "character", NA_character_),
        sha256 = wire_field(trees, "sha256", "character", NA_character_),
        release = wire_field(trees, "release", "character", NA_character_),
        collection = wire_field(trees, "collection", "character", NA_character_),
        member_count = proteins_count(trees, "member_count"),
        # NULL, not character(0): every gene in the file was kept, which is not "restricted to none".
        restricted_to_gene_sets = if (is.list(restricted)) wire_strings(restricted) else NULL
      ),
      gene_sets = proteins_species_frame(data[["gene_sets"]], c(
        "source_file_name", "sha256", "release", "genome_build", "genebuild_last_updated", "gene_count"
      )),
      row_count = proteins_count(data, "row_count"),
      outcome_counts = outcomes,
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data),
      proteins_included = isTRUE(data[["proteins_included"]]),
      proteins = proteins_parse_aligned_proteins(data[["proteins"]]),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_residue_mappings"
  )
}

# ---------------------------------------------------------------- the public surface

#' Map residues across species through Ensembl Compara's gene-tree alignments
#'
#' Answers **which mouse residue is human p53's Ser15?** Ensembl Compara has already aligned every
#' protein of each gene tree, across all its species, to build the tree. A residue's column in that
#' alignment is the published answer. This wraps mzLib's `ComparaGeneTreeAlignment` (mzLib #1465):
#' `Load` reads Compara's peptide alignment of every gene tree, keeping the canonical proteins of
#' your species' genes, and `MapResidue` answers each query: which residue of the target protein
#' sits in the same alignment column as residue `position` of the source protein. Local files only;
#' nothing is downloaded.
#'
#' Every query gets exactly one outcome, mzLib's own name for it: `"Aligned"` (with the column, the
#' target position and its residue), `"GapInTarget"` (no target residue in that column - never the
#' nearest one), `"DifferentAlignments"` (the proteins are in different gene trees), or
#' `"SourceNotInAlignment"` / `"TargetNotInAlignment"` (that protein was not kept: only each gene's
#' canonical protein is aligned).
#'
#' **Positions are Ensembl's.** To start from a UniProt position, align the UniProt sequence to the
#' canonical Ensembl protein with [proteins_align()] first, and align the target's Ensembl protein to
#' its UniProt sequence after. `proteins = TRUE` returns the Ensembl sequences for that.
#'
#' **A position outside the source protein fails the whole call.** mzLib's `MapResidue` throws for it
#' rather than returning an outcome, and mzLibR passes positions through unchecked, so keep them
#' within each source protein: with `proteins = TRUE`, a sequence's length is its last valid
#' position.
#'
#' @param queries The `(source protein, position, target protein)` queries, one result row each, in
#'   this order, duplicates included: a data.frame with three columns, `source_protein_id`,
#'   `position` and `target_protein_id`, or a list of `list(source_protein_id, position,
#'   target_protein_id)`. The ids are Ensembl protein stable ids; the position is a whole number, the
#'   1-based residue in the source protein. They travel on stdin after the species lines.
#' @param alignment Compara's peptide alignments,
#'   `Compara.<release>.protein_<collection>.aa.fasta.gz` (908 MB for `protein_default` in release
#'   116), from `emf/ensembl-compara/homologies/`. Required.
#' @param gene_trees Compara's gene-tree file of the same release,
#'   `<collection>.GeneTree_content.default.e<release>.txt.gz`. It names each gene's tree and its
#'   canonical protein. Required.
#' @param gtfs A named character vector from species to its Ensembl primary-assembly GTF, as for
#'   [proteins_build_orthology()]: `c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf.gz")`. Only the
#'   genes in these are kept, and so only their canonical proteins. `NULL` (the default) with no
#'   `gene_sets`: mzLib's default keeps every gene in the gene-tree file - millions for the
#'   vertebrates - and a caveat says so. The answers for your species' proteins are the same.
#' @param gene_sets Instead of a GTF, a compact gene table mzLib's `EnsemblGeneSetWriter` made from
#'   one, by species in the same way. A species may be in `gtfs` or `gene_sets`, not both.
#' @param proteins Also return the kept proteins with their ungapped sequences, as `proteins`.
#'   `FALSE` (the default): `proteins` is `NULL`.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait: for human, mouse and rat at
#'   release 116 the whole call takes about a minute, nearly all of it reading the alignment.
#'
#' @return An `mzlibr_residue_mappings`.
#'
#'   `records` is a data.frame with one row per query, in the order given: `source_protein_id`,
#'   `position` (the 1-based residue in the source, as given) and `target_protein_id`, as given;
#'   `outcome`, one of the five outcomes above; `column`, the source residue's 1-based column in the
#'   gene tree's alignment, `NA` when the two proteins share no alignment; `target_position`, the
#'   1-based position of the target's residue in that column, in the target's ungapped sequence,
#'   `NA` unless `"Aligned"`; and `target_residue`, the target's residue there, one letter, `NA`
#'   unless `"Aligned"`. It may be a different amino acid from the source's: a column pairs
#'   homologous positions, not identical ones.
#'
#'   `row_count` counts the queries, one row each (queries). `outcome_counts` is a named vector,
#'   outcome to the queries with it (queries), every outcome present in mzLib's order, zeros
#'   included. `column_names` lists the table's columns.
#'
#'   `alignment` is a list describing the alignment file: `source_file_name`; `sha256`, of the
#'   file's bytes as read (the compressed bytes for a `.gz`); `release` and `collection`, read from
#'   the file name (`NA` when the name carries none); `alignment_count`, the alignments in the file,
#'   one per gene tree, kept proteins or not (alignments); `protein_count`, the canonical proteins
#'   kept (proteins); and `restricted_to_gene_tree_content`, the sha256 of the gene-tree file whose
#'   canonical proteins were kept, the same as `gene_trees$sha256`.
#'
#'   `gene_trees` is a list, as [proteins_build_orthology()] reports it: `source_file_name`,
#'   `sha256`, `release`, `collection`, `member_count`, the genes kept from the file (genes), and
#'   `restricted_to_gene_sets`, the gene sets' sha256s, or `NULL` when no gene set was given and
#'   every gene was kept.
#'
#'   `gene_sets` is a data.frame, one row per species given (none without `gtfs` or `gene_sets`):
#'   `species`, `source_file_name`, `sha256`, `release`, `genome_build` and
#'   `genebuild_last_updated` (`NA` when the GTF's name or headers carry none), and `gene_count`, in
#'   genes, every biotype.
#'
#'   `proteins_included` says whether `proteins = TRUE` was given. `proteins` is `NULL` without it;
#'   with it, a list: `row_count` (proteins kept), `column_names`, and `records`, a data.frame with
#'   one row per kept protein in order of protein id: `protein_id`, the Ensembl protein stable id;
#'   `alignment_index`, which alignment holds it (**1-based in R**, in file order; two proteins share
#'   a column only when they share this); `tree_id`, its gene's Compara gene tree; and `sequence`,
#'   its residues with the alignment's gaps removed, the sequence every position counts in.
#'   `setNames(m$proteins$records$sequence, m$proteins$records$protein_id)` gives the sequences by
#'   id, for [proteins_align()]. `caveats` names this call's traps: gene trees read unrestricted,
#'   files named for different releases or for none, queries naming a protein the alignment did not
#'   keep.
#'
#' @section What each outcome means:
#'
#' `"Aligned"`: the target has a residue in the source residue's column; `column`, `target_position`
#' and `target_residue` are set. `"GapInTarget"`: the target has a gap in that column, so no residue
#' corresponds; `column` is set, the target's are `NA`. `"DifferentAlignments"`: the two proteins
#' are in different gene trees, so this alignment says nothing about them. It does not say they are
#' unrelated. `"SourceNotInAlignment"` and `"TargetNotInAlignment"`: that protein was not kept - not
#' in the file, not its gene's canonical protein, or outside your gene sets. Another isoform is
#' never mapped through the canonical one.
#'
#' @section The traps:
#'
#' **A position past the end of the source protein fails the whole call** with mzLib's message (for
#' example "ENSP00000269305 has 393 residues."), as an `mzlib_usage_error`, and no row is returned,
#' not even for the queries that were fine. A source protein that was not kept is an outcome,
#' `"SourceNotInAlignment"`, checked before the position. **Use Ensembl ids and Ensembl
#' positions**: a UniProt accession is `"SourceNotInAlignment"`, never an error. **Use the alignment
#' and the gene trees of one release**: mzLib does not compare them, and a gene's canonical protein
#' and tree can change between releases. Keep Ensembl's file names, so the result can report the
#' release each one carries. **Sharing a gene tree is not orthology**: a tree holds paralogs too,
#' so ask [proteins_build_orthology()] which genes are orthologs, then map their residues here.
#' **Check `target_residue`** before carrying a modification across: a phosphosite on S is not a
#' site on A.
#'
#' @seealso [proteins_align()] for the UniProt legs, [proteins_build_orthology()]
#' @spec genes.map-residues
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # A real cut of Ensembl release 116: the TP53 and GAPDH gene trees' alignments, with the human,
#' # mouse and rat rows exactly as Ensembl publishes them.
#' gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
#'           mus_musculus = "Mus_musculus.GRCm39.116.gtf",
#'           rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
#' site <- proteins_map_residues(list(list("ENSP00000269305", 15, "ENSMUSP00000104298")),
#'   alignment = "Compara.116.protein_default.aa.fasta.gz",
#'   gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz", gtfs = gtfs)
#' site
#' # Human p53's Ser15 is mouse p53's Ser18:
#' site$records[, c("outcome", "column", "target_position", "target_residue")]
#' c(site$alignment$release, site$gene_trees$collection)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_map_residues <- function(queries, alignment, gene_trees, gtfs = NULL, gene_sets = NULL, proteins = FALSE,
                                  timeout = NULL) {
  # Every argument is checked before any process is started, the verb check included. A missing
  # alignment or gene-tree file is a usage error that says what to pass, not R's "argument is missing".
  request <- proteins_build_map_residues_request(
    if (missing(queries)) list() else queries,
    if (missing(alignment)) NULL else alignment,
    if (missing(gene_trees)) NULL else gene_trees,
    gtfs, gene_sets, proteins
  )
  bridge_require_verb("genes map-residues", PROTEINS_MAP_RESIDUES_SINCE)
  proteins_parse_residue_mappings(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Print residue mappings
#'
#' @param x A [proteins_map_residues()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_residue_mappings <- function(x, ...) {
  cat("<mzlibr_residue_mappings> ", format(x$row_count), " queries through ", x$alignment$source_file_name,
    " (Ensembl ", if (is.na(x$alignment$release)) "release unknown" else x$alignment$release, ", ",
    format(x$alignment$protein_count), " proteins kept)\n",
    sep = ""
  )
  shown <- x$outcome_counts[x$outcome_counts > 0]
  if (length(shown) > 0L) {
    cat("  ", paste0(names(shown), " ", shown, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$proteins)) {
    cat("  proteins: ", format(x$proteins$row_count), " kept, with their sequences\n", sep = "")
  }
  for (caveat in x$caveats) cat("  ! ", caveat, "\n", sep = "")
  invisible(x)
}
