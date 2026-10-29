# Orthologs across species: mzLib's Ensembl Compara orthology snapshot (mzLib #1381).
#
#   proteins_build_orthology()   which genes of one species are orthologous to which genes of
#                                another, kept exactly as Compara states it, and with pair_status
#                                every gene's status in one other species; out = writes it as Parquet
#   proteins_orthology_pair()    one species pair's counts, given in either order
#
# Ported from pyMzLib's `proteins.py` (build_orthology, OrthologySnapshot and its parts), which
# decided the verb, the wire fields and the caveats. What changes here is only the projection into
# R's idiom: the per-species maps (`gene_sets`, `dumps`) are data.frames with a `species` column,
# the pairs a data.frame, the gene-tree file and the written snapshot lists, and the species
# arguments named character vectors from Ensembl production name to path - the way sdrf_pool()
# names its documents by label.
#
# Local files only: nothing is downloaded and no default is invented. Every input is the caller's,
# the release included, and mzLib's refusals (a missing species' dump, a file of another release)
# reach the caller as mzLib wrote them.

# The pyMzLib release whose bridge first dispatches `genes orthology`.
PROTEINS_ORTHOLOGY_SINCE <- "0.6.0"

# Every status mzLib's OrthologySnapshot.StatusName can report, in mzLib's order. Every gene of a
# pair-status table gets exactly one; not_in_gene_set never occurs there.
PROTEINS_ORTHOLOGY_STATUSES <- c(
  "has_ortholog", "no_edge_in_shared_tree", "tree_lacks_target_species", "not_in_any_tree",
  "not_in_gene_set"
)

# ---------------------------------------------------------------- argument shaping

proteins_orthology_release <- function(release) {
  ok <- (is.character(release) || is.numeric(release)) && !is.logical(release) &&
    length(release) == 1L && !is.na(release)
  text <- if (!ok) {
    ""
  } else if (is.numeric(release)) {
    if (is.finite(release) && release == round(release) && release >= 0) sprintf("%.0f", release) else ""
  } else {
    trimws(release)
  }
  if (!nzchar(text) || grepl("\\s", text)) {
    stop(mzlib_usage_error(paste0(
      "release must be an Ensembl release, e.g. \"116\"; got ", paste(deparse(release), collapse = " "), "."
    )))
  }
  text
}

# The stdin species lines, `species<TAB>gtf|gene-set<TAB>path`, in the order given. genes
# orthology needs at least one species; genes map-residues takes none (mzLib keeps every gene).
proteins_species_lines <- function(gtfs, gene_sets, required = TRUE) {
  lines <- character(0)
  seen <- character(0)
  for (kind in c("gtf", "gene-set")) {
    mapping <- if (kind == "gtf") gtfs else gene_sets
    what <- if (kind == "gtf") "gtfs" else "gene_sets"
    if (is.null(mapping)) next
    if (!is.character(mapping) || is.null(names(mapping))) {
      stop(mzlib_usage_error(paste0(
        what, " must be a named character vector from Ensembl production name to path, e.g. ",
        "c(homo_sapiens = \"Homo_sapiens.GRCh38.116.gtf.gz\"); got ",
        paste(deparse(mapping), collapse = " "), "."
      )))
    }
    for (i in seq_along(mapping)) {
      species <- trimws(names(mapping)[[i]])
      if (is.na(species) || !nzchar(species)) {
        stop(mzlib_usage_error(paste0("Every species in ", what, " must be named; element ", i, " is not.")))
      }
      if (grepl("[\t\r\n,]", species)) {
        stop(mzlib_usage_error(paste0("A species name cannot hold a tab, line break or comma: '", species, "'.")))
      }
      if (species %in% seen) {
        stop(mzlib_usage_error(paste0(
          "'", species, "' is given twice (in gtfs and gene_sets, or twice in one); give each species one gene set."
        )))
      }
      seen <- c(seen, species)
      path <- proteins_paths(unname(mapping[[i]]), paste(what, species))
      lines <- c(lines, paste0(species, "\t", kind, "\t", path))
    }
  }
  if (length(lines) == 0L && required) {
    stop(mzlib_usage_error(paste0(
      "At least one species is required: gtfs = c(homo_sapiens = \"Homo_sapiens.GRCh38.116.gtf.gz\", ...) ",
      "(or gene_sets = for compact gene tables)."
    )))
  }
  lines
}

proteins_build_orthology_request <- function(release, gene_trees, dumps, gtfs, gene_sets,
                                             restrict_gene_trees, out, pair_status) {
  release_text <- proteins_orthology_release(release)
  species_lines <- proteins_species_lines(gtfs, gene_sets)
  dump_lines <- proteins_paths(dumps, "homology dump")
  if (length(dump_lines) == 0L) {
    stop(mzlib_usage_error(paste0(
      "dumps is empty. Give every species' Compara homology dump: each pair's relationships are ",
      "split between the pair's two dumps, so mzLib requires them all."
    )))
  }
  if (is.null(gene_trees)) {
    stop(mzlib_usage_error(
      "gene_trees is required: Compara's <collection>.GeneTree_content.default.e<release>.txt.gz."
    ))
  }
  args <- c(
    "genes", "orthology", "--release", release_text,
    "--gene-trees", proteins_single_path(gene_trees, "gene_trees")
  )
  if (!is.logical(restrict_gene_trees) || length(restrict_gene_trees) != 1L || is.na(restrict_gene_trees)) {
    stop(mzlib_usage_error("restrict_gene_trees must be TRUE or FALSE."))
  }
  if (restrict_gene_trees) {
    args <- c(args, "--restrict-gene-trees")
  }
  if (!is.null(out)) {
    args <- c(args, "--out", proteins_single_path(out, "out"))
  }
  if (!is.null(pair_status)) {
    if (!is.character(pair_status) || length(pair_status) != 2L) {
      stop(mzlib_usage_error(paste0(
        "pair_status must be c(from_species, to_species), e.g. c(\"homo_sapiens\", \"mus_musculus\"); got ",
        paste(deparse(pair_status), collapse = " "), "."
      )))
    }
    names <- trimws(unname(pair_status))
    if (anyNA(names) || any(!nzchar(names)) || any(grepl("[\t\r\n,]", names))) {
      stop(mzlib_usage_error(paste0(
        "Each species in pair_status must be a production name; got ", paste(deparse(pair_status), collapse = " "), "."
      )))
    }
    args <- c(args, "--pair-status", paste(names, collapse = ","))
  }
  stdin <- paste0(paste(c(species_lines, PROTEINS_SECTION, dump_lines), collapse = "\n"), "\n")
  list(args = args, stdin = stdin)
}

# ---------------------------------------------------------------- parsing

# A `{species: {...}}` map as a data.frame, one row per species in wire order, `species` first.
proteins_species_frame <- function(value, fields) {
  if (!is.list(value) || length(value) == 0L) {
    columns <- c(list(species = character(0)), rep(list(character(0)), length(fields)))
    names(columns) <- c("species", fields)
    return(as.data.frame(columns, stringsAsFactors = FALSE, optional = TRUE))
  }
  rows <- lapply(names(value), function(species) c(list(species = species), value[[species]]))
  frame <- wire_objects(rows, empty = c("species", fields))
  for (field in setdiff(c("species", fields), names(frame))) {
    frame[[field]] <- rep(NA, nrow(frame))
  }
  frame[, c("species", fields, setdiff(names(frame), c("species", fields))), drop = FALSE]
}

proteins_parse_orthology <- function(data) {
  trees <- if (is.list(data[["gene_trees"]])) data[["gene_trees"]] else list()
  restricted <- trees[["restricted_to_gene_sets"]]
  written <- data[["written"]]
  status <- data[["pair_status"]]
  structure(
    list(
      release = wire_field(data, "release", "character", NA_character_),
      species = wire_strings(data[["species"]]),
      format = wire_field(data, "format", "character", NA_character_),
      format_version = proteins_count(data, "format_version"),
      snapshot_id = wire_field(data, "snapshot_id", "character", NA_character_),
      identical_duplicates_dropped = proteins_count(data, "identical_duplicates_dropped"),
      gene_sets = proteins_species_frame(data[["gene_sets"]], c(
        "source_file_name", "sha256", "release", "genome_build", "genebuild_last_updated", "gene_count"
      )),
      gene_trees = list(
        source_file_name = wire_field(trees, "source_file_name", "character", NA_character_),
        sha256 = wire_field(trees, "sha256", "character", NA_character_),
        release = wire_field(trees, "release", "character", NA_character_),
        collection = wire_field(trees, "collection", "character", NA_character_),
        member_count = proteins_count(trees, "member_count"),
        # NULL, not character(0): every gene in the file was kept, which is not "restricted to none".
        restricted_to_gene_sets = if (is.list(restricted)) wire_strings(restricted) else NULL
      ),
      dumps = proteins_species_frame(data[["dumps"]], c(
        "source_file_name", "sha256", "release", "rows_read", "rows_kept"
      )),
      pairs = wire_objects(data[["pairs"]], empty = c(
        "species_a", "species_b", "row_count", "ortholog_count", "paralog_count", "homoeolog_count"
      )),
      written = if (is.list(written)) {
        directory <- wire_field(written, "directory", "character", NA_character_)
        list(
          directory = directory,
          manifest = file.path(directory, "manifest.json"),
          files = wire_objects(written[["files"]], empty = c("path", "rows", "bytes", "sha256"))
        )
      } else {
        NULL
      },
      pair_status = if (is.list(status)) {
        list(
          from_species = wire_field(status, "from_species", "character", NA_character_),
          to_species = wire_field(status, "to_species", "character", NA_character_),
          row_count = proteins_count(status, "row_count"),
          status_counts = proteins_parse_counts(status[["status_counts"]]),
          column_names = wire_strings(status[["column_names"]]),
          records = proteins_parse_table(status)
        )
      } else {
        NULL
      },
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_orthology_snapshot"
  )
}

# ---------------------------------------------------------------- the public surface

#' Build an Ensembl Compara orthology snapshot for your species, and write it as Parquet
#'
#' Answers **which mouse genes are orthologous to this human gene, and if none, why not?** Wraps
#' mzLib's `OrthologySnapshot.Build` and `OrthologySnapshotWriter.WriteAsync` (mzLib #1381): which
#' genes of one species are orthologous to which genes of another, kept exactly as Compara states
#' it. One-to-many and many-to-many orthologies are several rows, paralogs are kept, and nothing is
#' picked, merged or collapsed. Local files only; nothing is downloaded.
#'
#' **Every species' homology dump is required, and mzLib refuses a call without one.** Compara puts
#' each of a pair's relationships in one of the pair's two dumps, arbitrarily, and in practice they
#' partition: in release 116 every human-mouse orthology is in the *mouse* dump and none is in the
#' human one. Reading the dump named after your species gives a complete-looking, badly incomplete
#' answer.
#'
#' **Orthology is pairwise and not transitive.** A~B and B~C do not make A~C, so never chain pairs
#' into groups. The written `views.sql` has a `species_set3` macro that forms a tuple only when
#' every two members are orthologs.
#'
#' The relationships are in the Parquet files `out` writes, which mzLib cannot read back and mzLibR
#' does not read (it has no dependencies): query them with the `arrow` or `duckdb` packages, as the
#' proteins article shows. Without `out`, the call checks that the inputs fit together and, with
#' `pair_status`, says gene by gene how far Compara's inference got in one other species.
#'
#' @param release The Ensembl release of the Compara inputs, e.g. `"116"` (a whole number is sent
#'   as its digits). Required: mzLib holds the gene trees and every dump to it when their file
#'   names carry a release. A GTF's own number is not compared with it (Ensembl publishes some
#'   species with their own annotation release).
#' @param gene_trees Compara's gene-tree file, `<collection>.GeneTree_content.default.e<release>.txt.gz`.
#'   It tells a gene whose tree holds no gene of the other species from one whose tree does but got
#'   no ortholog call.
#' @param dumps One Compara homology dump per species
#'   (`homologies/<species>/Compara.<release>.protein_default.homologies.tsv.gz`), in any order: each
#'   dump's species is read from its rows. Sent on stdin after the species lines.
#' @param gtfs A named character vector from species to its Ensembl primary-assembly GTF, plain or
#'   `.gz`: `c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf.gz")`. The names are Ensembl production
#'   names, as the dumps' `species` column spells them, and they are the snapshot's species. Sent on
#'   stdin, one line per species.
#' @param gene_sets Instead of a GTF, a compact gene table mzLib's `EnsemblGeneSetWriter` made from
#'   one, by species in the same way (it keeps the GTF's provenance, so the snapshot is the same). A
#'   species may be in `gtfs` or `gene_sets`, not both.
#' @param restrict_gene_trees Keep only the gene-tree members in your gene sets. `FALSE` (the
#'   default) is mzLib's default and keeps every gene in the file - 4.2 million for the vertebrates
#'   of release 116. The snapshot is the same either way; only memory differs.
#' @param out A directory to write the snapshot to: Parquet tables (`genes/`, `members/`, `pairs/`),
#'   `views.sql` and `manifest.json`. It must not exist yet or be empty: mzLib writes a snapshot once
#'   and never updates one in place. `NULL` (the default): nothing is written.
#' @param pair_status `c(from_species, to_species)`: return every gene of `from_species` with its
#'   status in `to_species`. Two different species of the snapshot, in that order: a status is not
#'   symmetric. `NULL` (the default): not computed.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait: for human, mouse and rat at
#'   release 116, mzLib reports about 20 s to build and 40 s more to write.
#'
#' @return An `mzlibr_orthology_snapshot`. The relationships themselves are in the Parquet files
#'   `out` writes, not here: this says what went in, what came out, and, with `pair_status`, every
#'   gene's status in one other species.
#'
#'   `release` is the release you gave; `species` the snapshot's species, in sorted order; `format`
#'   is `"ensembl-orthology-snapshot"` and `format_version` its version, as `manifest.json` names
#'   them. `snapshot_id` is mzLib's id for the snapshot, a sha256 over the format, the release, the
#'   species and every input's sha256: the same inputs give the same id, written or not. Record it
#'   with any result drawn from the snapshot. `identical_duplicates_dropped` counts homology rows
#'   found in both of a pair's dumps, identical, and counted once, in rows: 0 in release 116, whose
#'   dumps partition.
#'
#'   `gene_sets` is a data.frame, one row per species: `species`, `source_file_name`, `sha256` (the
#'   value [proteins_resolve_genes()] reports as `gene_set_sha256`, so a resolution and a snapshot
#'   can be checked against one gene set), `release`, `genome_build` and `genebuild_last_updated`
#'   (`NA` when the GTF's name or headers carry none), and `gene_count`, in genes, every biotype.
#'
#'   `gene_trees` is a list: `source_file_name`, `sha256`, `release` (`NA` when the name carries
#'   none, so it was not checked), `collection` (e.g. `"vertebrates"`), `member_count`, the genes
#'   held from the file, and `restricted_to_gene_sets`, the gene sets' sha256s the members were
#'   restricted to, or `NULL` when every gene in the file was kept.
#'
#'   `dumps` is a data.frame, one row per species, the species read from the dump's rows:
#'   `source_file_name` (the same for every species, so it alone does not say whose dump it is),
#'   `sha256`, `release`, `rows_read`, every data row read and checked, in rows, and `rows_kept`, the
#'   rows whose two species are both snapshot species, in rows.
#'
#'   `pairs` is a data.frame, one row per species pair `(a, b)` with `a` sorting first, `(s, s)`
#'   included: `species_a`, `species_b`, and the relationships between them in rows - `row_count`,
#'   one row each, never merged or picked; `ortholog_count`, always 0 for `(s, s)`;
#'   `paralog_count`, always 0 between two species; and `homoeolog_count`. See
#'   [proteins_orthology_pair()].
#'
#'   `written` is `NULL` without `out`. With it, a list: `directory`, absolute; `manifest`, the path
#'   of `manifest.json`, which lists every input and file with its sha256; and `files`, a data.frame
#'   of every file written but the manifest, in order of path: `path` (relative to `directory`, with
#'   forward slashes), `rows` (rows in the Parquet table, `NA` for `views.sql`, 0 for a pair with no
#'   relationship, which is still written), `bytes` (the file's size in bytes) and `sha256`.
#'
#'   `pair_status` is `NULL` without `pair_status`. With it, a list: `from_species`, `to_species`,
#'   `row_count` (every gene of `from_species`, in genes), `status_counts` (a named vector, status to
#'   genes with it, every status present in mzLib's order, zeros included), `column_names`, and
#'   `records`, a data.frame with one row per gene in order of gene id: `gene_id`, `gene_name` (`NA`
#'   when the GTF gives none: a label, never a key), `gene_biotype` and `status`. `caveats` names this
#'   call's traps.
#'
#' @section What each status means:
#'
#' It says how far Compara's inference got, not what biology holds. `has_ortholog`: at least one
#' ortholog in `to_species`. `no_edge_in_shared_tree`: the gene's tree holds a `to_species` gene,
#' and Compara called no ortholog to it - the nearest the data comes to "no ortholog exists".
#' `tree_lacks_target_species`: the gene's tree holds no `to_species` gene, so nothing was
#' compared. `not_in_any_tree`: the gene is in the gene set but in no gene tree, so no call was
#' attempted. **Every biotype is a row**, so choose your denominator by filtering on
#' `gene_biotype`. One status's genes are
#' `ps$records$gene_id[ps$records$status == "not_in_any_tree"]`, with `ps <- snap$pair_status`.
#'
#' @section The traps:
#'
#' **The gene set and the dumps must be one release.** A homology naming a gene its gene set lacks
#' is refused as an `mzlib_bridge_error` (`InvalidDataException`): that is what a GTF and dumps from
#' different releases look like. Use the primary-assembly GTF, as for [proteins_resolve_genes()].
#' **Keep Ensembl's file names**: mzLib holds the gene trees and the dumps to `release` only when
#' their names carry it, and a caveat says which were not checked. **A directory that is not empty
#' is refused** as `out`: a snapshot is written once, and next release's goes in a new directory.
#'
#' @seealso [proteins_orthology_pair()], [proteins_resolve_genes()]
#' @spec genes.orthology
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own test world: four human genes, three mouse and two rat, with every human-mouse
#' # row in the mouse dump, as in release 116.
#' gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
#'           mus_musculus = "Mus_musculus.GRCm39.116.gtf",
#'           rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
#' dumps <- paste0("homologies/", names(gtfs), "/Compara.116.protein_default.homologies.tsv.gz")
#' snap <- proteins_build_orthology("116", gtfs = gtfs, dumps = dumps,
#'   gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
#'   pair_status = c("homo_sapiens", "mus_musculus"))
#' snap
#' snap$pair_status$records[, c("gene_id", "gene_biotype", "status")]
#' snap$dumps[, c("species", "rows_kept")]
#' proteins_orthology_pair(snap, "mus_musculus", "homo_sapiens")$row_count
#'
#' # With out = the snapshot is written, and its id is the same:
#' written <- proteins_build_orthology("116", gtfs = gtfs, dumps = dumps,
#'   gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
#'   restrict_gene_trees = TRUE, out = "compara-116")
#' identical(written$snapshot_id, snap$snapshot_id)
#' head(written$written$files[, c("path", "rows")], 7)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_build_orthology <- function(release, gene_trees, dumps, gtfs = NULL, gene_sets = NULL,
                                     restrict_gene_trees = FALSE, out = NULL, pair_status = NULL,
                                     timeout = NULL) {
  # The release, the gene trees and the dumps are required, with no default to fall back on: a
  # missing one is a usage error that says what to pass, not R's "argument is missing".
  request <- proteins_build_orthology_request(
    if (missing(release)) NULL else release, if (missing(gene_trees)) NULL else gene_trees,
    if (missing(dumps)) character(0) else dumps,
    gtfs, gene_sets, restrict_gene_trees, out, pair_status
  )
  bridge_require_verb("genes orthology", PROTEINS_ORTHOLOGY_SINCE)
  proteins_parse_orthology(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' One species pair's counts, in either order
#'
#' The row of a snapshot's `pairs` for species `a` and `b`, which may be given either way round:
#' side A is always the species that sorts first.
#'
#' @param snapshot A [proteins_build_orthology()] result.
#' @param a,b Two species of the snapshot, Ensembl production names; the same one twice for its
#'   paralogs.
#' @return A one-row data.frame: `species_a`, `species_b`, and `row_count`, `ortholog_count`,
#'   `paralog_count` and `homoeolog_count`, in rows.
#' @seealso [proteins_build_orthology()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' gtfs <- c(homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
#'           mus_musculus = "Mus_musculus.GRCm39.116.gtf",
#'           rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf")
#' snap <- proteins_build_orthology("116", gtfs = gtfs,
#'   dumps = paste0("homologies/", names(gtfs), "/Compara.116.protein_default.homologies.tsv.gz"),
#'   gene_trees = "vertebrates.GeneTree_content.default.e116.txt.gz",
#'   pair_status = c("homo_sapiens", "mus_musculus"))
#' proteins_orthology_pair(snap, "mus_musculus", "mus_musculus")$paralog_count
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_orthology_pair <- function(snapshot, a, b) {
  if (!inherits(snapshot, "mzlibr_orthology_snapshot")) {
    stop(mzlib_usage_error("snapshot must be a proteins_build_orthology() result."))
  }
  for (species in list(a, b)) {
    if (!is.character(species) || length(species) != 1L || is.na(species)) {
      stop(mzlib_usage_error("a and b must each be one species name, e.g. \"homo_sapiens\"."))
    }
  }
  ends <- sort(c(a, b), method = "radix")
  pairs <- snapshot$pairs
  hit <- pairs[pairs$species_a == ends[[1L]] & pairs$species_b == ends[[2L]], , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop(mzlib_usage_error(paste0(
      "('", a, "', '", b, "') is not a pair of this snapshot's species: ",
      paste(snapshot$species, collapse = ", "), "."
    )))
  }
  rownames(hit) <- NULL
  hit
}

#' Print an orthology snapshot
#'
#' @param x A [proteins_build_orthology()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_orthology_snapshot <- function(x, ...) {
  cat("<mzlibr_orthology_snapshot> Ensembl ", x$release, ", ", length(x$species), " species: ",
    paste(x$species, collapse = ", "), "\n",
    sep = ""
  )
  cat("  snapshot_id ", substr(x$snapshot_id, 1L, 12L), "... (record it with any result)\n", sep = "")
  between <- x$pairs[x$pairs$species_a != x$pairs$species_b, , drop = FALSE]
  if (nrow(between) > 0L) {
    cat("  orthologs: ", paste0(between$species_a, "~", between$species_b, " ", between$ortholog_count,
      collapse = ", "
    ), "\n", sep = "")
  }
  if (!is.null(x$pair_status)) {
    shown <- x$pair_status$status_counts[x$pair_status$status_counts > 0]
    cat("  ", x$pair_status$from_species, " in ", x$pair_status$to_species, ": ",
      paste0(names(shown), " ", shown, collapse = ", "), "\n",
      sep = ""
    )
  }
  if (!is.null(x$written)) {
    cat("  written: ", x$written$directory, " (", nrow(x$written$files), " files and manifest.json)\n", sep = "")
  }
  for (caveat in x$caveats) cat("  ! ", caveat, "\n", sep = "")
  invisible(x)
}
