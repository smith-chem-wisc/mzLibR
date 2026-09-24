# Protein databases: what each protein IS, which gene it belongs to, and whether a peptide
# identifies it.
#
# Three questions a search result cannot answer on its own, each answered by mzLib from the
# protein database you searched:
#
#   proteins_read()               what is this accession - organism, NCBI taxon, gene names,
#                                 length, mass, and on request GO terms and Ensembl gene links
#   proteins_resolve_genes()      which stable Ensembl gene is it, counted against a gene set the
#                                 caller supplies and pins, with one outcome per protein
#   proteins_classify_peptides()  does this peptide identify one protein, with I and L treated as
#                                 one residue
#
# Ported from pyMzLib's `proteins.py`, which decided the verbs, the wire fields and the caveats.
# What changes here is only the projection into R's idiom: data.frames, `NA` rather than `None`,
# list columns where a cell holds several values, `source_index` 1-based so that it indexes
# `files` directly, and an S3 class with a `print` method rather than a dataclass.
#
# All three take ONE database or MANY. Several are read in one bridge call, on stdin, with the
# thread count stated on the wire; there is no loop over databases here (the bridge's
# PARALLELISM.md). The wire picks the transport from the argument's shape, as pyMzLib does: one
# database travels as `--path` (plus `--contaminant`), several as `--paths-stdin` with a
# contaminant's line tagged. When a second list - accessions, or peptides - shares stdin with the
# paths, a line holding only `--` separates them (the bridge's BULK.md section 2).

# The tables proteins_read() can return. "proteins" alone is the default.
PROTEINS_TABLES <- c("proteins", "go_terms", "ensembl_genes")

# The line separating the database paths from a second list on stdin.
PROTEINS_SECTION <- "--"

# ---------------------------------------------------------------- argument shaping

# One path or several, as clean strings. NULL gives none.
proteins_paths <- function(value, what) {
  if (is.null(value)) {
    return(character(0))
  }
  if (!is.character(value)) {
    stop(mzlib_usage_error(paste0(
      "Every ", what, " must be a path given as a character string; got ", class(value)[1L], "."
    )))
  }
  paths <- trimws(unname(value))
  for (path in paths) {
    if (is.na(path) || !nzchar(path)) {
      stop(mzlib_usage_error(paste0("Every ", what, " path must be non-empty; got a blank or NA.")))
    }
    if (grepl("[\t\r\n]", path)) {
      stop(mzlib_usage_error(paste0("A ", what, " path contains a tab or newline: '", path, "'.")))
    }
    if (identical(path, PROTEINS_SECTION)) {
      stop(mzlib_usage_error(paste0("'", PROTEINS_SECTION, "' is not a usable ", what, " path.")))
    }
  }
  paths
}

# A list of accessions or peptides for stdin, one per line.
proteins_lines <- function(values, what) {
  if (!is.character(values)) {
    stop(mzlib_usage_error(paste0(
      what, " must be a character vector; got ", class(values)[1L], "."
    )))
  }
  lines <- trimws(unname(values))
  if (any(is.na(lines) | !nzchar(lines))) {
    stop(mzlib_usage_error(paste0(what, " contains a blank or NA entry.")))
  }
  bad <- grepl("[\t\r\n]", lines)
  if (any(bad)) {
    stop(mzlib_usage_error(paste0(
      "An entry in ", what, " contains a tab or newline: '", lines[bad][1L], "'."
    )))
  }
  lines
}

# The database options, and the stdin path lines when there are several databases.
proteins_database_args <- function(databases, contaminants) {
  targets <- proteins_paths(databases, "database")
  extra <- proteins_paths(contaminants, "contaminant database")
  everything <- c(targets, extra)
  if (length(everything) == 0L) {
    stop(mzlib_usage_error("At least one protein database is required, e.g. \"human.xml\"."))
  }
  keys <- normalizePath(everything, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") {
    keys <- tolower(keys)
  }
  if (anyDuplicated(keys) > 0L) {
    stop(mzlib_usage_error(paste0(
      "A database is listed twice: '", everything[anyDuplicated(keys)], "'."
    )))
  }

  if (length(everything) == 1L) {
    args <- c("--path", everything[[1L]])
    if (length(extra) == 1L) {
      args <- c(args, "--contaminant")
    }
    return(list(args = args, lines = NULL))
  }
  list(
    args = "--paths-stdin",
    lines = c(targets, if (length(extra) > 0L) paste0(extra, "\tcontaminant") else character(0))
  )
}

# The path lines and a second list, separated by a `--` line when both travel on stdin.
proteins_stdin <- function(path_lines, second) {
  if (is.null(path_lines) && is.null(second)) {
    return(NULL)
  }
  lines <- if (is.null(path_lines)) {
    second
  } else if (is.null(second)) {
    path_lines
  } else {
    c(path_lines, PROTEINS_SECTION, second)
  }
  paste0(paste(lines, collapse = "\n"), "\n")
}

proteins_build_read_request <- function(databases, tables, accessions, contaminants, sequences,
                                        threads, on_error) {
  if (!is.character(tables) || length(tables) == 0L || anyNA(tables) ||
    !all(tables %in% PROTEINS_TABLES)) {
    stop(mzlib_usage_error(paste0(
      "tables must be a non-empty selection from ", paste0("\"", PROTEINS_TABLES, "\"", collapse = ", "),
      "; got ", paste(deparse(tables), collapse = " "), "."
    )))
  }
  if (!is.logical(sequences) || length(sequences) != 1L || is.na(sequences)) {
    stop(mzlib_usage_error("sequences must be TRUE or FALSE."))
  }
  databases <- proteins_database_args(databases, contaminants)
  args <- c(
    "proteins", "read", databases$args,
    "--tables", paste(unique(tables), collapse = ","),
    "--threads", bulk_check_threads(threads),
    "--on-error", bulk_check_on_error(on_error)
  )
  if (isTRUE(sequences)) {
    args <- c(args, "--sequences")
  }
  filter <- NULL
  if (!is.null(accessions)) {
    filter <- proteins_lines(accessions, "accessions")
    if (length(filter) == 0L) {
      stop(mzlib_usage_error("accessions is empty; pass NULL to read every protein."))
    }
    args <- c(args, "--accessions-stdin")
  }
  list(args = args, stdin = proteins_stdin(databases$lines, filter))
}

proteins_single_path <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(mzlib_usage_error(paste0(name, " must be a single non-empty path, or NULL.")))
  }
  trimws(value)
}

proteins_build_resolve_request <- function(databases, gtf, gene_set, xref, contaminants, threads,
                                           on_error) {
  gtf <- proteins_single_path(gtf, "gtf")
  gene_set <- proteins_single_path(gene_set, "gene_set")
  if (is.null(gtf) == is.null(gene_set)) {
    stop(mzlib_usage_error(paste0(
      "Give exactly one of gtf (an Ensembl GTF, e.g. \"Homo_sapiens.GRCh38.116.gtf.gz\") or ",
      "gene_set (a compact table made from one). There is no default gene set: the release it ",
      "pins is part of the answer."
    )))
  }
  xref <- proteins_single_path(xref, "xref")
  databases <- proteins_database_args(databases, contaminants)
  args <- c(
    "genes", "resolve", databases$args,
    "--threads", bulk_check_threads(threads),
    "--on-error", bulk_check_on_error(on_error),
    if (!is.null(gtf)) c("--gtf", gtf) else c("--gene-set", gene_set),
    if (!is.null(xref)) c("--xref", xref)
  )
  list(args = args, stdin = proteins_stdin(databases$lines, NULL))
}

proteins_build_classify_request <- function(peptides, databases, contaminants, threads) {
  lines <- proteins_lines(peptides, "peptides")
  if (length(lines) == 0L) {
    stop(mzlib_usage_error("At least one peptide is required, e.g. \"PEPTIDEK\"."))
  }
  databases <- proteins_database_args(databases, contaminants)
  args <- c(
    "proteins", "classify-peptides", databases$args,
    "--threads", bulk_check_threads(threads)
  )
  list(args = args, stdin = proteins_stdin(databases$lines, lines))
}

# ---------------------------------------------------------------- parsing

# A wire table - `column_names` and a `columns` map - as a data.frame.
#
# A column whose cells are arrays (gene_names, accessions, evidence_codes, ...) is ALWAYS a list
# column of character vectors, empty ones included, so a cell with no values is character(0)
# rather than a row that silently vanished. `source_index` becomes 1-based.
proteins_parse_table <- function(block) {
  if (!is.list(block)) {
    return(NULL)
  }
  columns <- block[["columns"]]
  order <- wire_strings(block[["column_names"]])
  if (!is.list(columns)) {
    return(NULL)
  }
  if (length(order) == 0L) {
    order <- names(columns)
  }
  built <- lapply(order, function(name) {
    values <- columns[[name]]
    if (!is.list(values)) {
      return(logical(0))
    }
    if (any(vapply(values, is.list, logical(1L)))) {
      return(I(lapply(values, function(v) if (is.list(v)) wire_strings(v) else character(0))))
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
  table <- as.data.frame(built, stringsAsFactors = FALSE, optional = TRUE)
  if ("source_index" %in% names(table)) {
    table$source_index <- as.numeric(table$source_index) + 1
  }
  table
}

# A `{name: count}` map as a named numeric vector, in wire order.
proteins_parse_counts <- function(value) {
  if (!is.list(value) || length(value) == 0L) {
    return(proteins_named_numeric(character(0)))
  }
  counts <- vapply(value, function(v) if (wire_null(v)) NA_real_ else as.numeric(v), numeric(1L))
  names(counts) <- names(value)
  counts
}

proteins_named_numeric <- function(names) {
  out <- numeric(length(names))
  names(out) <- names
  out
}

# A flat wire object as a list, `null` fields as NA.
proteins_parse_object <- function(value) {
  if (!is.list(value) || is.null(names(value))) {
    return(NULL)
  }
  lapply(value, function(v) if (wire_null(v)) NA else v)
}

proteins_count <- function(data, name) {
  as.numeric(wire_field(data, name, "numeric", NA_real_))
}

proteins_parse_database <- function(data) {
  filter_count <- data[["accession_filter_count"]]
  not_found <- data[["accessions_not_found"]]
  structure(
    list(
      file_count = proteins_count(data, "file_count"),
      read_count = proteins_count(data, "read_count"),
      failed_count = proteins_count(data, "failed_count"),
      record_count = proteins_count(data, "record_count"),
      tables = wire_strings(data[["tables"]]),
      accession_filter_count = if (wire_null(filter_count)) NA_real_ else as.numeric(filter_count),
      # NULL (no filter was given) is a different answer from character(0) (every accession found).
      accessions_not_found = if (is.list(not_found)) wire_strings(not_found) else NULL,
      caveats = wire_strings(data[["caveats"]]),
      files = bulk_parse_files(data[["files"]]),
      column_names = wire_strings(data[["column_names"]]),
      proteins = proteins_parse_table(data),
      go_terms = proteins_parse_table(data[["go_terms"]]),
      ensembl_genes = proteins_parse_table(data[["ensembl_genes"]])
    ),
    class = "mzlibr_protein_database"
  )
}

proteins_parse_resolutions <- function(data) {
  structure(
    list(
      file_count = proteins_count(data, "file_count"),
      read_count = proteins_count(data, "read_count"),
      failed_count = proteins_count(data, "failed_count"),
      protein_count = proteins_count(data, "protein_count"),
      record_count = proteins_count(data, "record_count"),
      gene_set = proteins_parse_object(data[["gene_set"]]),
      xref = proteins_parse_object(data[["xref"]]),
      outcome_counts = proteins_parse_counts(data[["outcome_counts"]]),
      caveats = wire_strings(data[["caveats"]]),
      files = bulk_parse_files(data[["files"]]),
      column_names = wire_strings(data[["column_names"]]),
      resolutions = proteins_parse_table(data)
    ),
    class = "mzlibr_gene_resolutions"
  )
}

proteins_parse_classification <- function(data) {
  structure(
    list(
      file_count = proteins_count(data, "file_count"),
      peptide_count = proteins_count(data, "peptide_count"),
      target_protein_count = proteins_count(data, "target_protein_count"),
      decoy_proteins_ignored = proteins_count(data, "decoy_proteins_ignored"),
      i_and_l_equivalent = !isFALSE(data[["i_and_l_equivalent"]]),
      sharing_counts = proteins_parse_counts(data[["sharing_counts"]]),
      files = bulk_parse_files(data[["files"]]),
      column_names = wire_strings(data[["column_names"]]),
      peptides = proteins_parse_table(data)
    ),
    class = "mzlibr_peptide_classification"
  )
}

# ---------------------------------------------------------------- the public surface

#' Read protein databases: one row per protein, with GO terms and Ensembl genes on request
#'
#' Loads each database with mzLib's `ProteinDbLoader`, exactly as a search would but with no decoys
#' generated, and answers "what is this accession?": organism, NCBI taxon, gene names, length and
#' mass, plus - when you ask for them - its Gene Ontology terms and Ensembl gene links as long
#' tables.
#'
#' **A FASTA knows no GO terms and no Ensembl genes.** Its headers carry organism (`OS=`), taxon
#' (`OX=`) and gene name (`GN=`) and nothing else, so a FASTA contributes no `go_terms` or
#' `ensembl_genes` rows. That is the format's silence, not a biological absence, and `files` says
#' so in its `absent_fields` column. Read the UniProt XML of the same proteome for those.
#'
#' @param databases A UniProt XML (`.xml`) or FASTA (`.fasta`, `.fa`, `.faa`, `.fas`), each
#'   optionally gzipped, or a character vector of them. Read in the order given, in one bridge
#'   call.
#' @param tables Which tables to return, from `"proteins"`, `"go_terms"` and `"ensembl_genes"`.
#'   Default `"proteins"`. Ask for the other two when you want them: a whole human proteome has
#'   about twenty GO rows per protein.
#' @param accessions Keep only proteins whose accession is one of these - an **exact** match, so
#'   `"P04406-1"` does not find `"P04406"`. Applies to every table; misses are listed in
#'   `accessions_not_found`. Sent on stdin, so tens of thousands are fine. `NULL` keeps every
#'   protein.
#' @param contaminants Databases to load as contaminants; `is_contaminant` is `TRUE` on their rows.
#'   They come after `databases` in `files`.
#' @param sequences Add a `sequence` column to the `proteins` table.
#' @param threads Databases read at once. Default 1; `-1` means every core. A resource choice only:
#'   the result is identical at any value.
#' @param on_error `"fail"` (the default) raises on the first unreadable database, in input order;
#'   `"skip"` records the failure in that database's row of `files` and reads the rest.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait: a proteome XML takes a while.
#'
#' @return An `mzlibr_protein_database`. `file_count`, `read_count` and `failed_count` count
#'   database files; `record_count` counts the proteins that passed the accession filter;
#'   `accession_filter_count` counts the distinct accessions in the filter, and is `NA` with no
#'   filter; `accessions_not_found` is `NULL` with no filter. `tables` names what was returned,
#'   `caveats` the traps in the result as a whole, and `files` has one row per input (its
#'   `source_index`, `contaminant`, `protein_count`, `decoy_count`, `record_count`, `caveats`,
#'   `absent_fields` and any `error_*`).
#'
#'   `proteins` is a data.frame with one row per protein, or `NULL` when it was not requested:
#'   `length` in residues and `monoisotopic_mass` in Da (of the unmodified sequence as written,
#'   plus one water; `NA` when the sequence holds a letter with no defined mass). `gene_names` and
#'   `ensembl_gene_ids` are list columns. `go_terms` and `ensembl_genes` are data.frames, or `NULL`
#'   when not requested. **`source_index` is 1-based** in all of them, so
#'   `db$files[db$proteins$source_index, ]` finds each row's database.
#'
#' @section Masses are of the precursor as written:
#'
#' `monoisotopic_mass` is the unmodified sequence plus one water - the initiator methionine, signal
#' peptide and propeptide included - so it is not the mass of the mature protein.
#'
#' @section No decoys, no variants:
#'
#' A database is read as written. A decoy already in the file (an accession starting `DECOY`) is
#' kept and flagged `is_decoy`; none are generated. A UniProt XML that records a **genotype** still
#' has those variants applied by mzLib, which renames the accession (`P38936_C117Y`); the file's
#' `caveats` say how many.
#'
#' @seealso [proteins_resolve_genes()], [proteins_classify_peptides()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' db <- proteins_read(c("human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"),
#'   contaminants = "contaminants.fasta",
#'   tables = c("proteins", "go_terms", "ensembl_genes"))
#' db
#' db$proteins[, c("accession", "organism", "ncbi_taxonomy_id", "length", "monoisotopic_mass")]
#' head(db$go_terms[, c("accession", "go_id", "aspect", "term_name")])
#' # A FASTA says what it cannot know, rather than leaving the tables silently empty:
#' db$files[, c("path", "file_type", "contaminant")]
#' db$files$absent_fields[[2]]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec proteins.read
#' @export
proteins_read <- function(databases, tables = "proteins", accessions = NULL, contaminants = NULL,
                          sequences = FALSE, threads = 1, on_error = "fail", timeout = NULL) {
  request <- proteins_build_read_request(
    databases, tables, accessions, contaminants, sequences, threads, on_error
  )
  proteins_parse_database(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Resolve every protein to stable Ensembl gene ids, against a gene set you pin
#'
#' Wraps mzLib's `EnsemblGeneResolver`. The gene links come from the database itself (UniProt's
#' Ensembl `dbReference`s), and each is **counted against the gene set you pass**: a gene in the set
#' counts; a gene outside it (an ALT haplotype, a patch) is dropped from `n_genes` but counted in
#' `off_primary_genes`, and a protein with only such genes is `off_primary_only`. Every protein gets
#' exactly one `outcome`: `resolved`, `multi_gene`, `off_primary_only`, `not_in_source`,
#' `unrecognized_accession` or `contaminant_not_mapped`.
#'
#' **You supply the gene set, and you should pin it.** Gene ids, versions, symbols and the very set
#' of genes on the primary assembly change between Ensembl releases, so a resolution means something
#' only relative to one release. Use Ensembl's **primary-assembly** GTF
#' (`Species.Assembly.Release.gtf.gz`, not the `chr_patch_hapl_scaff` one), keep Ensembl's file name
#' so the release is recorded, and keep `gene_set$sha256` with your results. Nothing is downloaded
#' or defaulted here.
#'
#' @param databases The search database or databases, UniProt XML for real answers. A FASTA carries
#'   no gene links, so every FASTA protein is `not_in_source`; a caveat says so.
#' @param gtf An Ensembl GTF, plain or `.gz`. Only its `gene` rows are read. Give exactly one of
#'   `gtf` and `gene_set`.
#' @param gene_set Instead of `gtf`, a compact gene table written by mzLib's
#'   `EnsemblGeneSetWriter` (about 0.5 MB for human, against a 141 MB GTF). It carries the GTF's
#'   provenance, so rows are keyed exactly as against the GTF.
#' @param xref Optional: Ensembl's `Species.Assembly.Release.uniprot.tsv.gz`, for a second opinion.
#'   Each gene row then says whether Ensembl agrees, and a gene only Ensembl links gets its own row
#'   (`source == "ensembl_xref"`). Without it, `ensembl_xref_agrees` is `NA` - unknown, not false.
#' @param contaminants Databases to load as contaminants. Their proteins are
#'   `contaminant_not_mapped`: never mapped, and never silently dropped.
#' @param threads Databases read (and hashed) at once. Default 1; `-1` means every core. Same answer
#'   at any value.
#' @param on_error `"fail"` (the default) or `"skip"`, as for [proteins_read()].
#' @param timeout Seconds to allow, or `NULL` (the default) to wait: a gzipped GTF takes tens of
#'   seconds.
#'
#' @return An `mzlibr_gene_resolutions`. `file_count`, `read_count` and `failed_count` count
#'   database files; `protein_count` counts the proteins resolved; `record_count` counts the rows of
#'   `resolutions` - at least one per protein, one per gene for `multi_gene`, plus any
#'   `ensembl_xref` rows. `outcome_counts` is a named vector: outcome to number of proteins, every
#'   outcome present, zeros included. `gene_set` is a list (`source_file_name`, `sha256`,
#'   `release`, `genome_build`, `genebuild_last_updated`, `gene_count`), and `xref` one like it, or
#'   `NULL` without `xref`. `files` has one row per input, with its `search_database_sha256` (of
#'   the decompressed database).
#'
#'   `resolutions` is a data.frame whose columns after `source_index` (1-based) and `source_path`
#'   are mzLib's own `GeneResolutionTsv` schema. `n_genes` counts genes on the gene set's assembly;
#'   `off_primary_genes` counts linked genes outside it. `gene_id`, `versioned_gene_id` and
#'   `gene_biotype` are `NA` when the row has no gene; `isoform` is `NA` without a `-N` suffix.
#'
#' @section One row per gene, never a pick:
#'
#' A `multi_gene` protein has one row per gene, and a protein with no gene one outcome row. No cell
#' joins several genes, and nothing chooses among them for you.
#'
#' @seealso [proteins_read()], [proteins_classify_peptides()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' genes <- proteins_resolve_genes(c("human_subset.xml", "human_extra.fasta"),
#'   contaminants = "contaminants.fasta",
#'   gtf = "Homo_sapiens.GRCh38.116.gtf", xref = "Homo_sapiens.GRCh38.116.uniprot.tsv")
#' genes
#' genes$gene_set[c("release", "genome_build", "sha256")]
#' genes$outcome_counts
#' genes$resolutions[, c("accession", "outcome", "n_genes", "gene_id", "gene_symbol")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec genes.resolve
#' @export
proteins_resolve_genes <- function(databases, gtf = NULL, gene_set = NULL, xref = NULL,
                                   contaminants = NULL, threads = 1, on_error = "fail",
                                   timeout = NULL) {
  request <- proteins_build_resolve_request(
    databases, gtf, gene_set, xref, contaminants, threads, on_error
  )
  proteins_parse_resolutions(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Classify peptides by how widely they are shared across protein databases, with I = L
#'
#' Wraps mzLib's `PeptideUniquenessClassifier.Classify`. Each peptide is one of:
#'
#' `"Unique"` - every protein containing it has the **same sequence** (identical entries under two
#' accessions count as one; both are listed). `"SharedWithinGene"` - several distinct sequences, all
#' sharing a gene: isoforms, so the peptide supports the gene, not an isoform.
#' `"SharedAcrossGenes"` - sequences with no gene in common. `"NotInDatabase"` - no target protein
#' contains it.
#'
#' A protein contains a peptide when its sequence contains it **anywhere, whatever the protease** -
#' deliberately conservative, so a peptide called `"Unique"` cannot be explained by another entry
#' at a site the search's cleavage rules happened to skip. **I and L are the same residue**
#' throughout, because a mass spectrometer cannot tell them apart.
#'
#' @param peptides A character vector of unmodified base sequences, upper case (`"PEPTIDEK"`), one
#'   result row each in the same order, duplicates included. Sent on stdin. Strip modifications
#'   first: mzLib refuses a sequence with a bracketed modification, or in lower case
#'   (`"peptidek"`), rather than guess.
#' @param databases The search database or databases: UniProt XML or FASTA, optionally gzipped.
#' @param contaminants Contaminant databases, searched alongside: they are real sequences in the
#'   search space, so a peptide a contaminant shares with a target is shared.
#' @param threads Databases read at once (the classification itself is one pass). Default 1; `-1`
#'   means every core. Same answer at any value.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait.
#'
#' @return An `mzlibr_peptide_classification`. `file_count` counts the database files searched;
#'   `peptide_count` counts the peptides (rows); `target_protein_count` counts the proteins
#'   searched, contaminants included; `decoy_proteins_ignored` counts the decoy proteins, never
#'   searched. `sharing_counts` is a named vector: class to number of peptides, every class
#'   present. `i_and_l_equivalent` is always `TRUE`, carried so the rule travels with the result.
#'
#'   `peptides` is a data.frame with one row per peptide, in the order given: `peptide` (exactly as
#'   given, not I/L-folded), `sharing`, `accession_count` (the target proteins containing it), and
#'   the list columns `accessions` and `shared_gene_keys` (the gene keys common to all of them, such
#'   as `"ensembl:ENSG00000111640"`; empty for `"NotInDatabase"` and `"SharedAcrossGenes"`).
#'
#' @section Why there is no on_error:
#'
#' Every database is part of one search space. Skipping one that failed to load would report the
#' peptides it contains as unique, so any failure raises.
#'
#' @seealso [proteins_read()], [proteins_resolve_genes()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' peptides <- c("VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK")
#' calls <- proteins_classify_peptides(peptides, c("human_subset.xml", "human_extra.fasta"),
#'   contaminants = "contaminants.fasta")
#' calls
#' calls$peptides[, c("peptide", "sharing", "accession_count")]
#' # YLYEIAR is in human and bovine albumin, which share no gene:
#' calls$peptides$accessions[[4]]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec proteins.classify-peptides
#' @export
proteins_classify_peptides <- function(peptides, databases, contaminants = NULL, threads = 1,
                                       timeout = NULL) {
  request <- proteins_build_classify_request(peptides, databases, contaminants, threads)
  proteins_parse_classification(
    bridge_invoke(request$args, stdin = request$stdin, timeout = timeout)
  )
}

# ---------------------------------------------------------------- print methods

proteins_print_files <- function(x) {
  failed <- bulk_failed(x)
  for (i in seq_len(nrow(failed))) {
    cat("  ! ", basename(failed$path[[i]]), ": ", failed$error_message[[i]], "\n", sep = "")
  }
  for (caveat in x$caveats) {
    cat("  ! ", caveat, "\n", sep = "")
  }
}

#' Print a protein database read
#'
#' @param x A [proteins_read()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_protein_database <- function(x, ...) {
  cat("<mzlibr_protein_database> ", format(x$read_count), " of ", format(x$file_count),
    " database(s) read, ", format(x$record_count), " proteins\n",
    sep = ""
  )
  cat("  tables: ", paste(x$tables, collapse = ", "), "\n", sep = "")
  for (name in c("go_terms", "ensembl_genes")) {
    if (!is.null(x[[name]])) cat("  ", name, ": ", nrow(x[[name]]), " rows\n", sep = "")
  }
  if (!is.null(x$accessions_not_found) && length(x$accessions_not_found) > 0L) {
    cat("  ! accessions not found (matching is exact): ",
      paste(x$accessions_not_found, collapse = ", "), "\n",
      sep = ""
    )
  }
  absent <- x$files$absent_fields
  if (!is.null(absent)) {
    silent <- vapply(absent, function(a) length(a) > 0L, logical(1L))
    if (any(silent)) {
      cat("  ", sum(silent), " database(s) cannot say: ",
        paste(unique(unlist(absent[silent])), collapse = ", "), " (see files$absent_fields)\n",
        sep = ""
      )
    }
  }
  proteins_print_files(x)
  invisible(x)
}

#' Print a gene resolution
#'
#' @param x A [proteins_resolve_genes()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_gene_resolutions <- function(x, ...) {
  cat("<mzlibr_gene_resolutions> ", format(x$protein_count), " proteins, ",
    format(x$record_count), " rows\n",
    sep = ""
  )
  if (!is.null(x$gene_set)) {
    cat("  gene set: ", x$gene_set$source_file_name,
      if (!is.na(x$gene_set$release)) paste0(", Ensembl release ", x$gene_set$release) else ", no release in its name",
      "\n",
      sep = ""
    )
  }
  if (length(x$outcome_counts) > 0L) {
    shown <- x$outcome_counts[x$outcome_counts > 0]
    cat("  outcomes: ", paste0(names(shown), " ", shown, collapse = ", "), "\n", sep = "")
  }
  if (is.null(x$xref)) {
    cat("  no xref: ensembl_xref_agrees is NA (unknown, not false)\n")
  }
  proteins_print_files(x)
  invisible(x)
}

#' Print a peptide classification
#'
#' @param x A [proteins_classify_peptides()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_peptide_classification <- function(x, ...) {
  cat("<mzlibr_peptide_classification> ", format(x$peptide_count), " peptides against ",
    format(x$target_protein_count), " proteins in ", format(x$file_count), " database(s)\n",
    sep = ""
  )
  if (length(x$sharing_counts) > 0L) {
    cat("  ", paste0(names(x$sharing_counts), " ", x$sharing_counts, collapse = ", "), "\n", sep = "")
  }
  if (isTRUE(x$i_and_l_equivalent)) {
    cat("  I and L are one residue for matching\n")
  }
  proteins_print_files(x)
  invisible(x)
}
