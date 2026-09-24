# Protein databases: proteins_read(), proteins_resolve_genes(), proteins_classify_peptides().
#
# Ported from pyMzLib's test_proteins.py. Offline except for the live tests at the end. This
# layer's job is to shape arguments - one database or many, contaminants marked, a second list
# after a `--` line on stdin - and to keep "this format cannot say" apart from "there is nothing".
# The biology is mzLib's, and the bridge's C# suite covers it. The four payloads were recorded from
# the real bridge by pyMzLib and are shared verbatim.

proteins_fixture <- function(name) {
  mz$json_parse(paste(readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
}

proteins_human <- function() mz$proteins_parse_database(proteins_fixture("proteins_read_human.json"))
proteins_genes <- function() mz$proteins_parse_resolutions(proteins_fixture("genes_resolve_human.json"))
proteins_calls <- function() {
  mz$proteins_parse_classification(proteins_fixture("proteins_classify_peptides.json"))
}

# ---------------------------------------------------------------- proteins_read: the payload

test_that("a protein read reports organism and taxon per accession", {
  db <- proteins_human()
  expect_true(inherits(db, "mzlibr_protein_database"))
  taxa <- setNames(db$proteins$ncbi_taxonomy_id, db$proteins$accession)
  expect_identical(taxa[["P04406"]], "9606") # UniProt XML dbReference
  expect_identical(taxa[["Q9Z0X1"]], "10090") # FASTA OX=
  organisms <- setNames(db$proteins$organism, db$proteins$accession)
  expect_identical(organisms[["P02769"]], "Bos taurus")
})

test_that("every input is kept in order with its role, and source_index is 1-based", {
  db <- proteins_human()
  expect_identical(db$files$file_type, c("UniProtXml", "Fasta", "Fasta", "Fasta"))
  expect_identical(db$files$contaminant, c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(db$files$source_index, c(1, 2, 3, 4))
  expect_identical(db$proteins$source_index, sort(db$proteins$source_index))
  expect_identical(db$proteins$accession[db$proteins$is_contaminant], "P02769")
  # 1-based, so it indexes files directly.
  expect_identical(db$files$path[db$proteins$source_index], db$proteins$source_path)
})

test_that("a FASTA says GO and Ensembl are absent rather than empty", {
  db <- proteins_human()
  expect_true(all(c("go_terms", "ensembl_genes") %in% db$files$absent_fields[[2L]]))
  expect_true(any(grepl("FASTA carries no GO", db$files$caveats[[2L]], fixed = TRUE)))
  expect_identical(db$files$absent_fields[[1L]], character(0))
  expect_identical(unique(db$go_terms$source_index), 1, info = "only the XML contributes GO rows")
})

test_that("GO terms are a long table with aspect and evidence as list columns", {
  go <- proteins_human()$go_terms
  expect_true(is.data.frame(go))
  expect_identical(nrow(go), 47L)
  expect_true(all(startsWith(go$go_id, "GO:")))
  expect_true(all(go$aspect %in% c("BiologicalProcess", "CellularComponent", "MolecularFunction", "Unknown")))
  expect_true(is.list(go$evidence_codes))
  expect_true(is.character(go$evidence_codes[[1L]]))
})

test_that("a list cell with no values stays a row, as character(0)", {
  # The FASTA proteins have no Ensembl gene ids: an empty list cell, which a naive unlist() would
  # drop, leaving the column shorter than the table.
  db <- proteins_human()
  expect_identical(nrow(db$proteins), 8L)
  expect_true(is.list(db$proteins$ensembl_gene_ids))
  albumin <- which(db$proteins$accession == "P02768")
  expect_identical(db$proteins$ensembl_gene_ids[[albumin]], character(0))
  expect_identical(db$proteins$gene_names[[1L]][[1L]], "GAPDH")
})

test_that("Ensembl links carry stable and versioned ids, and unversioned is NA not 0", {
  ens <- proteins_human()$ensembl_genes
  gapdh <- ens[ens$accession == "P04406", ]
  expect_identical(unique(gapdh$gene_id), "ENSG00000111640")
  expect_identical(unique(gapdh$versioned_gene_id), "ENSG00000111640.15")
  psca <- ens[ens$accession == "O43653", ]
  expect_true(nrow(psca) > 0L)
  expect_true(all(is.na(psca$gene_version)))
})

test_that("mass and length are numbers", {
  db <- proteins_human()
  row <- db$proteins[db$proteins$accession == "P04406", ]
  expect_identical(row$length, 335)
  expect_true(abs(row$monoisotopic_mass - 36030.398) < 0.01)
})

test_that("a filtered read names the accessions it did not find", {
  db <- mz$proteins_parse_database(proteins_fixture("proteins_read_filtered.json"))
  expect_identical(db$proteins$accession, "O43653")
  expect_identical(db$accessions_not_found, "P04406-1")
  expect_identical(db$accession_filter_count, 2)
  expect_true(is.null(db$ensembl_genes))
})

test_that("with no filter, accessions_not_found is NULL and the count NA", {
  db <- proteins_human()
  expect_true(is.null(db$accessions_not_found))
  expect_true(is.na(db$accession_filter_count))
})

test_that("a skipped file carries its error", {
  db <- mz$proteins_parse_database(list(
    file_count = 1, read_count = 0, failed_count = 1, record_count = 0,
    files = list(list(source_index = 0, path = "x.fasta", error = list(
      kind = "usage", type = "usage", message = "Protein database not found: 'x.fasta'."
    )))
  ))
  expect_identical(db$files$error_kind, "usage")
  expect_true(grepl("x.fasta", db$files$error_message, fixed = TRUE))
  expect_true(is.null(db$proteins))
})

test_that("the print method names what a FASTA cannot say", {
  printed <- paste(capture.output(print(proteins_human())), collapse = "\n")
  expect_true(grepl("8 proteins", printed, fixed = TRUE))
  expect_true(grepl("go_terms, ensembl_genes", printed, fixed = TRUE))
})

# ---------------------------------------------------------------- proteins_read: what is sent

read_request <- function(databases, tables = "proteins", accessions = NULL, contaminants = NULL,
                         sequences = FALSE, threads = 1, on_error = "fail") {
  mz$proteins_build_read_request(databases, tables, accessions, contaminants, sequences, threads, on_error)
}

test_that("one database goes on argv", {
  request <- read_request("human.xml")
  expect_identical(request$args[1:4], c("proteins", "read", "--path", "human.xml"))
  expect_false("--contaminant" %in% request$args)
  expect_identical(request$args[which(request$args == "--tables") + 1L], "proteins")
  expect_identical(request$args[which(request$args == "--threads") + 1L], "1")
  expect_identical(request$args[which(request$args == "--on-error") + 1L], "fail")
  expect_true(is.null(request$stdin))
})

test_that("a lone contaminant database is marked", {
  request <- read_request(NULL, contaminants = "crap.fasta")
  expect_identical(request$args[3:5], c("--path", "crap.fasta", "--contaminant"))
})

test_that("many databases go on stdin with contaminants tagged", {
  request <- read_request(c("a.xml", "b.fasta"), contaminants = "c.fasta", threads = -1, on_error = "skip")
  expect_true("--paths-stdin" %in% request$args)
  expect_identical(request$stdin, "a.xml\nb.fasta\nc.fasta\tcontaminant\n")
  expect_identical(request$args[which(request$args == "--threads") + 1L], "-1")
  expect_identical(request$args[which(request$args == "--on-error") + 1L], "skip")
})

test_that("accessions follow the paths after a -- line", {
  request <- read_request(c("a.xml", "b.fasta"), accessions = c(" P04406 ", "Q9Z0X1"),
    tables = c("proteins", "go_terms", "ensembl_genes"), sequences = TRUE)
  expect_identical(request$stdin, "a.xml\nb.fasta\n--\nP04406\nQ9Z0X1\n")
  expect_true("--accessions-stdin" %in% request$args)
  expect_true("--sequences" %in% request$args)
  expect_identical(request$args[which(request$args == "--tables") + 1L], "proteins,go_terms,ensembl_genes")
})

test_that("with one database the accessions are all of stdin", {
  expect_identical(read_request("a.xml", accessions = "P04406")$stdin, "P04406\n")
})

test_that("bad read arguments fail before a process starts", {
  cases <- list(
    list(list(databases = character(0)), "At least one protein database"),
    list(list(databases = ""), "non-empty"),
    list(list(databases = c("a.xml", "a.xml")), "listed twice"),
    list(list(databases = "a\tb.xml"), "tab or newline"),
    list(list(databases = "--"), "not a usable"),
    list(list(databases = 42), "character string"),
    list(list(databases = "a.xml", tables = c("proteins", "peptides")), "tables must be"),
    list(list(databases = "a.xml", tables = character(0)), "tables must be"),
    list(list(databases = "a.xml", threads = 0), "threads must be"),
    list(list(databases = "a.xml", threads = TRUE), "threads must be"),
    list(list(databases = "a.xml", on_error = "ignore"), "on_error"),
    list(list(databases = "a.xml", accessions = character(0)), "accessions is empty"),
    list(list(databases = "a.xml", accessions = c("P04406", " ")), "blank"),
    list(list(databases = "a.xml", accessions = 42), "character vector"),
    list(list(databases = "a.xml", accessions = "P0\n4406"), "tab or newline"),
    list(list(databases = "a.xml", sequences = NA), "sequences must be")
  )
  for (case in cases) {
    expect_error(do.call(read_request, case[[1L]]), class = "mzlib_usage_error", contains = case[[2L]])
  }
})

test_that("the exported function refuses bad arguments before looking for a bridge", {
  old <- options(mzlibr.bridge = file.path(tempdir(), "no-such-bridge"))
  on.exit(options(old), add = TRUE)
  expect_error(proteins_read("a.xml", tables = "nope"), class = "mzlib_usage_error", contains = "tables")
  expect_error(proteins_resolve_genes("a.xml"), class = "mzlib_usage_error", contains = "exactly one of gtf")
  expect_error(proteins_classify_peptides(character(0), "a.xml"), class = "mzlib_usage_error")
})

# ---------------------------------------------------------------- proteins_resolve_genes

test_that("a resolution pins the release and every input by hash", {
  genes <- proteins_genes()
  expect_true(inherits(genes, "mzlibr_gene_resolutions"))
  expect_identical(genes$gene_set$release, "116")
  expect_identical(genes$gene_set$genome_build, "GRCh38.p14")
  expect_identical(nchar(genes$gene_set$sha256), 64L)
  expect_identical(genes$xref$release, "116")
  expect_true(all(nchar(genes$files$search_database_sha256) == 64L))
  expect_identical(unique(genes$resolutions$gene_set_sha256), genes$gene_set$sha256)
})

test_that("every outcome is counted and each protein has one", {
  genes <- proteins_genes()
  expect_identical(sort(names(genes$outcome_counts)), sort(c(
    "resolved", "multi_gene", "off_primary_only", "not_in_source", "unrecognized_accession",
    "contaminant_not_mapped"
  )))
  expect_identical(sum(genes$outcome_counts), genes$protein_count)
  rows <- genes$resolutions
  outcome <- setNames(rows$outcome, rows$accession)
  expect_identical(outcome[["P04406"]], "resolved")
  expect_identical(rows$gene_symbol[rows$accession == "P04406"], "GAPDH")
  expect_true(isTRUE(rows$ensembl_xref_agrees[rows$accession == "P04406"]))
  expect_identical(outcome[["O43653"]], "off_primary_only")
  expect_identical(outcome[["P02768"]], "not_in_source")
  expect_identical(outcome[["P02769"]], "contaminant_not_mapped")
  # No gene on the row: NA, which is "unknown", never FALSE.
  expect_true(is.na(rows$ensembl_xref_agrees[rows$accession == "P02768"]))
  expect_true(is.na(rows$isoform[[1L]]))
})

test_that("a resolution warns that FASTA proteins cannot resolve", {
  expect_true(any(grepl("came from a FASTA", proteins_genes()$caveats, fixed = TRUE)))
})

resolve_request <- function(databases = "a.xml", gtf = NULL, gene_set = NULL, xref = NULL,
                            contaminants = NULL, threads = 1, on_error = "fail") {
  mz$proteins_build_resolve_request(databases, gtf, gene_set, xref, contaminants, threads, on_error)
}

test_that("resolve sends the gene set and the xref", {
  request <- resolve_request(gtf = "g.gtf.gz", xref = "x.tsv.gz")
  expect_identical(request$args[1:2], c("genes", "resolve"))
  expect_identical(request$args[which(request$args == "--gtf") + 1L], "g.gtf.gz")
  expect_identical(request$args[which(request$args == "--xref") + 1L], "x.tsv.gz")
  expect_false("--gene-set" %in% request$args)
})

test_that("resolve accepts a compact gene set", {
  request <- resolve_request(c("a.xml", "b.xml"), gene_set = "genes.tsv")
  expect_identical(request$args[which(request$args == "--gene-set") + 1L], "genes.tsv")
  expect_identical(request$stdin, "a.xml\nb.xml\n")
})

test_that("bad resolve arguments fail before a process starts", {
  expect_error(resolve_request(), class = "mzlib_usage_error", contains = "exactly one of gtf")
  expect_error(resolve_request(gtf = "g.gtf", gene_set = "s.tsv"), class = "mzlib_usage_error",
    contains = "exactly one of gtf")
  expect_error(resolve_request(gtf = " "), class = "mzlib_usage_error", contains = "gtf must be")
  expect_error(resolve_request(gtf = "g.gtf", xref = ""), class = "mzlib_usage_error", contains = "xref must be")
})

# ---------------------------------------------------------------- proteins_classify_peptides

test_that("each sharing class has a real example", {
  calls <- proteins_calls()
  expect_true(inherits(calls, "mzlibr_peptide_classification"))
  sharing <- setNames(calls$peptides$sharing, calls$peptides$peptide)
  expect_identical(sharing[["VGVNGFGR"]], "Unique")
  expect_identical(sharing[["ALSEQINIFFDYSGR"]], "SharedWithinGene") # DYNC1I2 isoforms
  expect_identical(sharing[["YLYEIAR"]], "SharedAcrossGenes") # human and bovine albumin
  expect_identical(sharing[["AEFVEVTK"]], "Unique") # bovine only: a contaminant
  expect_identical(sharing[["PEPTIDEK"]], "NotInDatabase")
  expect_identical(sum(calls$sharing_counts), calls$peptide_count)
  expect_identical(calls$peptide_count, 6)
})

test_that("I and L are one residue, and the peptide is echoed as given", {
  calls <- proteins_calls()
  expect_identical(calls$peptides$peptide[[2L]], "LVLNGNPLTLFQER") # GAPDH has LVINGNPITIFQER
  expect_identical(calls$peptides$accessions[[2L]], "P04406")
  expect_true(calls$i_and_l_equivalent)
})

test_that("shared across genes lists both proteins and no common gene", {
  calls <- proteins_calls()
  expect_identical(calls$peptides$accessions[[4L]], c("P02768", "P02769"))
  expect_identical(calls$peptides$shared_gene_keys[[4L]], character(0))
  expect_identical(calls$peptides$accessions[[6L]], character(0))
})

test_that("peptides follow the paths on stdin, with no on-error", {
  request <- mz$proteins_build_classify_request(c("PEPTIDEK", " VGVNGFGR"), "a.xml", "c.fasta", 2)
  expect_identical(request$args[1:3], c("proteins", "classify-peptides", "--paths-stdin"))
  expect_identical(request$stdin, "a.xml\nc.fasta\tcontaminant\n--\nPEPTIDEK\nVGVNGFGR\n")
  expect_false("--on-error" %in% request$args, info = "there is no skip mode for one search space")
})

test_that("with one database the peptides are all of stdin", {
  request <- mz$proteins_build_classify_request("PEPTIDEK", "a.xml", NULL, 1)
  expect_identical(request$args[3:4], c("--path", "a.xml"))
  expect_identical(request$stdin, "PEPTIDEK\n")
})

test_that("bad peptides fail before a process starts", {
  expect_error(mz$proteins_build_classify_request(character(0), "a.xml", NULL, 1),
    class = "mzlib_usage_error", contains = "At least one peptide")
  expect_error(mz$proteins_build_classify_request(c("PEPTIDEK", ""), "a.xml", NULL, 1),
    class = "mzlib_usage_error", contains = "blank")
  expect_error(mz$proteins_build_classify_request(list("PEPTIDEK"), "a.xml", NULL, 1),
    class = "mzlib_usage_error", contains = "character vector")
})

# ---------------------------------------------------------------- live

# The databases pyMzLib recorded these fixtures from. The proteins verbs are new in the bridge
# pyMzLib 0.2.0 publishes, so a bridge that does not list them is skipped, not failed.
proteins_live <- function(verb) {
  skip_if_no_bridge()
  root <- Sys.getenv("MZLIB_PROTEINS_FILES")
  skip_if(!nzchar(root) || !file.exists(file.path(root, "human_subset.xml")),
    "MZLIB_PROTEINS_FILES does not hold pyMzLib's tests/fixtures/proteins")
  verbs <- mzlibr_bridge_version()$verbs
  skip_if(is.null(verbs) || !verb %in% verbs, paste0("the bridge does not dispatch '", verb, "'"))
  function(name) file.path(root, name)
}

test_that("live: proteins read matches the recording", {
  at <- proteins_live("proteins read")
  db <- proteins_read(at(c("human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta")),
    contaminants = at("contaminants.fasta"), tables = c("proteins", "go_terms", "ensembl_genes"))
  recorded <- proteins_human()
  expect_identical(db$proteins$accession, recorded$proteins$accession)
  expect_identical(db$proteins$ncbi_taxonomy_id, recorded$proteins$ncbi_taxonomy_id)
  expect_identical(nrow(db$go_terms), nrow(recorded$go_terms))
  expect_identical(db$files$absent_fields, recorded$files$absent_fields)
})

test_that("live: a filtered read names the accession it did not find", {
  at <- proteins_live("proteins read")
  db <- proteins_read(at("human_subset.xml"), accessions = c("O43653", "P04406-1"),
    tables = c("proteins", "go_terms"))
  expect_identical(db$proteins$accession, "O43653")
  expect_identical(db$accessions_not_found, "P04406-1")
})

test_that("live: genes resolve matches the recording", {
  at <- proteins_live("genes resolve")
  genes <- proteins_resolve_genes(at(c("human_subset.xml", "human_extra.fasta")),
    contaminants = at("contaminants.fasta"), gtf = at("Homo_sapiens.GRCh38.116.gtf"),
    xref = at("Homo_sapiens.GRCh38.116.uniprot.tsv"))
  recorded <- proteins_genes()
  expect_identical(genes$outcome_counts, recorded$outcome_counts)
  expect_identical(genes$resolutions$outcome, recorded$resolutions$outcome)
  # Not the sha256: a checkout that rewrote line endings changes the bytes, and so the hash.
  expect_identical(genes$gene_set$release, recorded$gene_set$release)
})

test_that("live: classify-peptides matches the recording", {
  at <- proteins_live("proteins classify-peptides")
  calls <- proteins_classify_peptides(
    c("VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"),
    at(c("human_subset.xml", "human_extra.fasta")), contaminants = at("contaminants.fasta"))
  recorded <- proteins_calls()
  expect_identical(calls$peptides$sharing, recorded$peptides$sharing)
  expect_identical(calls$peptides$accessions, recorded$peptides$accessions)
})
