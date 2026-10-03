# Gene Ontology for protein groups: proteins_annotate_go(), proteins_update_go().
#
# Ported from pyMzLib's test_proteins_go.py and test_proteins_go_live.py. The payloads were recorded
# from the real bridge by pyMzLib over mzLib's own PXD036557 MetaMorpheus 1.1.11 table and are shared
# verbatim. This layer shapes arguments and projects the table; the annotation is mzLib's.

go_fixture <- function(name) {
  parsed <- mz$json_parse(paste(readLines(fixture_path(name), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}
go_recorded <- function() mz$proteins_parse_go_annotations(go_fixture("proteins_annotate_go_pxd036557.json"))
go_rows <- function(go, group) go$annotations[go$annotations$protein_group == group, , drop = FALSE]

# ---------------------------------------------------------------- the payload

test_that("the statuses are mzLib's four", {
  expect_identical(mz$PROTEINS_ANNOTATION_STATUSES, c("annotated", "no_go_terms", "no_entry", "contaminant"))
})

test_that("every non-decoy group is annotated and the decoy skipped", {
  go <- go_recorded()
  expect_true(inherits(go, "mzlibr_go_annotations"))
  expect_identical(c(go$table_row_count, go$decoy_group_count, go$group_count), c(6, 1, 5))
  expect_identical(c(go$row_count, go$returned_count), c(563, 563))
  expect_identical(nrow(go$annotations), 563L)
  expect_false(go$truncated)
  expect_identical(unique(go$annotations$protein_group),
                   c("P68363", "P05141", "P0C0S5|Q71UI9", "P02769", "P63104"))
})

test_that("a contaminant is one term-less row that says why", {
  albumin <- go_rows(go_recorded(), "P02769")
  expect_identical(nrow(albumin), 1L)
  expect_identical(albumin$annotation_status, "contaminant")
  expect_true(is.na(albumin$go_id) && is.na(albumin$aspect))
  expect_true(is.na(albumin$inherited) && is.na(albumin$propagated))
  expect_identical(albumin$n_with, 0)
  expect_identical(albumin$accession_used[[1L]], character(0))
})

test_that("the union keeps both histones, and the consensus is a filter", {
  histones <- go_rows(go_recorded(), "P0C0S5|Q71UI9")
  expect_identical(nrow(histones), 106L)
  expect_identical(sum(histones$n_with == histones$n_members), 57L)
  only_one <- histones[histones$go_id %in% "GO:0000791", ]
  expect_identical(only_one$go_name, "euchromatin")
  expect_identical(only_one$accession_used[[1L]], "P0C0S5")
})

test_that("evidence is kept per member, as a named list", {
  histones <- go_rows(go_recorded(), "P0C0S5|Q71UI9")
  nucleosome <- histones[histones$go_id %in% "GO:0000786", ]
  expect_identical(nucleosome$evidence_by_member[[1L]], list(P0C0S5 = "ECO:0000353", Q71UI9 = "ECO:0000353"))
  albumin <- go_rows(go_recorded(), "P02769")
  expect_identical(length(albumin$evidence_by_member[[1L]]), 0L)
})

test_that("a direct annotation is not propagated, and ancestors are", {
  rows <- go_rows(go_recorded(), "P05141")
  inner <- rows[rows$go_id %in% "GO:0005743", ]
  expect_identical(c(inner$go_name, inner$aspect), c("mitochondrial inner membrane", "cellular_component"))
  expect_identical(c(inner$propagated, inner$inherited), c(FALSE, FALSE))
  expect_true(any(rows$propagated, na.rm = TRUE))
})

test_that("the columns are mzLib's GoAnnotationTsv schema, in order", {
  go <- go_recorded()
  expect_identical(go$column_names[1L], "protein_group")
  expect_identical(utils::tail(go$column_names, 3L), c("go_release", "go_obo_sha256", "annotation_db_sha256"))
  expect_identical(length(go$column_names), 19L)
  expect_identical(names(go$annotations), go$column_names)
  for (column in c("accession_used", "accession_direct", "accession_inherited", "evidence", "entrapment_members")) {
    expect_true(is.list(go$annotations[[column]]), info = column)
  }
})

test_that("provenance pins every input", {
  go <- go_recorded()
  expect_identical(go$go$source_file_name, "go-pxd036557.obo")
  expect_identical(go$go$release, "releases/2026-07-26")
  expect_identical(go$go$term_count, 412)
  expect_identical(go$annotation_database$file_type, "UniProtXml")
  expect_identical(go$annotation_database$protein_count, 5)
  expect_identical(go$header[["source_file_sha256"]], go$groups_file_sha256)
  expect_identical(go$header[["annotation_db_sha256"]], go$annotation_database$sha256)
  expect_identical(unique(go$annotations$go_obo_sha256), go$go$sha256)
})

test_that("the header counts groups, not rows", {
  header <- go_recorded()$header
  expect_true(is.character(header))
  expect_identical(header[["counter_q_value_max"]], "0.01")
  expect_identical(unname(header[c("status_annotated", "status_contaminant", "n_multi_member_groups")]),
                   c("4", "1", "1"))
})

test_that("categories join on go_id and name the subcategory", {
  go <- go_recorded()
  cats <- go$categories
  expect_identical(list(cats$map_name, cats$map_version, cats$anchor_count, cats$row_count),
                   list("organelle", "1", 9, 30))
  by_term <- cats$records[!duplicated(cats$records$go_id), ]
  expect_identical(by_term$subcategory[by_term$go_id == "GO:0005743"], "mitochondrion:inner_membrane")
  expect_true(is.na(by_term$subcategory[by_term$go_id == "GO:0005739"]))
  expect_true(all(cats$records$go_id %in% go$annotations$go_id))
})

test_that("nothing is written unless asked, and no ids dropped in a strict run", {
  go <- go_recorded()
  expect_true(is.null(go$written) && is.null(go$categories_written))
  expect_false(go$skip_unknown_go_ids)
  expect_identical(go$unresolved_go_ids, character(0))
  expect_identical(go$caveats, character(0))
})

test_that("out= returns the summary alone, with the whole table's count", {
  go <- mz$proteins_parse_go_annotations(go_fixture("proteins_annotate_go_out.json"))
  expect_identical(go$written$path, "go_annotations.tsv")
  expect_identical(go$written$row_count, 563)
  expect_identical(c(go$row_count, go$returned_count), c(563, 0))
  expect_true(go$truncated)
  expect_identical(nrow(go$annotations), 0L)
  expect_identical(names(go$annotations), go$column_names)
  expect_true(is.null(go$categories))
})

# ---------------------------------------------------------------- arguments

go_request <- function(...) {
  defaults <- list(groups = "g.tsv", database = "db.xml", go_obo = "go.obo", category_map = NULL,
                   skip_unknown_go_ids = FALSE, out = NULL, categories_out = NULL, limit = NULL, offset = 0)
  given <- list(...)
  defaults[names(given)] <- given
  do.call(mz$proteins_build_go_request, defaults)
}

test_that("the required inputs go on argv", {
  expect_identical(go_request(), c("proteins", "annotate-go", "--groups", "g.tsv", "--database", "db.xml",
                                   "--go-obo", "go.obo", "--offset", "0"))
})

test_that("the optional inputs follow", {
  args <- go_request(category_map = "map.tsv", skip_unknown_go_ids = TRUE, out = "go.tsv",
                     categories_out = "cats.tsv", limit = 0, offset = 5)
  for (pair in list(c("--category-map", "map.tsv"), c("--out", "go.tsv"), c("--categories-out", "cats.tsv"),
                    c("--limit", "0"), c("--offset", "5"))) {
    expect_identical(args[which(args == pair[1L]) + 1L], pair[2L])
  }
  expect_true("--skip-unknown-go-ids" %in% args)
})

test_that("bad arguments fail before a process starts", {
  expect_error(go_request(limit = -1), class = "mzlib_usage_error", contains = "limit must be a whole number")
  expect_error(go_request(offset = -2), class = "mzlib_usage_error", contains = "offset must be a whole number")
  expect_error(go_request(limit = TRUE), class = "mzlib_usage_error", contains = "limit must be a whole number")
  expect_error(go_request(limit = 1.5), class = "mzlib_usage_error", contains = "limit must be a whole number")
  expect_error(go_request(out = ""), class = "mzlib_usage_error", contains = "out must be a single non-empty")
  expect_error(go_request(category_map = 3), class = "mzlib_usage_error", contains = "category_map must be")
  expect_error(go_request(go_obo = NULL), class = "mzlib_usage_error", contains = "go_obo is required")
  expect_error(go_request(skip_unknown_go_ids = NA), class = "mzlib_usage_error", contains = "TRUE or FALSE")
})

test_that("update_go needs a path, and reports the release now on disk", {
  expect_error(proteins_update_go(""), class = "mzlib_usage_error", contains = "go_obo must be")
  update <- mz$proteins_parse_go_update(go_fixture("proteins_update_go.json"))
  expect_true(inherits(update, "mzlibr_go_update"))
  expect_identical(c(update$existed_before, update$changed), c(FALSE, TRUE))
  expect_true(is.na(update$previous_sha256))
  expect_identical(update$go$release, "releases/2026-07-26")
  expect_identical(update$go$term_count, 48340)
  expect_identical(update$url, "https://purl.obolibrary.org/obo/go.obo")
})

test_that("both verbs replay end to end, verb check included, and out= picks its recording", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  go <- proteins_annotate_go("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
                             go_obo = "go-pxd036557.obo", category_map = "organelle_map.tsv")
  expect_identical(go$group_count, 5)
  summary <- proteins_annotate_go("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
                                  go_obo = "go-pxd036557.obo", out = "go_annotations.tsv", limit = 0)
  expect_identical(summary$written$row_count, 563)
  expect_identical(proteins_update_go("go.obo")$go$term_count, 48340)
})

# ---------------------------------------------------------------- against a real mzLib

go_live_files <- function() {
  dir <- Sys.getenv("MZLIB_PROTEINS_FILES", "")
  files <- file.path(dir, c("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
                            "go-pxd036557.obo", "organelle_map.tsv"))
  skip_if(!nzchar(dir) || !all(file.exists(files)), "MZLIB_PROTEINS_FILES does not hold pyMzLib's tests/fixtures/proteins")
  files
}

test_that("LIVE: annotate-go matches the recording", {
  skip_unless_bridge_has("proteins annotate-go")
  files <- go_live_files()
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  live <- proteins_annotate_go(files[[1L]], files[[2L]], go_obo = files[[3L]], category_map = files[[4L]])
  expect_identical(c(live$group_count, live$row_count, live$go$term_count), c(5, 563, 412))
  expect_identical(unname(live$header[c("status_annotated", "status_contaminant")]), c("4", "1"))
  histones <- go_rows(live, "P0C0S5|Q71UI9")
  expect_identical(sum(histones$n_with == histones$n_members), 57L)
  expect_identical(live$categories$row_count, 30)
})

test_that("LIVE: out= holds the whole table while the wire holds none", {
  skip_unless_bridge_has("proteins annotate-go")
  files <- go_live_files()
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  out <- tempfile("mzlibr-go-", fileext = ".tsv")
  on.exit(unlink(out), add = TRUE)
  live <- proteins_annotate_go(files[[1L]], files[[2L]], go_obo = files[[3L]], out = out, limit = 0)
  lines <- readLines(out, warn = FALSE, encoding = "UTF-8")
  expect_identical(lines[[1L]], "#!go_annotation_format 1")
  expect_identical(sum(!startsWith(lines, "#!")), 564L) # the column header and 563 rows
  expect_identical(live$written$row_count, 563)
  expect_identical(c(live$returned_count, nrow(live$annotations)), c(0, 0))
  expect_true(live$truncated)
})

test_that("LIVE: a wrong extension is refused, and a missing go.obo is a usage error", {
  skip_unless_bridge_has("proteins annotate-go")
  files <- go_live_files()
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  expect_error(proteins_annotate_go(files[[1L]], files[[2L]], go_obo = files[[3L]],
                                    out = tempfile(fileext = ".csv")),
               class = "mzlib_usage_error", contains = ".tsv")
  expect_error(proteins_annotate_go(files[[1L]], files[[2L]], go_obo = tempfile(fileext = ".obo")),
               class = "mzlib_usage_error")
})
