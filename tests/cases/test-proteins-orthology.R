# proteins_build_orthology(), the projection of genes orthology, and proteins_orthology_pair().
#
# Ported from pyMzLib's test_proteins_orthology.py and the orthology cases of
# test_proteins_live.py. Both payloads are pyMzLib's recordings from the real bridge over mzLib's
# own three-species test world, which is copied byte for byte into tests/orthology/ (see its
# PROVENANCE.txt; .gitattributes keeps it byte-exact, because its sha256s are in the recordings,
# and it sits outside tests/fixtures/ so that no path in the built package passes 100 bytes):
#
#   genes_orthology_pair_status.json  human, mouse and rat, unrestricted gene trees, with
#                                     pair_status = c("homo_sapiens", "mus_musculus")
#   genes_orthology_written.json      the same inputs with restrict_gene_trees = TRUE and out =
#
# What this layer must get right: what it sends (only the options given; the species lines and the
# dumps on stdin, split by `--`), what it refuses before spawning a process, and that it projects
# mzLib's answer without repairing it: a null stays NA or NULL, every status is mzLib's word.

ORTHOLOGY_HUMAN <- "homo_sapiens"
ORTHOLOGY_MOUSE <- "mus_musculus"
ORTHOLOGY_RAT <- "rattus_norvegicus"
ORTHOLOGY_GTFS <- c(
  homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
  mus_musculus = "Mus_musculus.GRCm39.116.gtf",
  rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf"
)
ORTHOLOGY_DUMPS <- paste0("homologies/", names(ORTHOLOGY_GTFS), "/Compara.116.protein_default.homologies.tsv.gz")
ORTHOLOGY_TREES <- "vertebrates.GeneTree_content.default.e116.txt.gz"

orthology_call <- function(fixture, ...) {
  arguments <- utils::modifyList(
    list(release = "116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES, dumps = ORTHOLOGY_DUMPS),
    list(...)
  )
  recorded_bridge_call(fixture, function() do.call(proteins_build_orthology, arguments))
}
orthology_status <- function() mz$proteins_parse_orthology(recorded_data("genes_orthology_pair_status.json"))
orthology_written <- function() mz$proteins_parse_orthology(recorded_data("genes_orthology_written.json"))

# ---------------------------------------------------------------- what is sent

test_that("the species lines, then the dumps, go on stdin, and no option that was not given", {
  seen <- orthology_call("genes_orthology_pair_status.json")
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("genes", "orthology", "--release", "116", "--gene-trees", ORTHOLOGY_TREES))
  expect_identical(seen$stdin, paste0(
    "homo_sapiens\tgtf\tHomo_sapiens.GRCh38.116.gtf\n",
    "mus_musculus\tgtf\tMus_musculus.GRCm39.116.gtf\n",
    "rattus_norvegicus\tgtf\tRattus_norvegicus.GRCr8.116.gtf\n",
    "--\n", paste(ORTHOLOGY_DUMPS, collapse = "\n"), "\n"
  ))
  expect_true(is.null(seen$timeout))
})

test_that("every option given is sent in its wire spelling", {
  request <- mz$proteins_build_orthology_request("116", ORTHOLOGY_TREES, ORTHOLOGY_DUMPS, ORTHOLOGY_GTFS, NULL,
    TRUE, "compara-116", c(ORTHOLOGY_MOUSE, ORTHOLOGY_HUMAN))
  args <- request$args
  expect_true("--restrict-gene-trees" %in% args)
  expect_identical(args[[match("--out", args) + 1L]], "compara-116")
  # From, then to: a status is not symmetric.
  expect_identical(args[[match("--pair-status", args) + 1L]], "mus_musculus,homo_sapiens")
  seen <- orthology_call("genes_orthology_written.json", restrict_gene_trees = TRUE, out = "compara-116", timeout = 600)
  expect_identical(seen$timeout, 600)
})

test_that("gene sets travel as gene-set lines beside gtf lines", {
  request <- mz$proteins_build_orthology_request("116", ORTHOLOGY_TREES, ORTHOLOGY_DUMPS,
    c(homo_sapiens = "human.gtf.gz"), c(mus_musculus = "mouse.genes.tsv", rattus_norvegicus = "rat.genes.tsv"),
    FALSE, NULL, NULL)
  species <- strsplit(sub("\n--\n.*$", "", request$stdin), "\n", fixed = TRUE)[[1L]]
  expect_identical(species, c(
    "homo_sapiens\tgtf\thuman.gtf.gz",
    "mus_musculus\tgene-set\tmouse.genes.tsv",
    "rattus_norvegicus\tgene-set\trat.genes.tsv"
  ))
})

test_that("a whole-number release is sent as its digits", {
  for (release in list(116, 116L, " 116 ")) {
    request <- mz$proteins_build_orthology_request(release, ORTHOLOGY_TREES, ORTHOLOGY_DUMPS, ORTHOLOGY_GTFS, NULL,
      FALSE, NULL, NULL)
    expect_identical(request$args[3:4], c("--release", "116"))
  }
})

# ---------------------------------------------------------------- what is refused before the bridge

test_that("malformed arguments are usage errors before the bridge runs", {
  cases <- list(
    list(list(gtfs = NULL), "At least one species is required"),
    list(list(gtfs = character(0)), "named character vector"),
    list(list(gtfs = list(homo_sapiens = "x.gtf")), "named character vector"),
    list(list(gtfs = "x.gtf"), "named character vector"),
    list(list(gene_sets = c(homo_sapiens = "h.tsv")), "'homo_sapiens' is given twice"),
    list(list(gtfs = c("homo\tsapiens" = "x.gtf")), "cannot hold a tab"),
    list(list(gtfs = stats::setNames("x.gtf", "")), "must be named"),
    list(list(gtfs = c(homo_sapiens = "")), "non-empty"),
    list(list(dumps = character(0)), "dumps is empty"),
    list(list(dumps = "--"), "is not a usable homology dump path"),
    list(list(gene_trees = NULL), "gene_trees is required"),
    list(list(restrict_gene_trees = "yes"), "restrict_gene_trees must be TRUE or FALSE"),
    list(list(pair_status = "homo_sapiens,mus_musculus"), "pair_status must be"),
    list(list(pair_status = c(ORTHOLOGY_HUMAN, "mus,musculus")), "must be a production name"),
    list(list(out = ""), "out must be a single non-empty path"),
    list(list(release = NULL), "release must be an Ensembl release"),
    list(list(release = TRUE), "release must be an Ensembl release"),
    list(list(release = ""), "release must be an Ensembl release"),
    list(list(release = "11 6"), "release must be an Ensembl release"),
    list(list(release = 116.5), "release must be an Ensembl release")
  )
  for (case in cases) {
    arguments <- list(release = "116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES, dumps = ORTHOLOGY_DUMPS)
    arguments[names(case[[1L]])] <- case[[1L]]
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() {
        do.call(proteins_build_orthology, arguments)
      }),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
  # Leaving a required argument out says what to pass, not R's "argument is missing".
  expect_error(proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, dumps = ORTHOLOGY_DUMPS),
    class = "mzlib_usage_error", contains = "gene_trees is required")
  expect_error(proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES),
    class = "mzlib_usage_error", contains = "dumps is empty")
})

# ---------------------------------------------------------------- what is projected

test_that("the snapshot says what went in", {
  snap <- orthology_call("genes_orthology_pair_status.json", pair_status = c(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE))$result
  expect_true(inherits(snap, "mzlibr_orthology_snapshot"))
  expect_identical(snap$release, "116")
  expect_identical(snap$species, c(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE, ORTHOLOGY_RAT))
  expect_identical(list(snap$format, snap$format_version), list("ensembl-orthology-snapshot", 1))
  expect_identical(nchar(snap$snapshot_id), 64L)
  expect_identical(substr(snap$snapshot_id, 1L, 12L), "5487ebcabe48")
  expect_identical(snap$identical_duplicates_dropped, 0)
  expect_identical(snap$gene_sets$species, snap$species)
  expect_identical(snap$gene_sets$gene_count[snap$gene_sets$species == ORTHOLOGY_HUMAN], 4)
  # The test GTFs carry no #!genome-build: NA, not "".
  expect_true(all(is.na(snap$gene_sets$genome_build)))
  expect_identical(snap$gene_trees$collection, "vertebrates")
  expect_identical(snap$gene_trees$member_count, 9)
  expect_true(is.null(snap$gene_trees$restricted_to_gene_sets))
  expect_identical(stats::setNames(snap$dumps$rows_kept, snap$dumps$species),
                   c(homo_sapiens = 1, mus_musculus = 5, rattus_norvegicus = 1))
  expect_identical(unique(snap$dumps$source_file_name), "Compara.116.protein_default.homologies.tsv.gz")
  expect_true(is.null(snap$written))
})

test_that("every human-mouse row came from the mouse dump, and none was lost", {
  snap <- orthology_status()
  expect_identical(snap$dumps$rows_kept[snap$dumps$species == ORTHOLOGY_HUMAN], 1)
  expect_identical(proteins_orthology_pair(snap, ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE)$ortholog_count, 3)
  # Either order.
  expect_identical(proteins_orthology_pair(snap, ORTHOLOGY_MOUSE, ORTHOLOGY_HUMAN),
                   proteins_orthology_pair(snap, ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE))
  expect_identical(proteins_orthology_pair(snap, ORTHOLOGY_MOUSE, ORTHOLOGY_MOUSE)$paralog_count, 1)
  expect_identical(snap$pairs[, c("species_a", "species_b")], data.frame(
    species_a = c(ORTHOLOGY_HUMAN, ORTHOLOGY_HUMAN, ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE, ORTHOLOGY_MOUSE, ORTHOLOGY_RAT),
    species_b = c(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE, ORTHOLOGY_RAT, ORTHOLOGY_MOUSE, ORTHOLOGY_RAT, ORTHOLOGY_RAT),
    stringsAsFactors = FALSE
  ))
  expect_identical(snap$pairs$row_count, c(0, 3, 1, 1, 2, 0))
  expect_identical(snap$pairs$ortholog_count, c(0, 3, 1, 0, 2, 0))
})

test_that("a pair of a species not in the snapshot is a usage error", {
  expect_error(proteins_orthology_pair(orthology_status(), ORTHOLOGY_HUMAN, "danio_rerio"),
    class = "mzlib_usage_error", contains = "not a pair of this snapshot's species")
  expect_error(proteins_orthology_pair(list(), ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE), class = "mzlib_usage_error")
})

test_that("pair_status is mzLib's status for every gene", {
  status <- orthology_status()$pair_status
  expect_identical(list(status$from_species, status$to_species, status$row_count), list(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE, 4))
  expect_identical(status$column_names, c("gene_id", "gene_name", "gene_biotype", "status"))
  expect_identical(names(status$records), status$column_names)
  expect_identical(status$records$status,
                   c("has_ortholog", "has_ortholog", "tree_lacks_target_species", "not_in_any_tree"))
  # No gene_name in the GTF: NA, not "".
  expect_true(all(is.na(status$records$gene_name)))
  # Every status, in mzLib's order, zeros included.
  expect_identical(names(status$status_counts), mz$PROTEINS_ORTHOLOGY_STATUSES)
  expect_identical(status$status_counts[["not_in_gene_set"]], 0)
  expect_identical(status$records$gene_id[status$records$status == "has_ortholog"], c("ENSG00000000001", "ENSG00000000002"))
  expect_identical(status$records[4L, "gene_id"], "ENSG00000000004")
  expect_identical(status$records[4L, "gene_biotype"], "protein_coding")
})

test_that("the caveats name the shared dump name and the unrestricted trees", {
  caveats <- orthology_status()$caveats
  expect_true(any(grepl("share the file name Compara.116.protein_default.homologies.tsv.gz", caveats, fixed = TRUE)))
  expect_true(any(grepl("read unrestricted: 9 genes", caveats, fixed = TRUE)))
})

test_that("a written snapshot lists every file but the manifest", {
  snap <- orthology_call("genes_orthology_written.json", restrict_gene_trees = TRUE, out = "compara-116")$result
  expect_true(is.null(snap$pair_status))
  # The zebrafish gene is not kept.
  expect_identical(snap$gene_trees$member_count, 8)
  expect_identical(snap$gene_trees$restricted_to_gene_sets, sort(snap$gene_sets$sha256, method = "radix"))
  files <- snap$written$files
  rownames(files) <- files$path
  expect_identical(nrow(files), 13L)
  expect_false("manifest.json" %in% files$path)
  # Not a table, so NA rows; an empty pair is still written; ENSG00000000004 is in no tree.
  expect_true(is.na(files["views.sql", "rows"]))
  expect_identical(files["pairs/homo_sapiens__homo_sapiens.parquet", "rows"], 0)
  expect_identical(files["members/homo_sapiens.parquet", "rows"], 3)
  expect_true(all(nchar(files$sha256) == 64L & files$bytes > 0))
  expect_identical(files$rows[startsWith(files$path, "pairs/")], c(0, 3, 1, 1, 2, 0))
  expect_true(endsWith(snap$written$manifest, "manifest.json"))
  expect_true(startsWith(snap$written$manifest, snap$written$directory))
})

test_that("restricting the gene trees does not change the snapshot", {
  unrestricted <- orthology_status()
  restricted <- orthology_written()
  expect_identical(restricted$snapshot_id, unrestricted$snapshot_id)
  expect_identical(restricted$pairs, unrestricted$pairs)
  expect_false(any(grepl("unrestricted", restricted$caveats, fixed = TRUE)))
})

test_that("a payload without its optional blocks projects them to NULL", {
  data <- recorded_data("genes_orthology_written.json")
  data$written <- NA
  snap <- mz$proteins_parse_orthology(data)
  expect_true(is.null(snap$written))
  expect_true(is.null(snap$pair_status))
})

test_that("an older bridge is named, not spawned, for the orthology verb", {
  expect_error(
    recorded_bridge_call(function(args, stdin) stop("reached the verb"), function() {
      proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES, dumps = ORTHOLOGY_DUMPS)
    }, version = "readers_formats.json"),
    class = "mzlib_usage_error", contains = "pyMzLib 0.6.0"
  )
})

test_that("the orthology examples replay end to end, each recording only for its own call", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  snap <- proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES, dumps = ORTHOLOGY_DUMPS,
    pair_status = c(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE))
  expect_identical(snap$pair_status$row_count, 4)
  written <- proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES, dumps = ORTHOLOGY_DUMPS,
    restrict_gene_trees = TRUE, out = "compara-116")
  expect_identical(nrow(written$written$files), 13L)
  # The other direction was not recorded, and a dump left out is another call.
  expect_error(proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES,
    dumps = ORTHOLOGY_DUMPS, pair_status = c(ORTHOLOGY_MOUSE, ORTHOLOGY_HUMAN)),
    class = "mzlib_usage_error", contains = "replay bridge")
  expect_error(proteins_build_orthology("116", gtfs = ORTHOLOGY_GTFS, gene_trees = ORTHOLOGY_TREES,
    dumps = ORTHOLOGY_DUMPS[1:2], pair_status = c(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE)),
    class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- the committed inputs

orthology_inputs <- function() {
  root <- Filter(function(d) file.exists(file.path(d, "PROVENANCE.txt")), c("orthology", file.path("tests", "orthology")))
  if (length(root) == 0L) stop("tests/orthology/ not found")
  root <- root[[1L]]
  list(
    gtfs = stats::setNames(file.path(root, ORTHOLOGY_GTFS), names(ORTHOLOGY_GTFS)),
    dumps = file.path(root, ORTHOLOGY_DUMPS),
    gene_trees = file.path(root, ORTHOLOGY_TREES)
  )
}

test_that("the orthology inputs are byte-exact, as their recordings pin them", {
  # The snapshot_id is a sha256 over these bytes, so a checkout that rewrote line endings would
  # change it: tests/orthology/ is -text in .gitattributes, as in pyMzLib.
  inputs <- orthology_inputs()
  recorded <- orthology_status()
  expect_identical(file.size(inputs$gtfs[["homo_sapiens"]]), 400)
  # No carriage return anywhere: a CRLF checkout would have added one per line.
  for (gtf in inputs$gtfs) {
    expect_false(as.raw(13) %in% readBin(gtf, "raw", file.size(gtf)), info = gtf)
  }
  expect_identical(file.size(inputs$gene_trees), 120)
  expect_identical(unname(file.size(inputs$dumps)), c(195, 295, 199))
  expect_identical(recorded$gene_trees$source_file_name, basename(inputs$gene_trees))
})

# ---------------------------------------------------------------- against a real mzLib

orthology_live <- function() {
  skip_unless_bridge_has("genes orthology")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  invisible(NULL)
}
orthology_live_build <- function(...) {
  inputs <- orthology_inputs()
  arguments <- utils::modifyList(
    list(release = "116", gtfs = inputs$gtfs, gene_trees = inputs$gene_trees, dumps = inputs$dumps),
    list(...)
  )
  do.call(proteins_build_orthology, arguments)
}

test_that("LIVE: orthology matches its recording over the committed inputs", {
  orthology_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  live <- orthology_live_build(pair_status = c(ORTHOLOGY_HUMAN, ORTHOLOGY_MOUSE))
  recorded <- orthology_status()
  # The inputs are byte-exact, so the id is too.
  expect_identical(live$snapshot_id, recorded$snapshot_id)
  expect_identical(live$pair_status$records, recorded$pair_status$records)
  expect_identical(live$pairs$row_count, recorded$pairs$row_count)
  expect_identical(live$caveats, recorded$caveats)
})

test_that("LIVE: orthology writes Parquet once and refuses to overwrite", {
  orthology_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  out <- file.path(tempfile("orthology-"), "compara-116")
  on.exit(unlink(dirname(out), recursive = TRUE), add = TRUE)
  snap <- orthology_live_build(restrict_gene_trees = TRUE, out = out)
  expect_identical(normalizePath(snap$written$directory, winslash = "/"), normalizePath(out, winslash = "/"))
  expect_true(file.exists(snap$written$manifest))
  manifest <- mz$json_parse(paste(readLines(snap$written$manifest, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  expect_identical(manifest$snapshot_id, snap$snapshot_id)
  expect_identical(vapply(manifest$files, function(f) f$path, character(1L)), snap$written$files$path)
  expect_identical(unname(file.size(file.path(out, snap$written$files$path))), snap$written$files$bytes)
  expect_error(orthology_live_build(out = out), class = "mzlib_usage_error", contains = "is not empty")
})

test_that("LIVE: orthology without every dump is refused by mzLib", {
  orthology_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  dumps <- orthology_inputs()$dumps
  expect_error(orthology_live_build(dumps = dumps[c(1L, 3L)]), class = "mzlib_usage_error",
    contains = "No homology dump for mus_musculus")
})

test_that("LIVE: a pair status outside the species is refused before reading", {
  orthology_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  expect_error(orthology_live_build(pair_status = c(ORTHOLOGY_HUMAN, "danio_rerio")), class = "mzlib_usage_error",
    contains = "not one of the species on stdin")
})

test_that("LIVE: a malformed dump is mzLib's correctness failure", {
  orthology_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  bad <- file.path(tempfile("orthology-bad-"), "Compara.116.protein_default.homologies.tsv")
  dir.create(dirname(bad))
  on.exit(unlink(dirname(bad), recursive = TRUE), add = TRUE)
  writeLines("gene_stable_id\tspecies", bad)
  dumps <- orthology_inputs()$dumps
  condition <- expect_error(orthology_live_build(dumps = c(dumps[1:2], bad)), class = "mzlib_bridge_error")
  expect_identical(condition$error_type, "InvalidDataException")
  expect_true(grepl("unexpected columns", conditionMessage(condition), fixed = TRUE))
})
