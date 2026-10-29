# proteins_map_residues(), the projection of genes map-residues.
#
# Ported from pyMzLib's test_proteins_map_residues.py and test_proteins_map_residues_live.py. Both
# payloads are pyMzLib's recordings from the real bridge over a real cut of Ensembl release 116,
# copied byte for byte into tests/map_residues/ (see its PROVENANCE.txt; .gitattributes keeps it
# byte-exact, because the files' sha256s are in the recordings):
#
#   genes_map_residues_tp53.json   human, mouse and rat gene sets, nine queries through the TP53 and
#                                  GAPDH gene trees' alignments, one per outcome and more, proteins = TRUE
#   genes_map_residues_chain.json  the same inputs, one query (human p53 Ser15 to mouse), no proteins
#
# What this layer must get right: what it sends (the species lines, then `--`, then one query per
# line, and only the options it was given), what it refuses before spawning a process, and that it
# projects mzLib's answer without repairing it: an outcome is mzLib's own name, a gap stays NA.

MAP_GTFS <- c(
  homo_sapiens = "Homo_sapiens.GRCh38.116.gtf",
  mus_musculus = "Mus_musculus.GRCm39.116.gtf",
  rattus_norvegicus = "Rattus_norvegicus.GRCr8.116.gtf"
)
MAP_ALIGNMENT <- "Compara.116.protein_default.aa.fasta.gz"
MAP_TREES <- "vertebrates.GeneTree_content.default.e116.txt.gz"
MAP_TP53 <- "ENSP00000269305"
MAP_TRP53 <- "ENSMUSP00000104298"
MAP_RAT_TP53 <- "ENSRNOP00000074031"
MAP_QUERIES <- list(
  list(MAP_TP53, 15, MAP_TRP53), list(MAP_TP53, 20, MAP_TRP53), list(MAP_TP53, 46, MAP_TRP53),
  list(MAP_TP53, 392, MAP_TRP53), list(MAP_TRP53, 18, MAP_TP53), list(MAP_TP53, 15, MAP_RAT_TP53),
  list(MAP_TP53, 15, "ENSMUSP00000113213"), list("ENSP00000352610", 15, MAP_TRP53),
  list(MAP_TP53, 15, "ENSMUSP00000005371")
)

map_call <- function(fixture, queries = MAP_QUERIES, ...) {
  arguments <- utils::modifyList(
    list(queries = queries, alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS),
    list(...)
  )
  recorded_bridge_call(fixture, function() do.call(proteins_map_residues, arguments))
}
map_tp53 <- function() mz$proteins_parse_residue_mappings(recorded_data("genes_map_residues_tp53.json"))
map_chain <- function() mz$proteins_parse_residue_mappings(recorded_data("genes_map_residues_chain.json"))

# ---------------------------------------------------------------- what is sent

test_that("map_residues sends the species lines, then the queries, and no option it was not given", {
  seen <- map_call("genes_map_residues_chain.json", list(list(MAP_TP53, 15, MAP_TRP53)))
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("genes", "map-residues", "--alignment", MAP_ALIGNMENT, "--gene-trees", MAP_TREES))
  expect_identical(seen$stdin, paste0(
    "homo_sapiens\tgtf\tHomo_sapiens.GRCh38.116.gtf\n",
    "mus_musculus\tgtf\tMus_musculus.GRCm39.116.gtf\n",
    "rattus_norvegicus\tgtf\tRattus_norvegicus.GRCr8.116.gtf\n",
    "--\nENSP00000269305\t15\tENSMUSP00000104298\n"
  ))
  expect_true(is.null(seen$timeout))
})

test_that("map_residues sends --proteins and the timeout when given, and takes a data.frame of queries", {
  queries <- data.frame(source_protein_id = MAP_TP53, position = 15L, target_protein_id = MAP_TRP53,
    stringsAsFactors = FALSE)
  seen <- map_call("genes_map_residues_tp53.json", queries, proteins = TRUE, timeout = 600)
  expect_identical(seen$args[[length(seen$args)]], "--proteins")
  expect_identical(seen$timeout, 600)
  expect_true(endsWith(seen$stdin, "--\nENSP00000269305\t15\tENSMUSP00000104298\n"))
})

test_that("without gene sets no species line is sent, and mzLib keeps every gene", {
  seen <- map_call("genes_map_residues_chain.json", list(list(MAP_TP53, 15, MAP_TRP53)), gtfs = NULL)
  expect_identical(seen$stdin, "--\nENSP00000269305\t15\tENSMUSP00000104298\n")
})

test_that("a position is sent as given, never checked against a length", {
  # mzLib owns the rule (MapResidue throws for a position outside the source), so mzLibR sends 0, a
  # negative and a huge position unchanged.
  request <- mz$proteins_build_map_residues_request(
    list(list(MAP_TP53, 0, MAP_TRP53), list(MAP_TP53, -3L, MAP_TRP53), list(MAP_TP53, 1e6, MAP_TRP53)),
    MAP_ALIGNMENT, MAP_TREES, MAP_GTFS, NULL, FALSE
  )
  queries <- strsplit(strsplit(request$stdin, "--\n", fixed = TRUE)[[1L]][[2L]], "\n", fixed = TRUE)[[1L]]
  expect_identical(queries, paste(MAP_TP53, c("0", "-3", "1000000"), MAP_TRP53, sep = "\t"))
})

test_that("gene sets are sent as gene-set lines beside gtf lines", {
  request <- mz$proteins_build_map_residues_request(list(list(MAP_TP53, 15, MAP_TRP53)), MAP_ALIGNMENT, MAP_TREES,
    c(homo_sapiens = "human.gtf.gz"), c(mus_musculus = "mouse.genes.tsv"), FALSE)
  species <- strsplit(strsplit(request$stdin, "--\n", fixed = TRUE)[[1L]][[1L]], "\n", fixed = TRUE)[[1L]]
  expect_identical(species, c("homo_sapiens\tgtf\thuman.gtf.gz", "mus_musculus\tgene-set\tmouse.genes.tsv"))
})

test_that("malformed arguments are usage errors before the bridge runs", {
  never <- function(args, stdin) stop("the bridge must not be called")
  refused <- function(contains, queries = list(list(MAP_TP53, 15, MAP_TRP53)), ...) {
    arguments <- utils::modifyList(
      list(queries = queries, alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS),
      list(...)
    )
    expect_error(recorded_bridge_call(never, function() do.call(proteins_map_residues, arguments)),
      class = "mzlib_usage_error", contains = contains)
  }
  refused("queries is empty", list())
  refused("queries must be a data.frame", "ENSP00000269305\t15\tENSMUSP00000104298")
  refused("must be list(source_protein_id, position, target_protein_id)", list(list(MAP_TP53, 15)))
  refused("position must be a whole number", list(list(MAP_TP53, "15", MAP_TRP53)))
  refused("position must be a whole number", list(list(MAP_TP53, 15.5, MAP_TRP53)))
  refused("position must be a whole number", list(list(MAP_TP53, TRUE, MAP_TRP53)))
  refused("position must be a whole number", list(list(MAP_TP53, NA_real_, MAP_TRP53)))
  refused("source must be a protein id", list(list("", 15, MAP_TRP53)))
  refused("target must be a protein id", list(list(MAP_TP53, 15, NULL)))
  refused("contains a tab or line break", list(list(MAP_TP53, 15, "ENSMUSP\t1")))
  refused("three columns", data.frame(a = MAP_TP53, b = 15))
  refused("alignment is required", alignment = NULL)
  refused("gene_trees must be a single non-empty path", gene_trees = "")
  refused("proteins must be TRUE or FALSE", proteins = "yes")
  refused("is given twice", gene_sets = c(homo_sapiens = "h.tsv"))
  refused("gtfs must be a named character vector", gtfs = list(homo_sapiens = "x.gtf"))
  expect_error(proteins_map_residues(list(list(MAP_TP53, 15, MAP_TRP53)), alignment = MAP_ALIGNMENT),
    class = "mzlib_usage_error", contains = "gene_trees is required")
})

test_that("an older bridge is named, not spawned, for the map-residues verb", {
  expect_error(
    recorded_bridge_call(function(args, stdin) stop("reached the verb"), function() {
      proteins_map_residues(MAP_QUERIES, alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS)
    }, version = "readers_formats.json"),
    class = "mzlib_usage_error", contains = "pyMzLib 0.7.0"
  )
})

# ---------------------------------------------------------------- what is projected

test_that("every outcome is mzLib's name, and a gap is never the nearest residue", {
  mapped <- map_tp53()
  expect_true(inherits(mapped, "mzlibr_residue_mappings"))
  expect_identical(mapped$row_count, 9)
  expect_identical(mapped$column_names, c("source_protein_id", "position", "target_protein_id", "outcome", "column",
    "target_position", "target_residue"))
  r <- mapped$records
  expect_identical(names(r), mapped$column_names)
  expect_identical(r$outcome, c("Aligned", "Aligned", "GapInTarget", "Aligned", "Aligned", "Aligned",
    "DifferentAlignments", "SourceNotInAlignment", "TargetNotInAlignment"))
  expect_identical(names(mapped$outcome_counts), mz$PROTEINS_COLUMN_OUTCOMES, info = "every outcome, mzLib's order")
  expect_identical(unname(mapped$outcome_counts), c(5, 1, 1, 1, 1))
  # Human S15, S20 and S392 are mouse S18, S23 and S389; S46 faces a gap in mouse.
  expect_identical(r$target_position, c(18, 23, NA, 389, 15, 15, NA, NA, NA))
  expect_identical(r$target_residue, c("S", "S", NA, "S", "S", "S", NA, NA, NA))
  expect_identical(r$column, c(476, 547, 863, 3036, 476, 476, NA, NA, NA))
  gap <- r[3L, ]
  expect_identical(c(gap$source_protein_id, gap$target_protein_id, gap$outcome), c(MAP_TP53, MAP_TRP53, "GapInTarget"))
  expect_identical(c(gap$position, gap$column), c(46, 863))
})

test_that("the envelope says what went in", {
  mapped <- map_tp53()
  a <- mapped$alignment
  expect_identical(c(a$source_file_name, a$release, a$collection), c(MAP_ALIGNMENT, "116", "protein_default"))
  expect_identical(c(a$alignment_count, a$protein_count), c(2, 5))
  expect_identical(a$restricted_to_gene_tree_content, mapped$gene_trees$sha256)
  expect_identical(mapped$gene_trees$member_count, 5)
  expect_identical(mapped$gene_trees$restricted_to_gene_sets, sort(mapped$gene_sets$sha256, method = "radix"))
  expect_identical(mapped$gene_sets$species, names(MAP_GTFS))
  expect_identical(mapped$gene_sets$genome_build[[1L]], "GRCh38.p14")
  expect_identical(mapped$gene_sets$gene_count, c(2, 2, 1))
  expect_true(any(grepl("2 queries name a protein the alignment did not keep", mapped$caveats, fixed = TRUE)))
})

test_that("proteins hold the kept canonical proteins, with their ungapped sequences", {
  mapped <- map_tp53()
  expect_true(mapped$proteins_included)
  kept <- mapped$proteins
  expect_identical(kept$row_count, 5)
  expect_identical(kept$column_names, c("protein_id", "alignment_index", "tree_id", "sequence"))
  expect_identical(kept$records$protein_id, sort(kept$records$protein_id, method = "radix"), info = "mzLib's ordinal order")
  sequences <- stats::setNames(kept$records$sequence, kept$records$protein_id)
  expect_identical(nchar(unname(sequences[c(MAP_TP53, MAP_TRP53)])), c(393L, 390L))
  expect_false(any(grepl("-", sequences, fixed = TRUE)))
  expect_identical(c(substr(sequences[[MAP_TP53]], 15, 15), substr(sequences[[MAP_TRP53]], 18, 18)), c("S", "S"))
  trees <- stats::setNames(kept$records$tree_id, kept$records$protein_id)
  expect_identical(unname(trees[c(MAP_TP53, MAP_TRP53)]), rep("ENSGT00950000183153", 2L))
  expect_false("ENSP00000352610" %in% names(sequences), info = "a TP53 isoform Compara did not align")
  # 1-based in R, like every position: the wire's 0 and 1 are 1 and 2.
  wire <- vapply(recorded_data("genes_map_residues_tp53.json")$proteins$columns$alignment_index, as.numeric, numeric(1L))
  expect_identical(kept$records$alignment_index, wire + 1)
  expect_identical(sort(unique(kept$records$alignment_index)), c(1, 2))
})

test_that("without proteins the block is NULL", {
  site <- map_chain()
  expect_false(site$proteins_included)
  expect_true(is.null(site$proteins) && "proteins" %in% names(site))
  expect_identical(site$records$outcome, "Aligned")
  expect_identical(c(site$records$position, site$records$column, site$records$target_position), c(15, 476, 18))
  expect_identical(site$records$target_residue, "S")
  expect_identical(site$caveats, character(0))
})

test_that("a payload without gene sets projects an empty gene_sets table and NULL restriction", {
  data <- recorded_data("genes_map_residues_chain.json")
  data$gene_sets <- structure(list(), names = character(0))
  data$gene_trees$restricted_to_gene_sets <- NA
  mapped <- mz$proteins_parse_residue_mappings(data)
  expect_identical(nrow(mapped$gene_sets), 0L)
  expect_true("species" %in% names(mapped$gene_sets))
  expect_true(is.null(mapped$gene_trees$restricted_to_gene_sets))
})

test_that("the map-residues examples replay end to end, each recording only for its own call", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  mapped <- proteins_map_residues(MAP_QUERIES, alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS,
    proteins = TRUE)
  expect_identical(mapped$proteins$row_count, 5)
  site <- proteins_map_residues(MAP_QUERIES[1L], alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS)
  expect_identical(site$records$target_position, 18)
  # Another query, another --proteins choice, or another alignment file is another call.
  expect_error(proteins_map_residues(MAP_QUERIES[2L], alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS),
    class = "mzlib_usage_error", contains = "replay bridge")
  expect_error(proteins_map_residues(MAP_QUERIES[1L], alignment = MAP_ALIGNMENT, gene_trees = MAP_TREES, gtfs = MAP_GTFS,
    proteins = TRUE), class = "mzlib_usage_error", contains = "replay bridge")
  expect_error(proteins_map_residues(MAP_QUERIES[1L], alignment = "Compara.115.protein_default.aa.fasta.gz",
    gene_trees = MAP_TREES, gtfs = MAP_GTFS), class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- the committed inputs

map_inputs <- function() {
  root <- Filter(function(d) file.exists(file.path(d, "PROVENANCE.txt")), c("map_residues", file.path("tests", "map_residues")))
  if (length(root) == 0L) stop("tests/map_residues/ not found")
  root <- root[[1L]]
  list(
    gtfs = stats::setNames(file.path(root, MAP_GTFS), names(MAP_GTFS)),
    alignment = file.path(root, MAP_ALIGNMENT),
    gene_trees = file.path(root, MAP_TREES)
  )
}

test_that("the map-residues inputs are byte-exact, as their recordings pin them", {
  # The recordings carry each file's sha256, so a checkout that rewrote line endings would change
  # them: tests/map_residues/ is -text in .gitattributes, as in pyMzLib.
  inputs <- map_inputs()
  for (gtf in inputs$gtfs) {
    expect_false(as.raw(13) %in% readBin(gtf, "raw", file.size(gtf)), info = gtf)
  }
  expect_identical(unname(file.size(inputs$gtfs)), c(488, 495, 300))
  expect_identical(c(file.size(inputs$alignment), file.size(inputs$gene_trees)), c(2035, 647))
})

# ---------------------------------------------------------------- against a real mzLib

map_live <- function() {
  skip_unless_bridge_has("genes map-residues")
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  invisible(NULL)
}
map_live_call <- function(queries, ...) {
  inputs <- map_inputs()
  arguments <- utils::modifyList(
    list(queries = queries, alignment = inputs$alignment, gene_trees = inputs$gene_trees, gtfs = inputs$gtfs),
    list(...)
  )
  do.call(proteins_map_residues, arguments)
}

test_that("LIVE: each map-residues recording is what the bridge answers now", {
  map_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  for (case in list(list("genes_map_residues_tp53.json", TRUE), list("genes_map_residues_chain.json", FALSE))) {
    recorded <- mz$proteins_parse_residue_mappings(recorded_data(case[[1L]]))
    queries <- lapply(seq_len(nrow(recorded$records)), function(i) {
      list(recorded$records$source_protein_id[[i]], recorded$records$position[[i]], recorded$records$target_protein_id[[i]])
    })
    live <- map_live_call(queries, proteins = case[[2L]])
    expect_identical(live$records, recorded$records, info = case[[1L]])
    expect_identical(live$outcome_counts, recorded$outcome_counts, info = case[[1L]])
    # Byte-exact inputs, so the same sha256s and counts.
    expect_identical(live$alignment, recorded$alignment, info = case[[1L]])
    expect_identical(live$gene_trees, recorded$gene_trees, info = case[[1L]])
    expect_identical(live$caveats, recorded$caveats, info = case[[1L]])
    expect_identical(live$proteins$records, recorded$proteins$records, info = case[[1L]])
  }
})

test_that("LIVE: a position past the source fails the whole call with mzLib's message", {
  map_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  # mzLib's MapResidue throws for it (the bridge checks no length), so the first query's answer is
  # lost too: no row comes back.
  condition <- expect_error(map_live_call(list(list(MAP_TP53, 15, MAP_TRP53), list(MAP_TP53, 394, MAP_TRP53))),
    class = "mzlib_usage_error",
    contains = c("Query 2 (ENSP00000269305, position 394, ENSMUSP00000104298)", "ENSP00000269305 has 393 residues.",
      "no query of this call was answered"))
  expect_error(map_live_call(list(list(MAP_TP53, 0, MAP_TRP53))), class = "mzlib_usage_error", contains = "has 393 residues")
})

test_that("LIVE: a source not kept is an outcome before its position is looked at", {
  map_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  mapped <- map_live_call(list(list("ENSP00000352610", 100000, MAP_TRP53)))
  expect_identical(mapped$records$outcome, "SourceNotInAlignment")
})

test_that("LIVE: without gene sets mzLib's default keeps every gene, and a caveat says so", {
  map_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  mapped <- map_live_call(list(list(MAP_TP53, 15, MAP_TRP53)), gtfs = NULL)
  expect_identical(nrow(mapped$gene_sets), 0L)
  expect_true(is.null(mapped$gene_trees$restricted_to_gene_sets))
  expect_identical(mapped$records$target_position, 18, info = "the same answer for the same proteins")
  expect_true(any(grepl("No gene sets were given, so mzLib kept every gene", mapped$caveats, fixed = TRUE)))
})

test_that("LIVE: a ragged alignment is mzLib's correctness failure", {
  map_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  bad <- file.path(tempfile("map-bad-"), MAP_ALIGNMENT)
  dir.create(dirname(bad))
  on.exit(unlink(dirname(bad), recursive = TRUE), add = TRUE)
  connection <- gzfile(bad, "wb")
  writeLines(c(">ENSDARP00000000001", "PEPT", ">ENSDARP00000000002", "PEP", "//"), connection, sep = "\n")
  close(connection)
  condition <- expect_error(map_live_call(list(list(MAP_TP53, 15, MAP_TRP53)), alignment = bad), class = "mzlib_bridge_error")
  expect_identical(condition$error_type, "InvalidDataException")
  expect_true(grepl("is 3 columns wide, the alignment 4", conditionMessage(condition), fixed = TRUE))
})

test_that("LIVE: a missing alignment is a usage error before any input is read", {
  map_live()
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  expect_error(map_live_call(list(list(MAP_TP53, 15, MAP_TRP53)), alignment = file.path(tempdir(), "absent.aa.fasta.gz")),
    class = "mzlib_usage_error", contains = "Compara gene-tree alignment not found")
})
