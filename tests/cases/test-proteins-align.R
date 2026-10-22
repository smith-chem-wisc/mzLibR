# proteins_align() and its residue maps, and proteins_parse_accessions().
#
# Ported from pyMzLib's test_proteins_align.py and the align and parse-accession cases of
# test_proteins_live.py. Every payload is pyMzLib's recording from the real bridge:
#
#   proteins_align_deletion.json   mzLib's own deletion case (QRQ missing from the target)
#   proteins_align_many.json       that pair and mzLib's extra-N-terminus case, in one call
#   proteins_align_charged.json    the extra N-terminus with end_gaps = "charged"
#   proteins_align_no_map.json     the two pairs with residue_map = FALSE
#   proteins_align_dync1i2.json    DYNC1I2 (Q13409) against its isoforms 2B and 2C
#   proteins_align_albumin.json    human albumin (P02768) against bovine (P02769)
#   proteins_parse_accession.json  nine accessions, recognized and not
#
# What this layer must get right: it sends only the options it was given, sends the sequences and
# accessions unchanged, makes pair_index 1-based, and projects mzLib's answer residue for residue
# without repairing it - a residue facing a gap is NA, never a neighbour.

ALIGN_DELETION <- c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")
ALIGN_OVERHANG <- c("MLLLPEPTIDE", "PEPTIDE")
PARSE_ACCESSIONS <- c(
  "P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P", "P04637_S70N_A80T", "P04637_2",
  "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"
)

# The sequences pyMzLib read from its fixture databases with sequences = TRUE, recorded.
align_sequences <- function() {
  db <- mz$proteins_parse_database(recorded_data("proteins_read_sequences.json"))
  stats::setNames(db$proteins$sequence, db$proteins$accession)
}
align_call <- function(fixture, ...) {
  arguments <- list(...)
  recorded_bridge_call(fixture, function() do.call(proteins_align, arguments))
}

# ---------------------------------------------------------------- what is sent

test_that("align sends the pair on stdin and no option it was not given", {
  seen <- align_call("proteins_align_deletion.json", list(ALIGN_DELETION))
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("proteins", "align"))
  expect_identical(seen$stdin, "MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\n")
  expect_true(is.null(seen$timeout))
})

test_that("each option given is sent, in its wire spelling, and nothing else", {
  cases <- list(
    list(list(gap_open = 5), c("--gap-open", "5")),
    list(list(gap_extend = 2), c("--gap-extend", "2")),
    list(list(max_cells = 1000), c("--max-cells", "1000")),
    list(list(max_cells = 2147483591), c("--max-cells", "2147483591")),
    list(list(end_gaps = "charged"), c("--end-gaps", "charged")),
    list(list(gap_open = 0, end_gaps = "free"), c("--gap-open", "0", "--end-gaps", "free")),
    list(list(residue_map = FALSE), "--no-residue-map")
  )
  for (case in cases) {
    request <- do.call(mz$proteins_build_align_request, c(
      list(pairs = list(ALIGN_DELETION)),
      utils::modifyList(list(gap_open = NULL, gap_extend = NULL, end_gaps = NULL, max_cells = NULL,
                             residue_map = TRUE), case[[1L]])
    ))
    expect_identical(request$args, c("proteins", "align", case[[2L]]))
  }
})

test_that("many pairs are one call, one stdin line each, in order", {
  seen <- align_call("proteins_align_many.json", list(ALIGN_DELETION, ALIGN_OVERHANG))
  expect_identical(seen$stdin, "MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\nMLLLPEPTIDE\tPEPTIDE\n")
  expect_identical(length(Filter(function(c) !identical(c$args, "version"), seen$calls)), 1L)
})

test_that("an empty sequence on one side is sent as it is", {
  request <- mz$proteins_build_align_request(list(c("", "PEP")), NULL, NULL, NULL, NULL, TRUE)
  expect_identical(request$stdin, "\tPEP\n")
})

test_that("bad pairs are refused before the bridge", {
  cases <- list(
    list(c("PEPTIDE", "PEPTIDE"), "list of c(source, target)"),
    list("PEPTIDE", "list of c(source, target)"),
    list(data.frame(source = "A", target = "B"), "Map(c, sources, targets)"),
    list(list("PEPTIDE"), "Pair 1 must be c(source, target)"),
    list(list(c("A", "B", "C")), "Pair 1 must be c(source, target)"),
    list(list(c("PEPTIDE", NA)), "Pair 1 must be c(source, target)"),
    list(list(c(1, 2)), "Pair 1 must be c(source, target)"),
    list(list(c("PEP\tTIDE", "PEPTIDE")), "tab or line break"),
    list(list(c("PEPTIDE", "PEP\nTIDE")), "tab or line break"),
    list(list(c("PEPTIDE", "PEPTIDE"), c("", " ")), "Pair 2 has two empty sequences"),
    list(list(), "At least one pair")
  )
  for (case in cases) {
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() proteins_align(case[[1L]])),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

test_that("bad options are refused before the bridge", {
  cases <- list(
    list(list(gap_open = TRUE), "gap_open must be a whole number"),
    list(list(gap_extend = 1.5), "gap_extend must be a whole number"),
    list(list(max_cells = "1000"), "max_cells must be a whole number"),
    list(list(gap_open = NA_real_), "gap_open must be a whole number"),
    list(list(end_gaps = "Free"), "end_gaps must be"),
    list(list(end_gaps = TRUE), "end_gaps must be"),
    list(list(residue_map = NA), "residue_map must be TRUE or FALSE")
  )
  for (case in cases) {
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() {
        do.call(proteins_align, c(list(list(ALIGN_DELETION)), case[[1L]]))
      }),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

# ---------------------------------------------------------------- what is projected

test_that("align projects mzLib's counts and rows", {
  aln <- align_call("proteins_align_deletion.json", list(ALIGN_DELETION))$result
  expect_true(inherits(aln, "mzlibr_alignments"))
  expect_identical(aln$aligner_id, "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free")
  row <- aln$records
  # 39 + 55 - 14: mzLib's own figure for this case.
  expect_identical(c(row$score, row$identical, row$substituted), c(80, 19, 0))
  expect_identical(c(row$source_residues_facing_gap, row$target_residues_facing_gap), c(3, 0))
  expect_identical(c(row$source_length, row$target_length), c(22, 19))
  expect_identical(row$aligned_target, "MKTAYIAK---ISFVKSHFSRQ")
  # 1-based in R, so it indexes the pairs given.
  expect_identical(row$pair_index, 1)
  expect_identical(length(aln$caveats), 7L)
})

test_that("a residue facing a gap maps to NA, never a neighbour", {
  aln <- align_call("proteins_align_deletion.json", list(ALIGN_DELETION))$result
  expect_identical(proteins_map_source_to_target(aln, 7:13), c(7, 8, NA, NA, NA, 9, 10))
  expect_identical(proteins_map_target_to_source(aln, 9), 12)
  one <- proteins_alignments(aln)[[1L]]
  expect_identical(c(length(one$source_to_target), length(one$target_to_source)), c(22L, 19L))
  expect_false(anyNA(proteins_map_target_to_source(one, 1:19)))
})

test_that("a position outside the sequence is refused", {
  aln <- align_call("proteins_align_deletion.json", list(ALIGN_DELETION))$result
  for (position in list(0, 23, -1, TRUE, "1", 1.5, NA_real_, numeric(0))) {
    expect_error(proteins_map_source_to_target(aln, position), class = "mzlib_usage_error", contains = "1-based")
  }
})

test_that("charged end gaps come back in the envelope", {
  aln <- align_call("proteins_align_charged.json", list(ALIGN_OVERHANG), end_gaps = "charged")$result
  expect_identical(aln$records$score, 39 - (11 + 4))
  expect_true(endsWith(aln$aligner_id, "end-gaps=charged"))
  expect_identical(aln$end_gaps, "charged")
})

test_that("many pairs keep input order, with one alignment each", {
  alns <- align_call("proteins_align_many.json", list(ALIGN_DELETION, ALIGN_OVERHANG))$result
  expect_identical(alns$record_count, 2)
  expect_identical(list(alns$matrix, alns$gap_open, alns$gap_extend, alns$end_gaps, alns$max_cells),
                   list("BLOSUM62", 11, 1, "free", 2e8))
  expect_identical(alns$records$pair_index, c(1, 2))
  expect_identical(nrow(alns$residue_map), 41L + 18L)
  expect_identical(sort(unique(alns$residue_map$pair_index)), c(1, 2))
  pairs <- proteins_alignments(alns)
  expect_identical(vapply(pairs, function(a) a$pair_index, numeric(1)), c(1, 2))
  # The free MLLL overhang faces a gap.
  expect_identical(proteins_map_source_to_target(pairs[[2L]], 1:5), c(NA, NA, NA, NA, 1))
  expect_identical(alns$records$aligned_target[[2L]], "----PEPTIDE")
  # A many-pair result is not one alignment: pick one.
  expect_error(proteins_map_source_to_target(alns, 1), class = "mzlib_usage_error", contains = "proteins_alignments")
})

test_that("residue_map = FALSE sends the flag and leaves no map", {
  seen <- align_call("proteins_align_no_map.json", list(ALIGN_DELETION, ALIGN_OVERHANG), residue_map = FALSE)
  expect_identical(seen$args, c("proteins", "align", "--no-residue-map"))
  alns <- seen$result
  expect_true(is.null(alns$residue_map))
  expect_identical(alns$records$score, c(80, 39))
  first <- proteins_alignments(alns)[[1L]]
  expect_true(is.null(first$source_to_target))
  expect_error(proteins_map_source_to_target(first, 1), class = "mzlib_usage_error", contains = "no residue map")
})

test_that("real isoforms shift every position after the missing exon", {
  seq <- align_sequences()
  alns <- align_call("proteins_align_dync1i2.json",
                     list(c(seq[["Q13409"]], seq[["Q13409-2"]]), c(seq[["Q13409"]], seq[["Q13409-3"]])))$result
  pairs <- proteins_alignments(alns)
  two_b <- pairs[[1L]]
  two_c <- pairs[[2L]]
  expect_identical(c(two_b$identical, two_b$source_residues_facing_gap), c(632, 6))
  expect_identical(proteins_map_source_to_target(two_b, 75), 75)
  expect_identical(proteins_map_source_to_target(two_b, 76:81), rep(NA_real_, 6))
  # V76 and V82 could each pair with the isoform's V76; mzLib's tie-break pairs V82.
  expect_identical(substring(seq[["Q13409"]], c(76, 82), c(76, 82)), c("V", "V"))
  expect_identical(proteins_map_source_to_target(two_b, 82), 76)
  expect_identical(c(two_c$identical, two_c$source_residues_facing_gap), c(612, 26))
  expect_identical(proteins_map_source_to_target(two_c, 638), 612)
  # The article's table, row by row.
  expect_identical(proteins_map_source_to_target(two_b, c(75, 76, 82, 100, 120, 133)), c(75, NA, 76, 94, 114, 127))
  expect_identical(proteins_map_source_to_target(two_c, c(75, 76, 82, 100, 120, 133)), c(75, NA, 76, 94, NA, 107))
})

test_that("real orthologs' substituted residues keep their own letters", {
  seq <- align_sequences()
  aln <- align_call("proteins_align_albumin.json", list(c(seq[["P02768"]], seq[["P02769"]])))$result
  row <- aln$records
  expect_identical(c(row$identical, row$substituted, row$source_residues_facing_gap), c(465, 142, 2))
  # Mature Cys34, conserved.
  expect_identical(proteins_map_source_to_target(aln, 58), 58)
  expect_identical(substr(seq[["P02769"]], 58, 58), "C")
  expect_identical(proteins_map_source_to_target(aln, 26), 26)
  expect_identical(c(substr(seq[["P02768"]], 26, 26), substr(seq[["P02769"]], 26, 26)), c("A", "T"))
  map <- aln$residue_map[aln$residue_map$sequence == "source", ]
  expect_identical(as.numeric(sum(!is.na(map$mapped_residue) & map$mapped_residue != map$residue)), row$substituted)
})

test_that("the alignments replay end to end, as the help page runs them", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(proteins_align(list(ALIGN_DELETION))$records$score, 80)
  expect_identical(proteins_align(list(ALIGN_OVERHANG), end_gaps = "charged")$records$score, 24)
  # The charged recording does not answer a call that leaves end_gaps to mzLib.
  expect_error(proteins_align(list(ALIGN_OVERHANG)), class = "mzlib_usage_error", contains = "replay bridge")
  # Nor does any recording with a residue map answer residue_map = FALSE.
  expect_error(proteins_align(list(ALIGN_DELETION), residue_map = FALSE), class = "mzlib_usage_error",
    contains = "other --no-residue-map choice")
})

# ---------------------------------------------------------------- parse_accessions

test_that("parse_accessions sends each accession unchanged", {
  seen <- recorded_bridge_call("proteins_parse_accession.json", function() {
    proteins_parse_accessions(c(" P04637", "P04637 ", "p04637", "P04637|P02768"))
  })
  expect_identical(seen$calls[[1L]]$args, "version")
  expect_identical(seen$args, c("proteins", "parse-accession"))
  # Not trimmed, not case-folded, not split.
  expect_identical(seen$stdin, " P04637\nP04637 \np04637\nP04637|P02768\n")
})

test_that("parse_accessions projects mzLib's reading", {
  parsed <- recorded_bridge_call("proteins_parse_accession.json", function() proteins_parse_accessions(PARSE_ACCESSIONS))$result
  expect_true(inherits(parsed, "mzlibr_parsed_accessions"))
  rows <- parsed$records
  rownames(rows) <- rows$input
  expect_identical(parsed$record_count, 9)
  expect_identical(parsed$namespace_counts, c(unrecognized = 4, uniprot = 4, refseq = 1))
  expect_identical(sort(names(parsed$namespace_counts)), sort(mz$PROTEINS_NAMESPACES))
  expect_identical(parsed$with_applied_variants_count, 2)
  expect_identical(c(rows["Q13409-2", "entry_accession"], format(rows["Q13409-2", "isoform"])), c("Q13409", "2"))
  refseq <- rows["NP_000537.3_R72P", ]
  expect_identical(list(refseq$entry_accession, refseq$version, refseq$parent_accession, refseq$applied_variants),
                   list("NP_000537", 3, "NP_000537.3", "R72P"))
  expect_true(is.na(refseq$isoform))
  expect_identical(rows["P04637_S70N_A80T", "applied_variants"], "S70N_A80T")
  expect_true(rows["P04637_S70N_A80T", "has_applied_variants"])
  for (unrecognized in c("P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768")) {
    expect_identical(rows[unrecognized, "namespace"], "unrecognized", info = unrecognized)
    # Kept verbatim, never repaired.
    expect_identical(rows[unrecognized, "entry_accession"], unrecognized, info = unrecognized)
    expect_true(is.na(rows[unrecognized, "applied_variants"]), info = unrecognized)
  }
  expect_identical(setNames(rows$entry_accession, rows$input)[["P04637-2"]], "P04637")
  expect_identical(length(parsed$caveats), 7L)
})

test_that("bad accession lists are refused before the bridge", {
  cases <- list(
    list(character(0), "At least one accession"),
    list(c("P04637", ""), "blank or NA"),
    list(c("P04637", "   "), "blank or NA"),
    list(c("P04637", NA), "blank or NA"),
    list("P04637\nP02768", "line break"),
    list("P04637\r", "line break"),
    list(list("P04637"), "character vector"),
    list(4637, "character vector")
  )
  for (case in cases) {
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the bridge"), function() proteins_parse_accessions(case[[1L]])),
      class = "mzlib_usage_error", contains = case[[2L]]
    )
  }
})

test_that("an older bridge is named, not spawned, for both verbs", {
  for (call in list(function() proteins_align(list(ALIGN_DELETION)), function() proteins_parse_accessions("P04637"))) {
    expect_error(
      recorded_bridge_call(function(args, stdin) stop("reached the verb"), call, version = "readers_formats.json"),
      class = "mzlib_usage_error", contains = "pyMzLib 0.6.0"
    )
  }
})

test_that("accessions replay end to end, and other inputs are not answered", {
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  expect_identical(proteins_parse_accessions(PARSE_ACCESSIONS)$record_count, 9)
  expect_error(proteins_parse_accessions("P04637"), class = "mzlib_usage_error", contains = "replay bridge")
})

# ---------------------------------------------------------------- against a real mzLib

align_live <- function(verb) {
  skip_unless_bridge_has(verb)
  options(mzlibr.bridge = Sys.getenv("MZLIB_BRIDGE"))
  invisible(NULL)
}

test_that("LIVE: align leaves every option to mzLib and says what it used", {
  align_live("proteins align")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  aln <- proteins_align(list(ALIGN_OVERHANG))
  expect_identical(list(aln$matrix, aln$gap_open, aln$gap_extend, aln$end_gaps, aln$max_cells),
                   list("BLOSUM62", 11, 1, "free", 2e8))
  expect_identical(aln$aligner_id, "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free")
  expect_identical(aln$records$score, 39)
})

test_that("LIVE: canonical DYNC1I2 maps onto its isoforms as recorded", {
  align_live("proteins align")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  seq <- align_sequences()
  live <- proteins_align(list(c(seq[["Q13409"]], seq[["Q13409-2"]]), c(seq[["Q13409"]], seq[["Q13409-3"]])))
  recorded <- mz$proteins_parse_alignments(recorded_data("proteins_align_dync1i2.json"))
  expect_identical(live$records, recorded$records)
  expect_identical(live$residue_map, recorded$residue_map)
  expect_identical(proteins_map_source_to_target(proteins_alignments(live)[[1L]], 82), 76)
})

test_that("LIVE: align options reach mzLib", {
  align_live("proteins align")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  aln <- proteins_align(list(ALIGN_OVERHANG), gap_open = 5, gap_extend = 2, end_gaps = "charged")
  expect_identical(aln$aligner_id, "global-affine;BLOSUM62;open=5;extend=2;end-gaps=charged")
  expect_identical(aln$records$score, 39 - (5 + 4 * 2))
})

test_that("LIVE: mzLib's refusals are usage errors", {
  align_live("proteins align")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  cases <- list(
    list(list(c("PEP-TIDE", "PEPTIDE")), list(), "gap character"),
    list(list(c("PEPTIDE", "PEPTIDE"), c("AAAA", "AAAC")), list(max_cells = 1), "Pair 1 was refused"),
    list(list(c("PEPTIDE", "PEPTIDE")), list(gap_open = -1), "--gap-open was refused by mzLib"),
    list(list(c("PEPTIDE", "PEPTIDE")), list(gap_extend = 0), "--gap-extend was refused by mzLib")
  )
  for (case in cases) {
    expect_error(do.call(proteins_align, c(list(case[[1L]]), case[[2L]])), class = "mzlib_usage_error", contains = case[[3L]])
  }
})

test_that("LIVE: parse_accessions never repairs what it is given", {
  align_live("proteins parse-accession")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  given <- c("P04637-2", " P04637", "P04637 ", "NP_000537.3_R72P", "P04637_2")
  parsed <- proteins_parse_accessions(given)
  expect_identical(parsed$records$input, given)
  expect_identical(parsed$records$namespace, c("uniprot", "unrecognized", "unrecognized", "refseq", "unrecognized"))
  expect_identical(parsed$records$applied_variants, c(NA, NA, NA, "R72P", NA))
  recorded <- mz$proteins_parse_accessions_result(recorded_data("proteins_parse_accession.json"))
  expect_identical(proteins_parse_accessions(PARSE_ACCESSIONS)$records, recorded$records)
})

test_that("LIVE: parse_accessions and resolve_genes share one namespace vocabulary", {
  align_live("proteins parse-accession")
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)
  root <- Sys.getenv("MZLIB_PROTEINS_FILES")
  skip_if(!nzchar(root) || !file.exists(file.path(root, "human_subset.xml")),
    "MZLIB_PROTEINS_FILES does not hold pyMzLib's tests/fixtures/proteins")
  genes <- proteins_resolve_genes(file.path(root, "human_subset.xml"),
    gtf = file.path(root, "Homo_sapiens.GRCh38.116.gtf"))
  parsed <- proteins_parse_accessions(genes$resolutions$accession)
  expect_identical(parsed$records$namespace, genes$resolutions$namespace)
  expect_identical(parsed$records$entry_accession, genes$resolutions$entry_accession)
  expect_identical(parsed$records$isoform, genes$resolutions$isoform)
})
