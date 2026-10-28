# Two questions about protein sequences and their names, each answered by mzLib:
#
#   proteins_align()             where is residue 120 of this protein in its isoform, or in the
#                                bovine protein? mzLib's PairwiseAligner (#1428), and every
#                                residue's counterpart in the other sequence
#   proteins_parse_accessions()  which entry, isoform and variants does P04637-2_S70N or
#                                NP_000537.3_R72P name? VariantApplication.ParseAccession (#1382)
#
# Ported from pyMzLib's `proteins.py` (align_many, align, Alignment, Alignments, parse_accessions,
# ParsedAccessions), which decided the verbs, the wire fields and the caveats. What changes here is
# only the projection into R's idiom. The pairs are a list of `c(source, target)` vectors, as the
# spec spells mzLibR's call, so one function serves one pair and many; `pair_index` is 1-based, so
# it indexes that list; the tables are data.frames; and `NA` stands for a residue facing a gap,
# where pyMzLib has `None`.
#
# Faithful projection, never a repair: the sequences and the accessions reach mzLib exactly as
# given, every option left out is not sent (mzLib's PairwiseAligner defaults apply, and the result
# says what was used), and nothing here maps, scores or splits anything itself.

# The pyMzLib release whose bridge first dispatches `proteins align` and `proteins parse-accession`.
PROTEINS_ALIGN_SINCE <- "0.6.0"

# The end-gap rules, as mzLib's aligner id spells them: "free" (mzLib's default) or "charged".
PROTEINS_END_GAPS <- c("free", "charged")

# Every `namespace` mzLib's AccessionNamespace can report, lower case, as proteins_resolve_genes()
# writes it.
PROTEINS_NAMESPACES <- c("uniprot", "refseq", "unrecognized")

# ---------------------------------------------------------------- argument shaping

# The pairs as stdin lines, `source<TAB>target`, unchanged, refusing what cannot travel.
proteins_pair_lines <- function(pairs) {
  if (is.data.frame(pairs) || !is.list(pairs)) {
    stop(mzlib_usage_error(paste0(
      "pairs must be a list of c(source, target) character vectors, e.g. ",
      "list(c(\"MKTAYIAK\", \"MKTSYIAK\")); got ", class(pairs)[1L], ". For two columns of ",
      "sequences, Map(c, sources, targets) makes that list."
    )))
  }
  if (length(pairs) == 0L) {
    stop(mzlib_usage_error("At least one pair is required, e.g. list(c(\"MKTAYIAK\", \"MKTSYIAK\"))."))
  }
  lines <- character(length(pairs))
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    if (!is.character(pair) || length(pair) != 2L || anyNA(pair)) {
      stop(mzlib_usage_error(paste0(
        "Pair ", i, " must be c(source, target), two character strings; got ",
        paste(deparse(pair), collapse = " "), "."
      )))
    }
    # One pair is one stdin line, its two sequences split on the tab, so neither may hold a tab or
    # a line break: either would shift or split every later pair.
    for (side in 1:2) {
      if (grepl("[\t\r\n]", pair[[side]])) {
        stop(mzlib_usage_error(paste0(
          "Pair ", i, "'s ", c("source", "target")[side], " contains a tab or line break."
        )))
      }
    }
    # The bridge skips a blank line, which would put every later result on the wrong pair.
    if (!nzchar(trimws(paste0(pair[[1L]], pair[[2L]])))) {
      stop(mzlib_usage_error(paste0("Pair ", i, " has two empty sequences; there is nothing to align.")))
    }
    lines[[i]] <- paste0(pair[[1L]], "\t", pair[[2L]])
  }
  lines
}

# A whole number, or NULL (not sent). The range is mzLib's to check.
proteins_whole_option <- function(value, name) {
  if (is.null(value)) {
    return(character(0))
  }
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value) ||
    value != round(value)) {
    stop(mzlib_usage_error(paste0(
      name, " must be a whole number or NULL; got ", paste(deparse(value), collapse = " "), "."
    )))
  }
  # sprintf, not formatC(format = "d"), which turns a number past an R integer into "NA".
  c(paste0("--", gsub("_", "-", name, fixed = TRUE)), sprintf("%.0f", value))
}

# The aligner options that were given, and only those: mzLib's defaults apply to the rest.
proteins_aligner_args <- function(gap_open, gap_extend, end_gaps, max_cells) {
  args <- c(
    proteins_whole_option(gap_open, "gap_open"),
    proteins_whole_option(gap_extend, "gap_extend"),
    proteins_whole_option(max_cells, "max_cells")
  )
  if (!is.null(end_gaps)) {
    if (!is.character(end_gaps) || length(end_gaps) != 1L || is.na(end_gaps) ||
      !end_gaps %in% PROTEINS_END_GAPS) {
      stop(mzlib_usage_error(paste0(
        "end_gaps must be \"free\" or \"charged\", or NULL for mzLib's default; got ",
        paste(deparse(end_gaps), collapse = " "), "."
      )))
    }
    args <- c(args, "--end-gaps", end_gaps)
  }
  args
}

proteins_build_align_request <- function(pairs, gap_open, gap_extend, end_gaps, max_cells, residue_map) {
  lines <- proteins_pair_lines(pairs)
  if (!is.logical(residue_map) || length(residue_map) != 1L || is.na(residue_map)) {
    stop(mzlib_usage_error("residue_map must be TRUE or FALSE."))
  }
  args <- c("proteins", "align", proteins_aligner_args(gap_open, gap_extend, end_gaps, max_cells))
  if (!residue_map) {
    args <- c(args, "--no-residue-map")
  }
  list(args = args, stdin = paste0(paste(lines, collapse = "\n"), "\n"))
}

# The accessions as stdin lines, unchanged: mzLib parses and never repairs, so neither does this.
#
# Unlike proteins_lines(), nothing is trimmed: a space around an accession is part of what mzLib
# is asked to read, and it reads it as unrecognized. That is the answer, not something to fix.
proteins_accession_lines <- function(accessions) {
  if (!is.character(accessions)) {
    stop(mzlib_usage_error(paste0(
      "accessions must be a character vector, e.g. c(\"P04637-2\", \"NP_000537.3_R72P\"); got ",
      class(accessions)[1L], "."
    )))
  }
  if (length(accessions) == 0L) {
    stop(mzlib_usage_error("At least one accession is required, e.g. \"P04637-2\"."))
  }
  values <- unname(accessions)
  # The bridge skips a blank line and splits on a line break, so either would put every later row
  # on the wrong accession. Refused here rather than silently re-aligned.
  if (anyNA(values) || any(!nzchar(trimws(values)))) {
    stop(mzlib_usage_error("accessions contains a blank or NA entry; every entry is one result row."))
  }
  broken <- grepl("[\r\n]", values)
  if (any(broken)) {
    stop(mzlib_usage_error(paste0("An accession contains a line break: ", deparse(values[broken][1L]), ".")))
  }
  values
}

# ---------------------------------------------------------------- parsing

proteins_parse_alignments <- function(data) {
  records <- proteins_parse_table(data)
  if (!is.null(records) && "pair_index" %in% names(records)) {
    records$pair_index <- as.numeric(records$pair_index) + 1
  }
  residue_map <- proteins_parse_table(data[["residue_map"]])
  if (!is.null(residue_map) && "pair_index" %in% names(residue_map)) {
    residue_map$pair_index <- as.numeric(residue_map$pair_index) + 1
  }
  structure(
    list(
      aligner_id = wire_field(data, "aligner_id", "character", NA_character_),
      matrix = wire_field(data, "matrix", "character", NA_character_),
      gap_open = proteins_count(data, "gap_open"),
      gap_extend = proteins_count(data, "gap_extend"),
      end_gaps = wire_field(data, "end_gaps", "character", NA_character_),
      max_cells = proteins_count(data, "max_cells"),
      record_count = proteins_count(data, "record_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = records,
      residue_map = residue_map,
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_alignments"
  )
}

proteins_parse_accessions_result <- function(data) {
  structure(
    list(
      record_count = proteins_count(data, "record_count"),
      namespace_counts = proteins_parse_counts(data[["namespace_counts"]]),
      with_applied_variants_count = proteins_count(data, "with_applied_variants_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = proteins_parse_table(data),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_parsed_accessions"
  )
}

# ---------------------------------------------------------------- the public surface

#' Align protein sequences with mzLib, and map every residue onto the other sequence
#'
#' Answers **where is residue 120 of this protein in its isoform, its ortholog, or last year's
#' entry?** A position means something only in the sequence it was counted in, and counting the
#' same number along the other sequence is wrong after the first insertion or deletion. This aligns
#' each pair with mzLib's `PairwiseAligner.Align` (mzLib #1428) and maps every residue to its
#' counterpart in the other sequence, or to `NA` when it has none.
#'
#' It is a **global** alignment with affine gaps (Gotoh's form of Needleman-Wunsch), scored with
#' **BLOSUM62**. All the pairs are aligned in one bridge call, in order, so a few thousand isoform
#' pairs cost one process start.
#'
#' **Every option you leave out is mzLib's to set.** mzLibR sends nothing for it, and the result
#' says what the aligner used (`gap_open`, `gap_extend`, `end_gaps`, `max_cells`). mzLib's defaults
#' in the pinned build are gap open 11, gap extend 1, free end gaps and a 200,000,000-cell limit:
#' BLASTP's gap costs for BLOSUM62, and EMBOSS needle's free end gaps.
#'
#' @param pairs A list of `c(source, target)` character vectors of one-letter sequences, one result
#'   row each, in order: `list(c(canonical, isoform))` for one pair. `Map(c, sources, targets)`
#'   builds the list from two vectors. Sent to mzLib exactly as given (case is kept, and compared
#'   case-insensitively). Either sequence may be empty, not both. mzLib refuses `-` and `.`; it
#'   scores any other letter BLOSUM62 does not name (`U`, `O`, a digit, a space) as `X`, so strip
#'   modifications first.
#' @param gap_open Cost of opening a gap, in BLOSUM62's half-bits, charged once per gap: a gap of
#'   `k` residues costs `gap_open + k * gap_extend`. `NULL` (the default) sends nothing, and mzLib's
#'   default applies: 11 half-bits. mzLib refuses a negative value.
#' @param gap_extend Cost of each residue in a gap, the first included, in half-bits. `NULL` (the
#'   default): mzLib's, 1 half-bit. mzLib refuses less than 1.
#' @param end_gaps `"free"` - a gap at either end of either sequence costs nothing, so a signal
#'   peptide or a longer tail on one side is not penalised - or `"charged"`. `NULL` (the default):
#'   mzLib's, `"free"`.
#' @param max_cells The largest table mzLib will build, `(nchar(source) + 1) * (nchar(target) + 1)`
#'   cells of one byte each. A larger pair is refused rather than allowed to exhaust memory. `NULL`
#'   (the default): mzLib's, 200,000,000 cells (about 200 MB). Two titin isoforms need about 1.2
#'   billion; mzLib accepts up to 2,147,483,591.
#' @param residue_map `TRUE` (the default) returns `residue_map`, one row per residue of both
#'   sequences. `FALSE` sends `--no-residue-map` and returns only the pairs table, for a large batch
#'   where you want the scores and counts.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait.
#'
#' @return An `mzlibr_alignments`. `aligner_id` is mzLib's name for the method and every parameter,
#'   e.g. `"global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"`: store it with any position you
#'   map, so the mapping says how it was made. `matrix` is always `"BLOSUM62"`. `gap_open` and
#'   `gap_extend` are the gap costs used, in BLOSUM62's half-bits; `end_gaps` the end-gap rule used;
#'   `max_cells` the table limit used, in cells: yours, or mzLib's default when you gave none.
#'   `record_count` counts the pairs aligned. `caveats` says what the score and the mapping do and do
#'   not mean.
#'
#'   `records` is a data.frame with one row per pair, in input order: `pair_index` (**1-based in R**,
#'   so it indexes `pairs`); `source_length` and `target_length` in residues; `score`, the BLOSUM62
#'   score in half-bits (see the traps); `identical` and `substituted`, the paired residues that are
#'   the same letter and that differ, in residues; `source_residues_facing_gap` and
#'   `target_residues_facing_gap`, the residues with no counterpart, in residues; and
#'   `aligned_source` and `aligned_target`, each sequence with `-` wherever it faces a residue of the
#'   other across a gap, the same length, to print one above the other.
#'
#'   `residue_map` is a data.frame with one row per residue of each sequence, pair by pair, each
#'   pair's source residues first and then its target residues: `pair_index` (1-based); `sequence`,
#'   `"source"` or `"target"`; `position`, 1-based in that sequence, as in mzLib and UniProt;
#'   `residue`, as given; `mapped_position`, the 1-based position of its counterpart in the other
#'   sequence, **`NA` when it faces a gap, never the nearest residue**; and `mapped_residue`, the
#'   counterpart's letter (`NA` facing a gap), a different letter from `residue` on a substitution.
#'   It is `NULL` with `residue_map = FALSE`.
#'
#' @section Reading one position:
#'
#' [proteins_map_source_to_target()] and [proteins_map_target_to_source()] look positions up for one
#' pair; [proteins_alignments()] splits a result into one `mzlibr_alignment` per pair for them. They
#' read the residue map mzLib's `TryMapSourceToTarget` and `TryMapTargetToSource` filled, and compute
#' nothing.
#'
#' @section The traps:
#'
#' **`score` is not an identity, a bit score or an E-value.** It is a raw sum of BLOSUM62 scores in
#' half-bits, less the gaps, and it grows with length. Compare scores only between alignments with
#' the same `aligner_id`. mzLibR computes no percent identity, because there are several
#' definitions: divide `identical` by the length you mean.
#'
#' **It is a global alignment.** Every residue of both sequences is placed, paired or facing a gap.
#' With free end gaps a short sequence is placed on the best-scoring stretch of a long one, even
#' when it has no real match there.
#'
#' **One optimal alignment, chosen the same way every time.** Where two alignments score the same,
#' mzLib's traceback picks one, always the same one. DYNC1I2 isoform 2B is a real case: canonical V76
#' and V82 are both V, either could pair with the isoform's V76, and mzLib pairs V82.
#'
#' **Always check the residue at a mapped position** before carrying a modification across: a
#' phosphosite on S is not a site on A. `mapped_residue` sits beside `residue` for that.
#'
#' @seealso [proteins_map_source_to_target()], [proteins_read()] with `sequences = TRUE` for the
#'   sequences of a database
#' @spec proteins.align
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # mzLib's own deletion case: QRQ (source 9 to 11) is missing from the target.
#' aln <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")))
#' aln
#' writeLines(c(aln$records$aligned_source, aln$records$aligned_target))
#' aln$records[, c("score", "identical", "source_residues_facing_gap")]
#' proteins_map_source_to_target(aln, 8:12)
#'
#' # Two pairs in one call; the free overhang MLLL faces a gap.
#' alns <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"),
#'                             c("MLLLPEPTIDE", "PEPTIDE")))
#' alns$records[, c("pair_index", "score", "aligned_target")]
#' proteins_map_source_to_target(proteins_alignments(alns)[[2]], 1:5)
#'
#' # Charging the end gaps changes the overhang's score, and the aligner id says so:
#' charged <- proteins_align(list(c("MLLLPEPTIDE", "PEPTIDE")), end_gaps = "charged")
#' charged$records$score
#' charged$aligner_id
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_align <- function(pairs, gap_open = NULL, gap_extend = NULL, end_gaps = NULL, max_cells = NULL,
                           residue_map = TRUE, timeout = NULL) {
  request <- proteins_build_align_request(pairs, gap_open, gap_extend, end_gaps, max_cells, residue_map)
  bridge_require_verb("proteins align", PROTEINS_ALIGN_SINCE)
  proteins_parse_alignments(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' One alignment per pair, each with its residue map
#'
#' Splits a [proteins_align()] result into one `mzlibr_alignment` per pair, in input order, for
#' [proteins_map_source_to_target()] and [proteins_map_target_to_source()]. Nothing is computed:
#' each carries its row of `records` and its pair's slice of `residue_map`.
#'
#' @param alignments A [proteins_align()] result.
#' @return A list with one `mzlibr_alignment` per pair. Each holds `pair_index` (1-based),
#'   `aligner_id`, `score` (half-bits), `identical`, `substituted`, `source_residues_facing_gap`
#'   and `target_residues_facing_gap` (residues), `source_length` and `target_length` (residues),
#'   `aligned_source`, `aligned_target`, `caveats`, and `source_to_target` and `target_to_source`:
#'   for each 1-based position `p` of one sequence, element `p` is the paired position in the other,
#'   or `NA` when that residue faces a gap. Those two are `NULL` when the call asked for no residue
#'   map.
#' @seealso [proteins_align()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' alns <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"),
#'                             c("MLLLPEPTIDE", "PEPTIDE")))
#' each <- proteins_alignments(alns)
#' vapply(each, function(a) proteins_map_target_to_source(a, 1), numeric(1))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_alignments <- function(alignments) {
  if (!inherits(alignments, "mzlibr_alignments")) {
    stop(mzlib_usage_error("alignments must be a proteins_align() result."))
  }
  records <- alignments$records
  map <- alignments$residue_map
  lapply(seq_len(nrow(records)), function(i) {
    row <- records[i, , drop = FALSE]
    table_for <- function(side, length) {
      if (is.null(map)) {
        return(NULL)
      }
      rows <- map[map$pair_index == row$pair_index & map$sequence == side, , drop = FALSE]
      out <- rep(NA_real_, length)
      out[rows$position] <- as.numeric(rows$mapped_position)
      out
    }
    structure(
      list(
        pair_index = row$pair_index,
        aligner_id = alignments$aligner_id,
        score = row$score,
        identical = row$identical,
        substituted = row$substituted,
        source_residues_facing_gap = row$source_residues_facing_gap,
        target_residues_facing_gap = row$target_residues_facing_gap,
        source_length = row$source_length,
        target_length = row$target_length,
        aligned_source = row$aligned_source,
        aligned_target = row$aligned_target,
        source_to_target = table_for("source", row$source_length),
        target_to_source = table_for("target", row$target_length),
        caveats = alignments$caveats
      ),
      class = "mzlibr_alignment"
    )
  })
}

# One alignment, from an mzlibr_alignment or a one-pair mzlibr_alignments.
proteins_one_alignment <- function(alignment) {
  if (inherits(alignment, "mzlibr_alignment")) {
    return(alignment)
  }
  if (inherits(alignment, "mzlibr_alignments")) {
    if (nrow(alignment$records) != 1L) {
      stop(mzlib_usage_error(paste0(
        "This result holds ", nrow(alignment$records), " pairs; pick one with ",
        "proteins_alignments(alignment)[[i]]."
      )))
    }
    return(proteins_alignments(alignment)[[1L]])
  }
  stop(mzlib_usage_error("alignment must be one pair's alignment: proteins_align() of one pair, or an element of proteins_alignments()."))
}

proteins_map_positions <- function(alignment, position, side) {
  one <- proteins_one_alignment(alignment)
  table <- one[[paste0(side, if (side == "source") "_to_target" else "_to_source")]]
  if (is.null(table)) {
    stop(mzlib_usage_error("This alignment has no residue map: it came from a call with residue_map = FALSE."))
  }
  length <- one[[paste0(side, "_length")]]
  if (!is.numeric(position) || length(position) == 0L || anyNA(position) ||
    any(position != round(position)) || any(position < 1 | position > length)) {
    stop(mzlib_usage_error(paste0(
      "position must be whole numbers from 1 to ", length, ", the ", side, "'s positions ",
      "(positions are 1-based, as in mzLib and UniProt); got ", paste(deparse(position), collapse = " "), "."
    )))
  }
  table[position]
}

#' The target position paired with each source position
#'
#' Looks positions up in the residue map mzLib's `TryMapSourceToTarget` filled. A residue facing a
#' gap has no counterpart, and its answer is `NA` - never the nearest residue.
#'
#' @param alignment One pair's alignment: [proteins_align()] of one pair, or an element of
#'   [proteins_alignments()].
#' @param position 1-based source positions, 1 to the source's length in residues; a vector is fine.
#' @return The 1-based target positions, one per `position`, `NA` where the source residue faces a
#'   gap.
#' @seealso [proteins_map_target_to_source()], [proteins_align()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' aln <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")))
#' proteins_map_source_to_target(aln, c(8, 9, 10, 11, 12))
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_map_source_to_target <- function(alignment, position) {
  proteins_map_positions(alignment, position, "source")
}

#' The source position paired with each target position
#'
#' The reverse of [proteins_map_source_to_target()], from mzLib's `TryMapTargetToSource`.
#'
#' @param alignment One pair's alignment: [proteins_align()] of one pair, or an element of
#'   [proteins_alignments()].
#' @param position 1-based target positions, 1 to the target's length in residues.
#' @return The 1-based source positions, one per `position`, `NA` where the target residue faces a
#'   gap.
#' @seealso [proteins_map_source_to_target()], [proteins_align()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' aln <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")))
#' proteins_map_target_to_source(aln, 9)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_map_target_to_source <- function(alignment, position) {
  proteins_map_positions(alignment, position, "target")
}

#' Read accessions into entry, isoform or version, and applied variants, with mzLib
#'
#' Answers **which entry, isoform and variants does `P04637-2_S70N` or `NP_000537.3_R72P` name?**
#' Wraps mzLib's `VariantApplication.ParseAccession` (mzLib #1382), and through it
#' `ProteinAccession.Parse`. It is the inverse of the names mzLib gives variant proteoforms:
#' `P04637_S70N` is entry `P04637` with `S70N` applied, and `NP_000537.3_R72P` is RefSeq
#' `NP_000537`, version 3, with `R72P`. Use it to link a search result's accession to anything you
#' keep per entry - genes, orthologs, GO terms - without losing the isoform or the variants.
#' Cutting at the first `_` would turn `NP_000537` into `NP`.
#'
#' **It parses and never repairs.** An accession outside UniProt's and RefSeq's protein grammars
#' comes back `"unrecognized"` and verbatim - lower case, a stray space, a FASTA header, a `DECOY_`
#' prefix, a `|`-joined protein group, and `P12345_2`, mzLib's loader counter for a *different*
#' entry whose accession collided. mzLibR sends each accession exactly as you give it.
#'
#' @param accessions A character vector, one accession per result row, in order, duplicates
#'   included. Sent on stdin, unchanged: not trimmed, not split. A blank or `NA` entry, or one with a
#'   line break, is refused, because it would shift every later row.
#' @param timeout Seconds to allow, or `NULL` (the default) to wait.
#'
#' @return An `mzlibr_parsed_accessions`. `record_count` counts the accessions parsed, duplicates
#'   included. `namespace_counts` is a named vector, namespace to the number of accessions in it,
#'   every namespace present (`uniprot`, `refseq`, `unrecognized`), zeros included: a large
#'   `unrecognized` count is worth reading, see `caveats`. `with_applied_variants_count` counts the
#'   accessions that name applied variants. `caveats` says what mzLib leaves unrecognized rather
#'   than repairs.
#'
#'   `records` is a data.frame with one row per accession, in the order given: `input`, exactly as
#'   given; `entry_accession`, the UniProt entry or unversioned RefSeq accession (`P04637` for
#'   `P04637-2_S70N`), equal to `input` when unrecognized - an **entry**, not a sequence;
#'   `isoform`, the UniProt `-n` number, `NA` without one (always for RefSeq) or past 2,147,483,647;
#'   `version`, the RefSeq `.n` number, `NA` without one (always for UniProt) or past
#'   2,147,483,647; `namespace`; `parent_accession`, the entry the variants were applied to,
#'   isoform or version kept (`P04637-2` for `P04637-2_Q5*`), equal to `input` without variants;
#'   `applied_variants`, the suffix as written (`S70N_A80T`), `NA` without one; and
#'   `has_applied_variants`.
#'
#' @section Joining on it:
#'
#' `entry_accession` answers "same entry?", **not** "same sequence?": `P04637` and `P04637-2` share
#' it, but they are different sequences. Join per-entry data on it, and sequence-level data on
#' `input`. `namespace` is spelled as [proteins_resolve_genes()] spells it, so the two tables join.
#' `setNames(parsed$records$entry_accession, parsed$records$input)` maps each accession to its
#' entry, for grouping.
#'
#' @seealso [proteins_resolve_genes()], [proteins_read()]
#' @spec proteins.parse-accession
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' parsed <- proteins_parse_accessions(c("P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P",
#'   "P04637_S70N_A80T", "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"))
#' parsed
#' parsed$records[, c("input", "namespace", "entry_accession", "isoform", "version",
#'                    "applied_variants")]
#' parsed$namespace_counts
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @export
proteins_parse_accessions <- function(accessions, timeout = NULL) {
  lines <- proteins_accession_lines(accessions)
  bridge_require_verb("proteins parse-accession", PROTEINS_ALIGN_SINCE)
  proteins_parse_accessions_result(bridge_invoke(
    c("proteins", "parse-accession"),
    stdin = paste0(paste(lines, collapse = "\n"), "\n"), timeout = timeout
  ))
}

# ---------------------------------------------------------------- print methods

#' Print protein alignments
#'
#' @param x A [proteins_align()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_alignments <- function(x, ...) {
  cat("<mzlibr_alignments> ", format(x$record_count), " pair(s), ", x$aligner_id, "\n", sep = "")
  if (!is.null(x$records) && nrow(x$records) > 0L) {
    shown <- utils::head(x$records, 3L)
    for (i in seq_len(nrow(shown))) {
      cat("  pair ", format(shown$pair_index[[i]]), ": score ", format(shown$score[[i]]), " half-bits, ",
        format(shown$identical[[i]]), " identical, ", format(shown$substituted[[i]]), " substituted, ",
        format(shown$source_residues_facing_gap[[i]]), " + ", format(shown$target_residues_facing_gap[[i]]),
        " facing a gap\n",
        sep = ""
      )
    }
    if (nrow(x$records) > 3L) cat("  ... see $records\n")
  }
  if (is.null(x$residue_map)) {
    cat("  no residue map (residue_map = FALSE)\n")
  } else {
    cat("  residue map: ", nrow(x$residue_map), " rows; NA mapped_position faces a gap\n", sep = "")
  }
  invisible(x)
}

#' Print one pair's alignment
#'
#' @param x An element of [proteins_alignments()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_alignment <- function(x, ...) {
  cat("<mzlibr_alignment> pair ", format(x$pair_index), ", score ", format(x$score), " half-bits (",
    x$aligner_id, ")\n",
    sep = ""
  )
  cat("  ", x$aligned_source, "\n  ", x$aligned_target, "\n", sep = "")
  invisible(x)
}

#' Print parsed accessions
#'
#' @param x A [proteins_parse_accessions()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_parsed_accessions <- function(x, ...) {
  cat("<mzlibr_parsed_accessions> ", format(x$record_count), " accession(s), ",
    format(x$with_applied_variants_count), " with applied variants\n",
    sep = ""
  )
  if (length(x$namespace_counts) > 0L) {
    cat("  ", paste0(names(x$namespace_counts), " ", x$namespace_counts, collapse = ", "), "\n", sep = "")
  }
  unrecognized <- x$namespace_counts[["unrecognized"]]
  if (!is.null(unrecognized) && isTRUE(unrecognized > 0)) {
    cat("  ! unrecognized ones are kept verbatim, never repaired: see $caveats\n")
  }
  invisible(x)
}
