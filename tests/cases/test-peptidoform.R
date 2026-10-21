# Digestion, fragmentation, and the modification census.
#
# Two recordings of serum albumin, each trimmed to the same two peptides.
#
# peptidoform_P02768_small.json is pyMzLib's, byte for byte: max_modifications = 1 and nothing else
# given, so mzLib's defaults and the bridge's HCD - b and y ions.
#
# peptidoform_P02768_etd_small.json is mzLibR's own, recorded from the 0.5.0 bridge (mzLib
# d10017ce) on 2026-10-21 with the call the old defaults made, spelled out: protease = "trypsin|P",
# dissociation = "ETD", max_modifications = 1. pyMzLib re-recorded its fixture under the new
# defaults, which left no ETD digest to test against, and ETD carries two fragmentation traps at
# once: a zDot series numbered past the c series, and zDot counts that come in *below* their own
# maximum because of proline suppression. (mzLib #1114 removed the spurious ETD y series that used
# to be a third trap here.) Every one of those is asserted here rather than only described,
# because these are the numbers the documentation quotes. Its two peptides are identical, ion for
# ion, to those of the recording it replaces, made from the 0.4.0 bridge.

recorded_digest_data <- function(name = "peptidoform_P02768_small.json") {
  mz$json_parse(paste(
    readLines(fixture_path(name), warn = FALSE),
    collapse = "\n"
  ))
}

recorded_digest <- function() {
  mz$peptidoform_parse(recorded_digest_data())
}

recorded_etd_digest <- function() {
  mz$peptidoform_parse(recorded_digest_data("peptidoform_P02768_etd_small.json"))
}

# Every option NULL unless given, as peptidoform_fragments() itself defaults them.
build_args <- function(accession = "P02768", protease = NULL, dissociation = NULL,
                       terminus = NULL, modifications = TRUE, missed_cleavages = NULL,
                       min_length = NULL, max_length = NULL, max_modifications = NULL,
                       max_isoforms = NULL) {
  mz$peptidoform_build_args(
    accession, protease, dissociation, terminus, modifications,
    missed_cleavages, min_length, max_length, max_modifications, max_isoforms
  )
}

paired <- function(args, flag) args[which(args == flag) + 1L]

# ---------------------------------------------------------------- parsing

test_that("the digest carries the run's own description", {
  digest <- recorded_digest()
  expect_true(inherits(digest, "mzlibr_digest"))
  expect_identical(digest$accession, "P02768")
  expect_identical(digest$name, "ALBU_HUMAN")
  expect_identical(digest$organism, "Homo sapiens")
  expect_identical(digest$sequence_length, 609)
  # What the bridge used for a call that named neither: mzLib's protease and the bridge's HCD.
  expect_identical(digest$protease, "trypsin")
  expect_identical(digest$dissociation, "HCD")
  expect_identical(digest$terminus, "Both")
  expect_identical(digest$max_modifications, 1)
  expect_true(digest$modifications_applied)

  etd <- recorded_etd_digest()
  expect_identical(etd$protease, "trypsin|P")
  expect_identical(etd$dissociation, "ETD")
})

test_that("HCD, the default dissociation, emits b and y", {
  series <- mz$digest_fragments_by_series(recorded_digest())
  expect_identical(series$product_type, c("b", "y"))
  # 30 + 34 backbone cleavages on each series, for peptides of 31 and 35 residues.
  expect_identical(series$n, c(64L, 64L))
})

test_that("peptides, fragments and modifications are plain data.frames with no factors", {
  digest <- recorded_digest()
  for (frame in list(digest$peptides, digest$fragments, digest$modifications)) {
    expect_true(is.data.frame(frame))
  }
  expect_false(is.factor(digest$peptides$base_sequence))
  expect_false(is.factor(digest$fragments$product_type))
  expect_identical(nrow(digest$peptides), 2L)
})

test_that("the three frames join on peptide_index", {
  digest <- recorded_digest()
  expect_identical(digest$peptides$peptide_index, 1:2)
  expect_true(all(digest$fragments$peptide_index %in% digest$peptides$peptide_index))
  expect_true(all(digest$modifications$peptide_index %in% digest$peptides$peptide_index))

  # The point of the key: narrow the peptides, then take exactly those fragments.
  wanted <- digest$peptides[digest$peptides$length == 31, ]
  theirs <- digest$fragments[digest$fragments$peptide_index %in% wanted$peptide_index, ]
  expect_identical(nrow(theirs), 60L)
})

test_that("modifications carry NA for a residue modification's absent terminus", {
  # `null` on the wire. Read as R NULL it would have deleted the field and shifted the frame;
  # read as NA it sits in its column where a caller can see it is not a terminal modification.
  digest <- recorded_digest()
  expect_identical(nrow(digest$modifications), 1L)
  expect_identical(digest$modifications$id, "N6-succinyllysine on K")
  expect_identical(digest$modifications$one_based_residue, 35)
  expect_true(is.na(digest$modifications$terminus))
  expect_true(is.character(digest$modifications$terminus))
})

# ---------------------------------------------------------------- the fragmentation traps

test_that("ETD emits c and zDot, no y (mzLib #1114)", {
  # mzLib #1114 removed y from the ETD/ECD product sets - radical N-Ca cleavage yields c/z-dot
  # only, while b/y come from amide cleavage. Pins the post-#1114 contract, so a fixture or bridge
  # regressed below f6b0f0d1 fails here rather than resurrecting the retired "spurious y" note.
  series <- mz$digest_fragments_by_series(recorded_etd_digest())
  present <- series$product_type

  expect_false("y" %in% present)
  expect_false("b" %in% present)
  expect_identical(sort(present), c("c", "zDot"))
})

test_that("the c series runs 1..length-1", {
  digest <- recorded_etd_digest()
  for (index in digest$peptides$peptide_index) {
    peptide <- digest$peptides[digest$peptides$peptide_index == index, ]
    c_ions <- digest$fragments[
      digest$fragments$peptide_index == index & digest$fragments$product_type == "c",
    ]
    expect_identical(nrow(c_ions), as.integer(peptide$length - 1), info = paste("peptide", index))
    expect_identical(max(c_ions$fragment_number), peptide$length - 1)
  }
})

test_that("the zDot series is numbered up to length, not length-1", {
  # The extra ion numbered `length` is the whole peptide minus NH2 — the N-Ca cleavage at
  # residue 1. Correct and deliberate, not a defect: an mzLibRust finding was raised about this
  # and correctly dropped after investigation.
  digest <- recorded_etd_digest()
  for (index in digest$peptides$peptide_index) {
    peptide <- digest$peptides[digest$peptides$peptide_index == index, ]
    z_ions <- digest$fragments[
      digest$fragments$peptide_index == index & digest$fragments$product_type == "zDot",
    ]
    expect_identical(max(z_ions$fragment_number), peptide$length, info = paste("peptide", index))
  }
})

test_that("z-dot ions are suppressed at proline while the complementary c ions are not", {
  # smith-chem-wisc/mzLib#1110. The asymmetry is the finding: if suppression were right, the
  # complementary c ion would go too. Here the zDot series is numbered to `length` but *counts*
  # fewer than that, and the shortfall is exactly the number of prolines.
  digest <- recorded_etd_digest()

  for (index in digest$peptides$peptide_index) {
    peptide <- digest$peptides[digest$peptides$peptide_index == index, ]
    mine <- digest$fragments[digest$fragments$peptide_index == index, ]
    z_count <- sum(mine$product_type == "zDot")
    c_count <- sum(mine$product_type == "c")

    prolines <- length(gregexpr("P", peptide$base_sequence, fixed = TRUE)[[1]])
    expect_identical(z_count, as.integer(peptide$length - prolines),
      info = paste("peptide", index, "prolines", prolines)
    )
    # The complementary c ions are all still there — that is the defect.
    expect_identical(c_count, as.integer(peptide$length - 1),
      info = paste("peptide", index)
    )
    expect_true(z_count < c_count + 1)
  }
})

test_that("the c and z-dot ladders close on the precursor mass", {
  # The strongest internal check available offline: for every complementary pair,
  # c_k + z_(L-k) = M + 1.00782503, the hydrogen atom. It holds to eight decimal places on both
  # peptides, so a fragment ladder that was subtly wrong - an off-by-one in numbering, a lost
  # neutral loss, the wrong terminus - could not survive it.
  #
  # Note which constant appears here. The c/z complementarity is a real chemical relation and
  # takes the **hydrogen atom**; peptide m/z takes the **proton** (see peptide_mz). Both are in
  # this package, for different and correct reasons, which is exactly why getting them mixed up
  # is easy and invisible.
  digest <- recorded_etd_digest()

  for (index in digest$peptides$peptide_index) {
    peptide <- digest$peptides[digest$peptides$peptide_index == index, ]
    mine <- digest$fragments[digest$fragments$peptide_index == index, ]
    c_ions <- mine[mine$product_type == "c", ]
    z_ions <- mine[mine$product_type == "zDot", ]

    pairs <- 0L
    for (number in c_ions$fragment_number) {
      partner <- z_ions$neutral_mass[z_ions$fragment_number == peptide$length - number]
      if (length(partner) != 1L) {
        next
      }
      total <- c_ions$neutral_mass[c_ions$fragment_number == number] + partner
      expect_true(abs(total - peptide$monoisotopic_mass - 1.00782503) < 1e-6,
        info = paste("peptide", index, "c", number)
      )
      pairs <- pairs + 1L
    }
    # And there really were pairs to check, so the loop cannot pass by doing nothing.
    expect_true(pairs > 25L, info = paste("peptide", index, "pairs", pairs))
  }
})

test_that("modification positions index the peptide's own residues", {
  # Not the parent protein's. The succinyl-lysine sits at residue 35 of a 35-residue peptide,
  # which begins at residue 509 of albumin - so a position read against the protein would be
  # wrong by 508 and still look perfectly plausible.
  digest <- recorded_digest()
  modifications <- digest$modifications
  for (row in seq_len(nrow(modifications))) {
    peptide <- digest$peptides[
      digest$peptides$peptide_index == modifications$peptide_index[row],
    ]
    expect_true(modifications$one_based_residue[row] >= 1,
      info = paste("row", row)
    )
    expect_true(modifications$one_based_residue[row] <= peptide$length,
      info = paste(
        "residue", modifications$one_based_residue[row],
        "in a peptide of length", peptide$length
      )
    )
  }
})

test_that("fragments carry neutral_mass and no m/z", {
  # Deliberate, and it must stay that way: a c or z ion carries only the fixed charges within
  # its own span, and per-fragment charge accounting does not exist on this wire. An mz column
  # here would be a number with no defensible meaning.
  digest <- recorded_digest()
  expect_true("neutral_mass" %in% names(digest$fragments))
  expect_false("mz" %in% names(digest$fragments))
  expect_false("charge" %in% names(digest$fragments))
})

test_that("fragments_by_series is ordered and totals the fragment frame", {
  digest <- recorded_etd_digest()
  series <- mz$digest_fragments_by_series(digest)
  expect_identical(series$product_type, sort(series$product_type))
  expect_identical(sum(series$n), nrow(digest$fragments))
  expect_identical(series$n, c(64L, 63L))
})

# ---------------------------------------------------------------- peptidoforms vs sequences

test_that("peptides are peptidoforms and distinct sequences are counted separately", {
  digest <- recorded_digest()
  expect_identical(nrow(digest$peptides), 2L)
  expect_identical(mz$digest_distinct_base_sequences(digest), 2L)
})

test_that("modified peptides are the ones with a modification count", {
  digest <- recorded_digest()
  modified <- mz$digest_modified_peptides(digest)
  expect_identical(nrow(modified), 1L)
  expect_true(all(modified$modification_count > 0))
})

test_that("this digest did not hit the isoform cap", {
  expect_false(mz$digest_truncated(recorded_digest()))
})

test_that("truncation is reported when peptides hit the cap", {
  digest <- recorded_digest()
  digest$peptides_at_isoform_cap <- 12
  expect_true(mz$digest_truncated(digest))
})

# ---------------------------------------------------------------- m/z

test_that("peptide m/z uses the proton mass, not the hydrogen atom", {
  digest <- recorded_digest()
  peptide <- digest$peptides[1, ]
  expected <- (peptide$monoisotopic_mass + 2 * 1.00727646677) / 2
  expect_equal(mz$peptide_mz(peptide, 2), expected)

  # The wrong convention is 1.1 ppm high at m/z 500 — a match versus a miss on an Orbitrap, and
  # invisible in the answer. This asserts the two are genuinely different at this mass.
  hydrogen <- (peptide$monoisotopic_mass + 2 * 1.007825) / 2
  expect_true(abs(hydrogen - expected) > 1e-4)
})

test_that("fixed charges are not double-counted", {
  # A mass that already carries a charge needs fewer protons added, not the same number. Adding
  # a full complement would put a 2+ trimethylated peptide half a Thomson high, on the most
  # important histone modification there is.
  peptides <- data.frame(monoisotopic_mass = 1000, fixed_charges = 1, stringsAsFactors = FALSE)
  expect_equal(mz$peptide_mz(peptides, 2), (1000 + 1 * 1.00727646677) / 2)

  neutral <- data.frame(monoisotopic_mass = 1000, fixed_charges = 0, stringsAsFactors = FALSE)
  expect_equal(mz$peptide_mz(neutral, 2), (1000 + 2 * 1.00727646677) / 2)

  # A peptide whose fixed charge equals the requested charge needs no protonation at all.
  expect_equal(mz$peptide_mz(peptides, 1), 1000 / 1)
})

test_that("a charge below the peptide's fixed charge is refused", {
  peptides <- data.frame(monoisotopic_mass = 1000, fixed_charges = 2, stringsAsFactors = FALSE)
  expect_error(mz$peptide_mz(peptides, 1),
    class = "mzlib_usage_error", contains = "cannot be observed at that charge"
  )
})

test_that("an impossible charge is refused", {
  peptides <- data.frame(monoisotopic_mass = 1000, fixed_charges = 0, stringsAsFactors = FALSE)
  for (charge in list(0, -1, 1.5, NA_real_, "2", c(1, 2))) {
    expect_error(mz$peptide_mz(peptides, charge), class = "mzlib_usage_error")
  }
})

test_that("peptide_mz is vectorised over rows", {
  digest <- recorded_digest()
  values <- mz$peptide_mz(digest$peptides, 2)
  expect_identical(length(values), 2L)
  expect_true(all(values > 0))
})

# ---------------------------------------------------------------- the modification census

test_that("the census counts sites, applied and annotated separately", {
  # `annotated_modification_sites` is NOT `annotated_modifications_loaded`. A histone carries
  # K9me1, K9me2, K9me3 and K9ac at one residue — four modifications, one site. Conflating them
  # made H3.1 look as though 93 annotations had been dropped when all had loaded.
  census <- recorded_digest()$census
  expect_identical(census$sites, 14)
  expect_identical(census$applied, 14)
  expect_identical(census$annotated, 38)
})

test_that("the census names the feature types it could not load", {
  by_type <- recorded_digest()$census$by_type
  expect_identical(nrow(by_type), 2L)

  glyco <- by_type[by_type$type == "glycosylation site", ]
  expect_identical(glyco$count, 24)
  expect_false(glyco$loaded)

  modified <- by_type[by_type$type == "modified residue", ]
  expect_identical(modified$count, 14)
  expect_true(modified$loaded)
})

test_that("census_explain states the counts and the real reason", {
  explanation <- mz$census_explain(recorded_digest())
  expect_true(grepl("14 of 38", explanation, fixed = TRUE), info = explanation)
  expect_true(grepl("24 x glycosylation site", explanation, fixed = TRUE), info = explanation)
  expect_true(grepl("mzLib#1112", explanation, fixed = TRUE), info = explanation)
})

test_that("census_explain does not state a reason that is false", {
  # This is the test that exists because the feature built to prevent silent wrongness *was*
  # the wrong thing, twice in one day. It first claimed the excluded sites had no defined
  # composition — false for 22 of the 24, which UniProt's ptmlist gives C6H10O5 and 162.052823
  # — and then implied they should therefore have been loaded, which is also false.
  #
  # A disclosure that states a false reason is worse than no disclosure, so both claims are
  # asserted absent, and the true grounds are asserted present.
  explanation <- mz$census_explain(recorded_digest())
  expect_false(grepl("no defined composition", explanation, fixed = TRUE), info = explanation)
  expect_false(grepl("should have been loaded", explanation, fixed = TRUE), info = explanation)
  expect_true(grepl("labile", explanation, fixed = TRUE), info = explanation)
  expect_true(grepl("in vitro", explanation, fixed = TRUE), info = explanation)
  expect_true(grepl("Read the annotations", explanation, fixed = TRUE), info = explanation)
})

test_that("a census with nothing excluded says so plainly", {
  digest <- recorded_digest()
  digest$census$annotated <- 14
  digest$census$by_type$loaded <- TRUE
  expect_true(grepl("All 14 annotated modifications were applied",
    mz$census_explain(digest),
    fixed = TRUE
  ))
})

test_that("unresolved modifications are named, not silently dropped", {
  # On histone H3.1 seven N6-lactoyllysine sites vanished while the type summary still reported
  # "modified residue ... loaded". The name is the only way a caller can tell.
  digest <- recorded_digest()
  digest$census$unresolved <- "N6-lactoyllysine on K"
  explanation <- mz$census_explain(digest)
  expect_true(grepl("N6-lactoyllysine on K", explanation, fixed = TRUE))
  expect_true(grepl("absent from its own modification list", explanation, fixed = TRUE))
})

# ---------------------------------------------------------------- argument assembly

test_that("a call that names nothing sends only the accession: the defaults are mzLib's", {
  # The rule the siblings follow (pyMzLib #83): an option the caller leaves out is not sent, so
  # mzLib's DigestionParams default applies, and for dissociation the bridge's HCD until
  # mzLib#1488. Built from peptidoform_fragments()'s own formals, so a default that crept back
  # into the signature would put an option on the wire and fail here.
  defaults <- lapply(formals(peptidoform_fragments)[-1L], eval)
  args <- mz$peptidoform_build_args(
    "P02768", defaults$protease, defaults$dissociation, defaults$terminus, defaults$modifications,
    defaults$missed_cleavages, defaults$min_length, defaults$max_length,
    defaults$max_modifications, defaults$max_isoforms
  )
  expect_identical(args, c("peptidoform", "fragments", "--accession", "P02768"))
  for (name in c("protease", "dissociation", "terminus", "missed_cleavages", "min_length",
                 "max_length", "max_modifications", "max_isoforms")) {
    expect_true(is.null(defaults[[name]]), info = name)
  }
})

test_that("a default call is answered by the default recording, end to end", {
  # Through the real argument handling, transport and parsing, against the stand-in bridge. A call
  # that still sent the old "trypsin|P" and ETD would fit no recording and fail.
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  old <- mz$replay_bridge_start()
  on.exit(mz$replay_bridge_stop(old), add = TRUE)
  digest <- peptidoform_fragments("P02768")
  expect_identical(digest$protease, "trypsin")
  expect_identical(digest$dissociation, "HCD")
  expect_identical(nrow(digest$peptides), 310L)
  expect_identical(digest_distinct_base_sequences(digest), 202L)
  # Left to mzLib, so not stamped with a number mzLibR would have had to invent.
  expect_true(is.na(digest$min_length) && is.na(digest$max_length) && is.na(digest$missed_cleavages))
  expect_true(grepl("min_length mzLib's default", paste(capture.output(print(digest)), collapse = "\n"),
    fixed = TRUE
  ))
})

test_that("every option that is given is sent, in the wire's spelling", {
  args <- build_args(
    protease = "trypsin|P", dissociation = "ETD", terminus = "N", missed_cleavages = 0,
    min_length = 1, max_length = 30, max_modifications = 3, max_isoforms = 1024
  )
  expect_identical(paired(args, "--protease"), "trypsin|P")
  expect_identical(paired(args, "--dissociation"), "ETD")
  expect_identical(paired(args, "--terminus"), "N")
  expect_identical(paired(args, "--missed-cleavages"), "0")
  expect_identical(paired(args, "--min-length"), "1")
  expect_identical(paired(args, "--max-length"), "30")
  expect_identical(paired(args, "--max-mods"), "3")
  expect_identical(paired(args, "--max-isoforms"), "1024")
  expect_false("--no-modifications" %in% args)
  # Each alone, too: one given option does not drag the others onto the wire.
  expect_identical(build_args(protease = "trypsin|P"),
    c("peptidoform", "fragments", "--accession", "P02768", "--protease", "trypsin|P"))
})

test_that("the max-mods 1 recording answers only a call that sends max_modifications = 1", {
  # Its echo says max_modifications 1, but so would mzLib's default if it were 1: the echo of a value
  # given and of a default look alike. So the replay is told how it was made (pyMzLib #83).
  skip_if(!nzchar(mz$replay_dir()), "the package's replay recordings are not installed")
  answer <- function(...) mz$json_parse(mz$replay_answer(c("peptidoform", "fragments", ...)))
  one <- answer("--accession", "P02768", "--max-mods", "1")
  expect_true(isTRUE(one$ok))
  expect_identical(one$data$peptide_count, 2)
  bare <- answer("--accession", "P02768")
  expect_identical(bare$data$peptide_count, 310)
  expect_identical(mz$replay_recorded_with_mismatch(c("max-mods" = "1"), list(accession = "P02768")),
    "recorded with --max-mods 1, which the call does not send")
})

test_that("max_length is sent only when given, and 0 is the wire's unbounded", {
  expect_false("--max-length" %in% build_args())
  expect_identical(paired(build_args(max_length = 0), "--max-length"), "0")
  expect_identical(paired(build_args(max_length = 30), "--max-length"), "30")
})

test_that("a max_length below min_length is refused rather than returning nothing", {
  expect_error(build_args(min_length = 7, max_length = 5),
    class = "mzlib_usage_error", contains = "below min_length"
  )
  # Only against a min_length the caller gave: mzLib's own is not restated here to check against.
  expect_identical(paired(build_args(max_length = 5), "--max-length"), "5")
  expect_identical(paired(build_args(min_length = 7, max_length = 0), "--max-length"), "0")
})

test_that("the modifications flag is not inverted", {
  # The R argument is `modifications`, the wire flag is `--no-modifications`. Getting the
  # polarity wrong silently returns bare sequences to a caller who asked for annotated ones, or
  # the reverse - a wrong answer that looks like a valid digest, which is why it is worth a test
  # of its own rather than being left to the live canary.
  expect_false("--no-modifications" %in% build_args(modifications = TRUE))
  expect_true("--no-modifications" %in% build_args(modifications = FALSE))
  expect_error(build_args(modifications = NA), class = "mzlib_usage_error")
  expect_error(build_args(modifications = "yes"), class = "mzlib_usage_error")
})

test_that("terminus is checked against mzLib's set, case-insensitively", {
  # A typo here would silently halve a fragment list rather than failing.
  expect_identical(paired(build_args(terminus = "both"), "--terminus"), "Both")
  for (value in c("n", "C", "None")) {
    expect_true(is.character(build_args(terminus = value)))
  }
  expect_error(build_args(terminus = "Nterm"),
    class = "mzlib_usage_error", contains = "must be one of"
  )
})

test_that("a UniProt accession is validated and upper-cased", {
  expect_identical(mz$peptidoform_normalise_accession("  p02768 "), "P02768")
  expect_identical(mz$peptidoform_normalise_accession("A0A0B4J2D5"), "A0A0B4J2D5")
  for (accession in list("", "   ", "P0", "not an accession", NULL, NA_character_)) {
    expect_error(mz$peptidoform_normalise_accession(accession), class = "mzlib_usage_error")
  }
})

test_that("a flag-like protease is refused", {
  expect_error(build_args(protease = "--no-modifications"),
    class = "mzlib_usage_error", contains = "may not begin with"
  )
  expect_error(build_args(protease = "  "), class = "mzlib_usage_error")
})

test_that("counts must be whole numbers in range", {
  expect_error(build_args(min_length = 0), class = "mzlib_usage_error", contains = "at least 1")
  expect_error(build_args(missed_cleavages = -1), class = "mzlib_usage_error")
  expect_error(build_args(max_isoforms = 0), class = "mzlib_usage_error")
  expect_error(build_args(min_length = 7.5), class = "mzlib_usage_error")
  expect_error(build_args(max_modifications = "2"), class = "mzlib_usage_error")
})

test_that("a large isoform cap is not written in scientific notation", {
  expect_identical(paired(build_args(max_isoforms = 100000), "--max-isoforms"), "100000")
})

# ---------------------------------------------------------------- printing

test_that("printing a digest shows its identity and census", {
  # mzLib #1114 removed the spurious ETD y series, so the print method no longer warns about
  # "y ions with no b ions" - there are none. What it still shows is the run's identity and the
  # modification census.
  output <- paste(capture.output(print(recorded_etd_digest())), collapse = "\n")
  expect_true(grepl("ALBU_HUMAN", output, fixed = TRUE))
  expect_true(grepl("peptidoforms over", output, fixed = TRUE))
  expect_false(grepl("y ions with no b ions", output, fixed = TRUE))
  expect_true(grepl("14 of 38", output, fixed = TRUE))
})

test_that("printing says which length settings were given and which were mzLib's", {
  digest <- recorded_digest()
  digest$min_length <- 1
  digest$max_length <- NA_real_
  digest$missed_cleavages <- NA_real_
  output <- paste(capture.output(print(digest)), collapse = "\n")
  expect_true(grepl("min_length 1, missed cleavages mzLib's default", output, fixed = TRUE),
    info = output
  )
})

test_that("printing a truncated digest says the list is incomplete", {
  digest <- recorded_digest()
  digest$peptides_at_isoform_cap <- 12
  output <- paste(capture.output(print(digest)), collapse = "\n")
  expect_true(grepl("hit the isoform cap", output, fixed = TRUE))
})

# ---------------------------------------------------------------- against UniProt

test_that("LIVE: albumin digests to the numbers the documentation quotes", {
  # The strongest single check in the offline-plus-live suite: it exercises the whole chain and
  # pins the two figures ?peptidoform_fragments quotes for the peptidoform-versus-sequence
  # distinction.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  digest <- tryCatch(
    peptidoform_fragments("P02768", timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )

  expect_identical(digest$accession, "P02768")
  expect_identical(digest$protease, "trypsin")
  expect_identical(digest$dissociation, "HCD")
  expect_identical(digest_distinct_base_sequences(digest), 202L)
  expect_identical(nrow(digest$peptides), 310L)
  # And live confirmation of the census figures quoted in ?census_explain.
  expect_identical(digest$census$applied, 14)
  expect_identical(digest$census$annotated, 38)
})

test_that("LIVE: min_length is what hides the short peptides", {
  # 202 at mzLib's default, 257 at 1 — a fifth of the digest lives below the default.
  #
  # Under the old "trypsin|P" default these were 195 and 243. PLAN.md once said 254 for the
  # latter; 243 was confirmed by driving the bridge from the shell and counting distinct
  # base_sequence values with no R involved, which is the only way to get ground truth that
  # cannot agree with a bug in this package. 257 was re-measured the same way on 2026-10-21
  # (bridge 0.5.0, mzLib d10017ce).
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  digest <- tryCatch(
    peptidoform_fragments("P02768", min_length = 1, timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  expect_identical(digest_distinct_base_sequences(digest), 257L)
})

test_that("LIVE: an unknown accession is an error, not an empty digest", {
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  condition <- tryCatch(
    {
      peptidoform_fragments("Q9ZZZ9", timeout = 300)
      NULL
    },
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e))),
    mzlib_error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_error"))
})

test_that("LIVE: an unknown protease names the alternatives", {
  # A stuck user reads the error message and nothing else, so the message has to carry the list.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  condition <- tryCatch(
    {
      peptidoform_fragments("P02768", protease = "trypsen", timeout = 300)
      NULL
    },
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e))),
    mzlib_error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_error"))
  expect_true(grepl("trypsin", conditionMessage(condition), fixed = TRUE),
    info = conditionMessage(condition)
  )
})

test_that("LIVE: the isoform cap truncates and says so", {
  # A truncated peptide list and a genuinely short one look identical from the outside, so the
  # only defence is that digest_truncated() is honest.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  capped <- tryCatch(
    peptidoform_fragments("P02768", max_modifications = 3, max_isoforms = 2, timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  expect_true(digest_truncated(capped))
  expect_true(capped$peptides_at_isoform_cap > 0)
  expect_true(grepl("hit the isoform cap", paste(capture.output(print(capped)), collapse = "\n"),
    fixed = TRUE
  ))
})

test_that("LIVE: modifications change the answer substantially", {
  # The control that shows the annotations are doing real work: modifications multiply
  # PEPTIDOFORMS, because each backbone gains one row per modification placement. The parity
  # test both siblings have had all along (pyMzLib test_modifications_change_the_answer_
  # substantially, mzLibRust modifications_change_the_answer_substantially); mzLibR had no
  # equivalent, having only ever tested this seam through the bug below.
  #
  # What this deliberately does NOT assert is that the distinct BASE sequences differ. They do
  # not, and should not: digestion produces backbones and modification decorates them, so
  # toggling modifications cannot add or remove a backbone. A test asserting otherwise stood
  # here until 2026-08-09 and passed only because mzLib had a defect - `--no-modifications` also
  # discarded UniProt's feature table, taking the signal-peptide cleavage boundaries with it and
  # costing albumin two peptides (pyMzLib#8, since fixed). That test pinned a bug rather than a
  # behaviour, so the fix broke it, which is exactly what such a test does. Deleted rather than
  # inverted: "backbones are unchanged" is a domain tautology nobody would write unprompted, and
  # if it ever needs guarding it belongs in mzLib's own suite, not in three bindings' canaries.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  with_mods <- tryCatch(
    peptidoform_fragments("P02768", max_modifications = 1, timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  bare <- peptidoform_fragments("P02768", modifications = FALSE, timeout = 600)

  expect_identical(nrow(bare$modifications), 0L)
  expect_true(
    nrow(with_mods$peptides) > nrow(bare$peptides),
    info = paste(nrow(with_mods$peptides), "peptidoforms with mods vs", nrow(bare$peptides), "bare")
  )
  expect_true(nrow(digest_modified_peptides(with_mods)) > 0L)
  expect_identical(nrow(digest_modified_peptides(bare)), 0L)
})

test_that("LIVE: the trypsin naming inversion changes the peptide count", {
  # smith-chem-wisc/mzLib#1106. mzLib's "trypsin|P" applies the Keil rule; plain "trypsin", the
  # default, does not — the reverse of the MaxQuant and Mascot convention. 195 against 202 on
  # albumin.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  keil <- tryCatch(
    peptidoform_fragments("P02768", protease = "trypsin|P", timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  expect_identical(digest_distinct_base_sequences(keil), 195L)
  plain <- peptidoform_fragments("P02768", protease = "trypsin", timeout = 600)
  expect_identical(digest_distinct_base_sequences(plain), 202L)
})

test_that("LIVE: the old defaults' answer is one pair of arguments away", {
  # NEWS.md promises it: protease = "trypsin|P", dissociation = "ETD" gives back what a bare call
  # gave before mzLibR left the defaults to mzLib.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  old <- tryCatch(
    peptidoform_fragments("P02768", protease = "trypsin|P", dissociation = "ETD", timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  expect_identical(old$protease, "trypsin|P")
  expect_identical(old$dissociation, "ETD")
  expect_identical(nrow(old$peptides), 303L)
  expect_identical(digest_distinct_base_sequences(old), 195L)
  expect_identical(sort(digest_fragments_by_series(old)$product_type), c("c", "zDot"))
})
