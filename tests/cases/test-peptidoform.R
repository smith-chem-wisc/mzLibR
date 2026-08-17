# Digestion, fragmentation, and the modification census.
#
# The fixture is a two-peptide ETD digest of serum albumin, and it happens to contain two
# fragmentation traps at once: a zDot series numbered past the c series, and zDot counts that come
# in *below* their own maximum because of proline suppression. (mzLib #1114 removed the spurious
# ETD y series that used to be a third trap here.) Every
# one of those is asserted here rather than only described, because these are the numbers the
# documentation quotes.

recorded_digest_data <- function() {
  mz$json_parse(paste(
    readLines(fixture_path("peptidoform_P02768_small.json"), warn = FALSE),
    collapse = "\n"
  ))
}

recorded_digest <- function() {
  mz$peptidoform_parse(recorded_digest_data())
}

build_args <- function(accession = "P02768", protease = "trypsin|P", dissociation = "ETD",
                       terminus = "Both", modifications = TRUE, missed_cleavages = 2,
                       min_length = 7, max_length = NULL, max_modifications = 2,
                       max_isoforms = 1024) {
  mz$peptidoform_build_args(
    accession, protease, dissociation, terminus, modifications,
    missed_cleavages, min_length, max_length, max_modifications, max_isoforms
  )
}

# ---------------------------------------------------------------- parsing

test_that("the digest carries the run's own description", {
  digest <- recorded_digest()
  expect_true(inherits(digest, "mzlibr_digest"))
  expect_identical(digest$accession, "P02768")
  expect_identical(digest$name, "ALBU_HUMAN")
  expect_identical(digest$organism, "Homo sapiens")
  expect_identical(digest$sequence_length, 609)
  expect_identical(digest$protease, "trypsin|P")
  expect_identical(digest$dissociation, "ETD")
  expect_identical(digest$terminus, "Both")
  expect_true(digest$modifications_applied)
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
  series <- mz$digest_fragments_by_series(recorded_digest())
  present <- series$product_type

  expect_false("y" %in% present)
  expect_false("b" %in% present)
  expect_identical(sort(present), c("c", "zDot"))
})

test_that("the c series runs 1..length-1", {
  digest <- recorded_digest()
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
  digest <- recorded_digest()
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
  digest <- recorded_digest()

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
  digest <- recorded_digest()

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
  # which begins at residue 42 of albumin - so a position read against the protein would be
  # wrong by 41 and still look perfectly plausible.
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
  digest <- recorded_digest()
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

test_that("the default arguments match the siblings", {
  args <- build_args()
  expect_identical(args[1:2], c("peptidoform", "fragments"))
  paired <- function(flag) args[which(args == flag) + 1L]
  expect_identical(paired("--accession"), "P02768")
  expect_identical(paired("--protease"), "trypsin|P")
  expect_identical(paired("--dissociation"), "ETD")
  expect_identical(paired("--terminus"), "Both")
  expect_identical(paired("--missed-cleavages"), "2")
  expect_identical(paired("--min-length"), "7")
  expect_identical(paired("--max-mods"), "2")
  expect_identical(paired("--max-isoforms"), "1024")
  expect_false("--no-modifications" %in% args)
})

test_that("max_length NULL is the wire's unbounded", {
  expect_identical(build_args()[which(build_args() == "--max-length") + 1L], "0")
  args <- build_args(max_length = 30)
  expect_identical(args[which(args == "--max-length") + 1L], "30")
})

test_that("a max_length below min_length is refused rather than returning nothing", {
  expect_error(build_args(min_length = 7, max_length = 5),
    class = "mzlib_usage_error", contains = "below min_length"
  )
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
  expect_identical(build_args(terminus = "both")[which(build_args() == "--terminus") + 1L], "Both")
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
  args <- build_args(max_isoforms = 100000)
  expect_identical(args[which(args == "--max-isoforms") + 1L], "100000")
})

# ---------------------------------------------------------------- printing

test_that("printing a digest shows its identity and census", {
  # mzLib #1114 removed the spurious ETD y series, so the print method no longer warns about
  # "y ions with no b ions" - there are none. What it still shows is the run's identity and the
  # modification census.
  output <- paste(capture.output(print(recorded_digest())), collapse = "\n")
  expect_true(grepl("ALBU_HUMAN", output, fixed = TRUE))
  expect_true(grepl("peptidoforms over", output, fixed = TRUE))
  expect_false(grepl("y ions with no b ions", output, fixed = TRUE))
  expect_true(grepl("14 of 38", output, fixed = TRUE))
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
  expect_identical(digest_distinct_base_sequences(digest), 195L)
  expect_identical(nrow(digest$peptides), 303L)
  # And live confirmation of the census figures quoted in ?census_explain.
  expect_identical(digest$census$applied, 14)
  expect_identical(digest$census$annotated, 38)
})

test_that("LIVE: min_length is what hides the short peptides", {
  # 195 at the default of 7, 243 at 1 — a quarter of the digest lives below the default.
  #
  # PLAN.md said 254 here. It is 243, confirmed by driving the bridge from the shell and
  # counting distinct base_sequence values with no R involved, which is the only way to get
  # ground truth that cannot agree with a bug in this package. 254 matches no neighbouring
  # setting either: plain trypsin at min_length 1 gives 257, three missed cleavages gives 324.
  # It was simply wrong, and it had never shipped — neither pyMzLib nor mzLibRust quotes a
  # min_length figure at all, so mzLibR is the first binding to document this trap with a
  # number, and now with a verified one.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  digest <- tryCatch(
    peptidoform_fragments("P02768", min_length = 1, timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  expect_identical(digest_distinct_base_sequences(digest), 243L)
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
  # smith-chem-wisc/mzLib#1106. mzLib's "trypsin|P" applies the Keil rule; plain "trypsin" does
  # not — the reverse of the MaxQuant and Mascot convention. 195 against 202 on albumin.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  options(mzlibr.bridge = live_bridge)
  on.exit(options(mzlibr.bridge = NULL), add = TRUE)

  plain <- tryCatch(
    peptidoform_fragments("P02768", protease = "trypsin", timeout = 600),
    mzlib_service_unavailable = function(e) skip(paste("UniProt unavailable:", conditionMessage(e)))
  )
  expect_identical(digest_distinct_base_sequences(plain), 202L)
})
