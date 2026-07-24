# In-silico digestion and fragmentation, backed by mzLib.
#
# Fetch a UniProt entry, apply its annotated modifications, digest it, and fragment every
# resulting peptidoform. The result is an `mzlibr_digest`: a few scalars plus three tidy
# data.frames - `peptides`, `fragments`, `modifications` - which is the shape R wants and which
# joins on `peptide_index`.
#
# Most of the comments in this file are about what the numbers mean rather than how they are
# computed, because every documented failure in this module is somebody reading a correct number
# as an answer to a question it does not answer.

# The mass of a proton, in daltons.
#
# **Not the hydrogen atom, which is 1.007825.** The difference is 0.55 mDa - 1.1 ppm at m/z 500,
# which on an Orbitrap is a match versus a miss. Libraries differ on this and rarely say which
# they used, so mzLibR says.
PROTON_MASS <- 1.00727646677

# mzLib's fragmentation termini. A typo here would silently halve a fragment list, so it is
# checked rather than passed through.
PEPTIDOFORM_TERMINI <- c("N", "C", "Both", "None")

# ---------------------------------------------------------------- validation

# A UniProt accession, upper-cased.
#
# Looser than the PRIDE check on purpose: UniProt accessions come in a 6-character and a
# 10-character form and the grammar has been extended before, so this rejects only what is
# obviously not one.
peptidoform_normalise_accession <- function(accession) {
  if (!is.character(accession) || length(accession) != 1L || is.na(accession)) {
    stop(mzlib_usage_error(
      "A UniProt accession is required, as a single string, e.g. 'P02768'."
    ))
  }
  candidate <- toupper(trimws(accession))
  if (!nzchar(candidate)) {
    stop(mzlib_usage_error("A UniProt accession is required, e.g. 'P02768'."))
  }
  if (!grepl("^[A-Z0-9]{6,10}$", candidate)) {
    stop(mzlib_usage_error(paste0(
      "'", accession, "' is not a valid UniProt accession. Expected six or ten alphanumeric ",
      "characters, e.g. 'P02768' or 'A0A0B4J2D5'."
    )))
  }
  candidate
}

# A non-empty option value that the bridge's parser will not read as another flag.
peptidoform_normalise_option <- function(name, value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(mzlib_usage_error(paste0(name, " must be a single non-empty string.")))
  }
  trimmed <- trimws(value)
  if (startsWith(trimmed, "-")) {
    stop(mzlib_usage_error(paste0(
      name, " may not begin with '-'; got '", trimmed, "'. That would be read as another option."
    )))
  }
  trimmed
}

# A whole number in range, refused early rather than at the far end of a subprocess.
peptidoform_normalise_count <- function(name, value, minimum) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value) ||
    value != round(value)) {
    stop(mzlib_usage_error(paste0(name, " must be a single whole number; got ",
      paste(deparse(value), collapse = " "), "."
    )))
  }
  if (value < minimum) {
    stop(mzlib_usage_error(paste0(name, " must be at least ", minimum, "; got ", value, ".")))
  }
  format(value, scientific = FALSE)
}

# ---------------------------------------------------------------- argument assembly

peptidoform_build_args <- function(accession, protease, dissociation, terminus, modifications,
                                   missed_cleavages, min_length, max_length, max_modifications,
                                   max_isoforms) {
  canonical <- peptidoform_normalise_accession(accession)

  wanted_terminus <- peptidoform_normalise_option("terminus", terminus)
  match <- PEPTIDOFORM_TERMINI[tolower(PEPTIDOFORM_TERMINI) == tolower(wanted_terminus)]
  if (length(match) != 1L) {
    stop(mzlib_usage_error(paste0(
      "terminus must be one of ", paste(PEPTIDOFORM_TERMINI, collapse = ", "),
      "; got '", terminus, "'."
    )))
  }

  minimum <- peptidoform_normalise_count("min_length", min_length, 1)

  # `--max-length 0` is the wire's spelling of "unbounded", so NULL and 0 are the same request.
  # Anything else must not be below min_length, which would silently return nothing.
  if (is.null(max_length)) {
    maximum <- "0"
  } else {
    maximum <- peptidoform_normalise_count("max_length", max_length, 0)
    if (as.numeric(maximum) > 0 && as.numeric(maximum) < as.numeric(minimum)) {
      stop(mzlib_usage_error(paste0(
        "max_length (", maximum, ") is below min_length (", minimum,
        "), which can only return nothing. Pass max_length = NULL for no upper bound."
      )))
    }
  }

  if (!is.logical(modifications) || length(modifications) != 1L || is.na(modifications)) {
    stop(mzlib_usage_error("modifications must be TRUE or FALSE."))
  }

  args <- c(
    "peptidoform", "fragments",
    "--accession", canonical,
    "--protease", peptidoform_normalise_option("protease", protease),
    "--dissociation", peptidoform_normalise_option("dissociation", dissociation),
    "--terminus", match,
    "--missed-cleavages", peptidoform_normalise_count("missed_cleavages", missed_cleavages, 0),
    "--min-length", minimum,
    "--max-length", maximum,
    "--max-mods", peptidoform_normalise_count("max_modifications", max_modifications, 0),
    "--max-isoforms", peptidoform_normalise_count("max_isoforms", max_isoforms, 1)
  )

  if (!modifications) {
    args <- c(args, "--no-modifications")
  }
  args
}

# ---------------------------------------------------------------- parsing

# One row per peptidoform.
#
# `peptide_index` is added, not on the wire: it is the key that joins `fragments` and
# `modifications` back to their peptide. Subsetting `peptides` and then joining on it is how a
# caller narrows a fragment list without hand-tracking positions.
peptidoform_parse_peptides <- function(entries) {
  if (length(entries) == 0L) {
    return(data.frame(
      peptide_index = integer(0), base_sequence = character(0), full_sequence = character(0),
      monoisotopic_mass = numeric(0), length = numeric(0), one_based_start = numeric(0),
      one_based_end = numeric(0), missed_cleavages = numeric(0), modification_count = numeric(0),
      fixed_charges = numeric(0), stringsAsFactors = FALSE
    ))
  }

  numeric_field <- function(name) {
    vapply(entries, wire_field, numeric(1L), name, "numeric", NA_real_)
  }

  data.frame(
    peptide_index = seq_along(entries),
    base_sequence = vapply(entries, wire_field, character(1L), "base_sequence", "character", NA_character_),
    full_sequence = vapply(entries, wire_field, character(1L), "full_sequence", "character", NA_character_),
    monoisotopic_mass = numeric_field("monoisotopic_mass"),
    length = numeric_field("length"),
    one_based_start = numeric_field("one_based_start"),
    one_based_end = numeric_field("one_based_end"),
    missed_cleavages = numeric_field("missed_cleavages"),
    modification_count = numeric_field("modification_count"),
    fixed_charges = numeric_field("fixed_charges"),
    stringsAsFactors = FALSE
  )
}

# One row per fragment ion, across every peptide.
#
# Long rather than nested. A digest of a whole protein has tens of thousands of fragments, and
# every question worth asking of them - counts by series, the ladder for one peptide, the ions
# above some mass - is a subset or a table of a long frame.
peptidoform_parse_fragments <- function(entries) {
  parts <- lapply(seq_along(entries), function(index) {
    fragments <- entries[[index]][["fragments"]]
    if (!is.list(fragments) || length(fragments) == 0L) {
      return(NULL)
    }
    data.frame(
      peptide_index = index,
      product_type = vapply(fragments, wire_field, character(1L), "product_type", "character", NA_character_),
      fragment_number = vapply(fragments, wire_field, numeric(1L), "fragment_number", "numeric", NA_real_),
      neutral_mass = vapply(fragments, wire_field, numeric(1L), "neutral_mass", "numeric", NA_real_),
      neutral_loss = vapply(fragments, wire_field, numeric(1L), "neutral_loss", "numeric", NA_real_),
      residue_position = vapply(fragments, wire_field, numeric(1L), "residue_position", "numeric", NA_real_),
      stringsAsFactors = FALSE
    )
  })

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(data.frame(
      peptide_index = integer(0), product_type = character(0), fragment_number = numeric(0),
      neutral_mass = numeric(0), neutral_loss = numeric(0), residue_position = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, parts)
}

# One row per applied modification.
peptidoform_parse_modifications <- function(entries) {
  parts <- lapply(seq_along(entries), function(index) {
    modifications <- entries[[index]][["modifications"]]
    if (!is.list(modifications) || length(modifications) == 0L) {
      return(NULL)
    }
    data.frame(
      peptide_index = index,
      one_based_residue = vapply(modifications, wire_field, numeric(1L), "one_based_residue", "numeric", NA_real_),
      # `null` on the wire for a residue modification, which is most of them - the field only
      # carries a value for a terminal modification. NA and not "" so the two are distinct.
      terminus = vapply(modifications, wire_field, character(1L), "terminus", "character", NA_character_),
      id = vapply(modifications, wire_field, character(1L), "id", "character", NA_character_),
      mass = vapply(modifications, wire_field, numeric(1L), "mass", "numeric", NA_real_),
      formal_charge = vapply(modifications, wire_field, numeric(1L), "formal_charge", "numeric", NA_real_),
      stringsAsFactors = FALSE
    )
  })

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(data.frame(
      peptide_index = integer(0), one_based_residue = numeric(0), terminus = character(0),
      id = character(0), mass = numeric(0), formal_charge = numeric(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, parts)
}

peptidoform_parse_census <- function(data) {
  by_type_entries <- data[["uniprot_features_by_type"]]
  if (!is.list(by_type_entries)) {
    by_type_entries <- list()
  }
  by_type <- data.frame(
    type = vapply(by_type_entries, wire_field, character(1L), "type", "character", NA_character_),
    count = vapply(by_type_entries, wire_field, numeric(1L), "count", "numeric", NA_real_),
    loaded = vapply(
      by_type_entries,
      function(entry) isTRUE(entry[["loaded"]]),
      logical(1L)
    ),
    stringsAsFactors = FALSE
  )

  unresolved <- data[["unresolved_modifications"]]
  unresolved <- if (is.list(unresolved) && length(unresolved) > 0L) {
    vapply(unresolved, function(value) as.character(value)[1L], character(1L))
  } else {
    character(0)
  }

  structure(
    list(
      # Distinct residue positions carrying at least one modification. A histone lists K9me1,
      # K9me2, K9me3 and K9ac at one residue - four modifications, one site - so this is always
      # the smaller number and is **not** a modification count. Conflating the two made H3.1
      # look as though 93 annotations had been dropped when all had loaded.
      sites = as.numeric(data[["annotated_modification_sites"]]),
      applied = as.numeric(data[["annotated_modifications_loaded"]]),
      annotated = as.numeric(data[["uniprot_annotated_features"]]),
      unresolved = unresolved,
      by_type = by_type
    ),
    class = "mzlibr_census"
  )
}

peptidoform_parse <- function(data) {
  entries <- data[["peptides"]]
  if (!is.list(entries)) {
    entries <- list()
  }

  scalar <- function(name, missing) {
    value <- data[[name]]
    if (is.null(value) || (length(value) == 1L && is.logical(value) && is.na(value))) missing else value
  }

  structure(
    list(
      accession = as.character(scalar("accession", NA_character_)),
      name = as.character(scalar("name", NA_character_)),
      full_name = as.character(scalar("full_name", NA_character_)),
      organism = as.character(scalar("organism", NA_character_)),
      sequence_length = as.numeric(scalar("sequence_length", NA_real_)),
      protease = as.character(scalar("protease", NA_character_)),
      dissociation = as.character(scalar("dissociation", NA_character_)),
      terminus = as.character(scalar("terminus", NA_character_)),
      modifications_applied = isTRUE(scalar("modifications_applied", FALSE)),
      max_modifications = as.numeric(scalar("max_modifications", NA_real_)),
      max_modification_isoforms = as.numeric(scalar("max_modification_isoforms", NA_real_)),
      peptides_at_isoform_cap = as.numeric(scalar("peptides_at_isoform_cap", 0)),
      census = peptidoform_parse_census(data),
      peptides = peptidoform_parse_peptides(entries),
      fragments = peptidoform_parse_fragments(entries),
      modifications = peptidoform_parse_modifications(entries)
    ),
    class = "mzlibr_digest"
  )
}

# ---------------------------------------------------------------- the public surface

#' Digest a protein and fragment its peptidoforms
#'
#' Fetches a UniProt entry, applies its annotated modifications, digests it with the named
#' protease, and computes fragment ions for every resulting peptidoform.
#'
#' @param accession A UniProt accession, e.g. `"P02768"` (serum albumin).
#' @param protease The protease, in **mzLib's naming**.
#'
#'   **Read this if you are coming from MaxQuant or Mascot.** mzLib's `"trypsin|P"` *applies* the
#'   Keil rule - it does not cleave before proline - while plain `"trypsin"` cleaves everywhere.
#'   That is the **reverse** of the MaxQuant and Mascot convention, where the `/P` suffix means
#'   "do cleave before proline". On albumin the two give **195** and **202** peptides, so the
#'   mistake is quiet and small enough to survive review (smith-chem-wisc/mzLib#1106).
#' @param dissociation The dissociation type, e.g. `"ETD"`, `"HCD"`, `"CID"`, `"ECD"`.
#'
#'   **`"ETD"` and `"ECD"` return three series - `c`, `zDot` *and* `y` - not two.** No
#'   fragmentation mechanism produces `y` without `b`, and about **a third** of every ETD
#'   fragment list is those `y` ions. Use [digest_fragments_by_series()] and select the series
#'   you mean rather than trusting a total (smith-chem-wisc/mzLib#1109; a fix is proposed in
#'   mzLib PR #1114, so check whether it has merged).
#' @param modifications Whether to apply UniProt's annotated modifications.
#'
#'   `FALSE` does more than drop modifications: it also discards **proteolysis products**, so
#'   the peptide list itself changes. On albumin two signal-peptide peptides disappear
#'   (smith-chem-wisc/pyMzLib#8).
#' @param missed_cleavages Maximum missed cleavage sites per peptide.
#' @param min_length Shortest peptide to keep.
#'
#'   The default of 7 **silently discards** everything shorter. Albumin goes from **195**
#'   distinct sequences at `min_length = 7` to **243** at `min_length = 1` - a fifth of the
#'   digest lives below the default. If you are looking for a short peptide and not finding it,
#'   look here first.
#' @param max_length Longest peptide to keep, or `NULL` for no limit.
#' @param max_modifications Maximum modifications considered per peptide.
#' @param max_isoforms Maximum modification isoforms generated per peptide.
#'
#'   **This cap truncates silently.** A truncated result and a genuinely short one look
#'   identical from the outside - histone H3.1 at four modifications loses about **30%**. Check
#'   [digest_truncated()] before treating a peptide list as exhaustive.
#' @param terminus Which terminus to fragment: `"N"`, `"C"`, `"Both"` or `"None"`.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_digest`: a list with the scalars describing the run, a `census` (see
#'   [census_explain()]), and three data.frames that join on `peptide_index` -
#'   `peptides`, `fragments` and `modifications`.
#'
#'   **`peptides` holds peptidoforms, not distinct sequences.** One row per
#'   sequence-and-modification-placement, so albumin at two modifications is **303** rows over
#'   **195** distinct sequences. Both are legitimate answers to "how many peptides" and they are
#'   not interchangeable; quoting one for the other is a large error, not a rounding one. See
#'   [digest_distinct_base_sequences()].
#'
#' @seealso [digest_fragments_by_series()], [digest_truncated()], [peptide_mz()]
#' @export
peptidoform_fragments <- function(accession, protease = "trypsin|P", dissociation = "ETD",
                                  modifications = TRUE, missed_cleavages = 2, min_length = 7,
                                  max_length = NULL, max_modifications = 2, max_isoforms = 1024,
                                  terminus = "Both", timeout = 300) {
  args <- peptidoform_build_args(
    accession, protease, dissociation, terminus, modifications,
    missed_cleavages, min_length, max_length, max_modifications, max_isoforms
  )
  digest <- peptidoform_parse(bridge_invoke(args, timeout = timeout))

  # The bridge reports max_modifications and max_modification_isoforms back, but not the three
  # settings that decide which peptides exist at all. They are stamped on here so that a digest
  # sitting in a workspace, or printed six months later, says what produced it - min_length in
  # particular, since its default of 7 silently removes a fifth of the digest and is the first
  # thing to check when an expected peptide is missing.
  digest$min_length <- as.numeric(min_length)
  digest$max_length <- if (is.null(max_length)) NA_real_ else as.numeric(max_length)
  digest$missed_cleavages <- as.numeric(missed_cleavages)
  digest
}

#' Whether a digest hit the isoform cap, and is therefore incomplete
#'
#' A short answer and a truncated answer look identical from the outside. Check this before
#' treating a peptide list as exhaustive.
#'
#' @param digest An [peptidoform_fragments()] result.
#' @return `TRUE` if any peptide hit `max_isoforms`.
#' @export
digest_truncated <- function(digest) {
  stopifnot(inherits(digest, "mzlibr_digest"))
  isTRUE(digest$peptides_at_isoform_cap > 0)
}

#' How many distinct base sequences a digest produced
#'
#' `nrow(digest$peptides)` counts **peptidoforms** - one per sequence-and-modification-placement.
#' This counts distinct sequences. On albumin at two modifications the two are **303** and
#' **195**.
#'
#' @param digest An [peptidoform_fragments()] result.
#' @return A single count.
#' @export
digest_distinct_base_sequences <- function(digest) {
  stopifnot(inherits(digest, "mzlibr_digest"))
  length(unique(digest$peptides$base_sequence))
}

#' Only the peptidoforms carrying at least one modification
#'
#' @param digest An [peptidoform_fragments()] result.
#' @return The rows of `digest$peptides` with `modification_count > 0`.
#' @export
digest_modified_peptides <- function(digest) {
  stopifnot(inherits(digest, "mzlibr_digest"))
  digest$peptides[!is.na(digest$peptides$modification_count) &
    digest$peptides$modification_count > 0, ]
}

#' Fragment ions per product type
#'
#' **Prefer this to `nrow(digest$fragments)` whenever the ion series matter, which for ETD is
#' always.** A bare total folds in two things that are not comparable ions: the spurious `y`
#' series mzLib emits for ETD (about a third of the total), and one extra full-length `z-dot` per
#' peptide.
#'
#' The `zDot` series runs `1..length`, not `1..length-1`. The extra ion numbered `length` is the
#' whole peptide minus NH2 - the N-Ca cleavage at residue 1 - which is **correct and
#' deliberate**, not a defect. It is absent when the peptide starts with proline.
#'
#' Note also that `zDot` counts come in **below** `length` per peptide rather than above it,
#' because z-dot ions are suppressed N-terminal to proline while the complementary `c` ions are
#' not. On albumin that is 138 proline sites, 138 suppressed `z-dot` and **0** suppressed `c`
#' (smith-chem-wisc/mzLib#1110).
#'
#' @param digest An [peptidoform_fragments()] result.
#' @return A data.frame of `product_type` and `n`, ordered by product type.
#' @export
digest_fragments_by_series <- function(digest) {
  stopifnot(inherits(digest, "mzlibr_digest"))
  if (nrow(digest$fragments) == 0L) {
    return(data.frame(product_type = character(0), n = integer(0), stringsAsFactors = FALSE))
  }
  counts <- table(digest$fragments$product_type)
  data.frame(
    product_type = names(counts),
    n = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

#' The m/z of intact peptides at a given charge
#'
#' Two conventions are handled explicitly, because getting either wrong is invisible in the
#' answer.
#'
#' **The proton mass, not the hydrogen atom.** 1.00727646677 against 1.007825 - a difference of
#' 0.55 mDa, or 1.1 ppm at m/z 500, which on an Orbitrap is a match versus a miss.
#'
#' **Fixed charges are not double-counted.** A peptide whose modification leaves a permanently
#' charged residue - trimethyl-lysine gives a quaternary ammonium, and UniProt records the delta
#' as 43.054227, which is C3H7 *minus an electron* - already carries that charge in its
#' `monoisotopic_mass`. Only `charge - fixed_charges` protons are added. Adding a full
#' complement would put a 2+ trimethylated peptide half a Thomson high, on the most important
#' histone modification there is.
#'
#' A peptide with a fixed charge is observable at that charge with no protonation at all, which
#' is why `charge` may not be below it.
#'
#' Note that **fragments carry `neutral_mass` and deliberately have no m/z**: a c or z ion
#' carries only the fixed charges within its own span, and per-fragment charge accounting does
#' not exist on this wire.
#'
#' @param peptides Rows of a `digest$peptides` data.frame.
#' @param charge Total charge, a positive whole number.
#' @return A numeric vector of m/z, one per row.
#' @export
peptide_mz <- function(peptides, charge) {
  if (!is.data.frame(peptides) ||
    !all(c("monoisotopic_mass", "fixed_charges") %in% names(peptides))) {
    stop(mzlib_usage_error(
      "peptides must be a data.frame from peptidoform_fragments()$peptides."
    ))
  }
  if (!is.numeric(charge) || length(charge) != 1L || is.na(charge) || !is.finite(charge) ||
    charge != round(charge) || charge < 1) {
    stop(mzlib_usage_error(paste0(
      "charge must be a positive whole number; got ", paste(deparse(charge), collapse = " "), "."
    )))
  }

  too_high <- !is.na(peptides$fixed_charges) & peptides$fixed_charges > charge
  if (any(too_high)) {
    stop(mzlib_usage_error(paste0(
      sum(too_high), " peptide(s) already carry more than ", charge,
      " fixed charge(s) from their modifications, so they cannot be observed at that charge. ",
      "The highest is ", max(peptides$fixed_charges, na.rm = TRUE), "."
    )))
  }

  (peptides$monoisotopic_mass + (charge - peptides$fixed_charges) * PROTON_MASS) / charge
}

#' Annotated features that could not be used
#'
#' `annotated - applied`. See [census_explain()] for why they were excluded, and why that
#' exclusion is correct.
#'
#' @param digest An [peptidoform_fragments()] result, or its `census`.
#' @return A single count.
#' @export
census_excluded <- function(digest) {
  census <- if (inherits(digest, "mzlibr_digest")) digest$census else digest
  if (!inherits(census, "mzlibr_census")) {
    stop(mzlib_usage_error("Expected a peptidoform_fragments() result or its census."))
  }
  max(0, census$annotated - census$applied)
}

#' What UniProt annotated, and what mzLib could actually use
#'
#' A plain-language account, because the alternative is a number arriving with no indication
#' that a rule was ever applied. For serum albumin, 14 modifications are applied out of 38
#' annotated.
#'
#' **The exclusion is correct - do not try to defeat it.** mzLib loads only `modified residue`
#' and `lipid moiety-binding region` annotations. On albumin the census reports the 24 excluded
#' features under one feature *type*, `glycosylation site` - that is the label you will see in
#' `census$by_type`, and it is the only granularity the census has. At UniProt's finer
#' modification-*name* level, 22 of those 24 are specifically `N-linked (Glc) (glycation)
#' lysine`; the census does not surface that, so the "22" is something you confirm by reading
#' the UniProt entry, not a number this tool reports. Either way the exclusion is right:
#' glycation and glycosylation are labile, heterogeneous adducts, so assigning one an exact mass
#' and a clean fragment ladder would describe a species you cannot observe.
#'
#' **The defect is the silence, not the exclusion** (smith-chem-wisc/mzLib#1112). Note also that
#' mzLib reads no qualifiers for any feature type. It cannot tell an annotation marked
#' `; in vitro`, or one that exists only in a disease variant (albumin's Redhill and Casebrook),
#' from any other - different grounds for exclusion needing different judgements, which the
#' census cannot make for you. Read the UniProt entry before concluding anything about a
#' specific site.
#'
#' @param digest An [peptidoform_fragments()] result, or its `census`.
#' @return A single string.
#' @export
census_explain <- function(digest) {
  census <- if (inherits(digest, "mzlibr_digest")) digest$census else digest
  if (!inherits(census, "mzlibr_census")) {
    stop(mzlib_usage_error("Expected a peptidoform_fragments() result or its census."))
  }

  excluded <- max(0, census$annotated - census$applied)
  if (excluded == 0 && length(census$unresolved) == 0L) {
    return(paste0(
      "All ", census$annotated, " annotated modifications were applied, across ",
      census$sites, " residue positions."
    ))
  }

  sentences <- paste0(
    census$applied, " of ", census$annotated, " annotated modifications were applied, across ",
    census$sites, " residue positions."
  )

  excluded_types <- census$by_type[!census$by_type$loaded, ]
  if (nrow(excluded_types) > 0L) {
    sentences <- c(sentences, paste0(
      "Excluded by type: ",
      paste0(excluded_types$count, " x ", excluded_types$type, collapse = ", "),
      " - mzLib loads only 'modified residue' and 'lipid moiety-binding region' annotations, so ",
      "these were dropped on feature type alone. The exclusion is usually right: a glycation or ",
      "glycosylation annotation describes a labile, heterogeneous adduct, so assigning it one ",
      "exact mass and a clean fragment ladder would invent a species you cannot observe. But the ",
      "reason is not reported, and the qualifier is not read - some annotations are marked ",
      "'in vitro' and some exist only in disease variants, which are different grounds for ",
      "exclusion needing different judgements from you. Read the annotations on the UniProt ",
      "entry before concluding anything about a specific site; this census can only tell you ",
      "the count (smith-chem-wisc/mzLib#1112)."
    ))
  }

  if (length(census$unresolved) > 0L) {
    sentences <- c(sentences, paste0(
      "Could not be resolved to a mass: ", paste(census$unresolved, collapse = ", "),
      " - annotated by UniProt but absent from its own modification list, so they were dropped."
    ))
  }

  paste(sentences, collapse = " ")
}

#' Print a digest
#'
#' A compact summary, and the place several warnings actually reach someone: printing a digest
#' and reading a fragment total off it is exactly what a person about to count spurious ETD
#' \code{y} ions does.
#'
#' @param x An [peptidoform_fragments()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_digest <- function(x, ...) {
  cat("<mzlibr_digest> ", x$accession, " ", x$name, " (", x$organism, ")\n", sep = "")
  cat("  ", format(x$sequence_length), " residues, ", x$protease, ", ", x$dissociation,
    ", terminus ", x$terminus, "\n",
    sep = ""
  )
  if (!is.null(x$min_length)) {
    cat("  min_length ", format(x$min_length),
      if (is.na(x$max_length)) "" else paste0(", max_length ", format(x$max_length)),
      ", ", format(x$missed_cleavages), " missed cleavages",
      ", max ", format(x$max_modifications), " modifications\n",
      sep = ""
    )
  }
  cat("  ", nrow(x$peptides), " peptidoforms over ",
    digest_distinct_base_sequences(x), " distinct sequences\n",
    sep = ""
  )

  series <- digest_fragments_by_series(x)
  if (nrow(series) > 0L) {
    cat("  ", nrow(x$fragments), " fragments: ",
      paste0(series$product_type, "=", series$n, collapse = ", "), "\n",
      sep = ""
    )
    # Put the warning where the mistake is made. Someone printing a digest and reading off a
    # fragment total is exactly the person about to count spurious y ions as real ones.
    if (any(series$product_type == "y") && !any(series$product_type == "b")) {
      cat("  ! ", series$n[series$product_type == "y"],
        " y ions with no b ions - no fragmentation mechanism produces that.\n",
        "    See ?digest_fragments_by_series (smith-chem-wisc/mzLib#1109).\n",
        sep = ""
      )
    }
  }

  if (digest_truncated(x)) {
    cat("  ! ", format(x$peptides_at_isoform_cap),
      " peptides hit the isoform cap - this list is incomplete.\n",
      sep = ""
    )
  }

  if (x$census$annotated > x$census$applied) {
    cat("  ! ", format(x$census$applied), " of ", format(x$census$annotated),
      " annotated modifications applied. See ?census_explain.\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Print a modification census
#'
#' Prints what [census_explain()] returns.
#'
#' @param x The `census` of an [peptidoform_fragments()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_census <- function(x, ...) {
  cat(census_explain(x), "\n")
  invisible(x)
}
