# Predicting peptide properties: retention time, fragment intensities, CCS, detectability.
#
# mzLib ships clients for 37 published models on the Koina inference server (koina.wilhelmlab.org),
# across five families. This module calls them.
#
# Ported from pyMzLib's `prediction.py`, which decided the verbs, the wire fields and the caveats;
# what changes here is the projection into R's idiom - a data.frame rather than a map of arrays, a
# list column for the ragged fragment arrays, NA rather than None, S3 classes with print methods.
#
# Five things the design turns on, all of them mzLib's constraints rather than choices made here:
#
#  1. Koina is a public, shared, community-run GPU. The throttling defaults are mzLib's and are not
#     raised; `max_batches` and `throttle_ms` exist for a genuinely large job.
#  2. A prediction is an opinion, not a measurement. Nothing here is FDR-anything.
#  3. Retention time is not always in minutes - most of these models return iRT, a dimensionless
#     index, and only Chronologer_RT returns absolute minutes. The unit is a value, not prose.
#  4. Fragment arrays are ragged: Koina returns a fixed-width grid with -1 for impossible ions and
#     mzLib drops those, so each row's arrays are as long as THAT peptide's possible ions.
#  5. A peptide that cannot be predicted still gets a row, with NA and a warning, so predictions
#     line up with the peptides that were sent.

# ---------------------------------------------------------------- parsing

# What a model requires of one optional input parameter.
#
# A TRI-STATE, and mzLib expresses it as a nullable set whose emptiness means the opposite of what
# it looks like: null means "not applicable" and an EMPTY set means "required, any value". Reading
# the raw collection is how you conclude that Prosit_2020_intensity_CID accepts any collision
# energy - it accepts none, being fixed at NCE 35 - and that Prosit_2020_intensity_HCD accepts
# none, when it requires one.
prediction_parse_constraint <- function(payload) {
  if (!is.list(payload)) {
    return(structure(
      list(requirement = "not_applicable", values = NULL),
      class = "mzlibr_constraint"
    ))
  }

  values <- payload[["values"]]
  values <- if (is.list(values) && length(values) > 0L) unlist(values, use.names = FALSE) else NULL

  structure(
    list(
      requirement = as.character(
        wire_field(payload, "requirement", "character", "not_applicable")
      ),
      values = values
    ),
    class = "mzlibr_constraint"
  )
}

prediction_parse_models <- function(data) {
  entries <- data[["models"]]
  if (!is.list(entries) || length(entries) == 0L) {
    return(data.frame(
      model = character(0), family = character(0), verb = character(0),
      stringsAsFactors = FALSE
    ))
  }

  models <- data.frame(
    model = vapply(entries, wire_field, character(1L), "model", "character", NA_character_),
    family = vapply(entries, wire_field, character(1L), "family", "character", NA_character_),
    verb = vapply(entries, wire_field, character(1L), "verb", "character", NA_character_),
    type = vapply(entries, wire_field, character(1L), "type", "character", NA_character_),
    min_peptide_length = vapply(
      entries, wire_field, numeric(1L), "min_peptide_length", "numeric", NA_real_
    ),
    max_peptide_length = vapply(
      entries, wire_field, numeric(1L), "max_peptide_length", "numeric", NA_real_
    ),
    max_batch_size = vapply(
      entries, wire_field, numeric(1L), "max_batch_size", "numeric", NA_real_
    ),
    retention_time_unit = vapply(
      entries, wire_field, character(1L), "retention_time_unit", "character", NA_character_
    ),
    number_of_predicted_fragment_ions = vapply(
      entries, wire_field, numeric(1L),
      "number_of_predicted_fragment_ions", "numeric", NA_real_
    ),
    stringsAsFactors = FALSE
  )

  # List columns, because a filtered frame must keep them - the same reason readers' `views` is one.
  models$allowed_unimod_ids <- lapply(entries, function(entry) {
    ids <- entry[["allowed_unimod_ids"]]
    if (is.list(ids) && length(ids) > 0L) unlist(ids, use.names = FALSE) else integer(0)
  })
  for (parameter in c(
    "precursor_charge", "collision_energy", "instrument_type", "fragmentation_type"
  )) {
    models[[parameter]] <- lapply(entries, function(entry) {
      prediction_parse_constraint(entry[[parameter]])
    })
  }

  models
}

prediction_parse_predictions <- function(data) {
  caveats <- data[["caveats"]]
  caveats <- if (is.list(caveats) && length(caveats) > 0L) {
    vapply(caveats, function(c) as.character(c)[1L], character(1L), USE.NAMES = FALSE)
  } else {
    character(0)
  }

  column_names <- data[["column_names"]]
  column_names <- if (is.list(column_names) && length(column_names) > 0L) {
    vapply(column_names, function(n) as.character(n)[1L], character(1L), USE.NAMES = FALSE)
  } else {
    character(0)
  }

  structure(
    list(
      model = as.character(wire_field(data, "model", "character", NA_character_)),
      row_count = as.numeric(wire_field(data, "row_count", "numeric", NA_real_)),
      failed_row_count = as.numeric(
        wire_field(data, "failed_row_count", "numeric", NA_real_)
      ),
      retention_time_unit = as.character(
        wire_field(data, "retention_time_unit", "character", NA_character_)
      ),
      collisional_cross_section_unit = as.character(
        wire_field(data, "collisional_cross_section_unit", "character", NA_character_)
      ),
      intensity_scale = as.character(
        wire_field(data, "intensity_scale", "character", NA_character_)
      ),
      caveats = caveats,
      column_names = column_names,
      # Reuses the readers table parser, which already keeps a list column when a cell is itself
      # an array - exactly the shape the ragged fragment arrays need.
      predictions = readers_parse_records_table(data),
      output = readers_parse_output(data)
    ),
    class = "mzlibr_predictions"
  )
}

# ---------------------------------------------------------------- input

# The peptide table as the tab-separated text the bridge reads from stdin.
#
# On stdin rather than argv because a real prediction run is thousands of peptides and argv has a
# hard ceiling of roughly 32 KB - the same reason PRIDE's explicit file selection goes this way.
prediction_build_stdin <- function(peptides, columns) {
  if (is.character(peptides)) {
    peptides <- data.frame(sequence = peptides, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(peptides)) {
    stop(mzlib_usage_error(
      "peptides must be a character vector of sequences, or a data.frame of input columns."
    ))
  }
  if (nrow(peptides) == 0L) {
    stop(mzlib_usage_error("At least one peptide is required."))
  }

  # The crosslink family names its first column alpha_sequence; a caller who passed a bare
  # character vector still means "the sequence", so it is accepted under either name.
  if (!columns[1L] %in% names(peptides) && "sequence" %in% names(peptides)) {
    peptides[[columns[1L]]] <- peptides[["sequence"]]
  }
  if (!columns[1L] %in% names(peptides)) {
    stop(mzlib_usage_error(paste0(
      "peptides must have a '", columns[1L], "' column for this verb."
    )))
  }

  unknown <- setdiff(names(peptides), c(columns, "sequence"))
  if (length(unknown) > 0L) {
    stop(mzlib_usage_error(paste0(
      "Unknown input column(s): ", paste(unknown, collapse = ", "),
      ". This verb reads: ", paste(columns, collapse = ", "), "."
    )))
  }

  cells <- lapply(columns, function(column) {
    if (!column %in% names(peptides)) {
      return(rep("", nrow(peptides)))
    }
    values <- peptides[[column]]
    # An unset optional is an EMPTY cell, not the text "NA": the latter would reach the bridge as
    # a value and fail to parse as a number, reporting a confusing error about a column the caller
    # never set.
    ifelse(is.na(values), "", as.character(values))
  })

  rows <- do.call(paste, c(cells, sep = "\t"))
  paste0(paste(c(paste(columns, collapse = "\t"), rows), collapse = "\n"), "\n")
}

prediction_build_args <- function(verb, model, out, max_batches, throttle_ms) {
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(trimws(model))) {
    stop(mzlib_usage_error(
      "A model name is required. prediction_models() lists them with their constraints."
    ))
  }

  args <- c("predict", verb, "--model", trimws(model))

  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }

  # Deliberately not defaulted to anything faster than mzLib's own values: Koina is a shared
  # community server, and a binding that maximised throughput by default would be spending someone
  # else's GPU time.
  for (knob in list(
    list(name = "max-batches", value = max_batches, minimum = 1),
    list(name = "throttle-ms", value = throttle_ms, minimum = 0)
  )) {
    if (is.null(knob$value)) next
    if (!is.numeric(knob$value) || length(knob$value) != 1L || is.na(knob$value) ||
      knob$value != round(knob$value) || knob$value < knob$minimum) {
      stop(mzlib_usage_error(paste0(
        gsub("-", "_", knob$name), " must be a whole number of ", knob$minimum,
        " or more, or NULL for mzLib's default; got ",
        paste(deparse(knob$value), collapse = " "), "."
      )))
    }
    args <- c(args, paste0("--", knob$name), formatC(knob$value, format = "d"))
  }

  args
}

prediction_predict <- function(verb, model, peptides, columns, out, max_batches, throttle_ms,
                               timeout) {
  args <- prediction_build_args(verb, model, out, max_batches, throttle_ms)
  stdin <- prediction_build_stdin(peptides, columns)
  prediction_parse_predictions(bridge_invoke(args, stdin = stdin, timeout = timeout))
}

# ---------------------------------------------------------------- the public surface

#' Every Koina model mzLib can call
#'
#' Enumerated from mzLib itself rather than from a list maintained here, so it reflects the
#' installed version and cannot go stale. **This does not touch the network** - it describes what
#' mzLib can call, not what the server currently answers.
#'
#' @param family Restrict to one family, e.g. `"retention_time"`. `NULL` returns all 37.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return A data.frame with one row per model: `model`, `family`, `verb`, `type`, the peptide
#'   length bounds, `max_batch_size`, `retention_time_unit`,
#'   `number_of_predicted_fragment_ions`, an `allowed_unimod_ids` list column, and one
#'   `mzlibr_constraint` list column per optional input parameter.
#'
#' @section The constraints are more restrictive than they look:
#'
#' `max_peptide_length` is 30 for most models, which excludes a substantial fraction of real
#' tryptic peptides - and every one of them comes back as `NA` with a warning rather than an error.
#'
#' `allowed_unimod_ids` being **empty means the model accepts no modifications at all**, which is a
#' real answer rather than a missing one.
#'
#' @section Constraints are a tri-state, not a list:
#'
#' `precursor_charge`, `collision_energy`, `instrument_type` and `fragmentation_type` each hold an
#' `mzlibr_constraint` whose `requirement` is one of `"not_applicable"` (the model has no such
#' input - do not send it), `"any_value_required"` (you must send one, any value) or `"one_of"`
#' (you must send one of its `values`).
#'
#' mzLib encodes all three states in one nullable set, where `NULL` means *not applicable* and an
#' *empty* set means *required, any value*. Read straight, that makes `Prosit_2020_intensity_CID`
#' look permissive about collision energy - it accepts none, being fixed at NCE 35 - and
#' `Prosit_2020_intensity_HCD` look like it accepts none, when it requires one.
#'
#' @seealso [prediction_retention_time()], [prediction_fragments()]
#' @export
prediction_models <- function(family = NULL, timeout = 60) {
  args <- c("predict", "models")
  if (!is.null(family)) {
    if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !nzchar(trimws(family))) {
      stop(mzlib_usage_error("family must be a single non-empty name, or NULL for every family."))
    }
    args <- c(args, "--family", trimws(family))
  }

  prediction_parse_models(bridge_invoke(args, timeout = timeout))
}

#' Predict elution, one row per peptide
#'
#' @param model A retention-time model name from [prediction_models()], e.g. `"Prosit_2019_irt"`.
#' @param peptides A character vector of sequences in mzLib `FullSequence` notation, or a
#'   data.frame with a `sequence` column.
#' @param out Write a tab-separated table here and return only a summary.
#' @param max_batches In-flight requests. `NULL` uses mzLib's polite default.
#' @param throttle_ms Delay between request chunks. `NULL` uses mzLib's default.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_predictions`. `predictions` is a data.frame; `retention_time_unit` says what
#'   the numbers mean.
#'
#' @section Most of these models do not return minutes:
#'
#' They return **iRT** - indexed retention time, a dimensionless scale anchored to standard
#' peptides. Only `Chronologer_RT` returns absolute minutes. An iRT of 130 looks exactly like a
#' plausible 130-minute gradient, so check `retention_time_unit` before plotting anything against
#' a real gradient, and fit the iRT-to-minutes line on shared peptides first.
#'
#' @seealso [prediction_models()]
#' @export
prediction_retention_time <- function(model, peptides, out = NULL, max_batches = NULL,
                                      throttle_ms = NULL, timeout = NULL) {
  prediction_predict(
    "retention-time", model, peptides, "sequence", out, max_batches, throttle_ms, timeout
  )
}

#' Predict MS2 fragment m/z and relative intensity
#'
#' @param model A fragment-intensity model name, e.g. `"Prosit_2020_intensity_HCD"`.
#' @param peptides A character vector of sequences, or a data.frame whose columns may include
#'   `sequence`, `precursor_charge`, `collision_energy`, `instrument_type` and
#'   `fragmentation_type`. Check [prediction_models()] for which this model requires.
#' @param out Write a tab-separated table here; the arrays become `;`-joined lists.
#' @param max_batches In-flight requests.
#' @param throttle_ms Delay between request chunks.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_predictions` whose `fragment_annotations`, `fragment_mz` and
#'   `fragment_intensity` columns are **list columns**, one array per row.
#'
#' @section The arrays are ragged, and that is the correct answer:
#'
#' Koina returns a fixed-width grid with `-1` marking ions that cannot exist for a given peptide,
#' and mzLib drops those. So each row's three arrays are as long as **that peptide's** possible
#' ions - 28 for a short tryptic peptide from a model whose nominal count is 174 - and two rows in
#' one call will differ. They are index-aligned within a row, and are a list column precisely
#' because they cannot be a matrix.
#'
#' Intensities are **relative**, on Koina's own 0-1 scale, and are not comparable with a measured
#' intensity or between models.
#'
#' @seealso [prediction_models()]
#' @export
prediction_fragments <- function(model, peptides, out = NULL, max_batches = NULL,
                                 throttle_ms = NULL, timeout = NULL) {
  prediction_predict(
    "fragments", model, peptides,
    c(
      "sequence", "precursor_charge", "collision_energy",
      "instrument_type", "fragmentation_type"
    ),
    out, max_batches, throttle_ms, timeout
  )
}

#' Predict collisional cross-section, in square angstroms
#'
#' @param model A CCS model name, e.g. `"IM2Deep"`.
#' @param peptides A data.frame with `sequence` and `precursor_charge` columns. **CCS depends on
#'   charge**, so a bare character vector is rarely what you want.
#' @param out Write a tab-separated table here.
#' @param max_batches In-flight requests.
#' @param throttle_ms Delay between request chunks.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_predictions`. The unit is **square angstroms, never 1/K0** - converting to
#'   the reduced mobility a timsTOF reports needs drift-gas temperature and pressure, which mzLib
#'   does not carry, so no conversion is offered rather than a guessed one.
#'
#' @seealso [prediction_models()]
#' @export
prediction_ccs <- function(model, peptides, out = NULL, max_batches = NULL,
                           throttle_ms = NULL, timeout = NULL) {
  prediction_predict(
    "ccs", model, peptides, c("sequence", "precursor_charge"),
    out, max_batches, throttle_ms, timeout
  )
}

#' Predict flyability, as four class probabilities
#'
#' @param model A detectability model name, e.g. `"pfly_2024_fine_tuned"`.
#' @param peptides A character vector of sequences, or a data.frame with a `sequence` column.
#' @param out Write a tab-separated table here.
#' @param max_batches In-flight requests.
#' @param throttle_ms Delay between request chunks.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_predictions` with `not_detectable`, `low_detectability`,
#'   `intermediate_detectability` and `high_detectability` columns. They are a **distribution over
#'   classes** and sum to 1 per peptide - not an expected intensity, and not a probability of
#'   detection.
#'
#' @seealso [prediction_models()]
#' @export
prediction_detectability <- function(model, peptides, out = NULL, max_batches = NULL,
                                     throttle_ms = NULL, timeout = NULL) {
  prediction_predict(
    "detectability", model, peptides, "sequence", out, max_batches, throttle_ms, timeout
  )
}

#' Predict MS2 intensities for a crosslinked peptide pair
#'
#' @param model A crosslink model name, e.g. `"Prosit_2023_intensity_XL_CMS2"`.
#' @param pairs A data.frame with `alpha_sequence` and, for most models, `beta_sequence` columns,
#'   plus `precursor_charge` and `collision_energy`.
#' @param out Write a tab-separated table here.
#' @param max_batches In-flight requests.
#' @param throttle_ms Delay between request chunks.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_predictions`, with ragged fragment list columns as for
#'   [prediction_fragments()].
#'
#' @section This family takes a different sequence language:
#'
#' The other four functions accept mzLib's `FullSequence` notation and convert it; the crosslink
#' models **reject it** and require raw UNIMOD brackets - `K[UNIMOD:1896]`. That is mzLib's
#' constraint, not a choice made here, and it is repeated in the result's `caveats`.
#'
#' @seealso [prediction_models()]
#' @export
prediction_crosslink_fragments <- function(model, pairs, out = NULL, max_batches = NULL,
                                           throttle_ms = NULL, timeout = NULL) {
  prediction_predict(
    "crosslink-fragments", model, pairs,
    c("alpha_sequence", "beta_sequence", "precursor_charge", "collision_energy"),
    out, max_batches, throttle_ms, timeout
  )
}

#' Whether a model parameter should be sent at all
#'
#' @param x An `mzlibr_constraint`, from a [prediction_models()] row.
#' @return `TRUE` unless the model has no such input.
#' @export
constraint_applicable <- function(x) {
  if (!inherits(x, "mzlibr_constraint")) {
    stop(mzlib_usage_error("x must be an mzlibr_constraint from a prediction_models() row."))
  }
  !identical(x$requirement, "not_applicable")
}

# ---------------------------------------------------------------- printing

#' Print a model parameter constraint
#'
#' @param x An `mzlibr_constraint`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_constraint <- function(x, ...) {
  cat(switch(x$requirement,
    not_applicable = "<constraint> not applicable - this model has no such input\n",
    any_value_required = "<constraint> required; any value accepted\n",
    one_of = paste0(
      "<constraint> required; one of: ", paste(x$values, collapse = ", "), "\n"
    ),
    paste0("<constraint> ", x$requirement, "\n")
  ))
  invisible(x)
}

#' Print a prediction table
#'
#' @param x A prediction result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_predictions <- function(x, ...) {
  cat("<mzlibr_predictions> ", x$model, "\n", sep = "")
  cat("  ", format(x$row_count), " peptides", sep = "")
  if (!is.na(x$failed_row_count) && x$failed_row_count > 0) {
    # Not an error: too long, an untrained modification, or a missing required parameter are
    # normal outcomes, and the row survives so the alignment does.
    cat(", ", format(x$failed_row_count), " could not be predicted (see the warning column)",
      sep = ""
    )
  }
  cat("\n")
  for (unit in c("retention_time_unit", "collisional_cross_section_unit", "intensity_scale")) {
    if (!is.na(x[[unit]]) && nzchar(x[[unit]])) {
      cat("  ", unit, ": ", x[[unit]], "\n", sep = "")
    }
  }
  for (caveat in x$caveats) {
    cat("  ! ", caveat, "\n", sep = "")
  }
  if (!is.null(x$output)) {
    cat("  written to ", x$output$path, " (", x$output$format, ", ",
      format(x$output$row_count), " rows)\n",
      sep = ""
    )
  }
  invisible(x)
}
