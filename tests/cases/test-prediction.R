# Predicting peptide properties.
#
# The payloads are recorded from the real bridge against the real Koina server, so a wire-shape
# change shows up here as a parse failure rather than as a fixture that agrees with an R file and
# with nothing else.
#
# Two of these tests exist because of traps rather than features: `mzlibr_constraint` encodes a
# tri-state mzLib expresses as a nullable set whose emptiness means the opposite of what it looks
# like, and the fragment arrays are ragged, so a caller who expects a matrix gets a wrong answer
# rather than an error.

prediction_payload <- function(name) {
  mz$json_parse(paste(
    readLines(fixture_path(name), warn = FALSE),
    collapse = "\n"
  ))$data
}

recorded_models <- function() {
  mz$prediction_parse_models(prediction_payload("prediction_models.json"))
}

# ---------------------------------------------------------------- the catalogue

test_that("every model mzLib can call becomes a row", {
  models <- recorded_models()
  expect_identical(nrow(models), 37L)
  expect_true(is.data.frame(models))
  expect_false(is.factor(models$model))
})

test_that("the five families are all represented", {
  models <- recorded_models()
  expect_identical(
    sort(unique(models$family)),
    sort(c(
      "collisional_cross_section", "crosslink_intensity", "detectability",
      "fragment_intensity", "retention_time"
    ))
  )
})

test_that("every model names the verb that calls it", {
  models <- recorded_models()
  expect_true(all(nzchar(models$verb)))
  expect_false(any(is.na(models$verb)))
})

test_that("the retention-time unit is per model, not per family", {
  # The distinction that makes mzLib's bare IsIndexed boolean insufficient on the wire: one model
  # in the family returns absolute minutes and the rest return a dimensionless index.
  models <- recorded_models()
  unit <- function(name) models[models$model == name, ]$retention_time_unit

  expect_identical(unit("Prosit_2019_irt"), "indexed_retention_time")
  expect_identical(unit("Chronologer_RT"), "minutes")
  # ...and it is meaningless outside that family.
  expect_true(is.na(unit("IM2Deep")))
})

test_that("a constraint distinguishes not-applicable from required-any", {
  # The trap this type exists for. mzLib expresses both as a nullable set: NULL means "this model
  # has no such input", empty means "required, any value". Reading the raw collection makes CID
  # look permissive and HCD look impossible, which is backwards for both.
  models <- recorded_models()
  energy <- function(name) models[models$model == name, ]$collision_energy[[1]]

  expect_identical(energy("Prosit_2020_intensity_HCD")$requirement, "any_value_required")
  expect_true(mz$constraint_applicable(energy("Prosit_2020_intensity_HCD")))

  expect_identical(energy("Prosit_2020_intensity_CID")$requirement, "not_applicable")
  expect_false(mz$constraint_applicable(energy("Prosit_2020_intensity_CID")))
})

test_that("a constraint lists its values when it has them", {
  models <- recorded_models()
  altimeter <- models[models$model == "Altimeter_2024_intensities", ]$collision_energy[[1]]

  expect_identical(altimeter$requirement, "one_of")
  expect_identical(min(altimeter$values), 20)
  expect_identical(max(altimeter$values), 40)

  unispec <- models[models$model == "UniSpec", ]$instrument_type[[1]]
  expect_identical(unispec$requirement, "one_of")
  expect_true("LUMOS" %in% unispec$values)
})

test_that("a model that accepts no modifications says so with an empty list", {
  # An empty allowed-UNIMOD list is a real answer, not a missing one.
  models <- recorded_models()
  expect_identical(length(models[models$model == "pfly_2024_fine_tuned", ]$allowed_unimod_ids[[1]]), 0L)
  expect_true(length(models[models$model == "Prosit_2019_irt", ]$allowed_unimod_ids[[1]]) > 0L)
})

test_that("constraints survive subsetting the models frame", {
  # A list column for the same reason readers' views is one.
  models <- recorded_models()
  some <- models[models$family == "fragment_intensity", ]
  expect_true("collision_energy" %in% names(some))
  expect_true(all(vapply(some$collision_energy, inherits, logical(1L), "mzlibr_constraint")))
})

test_that("a blank family is refused before anything is spawned", {
  for (bad in list("", "   ", NA_character_, 42)) {
    expect_error(mz$prediction_models(bad), class = "mzlib_usage_error")
  }
})

test_that("printing a constraint explains the requirement rather than showing a set", {
  models <- recorded_models()
  for (name in c("Prosit_2020_intensity_HCD", "Prosit_2020_intensity_CID", "Altimeter_2024_intensities")) {
    output <- paste(
      capture.output(print(models[models$model == name, ]$collision_energy[[1]])),
      collapse = "\n"
    )
    expect_true(grepl("constraint", output, fixed = TRUE), info = output)
  }
})

# ---------------------------------------------------------------- predictions

test_that("a prediction table parses into a data.frame", {
  result <- mz$prediction_parse_predictions(prediction_payload("prediction_rt_irt.json"))

  expect_true(inherits(result, "mzlibr_predictions"))
  expect_identical(result$model, "Prosit_2019_irt")
  expect_identical(result$row_count, 2)
  expect_identical(nrow(result$predictions), 2L)
})

test_that("an iRT model is not reported as minutes", {
  result <- mz$prediction_parse_predictions(prediction_payload("prediction_rt_irt.json"))

  expect_identical(result$retention_time_unit, "indexed_retention_time")
  expect_true(any(grepl("iRT", result$caveats, fixed = TRUE)))
})

test_that("fragment arrays are ragged and index-aligned within a row", {
  # Koina returns a fixed-width grid with -1 for impossible ions and mzLib drops those, so each
  # row is as long as ITS peptide's possible ions. A list column is what R has for this; expecting
  # a matrix gives a wrong answer rather than an error.
  result <- mz$prediction_parse_predictions(prediction_payload("prediction_fragments.json"))

  lengths_seen <- vapply(result$predictions$fragment_mz, length, integer(1L))
  expect_true(length(unique(lengths_seen)) > 1L)
  # The model's published count is 174; a short tryptic peptide gets a fraction of it.
  expect_true(all(lengths_seen < 174L))
  expect_identical(
    vapply(result$predictions$fragment_intensity, length, integer(1L)),
    lengths_seen
  )
  expect_identical(
    vapply(result$predictions$fragment_annotations, length, integer(1L)),
    lengths_seen
  )
})

test_that("intensities are declared relative", {
  result <- mz$prediction_parse_predictions(prediction_payload("prediction_fragments.json"))
  expect_identical(result$intensity_scale, "relative")
  expect_true(any(grepl("RELATIVE", result$caveats, fixed = TRUE)))
})

test_that("a peptide that cannot be predicted still gets a row", {
  # Prosit_2020_intensity_HCD requires a collision energy. Omitting it must not lose the row, or
  # predictions would no longer line up with the peptides that were sent.
  result <- mz$prediction_parse_predictions(prediction_payload("prediction_fragments_warned.json"))

  expect_identical(result$row_count, 1)
  expect_identical(result$failed_row_count, 1)
  expect_true(grepl("CollisionEnergy", result$predictions$warning[1], fixed = TRUE))
})

test_that("the four detectability classes sum to one", {
  result <- mz$prediction_parse_predictions(prediction_payload("prediction_detectability.json"))

  total <- sum(
    result$predictions$not_detectable,
    result$predictions$low_detectability,
    result$predictions$intermediate_detectability,
    result$predictions$high_detectability
  )
  expect_true(abs(total - 1) < 1e-5)
})

test_that("printing a prediction shows the unit and the failures", {
  warned <- mz$prediction_parse_predictions(prediction_payload("prediction_fragments_warned.json"))
  output <- paste(capture.output(print(warned)), collapse = "\n")
  expect_true(grepl("could not be predicted", output, fixed = TRUE), info = output)

  irt <- mz$prediction_parse_predictions(prediction_payload("prediction_rt_irt.json"))
  output <- paste(capture.output(print(irt)), collapse = "\n")
  expect_true(grepl("retention_time_unit: indexed_retention_time", output, fixed = TRUE))
})

# ---------------------------------------------------------------- input assembly

test_that("a character vector of sequences becomes a one-column table", {
  stdin <- mz$prediction_build_stdin(c("PEPTIDEK", "ELVISLIVESK"), "sequence")
  expect_identical(stdin, "sequence\nPEPTIDEK\nELVISLIVESK\n")
})

test_that("an unset optional column is an empty cell, not the text NA", {
  # "NA" would reach the bridge as a value and fail to parse as a number, reporting a confusing
  # error about a column the caller never set.
  stdin <- mz$prediction_build_stdin(
    data.frame(sequence = "PEPTIDEK", precursor_charge = 2, stringsAsFactors = FALSE),
    c("sequence", "precursor_charge", "collision_energy")
  )
  expect_identical(stdin, "sequence\tprecursor_charge\tcollision_energy\nPEPTIDEK\t2\t\n")
})

test_that("an explicit NA is also an empty cell", {
  stdin <- mz$prediction_build_stdin(
    data.frame(sequence = "PEPTIDEK", collision_energy = NA_real_, stringsAsFactors = FALSE),
    c("sequence", "collision_energy")
  )
  expect_false(grepl("NA", stdin, fixed = TRUE))
})

test_that("a bare sequence vector still fills the crosslink alpha column", {
  # The crosslink family names its first column alpha_sequence, and a caller passing sequences
  # still means "the sequence".
  stdin <- mz$prediction_build_stdin("PEPTIDEK", c("alpha_sequence", "beta_sequence"))
  expect_identical(stdin, "alpha_sequence\tbeta_sequence\nPEPTIDEK\t\n")
})

test_that("an unknown input column is refused rather than ignored", {
  # A caller who wrote 'charge' where the contract says 'precursor_charge' asked for something
  # specific; predicting at a default charge would hand back a plausible wrong answer.
  expect_error(
    mz$prediction_build_stdin(
      data.frame(sequence = "PEPTIDEK", charge = 2, stringsAsFactors = FALSE),
      c("sequence", "precursor_charge")
    ),
    class = "mzlib_usage_error", contains = "Unknown input column"
  )
})

test_that("an empty peptide set is refused before anything is spawned", {
  expect_error(
    mz$prediction_build_stdin(
      data.frame(sequence = character(0), stringsAsFactors = FALSE), "sequence"
    ),
    class = "mzlib_usage_error", contains = "At least one peptide"
  )
})

test_that("a blank model name is refused with a pointer to the catalogue", {
  for (bad in list("", "   ", NA_character_, 42)) {
    expect_error(
      mz$prediction_build_args("retention-time", bad, NULL, NULL, NULL),
      class = "mzlib_usage_error", contains = "prediction_models()"
    )
  }
})

test_that("the politeness defaults are left alone unless asked", {
  # Koina is a shared community server. A binding that raised the throughput defaults would be
  # spending someone else's GPU time without being asked.
  args <- mz$prediction_build_args("retention-time", "Prosit_2019_irt", NULL, NULL, NULL)
  expect_identical(args, c("predict", "retention-time", "--model", "Prosit_2019_irt"))
})

test_that("a politeness knob below its floor is refused", {
  expect_error(
    mz$prediction_build_args("retention-time", "m", NULL, 0, NULL),
    class = "mzlib_usage_error", contains = "max_batches"
  )
  expect_error(
    mz$prediction_build_args("retention-time", "m", NULL, NULL, -1),
    class = "mzlib_usage_error", contains = "throttle_ms"
  )
  expect_error(
    mz$prediction_build_args("retention-time", "m", NULL, 1.5, NULL),
    class = "mzlib_usage_error", contains = "whole number"
  )
})

test_that("the politeness knobs are sent when asked", {
  args <- mz$prediction_build_args("fragments", "m", NULL, 200, 250)
  expect_true(all(c("--max-batches", "200", "--throttle-ms", "250") %in% args))
})
