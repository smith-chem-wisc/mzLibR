# Differential abundance, set up by mzLib: the experiment's design, and every peptide in every sample.
#
#   factor_spec(), covariate_spec()  the terms of a model, as mzLib's FactorSpec and CovariateSpec
#   diff_design()                    AnalysisDesign.ByStratum, DefaultContrasts, Estimability: the
#                                    design matrix per stratum, and the files stats_fit() reads
#   diff_observations()              FlashLfqObservations.FromPeptideTable and
#                                    SharedPeptideNormalization.Measure: every peptide in every
#                                    biological sample, as log2 intensity or NA with the reason
#
# Ported from pyMzLib's `diff.py` (mzLib #1450, #1455). Every number is mzLib's: this file writes the
# model or the runs as lines on stdin, and projects each result into R objects. A matrix is an R
# matrix with its sample and column names; a table is a data.frame, `records`, whose array cells
# (a contrast's weights, a peptide's protein groups) are list columns.

# The first pyMzLib whose bridge has the diff verbs (their specs' since.pymzlib).
DIFF_SINCE <- "0.6.0"

# Every value of diff_design()'s `estimability` column, as mzLib's AnalysisDesign writes them. Only
# "estimable" rows have `fitted_weights`.
DIFF_ESTIMABILITIES <- c("estimable", "not_estimable:confounded", "not_estimable:rank_deficient")

# Every value of diff_observations()' `state` column, mzLib's ObservationState from most to least
# informative. Only "quantified" and "mbr_transferred" rows can carry a value.
DIFF_OBSERVATION_STATES <- c(
  "quantified", "mbr_transferred", "ambiguous_peak", "identified_not_quantified", "not_detected"
)

# The values diff_observations() takes as `basis`: mzLib's QuantBasis machine names. There is no
# default; you choose.
DIFF_QUANT_BASES <- c("msms_only", "mbr_kept")

# ---------------------------------------------------------------- shared parsing

# A wire table as a data.frame, where a column whose cells are arrays becomes a list column of
# typed vectors, an empty array as a zero-length vector and a `null` element as NA. A `null` where
# an array belongs is NULL in its cell (diff design's fitted_weights for a contrast that is not
# estimable). Never shortened. `arrays` names the spec's array columns and their element type
# ("numeric" or "character"), so a column stays a list column even when every cell is null or
# holds one value; any other column with an array cell is typed by its values.
wire_array_table <- function(data, arrays = character(0)) {
  columns <- data[["columns"]]
  order <- wire_strings(data[["column_names"]])
  if (length(order) == 0L && is.list(columns)) {
    order <- names(columns)
  }
  built <- lapply(order, function(name) {
    values <- if (is.list(columns)) columns[[name]] else NULL
    if (!is.list(values) || length(values) == 0L) {
      return(if (name %in% names(arrays)) I(list()) else logical(0))
    }
    if (name %in% names(arrays) || any(vapply(values, is.list, logical(1L)))) {
      flat <- unlist(lapply(values, function(v) if (is.list(v)) v else list()), recursive = FALSE)
      numeric_column <- if (name %in% names(arrays)) {
        identical(arrays[[name]], "numeric")
      } else {
        all(vapply(flat, function(e) wire_null(e) || is.numeric(e), logical(1L)))
      }
      return(I(lapply(values, function(v) {
        if (!is.list(v)) {
          return(NULL)
        }
        if (numeric_column) {
          vapply(v, function(e) if (wire_null(e)) NA_real_ else as.numeric(e), numeric(1L))
        } else {
          vapply(v, function(e) if (wire_null(e)) NA_character_ else as.character(e), character(1L))
        }
      })))
    }
    unlist(lapply(values, function(v) if (wire_null(v)) NA else v), use.names = FALSE)
  })
  names(built) <- order
  lengths_seen <- vapply(built, length, integer(1L))
  if (length(unique(lengths_seen)) > 1L) {
    stop(mzlib_protocol_error(paste0(
      "The table's columns have different lengths (",
      paste(paste0(order, "=", lengths_seen), collapse = ", "), "), so they cannot form a table."
    )))
  }
  as.data.frame(built, stringsAsFactors = FALSE, optional = TRUE)
}

# A `{name: number}` map whose values may be null, as a named numeric vector (null as NA).
wire_named_numbers <- function(value) {
  if (!is.list(value) || length(value) == 0L) {
    return(structure(numeric(0), names = character(0)))
  }
  out <- vapply(value, function(v) if (wire_null(v)) NA_real_ else as.numeric(v), numeric(1L))
  names(out) <- names(value)
  out
}

# A wire matrix (an array of rows) as an R matrix with these row and column names.
wire_matrix <- function(rows, row_names, column_names) {
  values <- if (is.list(rows)) lapply(rows, function(r) vapply(r, function(e) if (wire_null(e)) NA_real_ else as.numeric(e), numeric(1L))) else list()
  out <- matrix(unlist(values), nrow = length(values), ncol = length(column_names), byrow = TRUE)
  if (nrow(out) == length(row_names)) rownames(out) <- row_names
  colnames(out) <- column_names
  out
}

# ---------------------------------------------------------------- the model's terms

diff_text <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(mzlib_usage_error(paste0(what, " must be a non-empty string; got ", paste(deparse(value), collapse = " "), ".")))
  }
  if (grepl("[\t\r\n]", value)) {
    stop(mzlib_usage_error(paste0(what, " contains a tab or a line break, which separate the fields: '", value, "'.")))
  }
  trimws(value)
}

diff_optional_text <- function(value, what) {
  if (is.null(value)) "" else diff_text(value, what)
}

diff_finite <- function(value, what) {
  if (is.null(value)) {
    return("")
  }
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(mzlib_usage_error(paste0(what, " must be a finite number; got ", paste(deparse(value), collapse = " "), ".")))
  }
  stats_number_text(value)
}

diff_whole <- function(value, what) {
  if (is.null(value) || (length(value) == 1L && is.na(value))) {
    return("")
  }
  if (!is.numeric(value) || length(value) != 1L || value != round(value)) {
    stop(mzlib_usage_error(paste0(what, " must be a whole number; got ", paste(deparse(value), collapse = " "), ".")))
  }
  formatC(value, format = "d")
}

# A stdin line, with trailing fields left out when they were not given.
diff_line <- function(...) {
  fields <- c(...)
  while (length(fields) > 0L && !nzchar(fields[length(fields)])) {
    fields <- fields[-length(fields)]
  }
  paste(fields, collapse = "\t")
}

#' A condition factor of a model, as mzLib's FactorSpec
#'
#' Names a sample-sheet column for [diff_design()] to treat-code against a reference level.
#'
#' @param name The sample sheet column holding each sample's level.
#' @param reference_level The level every other level is compared against, or `NULL` to compare
#'   every pair of levels (mzLib then warns, because the comparisons are many).
#'
#' @return An `mzlibr_factor_spec`: a list of `name` and `reference_level` (`NA` when none).
#'
#' @seealso [diff_design()], [covariate_spec()]
#' @examples
#' factor_spec("malat_load", "500ng")
#' factor_spec("condition")
#' @export
factor_spec <- function(name, reference_level = NULL) {
  structure(
    list(
      name = diff_text(name, "A factor name"),
      reference_level = if (is.null(reference_level)) NA_character_ else diff_text(reference_level, "A reference level")
    ),
    class = "mzlibr_factor_spec"
  )
}

#' A numeric covariate of a model, as mzLib's CovariateSpec
#'
#' Names a sample-sheet column for [diff_design()] to enter as `(value - centre) / scale`. Leave a
#' field `NULL` and nothing is sent for it, so mzLib's own default applies.
#'
#' @param name The sample sheet column holding each sample's value, in `unit`.
#' @param unit The stored unit, e.g. `"years"`.
#' @param scale Stored units per scaled unit, e.g. `10` years per decade. `NULL`: mzLib's 1.
#' @param scaled_unit The unit a slope is per, e.g. `"decade"`. `NULL`: `scale` stored units
#'   (`"years"` at scale 1, `"10 years"` at scale 10).
#' @param centre The value subtracted first, in `unit`. `NULL`: the stratum's mean.
#'
#' @return An `mzlibr_covariate_spec`: a list of `name`, `unit`, `scale`, `scaled_unit` and
#'   `centre`, each of the last three `NA` when not given (mzLib's default applies).
#'
#' @seealso [diff_design()], [factor_spec()]
#' @examples
#' covariate_spec("age", "years", scale = 10, scaled_unit = "decade", centre = 50)
#' covariate_spec("bmi", "kg/m2")
#' @export
covariate_spec <- function(name, unit, scale = NULL, scaled_unit = NULL, centre = NULL) {
  name <- diff_text(name, "A covariate name")
  unit <- diff_text(unit, paste0("The unit of covariate '", name, "'"))
  diff_finite(scale, paste0("The scale of covariate '", name, "'"))
  diff_finite(centre, paste0("The centre of covariate '", name, "'"))
  structure(
    list(
      name = name,
      unit = unit,
      scale = if (is.null(scale)) NA_real_ else as.numeric(scale),
      scaled_unit = if (is.null(scaled_unit)) NA_character_ else diff_text(scaled_unit, paste0("The scaled unit of covariate '", name, "'")),
      centre = if (is.null(centre)) NA_real_ else as.numeric(centre)
    ),
    class = "mzlibr_covariate_spec"
  )
}

# One element or a list of them, as a list. A bare spec is one term, not its fields.
diff_terms <- function(value, class) {
  if (is.null(value)) {
    return(list())
  }
  if (inherits(value, class)) {
    return(list(value))
  }
  if (!is.list(value) && !is.character(value)) {
    stop(mzlib_usage_error(paste0("Expected a list of terms; got ", class(value)[1L], ".")))
  }
  as.list(value)
}

# The model as the bridge's stdin lines: one term per line, its fields tab-separated, an empty
# optional field left out so that mzLib's default applies.
diff_model_stdin <- function(factors, covariates, interactions, strata, batch) {
  lines <- character(0)
  for (factor in diff_terms(factors, "mzlibr_factor_spec")) {
    if (is.character(factor) && length(factor) == 1L) {
      factor <- factor_spec(factor)
    }
    if (!inherits(factor, "mzlibr_factor_spec")) {
      stop(mzlib_usage_error("Each factor is a factor_spec() or a column name."))
    }
    lines <- c(lines, diff_line(
      "factor", diff_text(factor$name, "A factor name"),
      if (is.na(factor$reference_level)) "" else diff_text(factor$reference_level, "A reference level")
    ))
  }
  for (covariate in diff_terms(covariates, "mzlibr_covariate_spec")) {
    if (!inherits(covariate, "mzlibr_covariate_spec")) {
      stop(mzlib_usage_error("Each covariate is a covariate_spec()."))
    }
    lines <- c(lines, diff_line(
      "covariate", covariate$name, covariate$unit,
      if (is.na(covariate$scale)) "" else diff_finite(covariate$scale, paste0("The scale of covariate '", covariate$name, "'")),
      if (is.na(covariate$scaled_unit)) "" else covariate$scaled_unit,
      if (is.na(covariate$centre)) "" else diff_finite(covariate$centre, paste0("The centre of covariate '", covariate$name, "'"))
    ))
  }
  if (!is.null(interactions) && !is.list(interactions)) {
    stop(mzlib_usage_error("interactions takes a list of pairs of factor names, e.g. list(c(\"condition\", \"sex\"))."))
  }
  for (pair in if (is.null(interactions)) list() else interactions) {
    if (!is.character(pair) || length(pair) != 2L) {
      stop(mzlib_usage_error("Each interaction is a pair of factor names, e.g. c(\"condition\", \"sex\")."))
    }
    lines <- c(lines, diff_line(
      "interaction", diff_text(pair[[1L]], "An interaction's factor"), diff_text(pair[[2L]], "An interaction's factor")
    ))
  }
  if (!is.null(strata)) {
    if (!is.character(strata)) {
      stop(mzlib_usage_error("strata takes a character vector of sample sheet column names."))
    }
    for (column in strata) {
      lines <- c(lines, diff_line("stratum", diff_text(column, "A stratum column")))
    }
  }
  if (!is.null(batch)) {
    lines <- c(lines, diff_line("batch", diff_text(batch, "The batch column")))
  }
  # No trailing newline: system2(input =) ends the last line itself. No model sends nothing.
  if (length(lines) == 0L) NULL else paste(lines, collapse = "\n")
}

# The runs as the bridge's stdin lines: file, sample, then fraction and technical replicate when
# given. `runs` is a data.frame with `file` and `sample` columns (or the result's own `file_name`
# and `sample_id`), and optionally `fraction` and `technical_replicate`.
diff_runs_stdin <- function(runs) {
  if (!is.data.frame(runs)) {
    stop(mzlib_usage_error(paste0(
      "runs must be a data.frame with file and sample columns, and optionally fraction and ",
      "technical_replicate, one row per run; got ", class(runs)[1L], "."
    )))
  }
  file <- if ("file" %in% names(runs)) runs$file else runs$file_name
  sample <- if ("sample" %in% names(runs)) runs$sample else runs$sample_id
  if (is.null(file) || is.null(sample)) {
    stop(mzlib_usage_error("runs needs a file and a sample column: which biological sample each run belongs to."))
  }
  if (nrow(runs) == 0L) {
    stop(mzlib_usage_error("diff_observations() needs the runs: one row per file of the peptide table."))
  }
  lines <- vapply(seq_len(nrow(runs)), function(i) {
    diff_line(
      diff_text(as.character(file[[i]]), "A run's file name"),
      diff_text(as.character(sample[[i]]), "A run's sample"),
      diff_whole(if ("fraction" %in% names(runs)) runs$fraction[[i]] else NULL, "A run's fraction"),
      diff_whole(if ("technical_replicate" %in% names(runs)) runs$technical_replicate[[i]] else NULL, "A run's technical replicate")
    )
  }, character(1L))
  paste(lines, collapse = "\n")
}

diff_build_design_request <- function(samples, factors, covariates, interactions, strata, batch, out) {
  args <- c("diff", "design", "--samples", readers_normalise_path(samples, "samples"))
  stdin <- diff_model_stdin(factors, covariates, interactions, strata, batch)
  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  list(args = args, stdin = stdin)
}

diff_build_observations_request <- function(peptides, runs, basis, protein_groups, samples, normalize,
                                            min_peptides, out, responses_out) {
  if (!is.character(basis) || length(basis) != 1L || !basis %in% DIFF_QUANT_BASES) {
    stop(mzlib_usage_error(paste0(
      "basis must be one of ", paste(DIFF_QUANT_BASES, collapse = ", "), "; got ",
      paste(deparse(basis), collapse = " "), ". mzLib has no default: whether match-between-runs ",
      "values enter is your choice."
    )))
  }
  stats_flag(normalize, "normalize")
  args <- c("diff", "observations", "--peptides", readers_normalise_path(peptides, "peptides"), "--basis", basis)
  stdin <- diff_runs_stdin(runs)
  if (!is.null(protein_groups)) {
    args <- c(args, "--protein-groups", readers_normalise_path(protein_groups, "protein_groups"))
  }
  if (!is.null(samples)) {
    args <- c(args, "--samples", readers_normalise_path(samples, "samples"))
  }
  if (isTRUE(normalize)) {
    args <- c(args, "--normalize")
  }
  if (!is.null(min_peptides)) {
    args <- c(args, "--min-peptides", stats_whole(min_peptides, "min_peptides"))
  }
  if (!is.null(out)) {
    args <- c(args, "--out", readers_normalise_path(out, "out"))
  }
  if (!is.null(responses_out)) {
    args <- c(args, "--responses-out", readers_normalise_path(responses_out, "responses_out"))
  }
  list(args = args, stdin = stdin)
}

# ---------------------------------------------------------------- parsing

diff_parse_stratum <- function(block) {
  sample_ids <- wire_strings(block[["sample_ids"]])
  column_names <- wire_strings(block[["column_names"]])
  fitted_names <- wire_strings(block[["fitted_column_names"]])
  not_run <- wire_objects(block[["not_run"]], empty = c("label", "reason"))
  structure(
    list(
      stratum = wire_field(block, "stratum", "character", NA_character_),
      sample_count = stats_count(block, "sample_count"),
      sample_ids = sample_ids,
      column_names = column_names,
      fitted_column_names = fitted_names,
      rank = stats_count(block, "rank"),
      aliased_terms = wire_strings(block[["aliased_terms"]]),
      warnings = wire_strings(block[["warnings"]]),
      not_run = not_run,
      covariate_centres = wire_named_numbers(block[["covariate_centres"]]),
      matrix = wire_matrix(block[["matrix"]], sample_ids, column_names),
      fitted_matrix = wire_matrix(block[["fitted_matrix"]], sample_ids, fitted_names),
      design_file = wire_field(block, "design_file", "character", NA_character_),
      contrasts_file = wire_field(block, "contrasts_file", "character", NA_character_)
    ),
    class = "mzlibr_design_stratum"
  )
}

diff_parse_design <- function(data) {
  factors <- lapply(if (is.list(data[["factors"]])) data[["factors"]] else list(), function(f) {
    structure(list(
      name = wire_field(f, "name", "character", NA_character_),
      reference_level = wire_field(f, "reference_level", "character", NA_character_)
    ), class = "mzlibr_factor_spec")
  })
  covariates <- lapply(if (is.list(data[["covariates"]])) data[["covariates"]] else list(), function(c) {
    structure(list(
      name = wire_field(c, "name", "character", NA_character_),
      unit = wire_field(c, "unit", "character", NA_character_),
      scale = wire_field(c, "scale", "numeric", NA_real_),
      scaled_unit = wire_field(c, "scaled_unit", "character", NA_character_),
      centre = wire_field(c, "centre", "numeric", NA_real_)
    ), class = "mzlibr_covariate_spec")
  })
  interactions <- lapply(if (is.list(data[["interactions"]])) data[["interactions"]] else list(), function(p) {
    c(wire_field(p, "first", "character", NA_character_), wire_field(p, "second", "character", NA_character_))
  })
  strata <- lapply(if (is.list(data[["strata"]])) data[["strata"]] else list(), diff_parse_stratum)
  names(strata) <- vapply(strata, function(s) s$stratum, character(1L))
  written <- data[["written"]]
  structure(
    list(
      samples_file = wire_field(data, "samples_file", "character", NA_character_),
      sample_count = stats_count(data, "sample_count"),
      factors = factors,
      covariates = covariates,
      interactions = interactions,
      stratum_columns = wire_strings(data[["stratum_columns"]]),
      batch_column = wire_field(data, "batch_column", "character", NA_character_),
      unused_columns = wire_strings(data[["unused_columns"]]),
      stratum_count = stats_count(data, "stratum_count"),
      strata = strata,
      estimability_counts = proteins_parse_counts(data[["estimability_counts"]]),
      written = if (is.list(written)) {
        list(
          directory = wire_field(written, "directory", "character", NA_character_),
          file_count = stats_count(written, "file_count")
        )
      } else {
        NULL
      },
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data, c(weights = "numeric", fitted_weights = "numeric")),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_diff_design"
  )
}

diff_parse_observations <- function(data) {
  runs <- wire_objects(data[["runs"]], empty = c("file_name", "sample_id", "fraction", "technical_replicate"))
  groups <- wire_objects(data[["protein_groups"]],
    empty = c("name", "q_value", "is_decoy", "is_contaminant", "is_entrapment", "genes")
  )
  normalization <- data[["normalization"]]
  normalization <- if (is.list(normalization)) normalization else list()
  written <- data[["written"]]
  responses <- data[["responses_written"]]
  minimum <- data[["min_peptides"]]
  structure(
    list(
      peptides_file = wire_field(data, "peptides_file", "character", NA_character_),
      protein_groups_file = wire_field(data, "protein_groups_file", "character", NA_character_),
      samples_file = wire_field(data, "samples_file", "character", NA_character_),
      basis = wire_field(data, "basis", "character", NA_character_),
      normalize = isTRUE(data[["normalize"]]),
      min_peptides = if (wire_null(minimum)) NA_real_ else as.numeric(minimum),
      run_count = stats_count(data, "run_count"),
      runs = runs,
      sample_count = stats_count(data, "sample_count"),
      samples = wire_strings(data[["samples"]]),
      samples_without_values = wire_strings(data[["samples_without_values"]]),
      peptide_count = stats_count(data, "peptide_count"),
      protein_group_count = stats_count(data, "protein_group_count"),
      protein_groups = groups,
      state_counts = proteins_parse_counts(data[["state_counts"]]),
      normalization = list(
        setting = wire_field(normalization, "setting", "character", NA_character_),
        reference_set_size = stats_count(normalization, "reference_set_size"),
        minimum_reference_set_size = stats_count(normalization, "minimum_reference_set_size"),
        shifts = wire_named_numbers(normalization[["shifts"]]),
        warnings = wire_strings(normalization[["warnings"]])
      ),
      warnings = wire_strings(data[["warnings"]]),
      written = if (is.list(written)) {
        list(
          path = wire_field(written, "path", "character", NA_character_),
          format = wire_field(written, "format", "character", NA_character_),
          row_count = stats_count(written, "row_count")
        )
      } else {
        NULL
      },
      responses_written = if (is.list(responses)) {
        list(
          path = wire_field(responses, "path", "character", NA_character_),
          format = wire_field(responses, "format", "character", NA_character_),
          row_count = stats_count(responses, "row_count"),
          sample_count = stats_count(responses, "sample_count")
        )
      } else {
        NULL
      },
      row_count = stats_count(data, "row_count"),
      column_names = wire_strings(data[["column_names"]]),
      records = wire_array_table(data, c(protein_groups = "character")),
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_observations"
  )
}

# ---------------------------------------------------------------- the public surface

#' Build each stratum's design matrix and default contrasts from a sample sheet, with mzLib
#'
#' mzLib's `AnalysisDesign.ByStratum` splits the samples by stratum, and in each stratum builds one
#' additive model: an intercept, each factor treatment-coded against its reference level, batch,
#' each covariate as `(value - centre) / scale`, and any interactions. The matrix equals R's
#' `model.matrix()` for the same formula. It then lists `DefaultContrasts` - each level against its
#' reference, every pair when there is no reference, each covariate's slope - and judges whether
#' each can be estimated at all from these samples.
#'
#' @param samples The sample sheet: a TSV with a header row and one row per biological sample. First
#'   column the sample id; the other columns are available by header. A column the model does not
#'   name is ignored (`unused_columns`).
#' @param factors Condition factors, in column order: a list of [factor_spec()], or a character
#'   vector of column names (no reference: every pair of levels is compared), or a mix in a list.
#' @param covariates Numeric covariates, in column order: a list of [covariate_spec()], each with
#'   its stored unit.
#' @param interactions Two-factor interactions, as a list of pairs of factor names, e.g.
#'   `list(c("condition", "sex"))`.
#' @param strata Columns whose values split the samples into strata, each designed on its own (a
#'   character vector).
#' @param batch The column holding each sample's batch, fitted as a fixed block. `NULL`: no batch.
#' @param out A directory to write `design_<i>.tsv` (for [stats_fit()]'s `design`) and
#'   `contrasts_<i>.tsv` (for its `contrasts`) into, per stratum. `NULL`: nothing is written.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_diff_design`:
#'
#'   - `records`, a data.frame with one row per (stratum, default contrast), strata in order:
#'     `stratum`, the stratum key; `contrast_id`, `c1`, `c2`, ... numbered **within the stratum**;
#'     `label`, `factor=level vs factor=reference`, or a slope's column; `numerator` and
#'     `denominator`, `factor=level` (`NA` for a slope); `covariate`, `covariate_unit` (the unit a
#'     slope is per, e.g. `decade`) and `covariate_scale` (e.g. `years/10, centred at 50 years`),
#'     `NA` when not a slope; `estimability`, `"estimable"`, `"not_estimable:confounded"` or
#'     `"not_estimable:rank_deficient"`; and two list columns, `weights` (contrast weight per design
#'     column, unitless), one per entry of the stratum's `column_names`, and `fitted_weights`
#'     (contrast weight per fitted column, unitless), one per entry of its `fitted_column_names` -
#'     the contrasts file's row - and `NULL` when the contrast is not estimable.
#'   - `strata`, a list of `mzlibr_design_stratum`, named by stratum key (`"all"` for one stratum,
#'     else the stratum columns sorted by name as `name=value` joined by `;`), each with:
#'     `sample_count` samples; `sample_ids`, the matrix rows, ordinal by sample id; `column_names`,
#'     mzLib's columns (`intercept`, `factor=level` per non-reference level, `batch=level`,
#'     `name (per unit)` for a covariate, `a=level:b=level` for an interaction); `fitted_column_names`,
#'     the columns a fit takes (the design file's coefficient columns); `rank` (columns);
#'     `aliased_terms`, terms that share a direction with another; `warnings`; `not_run`, a
#'     data.frame of the contrasts this stratum cannot run (`label`, `reason`); `covariate_centres`
#'     (the covariate's stored unit), the centre subtracted from each covariate; `matrix` and
#'     `fitted_matrix`, R matrices with the sample ids as row names, in coding (0/1 indicators;
#'     covariate columns in scaled units), equal to R's `model.matrix()`; and `design_file` and
#'     `contrasts_file`, the files written for this stratum, `NA` without `out` (and
#'     `contrasts_file` `NA` when nothing here is estimable).
#'   - `samples_file`, as an absolute path; `sample_count` samples read, across every stratum;
#'     `factors`, `covariates` and `interactions` as given (a covariate field `NA` was not given, so
#'     mzLib's default applied); `stratum_columns`; `batch_column` (`NA` when batch is not in the
#'     model); `unused_columns`; `stratum_count` strata; `estimability_counts`, default contrasts
#'     (contrasts) per verdict across strata; `written`, `directory` and `file_count` (files), or
#'     `NULL` without `out`; `row_count` rows; `column_names`.
#'   - `caveats`: how the matrix is coded and what the contrast ids mean, plus counts of contrasts
#'     that cannot be estimated or run.
#'
#' @section When the samples cannot answer:
#'
#' A comparison the samples cannot separate - a condition run in only one batch - is
#' `"not_estimable:confounded"`, and one that needs a column no sample informs (an empty interaction
#' cell) is `"not_estimable:rank_deficient"`. mzLib says so instead of returning a number. The design
#' file still holds a design that can be fitted: `fitted_column_names` leaves out the later of two
#' aliased columns. Contrast ids are numbered within each stratum, so match strata by `label`.
#'
#' @seealso [stats_fit()], which reads the files `out` writes; [diff_observations()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # The MALAT1 dilution series stats_fit() analyses, described by each run's load rather than by a
#' # hand-made 0/1 design:
#' samples <- system.file("extdata", "diff", "malat_dilution_samples.tsv", package = "mzLibR")
#' design <- diff_design(samples, factors = list(factor_spec("malat_load", "500ng")),
#'   out = file.path(tempdir(), "malat_design"))
#' design
#' stratum <- design$strata[["all"]]
#' stratum$fitted_column_names
#' design$records[, c("contrast_id", "label", "estimability")]
#' design$records$fitted_weights
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec diff.design
#' @export
diff_design <- function(samples, factors = NULL, covariates = NULL, interactions = NULL, strata = NULL,
                        batch = NULL, out = NULL, timeout = 60) {
  request <- diff_build_design_request(samples, factors, covariates, interactions, strata, batch, out)
  bridge_require_verb("diff design", DIFF_SINCE)
  diff_parse_design(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

#' Every peptide in every biological sample, as log2 intensity or NA with the reason, by mzLib
#'
#' mzLib's `FlashLfqObservations.FromPeptideTable` reads a stored FlashLFQ or MetaMorpheus peptide
#' table, maps each run's detection type to a state, sums each technical replicate's fractions on the
#' linear scale and averages a sample's replicates in log2. Then `SharedPeptideNormalization.Measure`
#' measures how far each sample sits from the others on the peptides they share; `normalize = TRUE`
#' subtracts it.
#'
#' @param peptides A FlashLFQ `QuantifiedPeptides.tsv` or a MetaMorpheus `AllQuantifiedPeptides.tsv`.
#' @param runs Every run of the table, as a data.frame with `file` and `sample` columns, and
#'   optionally `fraction` and `technical_replicate` (whole numbers; `NA` or absent sends none, so
#'   mzLib's 0). `file` is the run's name as the `Intensity_` columns write it, without extension; a
#'   column with `-calib` after it matches too. The result's own `runs` is accepted as well.
#' @param basis `"msms_only"` (no match-between-runs value enters) or `"mbr_kept"`. Required: mzLib
#'   has no default.
#' @param protein_groups A MetaMorpheus `AllQuantifiedProteinGroups.tsv`, for each group's q-value
#'   and decoy, contaminant and entrapment flags. `NULL`: none described.
#' @param samples A TSV whose first column (after a header row) lists the samples to keep, e.g. one
#'   stratum's design file from [diff_design()]. `NULL`: every sample.
#' @param normalize Subtract each sample's shift. The shifts are measured and reported either way.
#' @param min_peptides The floor on the normalization's reference set (peptides). `NULL`: not sent,
#'   so mzLib's own floor, 100 at mzLib `d10017ce`.
#' @param out A `.tsv` path to write the long table to instead of returning its rows.
#' @param responses_out A `.tsv` path to also write the peptide-by-sample log2 table that
#'   [stats_fit()] reads as `responses`.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_observations`:
#'
#'   - `records`, a data.frame with one row per (peptide, sample), peptides ordinal by full
#'     sequence: `full_sequence`, the engine's modified sequence, the key; `base_sequence`;
#'     `protein_groups`, a list column of every group the peptide table names (`UNDEFINED` = in no
#'     surviving group); `unique_protein_group`, the one group it is unique to, `NA` when shared,
#'     `UNDEFINED` or none; `sample`; `log2_intensity` (log2 of FlashLFQ intensity units), fractions
#'     summed, replicates averaged in log2, less the shift with `normalize = TRUE`, and `NA` for
#'     **no value** (see `state`), never 0; and `state`, one of `"quantified"`, `"mbr_transferred"`,
#'     `"ambiguous_peak"`, `"identified_not_quantified"` or `"not_detected"`. Empty with `out`.
#'   - `peptides_file`, `protein_groups_file` and `samples_file` as absolute paths, the last two
#'     `NA` when not given; `basis` and `normalize` as given; `min_peptides` (peptides) as given, `NA`
#'     when mzLib's floor applied; `run_count` runs; `runs`, a data.frame of the runs as mzLib holds
#'     them (`file_name`, `sample_id`, `fraction`, `technical_replicate`), ordinal by file name;
#'     `sample_count` samples; `samples`; `samples_without_values`; `peptide_count` peptides;
#'     `protein_group_count` protein groups; `protein_groups`, a data.frame of those groups (`name`,
#'     `q_value` a fraction 0 to 1 or `NA` if not finite, `is_decoy`, `is_contaminant`,
#'     `is_entrapment`, `genes` or `NA`); `state_counts`, rows (cells) per state; `warnings`;
#'     `written` (`path`, `format`, `row_count`) or `NULL`; `responses_written` (`path`, `format`,
#'     `row_count` peptides, `sample_count` samples) or `NULL`; `row_count` rows, `peptide_count` x
#'     `sample_count`; `column_names`.
#'   - `normalization`, mzLib's `NormalizationResult`, measured whether or not it was applied:
#'     `setting`, `"shared_peptide_median"` (the peptides with a value in every sample),
#'     `"shared_peptide_median_half"` (in at least half; used when too few are in all) or `"none"`
#'     (too few either way: no sample shifted); `reference_set_size` (peptides) in the set used;
#'     `minimum_reference_set_size` (peptides), the floor used; `shifts`, each sample's shift in log2
#'     units, `NA` for a sample with no shift; and `warnings`.
#'   - `caveats`: how values are combined, what the basis does, what mzLib cannot check about the
#'     table, and what the normalization assumes.
#'
#' @section The traps:
#'
#' Median normalization assumes most peptides do not change between samples. When most of them do -
#' a table of the spiked organism alone, a loading series, a pull-down - normalizing deletes the
#' result: read `normalization$shifts` before you subtract them. The peptide table does not record
#' the design FlashLFQ quantified with, so give the runs that design. mzLib does not yet roll
#' peptides up to proteins or choose a test's rows, so neither does this.
#'
#' @seealso [diff_design()], [stats_fit()], [flashlfq_quantify()]
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' # Two runs of FlashLFQ's E. coli spike-in vignette, which their names say were spiked at 1x and
#' # at 2x:
#' peptides <- system.file("extdata", "diff", "EcoliSpikeIn_QuantifiedPeptides.tsv",
#'   package = "mzLibR")
#' runs <- data.frame(
#'   file = c("09-04-18_EcoliSpikeInSingleShot1x", "09-04-18_EcoliSpikeInSingleShot2x"),
#'   sample = c("spike_1x", "spike_2x"))
#' table <- diff_observations(peptides, runs, basis = "msms_only")
#' table
#' table$state_counts
#' table$normalization[c("setting", "reference_set_size")]
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec diff.observations
#' @export
diff_observations <- function(peptides, runs, basis, protein_groups = NULL, samples = NULL, normalize = FALSE,
                              min_peptides = NULL, out = NULL, responses_out = NULL, timeout = 600) {
  if (missing(basis)) {
    stop(mzlib_usage_error(paste0(
      "basis is required: one of ", paste(DIFF_QUANT_BASES, collapse = ", "), ". mzLib has no default: ",
      "whether match-between-runs values enter is your choice."
    )))
  }
  request <- diff_build_observations_request(peptides, runs, basis, protein_groups, samples, normalize,
    min_peptides, out, responses_out
  )
  bridge_require_verb("diff observations", DIFF_SINCE)
  diff_parse_observations(bridge_invoke(request$args, stdin = request$stdin, timeout = timeout))
}

# ---------------------------------------------------------------- printing

#' Print a factor of a model
#'
#' @param x A [factor_spec()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_factor_spec <- function(x, ...) {
  cat("<mzlibr_factor_spec> ", x$name,
    if (is.na(x$reference_level)) " (no reference: every pair of levels)" else paste0(" against ", x$reference_level),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' Print a covariate of a model
#'
#' @param x A [covariate_spec()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_covariate_spec <- function(x, ...) {
  scale <- if (is.na(x$scale)) "" else paste0(", scale ", format(x$scale))
  per <- if (is.na(x$scaled_unit)) "" else paste0(" per ", x$scaled_unit)
  centre <- if (is.na(x$centre)) ", centred at the stratum's mean" else paste0(", centred at ", format(x$centre))
  cat("<mzlibr_covariate_spec> ", x$name, " in ", x$unit, scale, per, centre, "\n", sep = "")
  invisible(x)
}

#' Print a design
#'
#' @param x A [diff_design()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_diff_design <- function(x, ...) {
  cat("<mzlibr_diff_design> ", format(x$sample_count), " samples in ", format(x$stratum_count),
    " stratum(s), ", format(x$row_count), " default contrast(s)\n",
    sep = ""
  )
  if (length(x$estimability_counts) > 0L) {
    cat("  ", paste0(names(x$estimability_counts), " ", x$estimability_counts, collapse = ", "), "\n", sep = "")
  }
  for (stratum in x$strata) {
    cat("  ", stratum$stratum, ": ", format(stratum$sample_count), " samples, columns ",
      paste(stratum$fitted_column_names, collapse = ", "), "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Print one stratum's design
#'
#' @param x One of the `strata` of a [diff_design()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_design_stratum <- function(x, ...) {
  cat("<mzlibr_design_stratum> ", x$stratum, ": ", format(x$sample_count), " samples, rank ",
    format(x$rank), ", columns ", paste(x$fitted_column_names, collapse = ", "), "\n",
    sep = ""
  )
  for (warning in x$warnings) {
    cat("  ! ", warning, "\n", sep = "")
  }
  invisible(x)
}

#' Print an observation table
#'
#' @param x A [diff_observations()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_observations <- function(x, ...) {
  cat("<mzlibr_observations> ", format(x$peptide_count), " peptides x ", format(x$sample_count),
    " samples, ", x$basis, if (x$normalize) ", normalized" else "", "\n",
    sep = ""
  )
  if (length(x$state_counts) > 0L) {
    cat("  ", paste0(names(x$state_counts), " ", x$state_counts, collapse = ", "), "\n", sep = "")
  }
  cat("  normalization: ", x$normalization$setting, " on ", format(x$normalization$reference_set_size),
    " peptides\n",
    sep = ""
  )
  invisible(x)
}
