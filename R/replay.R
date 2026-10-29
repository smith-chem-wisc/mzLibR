# A stand-in bridge that answers from recorded responses instead of running mzLib.
#
# Every help-page example in this package runs under `R CMD check`, on machines with no bridge, no
# .NET and no network. It can, because the example first points mzLibR at this stand-in (in a
# `\dontshow{}` block, so the reader sees only the call they would make): an executable that is
# invoked exactly as the real bridge is, with the same arguments, and answers with the same
# envelope - read from a fixture recorded from the real bridge. The package cannot tell the
# difference, so an example exercises the real argument handling, transport and parsing, and
# prints what mzLib really returned.
#
# It is a port of pyMzLib's `tests/replay_bridge.py`, and keeps its one rule: **a recording
# answers only a call it fits.** An example cannot print output recorded for different arguments:
#
#   * `--path` must name the file the recording was made from (compared by file name), and any
#     other `--option value` whose snake_case name is a top-level scalar of the recording must equal
#     it;
#   * `--limit`/`--offset` must reproduce the recording's `returned_count` from its `record_count`
#     (or `row_count`), and a `--flag` with a `<flag>_included` field must match it;
#   * a `--paths-stdin` call fits only a bulk recording (per-input `files[]` entries carrying `path`,
#     and no top-level `path`), and a one-path call only a one-document recording;
#   * a recording with a `written` key fits an `--out` call only if `written` is set, and a call
#     without `--out` only if it is null;
#   * a recording listed in `REPLAY_RECORDED_WITH` fits only a call that sends the options it was
#     made with;
#   * `--link` (stats fit-latent-mixed) must give the numbers of the recording's `link` echo, and a
#     `true`/`false` option (stats fit-link --with-depth) must equal its boolean echo.
#
# No fit, or more than one, is answered as a usage error naming every recording and why it did not
# fit - so the example fails and says why, rather than printing something plausible.
#
# The recordings live in inst/replay/, staged there from tests/fixtures/ by scripts/stage-replay.R,
# with inst/replay/TABLE mapping each wire verb to the recordings that may answer it (the verb
# specs' `examples`, plus a short list for verbs whose spec names none).

REPLAY_TABLE_FILE <- "TABLE"

replay_dir <- function() {
  system.file("replay", package = "mzLibR")
}

# The verb -> recordings map, as a named list of character vectors.
replay_read_table <- function(dir) {
  path <- file.path(dir, REPLAY_TABLE_FILE)
  if (!file.exists(path)) {
    return(list())
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  verbs <- vapply(parts, function(p) p[1L], character(1L))
  files <- vapply(parts, function(p) p[2L], character(1L))
  split(files, factor(verbs, levels = unique(verbs)))
}

# The verb words, then `--name value` pairs, with a bare `--flag` recorded as TRUE.
replay_parse_argv <- function(argv) {
  i <- 1L
  verb <- character(0)
  while (i <= length(argv) && !startsWith(argv[i], "--")) {
    verb <- c(verb, argv[i])
    i <- i + 1L
  }
  options <- list()
  while (i <= length(argv)) {
    name <- substring(argv[i], 3L)
    if (i + 1L <= length(argv) && !startsWith(argv[i + 1L], "--")) {
      options[[name]] <- argv[i + 1L]
      i <- i + 2L
    } else {
      options[[name]] <- TRUE
      i <- i + 1L
    }
  }
  list(verb = paste(verb, collapse = " "), options = options)
}

# A file name from a path written on either kind of machine: the recordings were made on Windows.
replay_base <- function(path) {
  parts <- strsplit(as.character(path), "[/\\\\]")[[1L]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) "" else parts[length(parts)]
}

replay_scalar <- function(value) {
  !is.list(value) && length(value) == 1L && !is.na(value)
}

replay_show <- function(value) {
  if (is.logical(value)) tolower(as.character(value)) else format(value)
}

# Wire options a verb echoes under another name, so a recording made with one value cannot answer
# a call with another (peptidoform fragments echoes --max-mods as max_modifications).
REPLAY_ECHOED_AS <- c(
  "max-mods" = "max_modifications",
  "max-isoforms" = "max_modification_isoforms",
  "psms" = "psm_file",
  "peptides" = "peptides_file",
  "from" = "source_format",
  "to" = "target_format"
)

# Recordings made with an option the call must send as well. The echo alone cannot say so: a call
# that leaves the option out takes a default (mzLib's) that the replay does not know, and the echo
# of a value that was given looks the same as the echo of a default. So a call without --max-mods
# cannot be answered from the max-mods 1 albumin digest.
REPLAY_RECORDED_WITH <- list(
  "peptidoform_P02768_small.json" = c("max-mods" = "1")
)

# proteins align, library read and sdrf isobaric-design (pyMzLib 0.6.0): each echoes the option,
# but a call that leaves it out takes mzLib's default, which the echo cannot tell apart.
REPLAY_RECORDED_WITH[["proteins_align_charged.json"]] <- c("end-gaps" = "charged")
REPLAY_RECORDED_WITH[["library_read_fragments.json"]] <- c("fragments" = "TRUE")
REPLAY_RECORDED_WITH[["library_read_index_only.json"]] <- c("index-only-threshold" = "0")
REPLAY_RECORDED_WITH[["sdrf_isobaric_design_PXD008841.json"]] <- c("plex-file-name-pattern" = "TMT_?pool(?<plex>\\d+)")

# Why a recording made with `recorded_with` cannot answer a call with these options, or "".
replay_recorded_with_mismatch <- function(recorded_with, options) {
  for (name in names(recorded_with)) {
    sent <- options[[name]]
    if (is.null(sent) || !identical(as.character(sent), recorded_with[[name]])) {
      return(paste0("recorded with --", name, " ", recorded_with[[name]], ", which the call does not send"))
    }
  }
  ""
}

# Verbs whose input travels on stdin and is echoed in the recording, so a recording answers only a
# call that sent the same input. Read only for these: other verbs may be given a terminal's stdin.
REPLAY_STDIN_ECHO <- c("quant flashlfq", "peptidoform convert")

# The data verbs of pyMzLib 0.6.0, each with its own check of the call's stdin against the
# recording: the same accessions, the same pairs, the same searched files, the same species.
REPLAY_STDIN_CHECKS <- list(
  # The same lines, in the same order, as the recording's `input` column.
  "proteins parse-accession" = function(data, lines) {
    columns <- data[["columns"]]
    recorded <- if (is.list(columns)) wire_strings(columns[["input"]]) else character(0)
    if (!identical(lines, recorded)) {
      return(paste0("recorded for inputs ", paste(recorded, collapse = ", "), ", not ", paste(lines, collapse = ", ")))
    }
    ""
  },
  # The same pairs, in order, as the recording's aligned rows with their gaps removed (mzLib
  # refuses a '-' in an input sequence, so removing them is exact).
  "proteins align" = function(data, lines) {
    columns <- data[["columns"]]
    recorded <- if (is.list(columns)) {
      paste0(gsub("-", "", wire_strings(columns[["aligned_source"]]), fixed = TRUE), "\t",
        gsub("-", "", wire_strings(columns[["aligned_target"]]), fixed = TRUE))
    } else {
      character(0)
    }
    if (!identical(lines, recorded)) {
      return(paste0("recorded for pairs starting ", paste(substr(recorded, 1L, 12L), collapse = ", "), ", not these"))
    }
    ""
  },
  # The same searched files, in order, as `searched_files`, and as many acquired names as
  # `acquired_name_count`.
  "sdrf restrict" = function(data, lines) {
    fields <- strsplit(lines, "\t", fixed = TRUE)
    sent <- vapply(fields, function(f) trimws(f[[1L]]), character(1L))
    recorded <- wire_strings(data[["searched_files"]])
    if (!identical(sent, recorded)) {
      return(paste0("recorded for searched files ", paste(recorded, collapse = ", "), ", not ", paste(sent, collapse = ", ")))
    }
    acquired <- sum(vapply(fields, function(f) length(f) > 1L && nzchar(trimws(f[[2L]])), logical(1L)))
    if (!identical(as.numeric(acquired), as.numeric(data[["acquired_name_count"]]))) {
      return(paste0("recorded with ", format(data[["acquired_name_count"]]), " acquired names, not ", acquired))
    }
    ""
  },
  # The same species (the first field of each line before `--`) as `species`, and as many dump
  # paths after it as the recording has `dumps`.
  "genes orthology" = function(data, lines) {
    lines <- trimws(lines)
    cut <- match("--", lines)
    if (is.na(cut)) cut <- length(lines) + 1L
    species <- sort(vapply(strsplit(lines[seq_len(cut - 1L)], "\t", fixed = TRUE), function(f) trimws(f[[1L]]),
      character(1L)), method = "radix")
    recorded <- wire_strings(data[["species"]])
    if (!identical(species, recorded)) {
      return(paste0("recorded for species ", paste(recorded, collapse = ", "), ", not ", paste(species, collapse = ", ")))
    }
    dumps <- length(lines) - cut
    if (dumps != length(data[["dumps"]])) {
      return(paste0("recorded with ", length(data[["dumps"]]), " dumps, not ", dumps))
    }
    ""
  }
)

# pyMzLib 0.7.0's two data verbs.
#
# genes map-residues: the same species as the recording's `gene_sets`, and the same queries after
# `--`, in the same order, as its source, position and target columns.
REPLAY_STDIN_CHECKS[["genes map-residues"]] <- function(data, lines) {
  lines <- trimws(lines)
  cut <- match("--", lines)
  if (is.na(cut)) cut <- length(lines) + 1L
  species <- sort(vapply(strsplit(lines[seq_len(cut - 1L)], "\t", fixed = TRUE), function(f) trimws(f[[1L]]),
    character(1L)), method = "radix")
  recorded_species <- sort(names(data[["gene_sets"]]), method = "radix")
  if (is.null(recorded_species)) recorded_species <- character(0)
  if (!identical(species, recorded_species)) {
    return(paste0("recorded for species ", paste(recorded_species, collapse = ", "), ", not ",
      paste(species, collapse = ", ")))
  }
  sent <- if (cut < length(lines)) lines[(cut + 1L):length(lines)] else character(0)
  sent <- vapply(strsplit(sent, "\t", fixed = TRUE), function(f) paste(trimws(f), collapse = "\t"), character(1L))
  columns <- data[["columns"]]
  recorded <- if (is.list(columns)) {
    paste(wire_strings(columns[["source_protein_id"]]),
      vapply(columns[["position"]], function(p) sprintf("%.0f", as.numeric(p)), character(1L)),
      wire_strings(columns[["target_protein_id"]]), sep = "\t")
  } else {
    character(0)
  }
  if (!identical(unname(sent), recorded)) {
    return(paste0("recorded for queries ", paste(utils::head(recorded, 2L), collapse = "; "), "..., not these"))
  }
  ""
}

# library calibrate-rt: the same lines, in the same order, as the recording's rows: each line's
# group, role (anchor before `--`; to_irt or to_minutes after it, by which time is given) and the
# times it gives, compared as numbers.
REPLAY_STDIN_CHECKS[["library calibrate-rt"]] <- function(data, lines) {
  number <- function(text) if (nzchar(trimws(text))) suppressWarnings(as.numeric(text)) else NA_real_
  value <- function(v) if (wire_null(v)) NA_real_ else as.numeric(v)
  queries <- FALSE
  sent <- list()
  for (line in lines) {
    if (identical(trimws(line), "--")) {
      queries <- TRUE
      next
    }
    f <- c(strsplit(line, "\t", fixed = TRUE)[[1L]], "", "")
    rt <- number(f[[2L]])
    irt <- number(f[[3L]])
    role <- if (!queries) "anchor" else if (!is.na(rt)) "to_irt" else "to_minutes"
    sent[[length(sent) + 1L]] <- list(trimws(f[[1L]]), role, if (role != "to_minutes") rt else NA_real_,
      if (role != "to_irt") irt else NA_real_)
  }
  columns <- data[["columns"]]
  groups <- if (is.list(columns)) wire_strings(columns[["group"]]) else character(0)
  roles <- if (is.list(columns)) wire_strings(columns[["role"]]) else character(0)
  recorded <- lapply(seq_along(groups), function(i) {
    list(groups[[i]], roles[[i]],
      if (roles[[i]] != "to_minutes") value(columns[["rt_minutes"]][[i]]) else NA_real_,
      if (roles[[i]] != "to_irt") value(columns[["irt"]][[i]]) else NA_real_)
  })
  if (!isTRUE(all.equal(sent, recorded, tolerance = 1e-12))) {
    return(paste0("recorded for ", length(recorded), " other lines"))
  }
  ""
}
REPLAY_STDIN_ECHO <- c(REPLAY_STDIN_ECHO, names(REPLAY_STDIN_CHECKS))

# The call's stdin must be the recording's input. quant flashlfq: the same mzML runs, by file
# name. peptidoform convert: the same sequences, in the same order, as the recording's `input`
# column.
replay_stdin_mismatch <- function(data, stdin, verb = "quant flashlfq") {
  if (verb %in% REPLAY_MODEL_STDIN) {
    return(replay_model_stdin_mismatch(data, stdin, verb))
  }
  lines <- stdin[nzchar(trimws(stdin))]
  if (verb %in% names(REPLAY_STDIN_CHECKS)) {
    return(REPLAY_STDIN_CHECKS[[verb]](data, lines))
  }
  if (identical(verb, "peptidoform convert")) {
    columns <- data[["columns"]]
    recorded <- if (is.list(columns)) wire_strings(columns[["input"]]) else character(0)
    if (!identical(lines, recorded)) {
      return(paste0("recorded for sequences ", paste(recorded, collapse = ", "), ", not ",
        paste(lines, collapse = ", ")))
    }
    return("")
  }
  sent <- sort(unique(vapply(lines, function(l) replay_base(trimws(strsplit(l, "\t", fixed = TRUE)[[1L]][1L])),
    character(1L), USE.NAMES = FALSE)))
  files <- data[["spectra_files"]]
  recorded <- if (is.list(files)) {
    sort(unique(vapply(files, function(f) replay_base(as.character(f[["full_path"]])), character(1L))))
  } else {
    character(0)
  }
  if (!identical(sent, recorded)) {
    return(paste0("recorded for runs ", paste(recorded, collapse = ", "), ", not ", paste(sent, collapse = ", ")))
  }
  ""
}

# ---------------------------------------------------------------- stats, diff and ptmqtl

# The rules pyMzLib 0.6.0's replay_bridge.py added for these verbs, ported as they are there.

# Recordings made with an option a call must send too: a call that leaves it out takes mzLib's
# default, which the echo cannot tell apart from a value that was given.
REPLAY_RECORDED_WITH[c(
  "stats_fit_mixed_lmm_ml.json", "stats_fit_malat_marginal.json",
  "diff_observations_spike_in_normalized.json", "ptmqtl_pairs_batch_run_level.json",
  "ptmqtl_trait_effects_strict.json"
)] <- list(
  c(estimator = "ml"),
  c(estimator = "marginal_likelihood"),
  c("min-peptides" = "50", normalize = "TRUE"),
  c("remove-run-level" = "TRUE"),
  c("min-replicates" = "6")
)

# pyMzLib 0.7.0's recordings, by the same rule: ptmqtl's new options; library calibrate-rt, which
# echoes the options it applied (mzLib's defaults when absent); stats fit-latent-mixed's ML fit; and
# stats fit-link, which echoes the value it passed, mzLib's declared default when absent.
REPLAY_RECORDED_WITH[c(
  "ptmqtl_pairs_bounds.json", "ptmqtl_pairs_depth_covariate_run_level.json", "ptmqtl_trait_effects_log2.json",
  "library_calibrate_rt_pbmc_skip.json", "library_calibrate_rt_pbmc_linear.json",
  "stats_fit_latent_mixed_link_ml.json", "stats_fit_link_simulated.json"
)] <- list(
  c("include-bounds" = "TRUE"),
  c("remove-run-level" = "TRUE"),
  c(scale = "log2_odds"),
  c("on-error" = "skip"),
  c(kind = "linear", "minimum-anchors" = "10"),
  c(estimator = "ml"),
  c("with-depth" = "true")
)

# Their input travels on stdin and is echoed: the coefficients, the model, the runs, the tests.
REPLAY_MODEL_STDIN <- c(
  "stats combine-p", "stats fit", "stats fit-logistic", "stats fit-mixed", "stats fit-latent-mixed",
  "diff design", "diff observations"
)
REPLAY_STDIN_ECHO <- c(REPLAY_STDIN_ECHO, REPLAY_MODEL_STDIN)

# Optional input files a verb echoes as `<name>_file`, null when the option was not given. A
# recording made with one answers only a call that gives one, and the reverse, or a plain
# correlation could print a partial one.
REPLAY_OPTIONAL_FILES <- c(
  "given", "protein-groups", "samples", "run-level-groups", "run-covariate", "contrasts",
  "responses", "counts", "depth", "residual-variance"
)

# Options echoed inside an object rather than at the top level (stats fit --estimator as
# prior.estimator, --confidence as contrasts.confidence_level).
REPLAY_NESTED_ECHOES <- list(estimator = c("prior", "estimator"), confidence = c("contrasts", "confidence_level"))

# Optional outputs reported in their own field, null when not asked for (as `written` is for --out).
REPLAY_OPTIONAL_OUTPUTS <- c("responses-out" = "responses_written")

# Echoes that say whether an option was given when its value echo cannot (stats entrapment-fdp
# echoes mzLib's default thresholds when --thresholds was absent).
REPLAY_GIVEN_ECHOES <- c(thresholds_given = "thresholds")

# Why the presence or absence of an option rules this recording out, or "".
replay_presence_mismatch <- function(data, options) {
  for (name in names(REPLAY_NESTED_ECHOES)) {
    where <- REPLAY_NESTED_ECHOES[[name]]
    echo <- data[[where[1L]]]
    sent <- options[[name]]
    if (!is.null(sent) && is.list(echo) && where[2L] %in% names(echo) && replay_scalar(echo[[where[2L]]])) {
      recorded <- echo[[where[2L]]]
      same <- if (is.numeric(recorded)) {
        isTRUE(all.equal(suppressWarnings(as.numeric(sent)), as.numeric(recorded), tolerance = 1e-12))
      } else {
        identical(as.character(recorded), as.character(sent))
      }
      if (!same) {
        return(paste0("recorded with ", where[1L], ".", where[2L], "=", replay_show(recorded), ", not ", sent))
      }
    }
  }
  for (key in names(REPLAY_GIVEN_ECHOES)) {
    given <- data[[key]]
    if (is.logical(given) && length(given) == 1L && !is.na(given) &&
      !is.null(options[[REPLAY_GIVEN_ECHOES[[key]]]]) != given) {
      return(paste0("recorded ", if (given) "with" else "without", " --", REPLAY_GIVEN_ECHOES[[key]]))
    }
  }
  for (name in REPLAY_OPTIONAL_FILES) {
    key <- paste0(gsub("-", "_", name, fixed = TRUE), "_file")
    if (key %in% names(data)) {
      used <- !wire_null(data[[key]])
      if (!is.null(options[[name]]) != used) {
        return(paste0("recorded ", if (used) "with" else "without", " --", name))
      }
    }
  }
  for (name in names(REPLAY_OPTIONAL_OUTPUTS)) {
    key <- REPLAY_OPTIONAL_OUTPUTS[[name]]
    if (key %in% names(data)) {
      used <- !wire_null(data[[key]])
      if (!is.null(options[[name]]) != used) {
        return(paste0("recorded ", if (used) "with" else "without", " --", name))
      }
    }
  }
  ""
}

# The model a diff design call sent, in the shape the recording echoes it.
replay_diff_model <- function(stdin) {
  model <- list(factors = list(), covariates = list(), interactions = list(), stratum_columns = character(0),
    batch_column = NA_character_)
  number <- function(text) if (nzchar(trimws(text))) as.numeric(text) else NA_real_
  for (line in stdin[nzchar(trimws(stdin))]) {
    f <- trimws(strsplit(line, "\t", fixed = TRUE)[[1L]])
    f <- c(f, rep("", 6L - length(f)))
    if (f[1L] == "factor") {
      model$factors[[length(model$factors) + 1L]] <- c(f[2L], if (nzchar(f[3L])) f[3L] else NA_character_)
    } else if (f[1L] == "covariate") {
      model$covariates[[length(model$covariates) + 1L]] <- list(f[2L], f[3L], number(f[4L]),
        if (nzchar(f[5L])) f[5L] else NA_character_, number(f[6L]))
    } else if (f[1L] == "interaction") {
      model$interactions[[length(model$interactions) + 1L]] <- c(f[2L], f[3L])
    } else if (f[1L] == "stratum") {
      model$stratum_columns <- c(model$stratum_columns, f[2L])
    } else if (f[1L] == "batch") {
      model$batch_column <- f[2L]
    }
  }
  model
}

replay_recorded_diff_model <- function(data) {
  text <- function(v) if (wire_null(v)) NA_character_ else as.character(v)
  value <- function(v) if (wire_null(v)) NA_real_ else as.numeric(v)
  list(
    factors = lapply(if (is.list(data$factors)) data$factors else list(), function(x) c(text(x$name), text(x$reference_level))),
    covariates = lapply(if (is.list(data$covariates)) data$covariates else list(), function(x) {
      list(text(x$name), text(x$unit), value(x$scale), text(x$scaled_unit), value(x$centre))
    }),
    interactions = lapply(if (is.list(data$interactions)) data$interactions else list(), function(x) c(text(x$first), text(x$second))),
    stratum_columns = wire_strings(data$stratum_columns),
    batch_column = text(data$batch_column)
  )
}

# The call's stdin must be the recording's input. stats fit, fit-logistic and fit-mixed: the same
# coefficients, in the same order, as `tested_coefficients` (none, for a fit of contrasts alone).
# diff design: the same model. diff observations: the same runs (with --samples, among those sent).
# stats combine-p: the same features in first-appearance order, the same number of tests, and
# weights and effects given or not as `weighted` and `signed` say.
replay_model_stdin_mismatch <- function(data, stdin, verb) {
  lines <- stdin[nzchar(trimws(stdin))]
  if (verb %in% c("stats fit", "stats fit-logistic", "stats fit-mixed", "stats fit-latent-mixed")) {
    sent <- trimws(lines)
    recorded <- wire_strings(data$tested_coefficients)
    if (!identical(sent, recorded)) {
      return(paste0("recorded for coefficients ", paste(recorded, collapse = ", "), ", not ", paste(sent, collapse = ", ")))
    }
    return("")
  }
  if (identical(verb, "diff design")) {
    sent <- replay_diff_model(lines)
    recorded <- replay_recorded_diff_model(data)
    if (!isTRUE(all.equal(sent, recorded, tolerance = 1e-12))) {
      return("recorded for another model")
    }
    return("")
  }
  if (identical(verb, "diff observations")) {
    key <- function(file, sample, fraction, replicate) {
      paste(file, sample, if (nzchar(fraction)) as.numeric(fraction) else 0, if (nzchar(replicate)) as.numeric(replicate) else 0, sep = "\r")
    }
    sent <- unique(vapply(lines, function(line) {
      f <- c(trimws(strsplit(line, "\t", fixed = TRUE)[[1L]]), "", "")
      key(f[1L], f[2L], f[3L], f[4L])
    }, character(1L), USE.NAMES = FALSE))
    recorded <- vapply(if (is.list(data$runs)) data$runs else list(), function(r) {
      key(as.character(r$file_name), as.character(r$sample_id), format(r$fraction), format(r$technical_replicate))
    }, character(1L))
    fits <- if (wire_null(data$samples_file)) setequal(recorded, sent) else all(recorded %in% sent)
    if (!fits) {
      return("recorded for other runs")
    }
    return("")
  }
  # stats combine-p
  fields <- lapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1L]])
  sent <- unique(vapply(fields, function(f) trimws(f[1L]), character(1L)))
  columns <- data$columns
  recorded <- if (is.list(columns)) wire_strings(columns$feature) else character(0)
  if (!identical(sent, recorded)) {
    return(paste0("recorded for features ", paste(recorded, collapse = ", "), ", not ", paste(sent, collapse = ", ")))
  }
  if (!identical(as.numeric(length(lines)), as.numeric(data$line_count))) {
    return(paste0("recorded for ", format(data$line_count), " tests, not ", length(lines)))
  }
  # strsplit drops a trailing empty field, so a line "A\t0.1\t" has length 2: a blank weight.
  weighted <- any(vapply(fields, function(f) length(f) >= 3L && !tolower(trimws(f[3L])) %in% c("", "na", "nan"), logical(1L)))
  signed <- any(vapply(lines, function(line) lengths(regmatches(line, gregexpr("\t", line, fixed = TRUE))) == 3L, logical(1L)))
  if (!identical(c(weighted, signed), c(isTRUE(data$weighted), isTRUE(data$signed)))) {
    return(paste0("recorded with weighted=", replay_show(isTRUE(data$weighted)), " signed=",
      replay_show(isTRUE(data$signed)), ", not weighted=", replay_show(weighted), " signed=", replay_show(signed)))
  }
  ""
}

# stats fit-latent-mixed --link against the recording's link echo: `identity` is (0, 1), and each
# number given must equal the echo's at its place (intercept, slope, depth_slope, overdispersion).
replay_link_mismatch <- function(echo, sent) {
  given <- if (identical(tolower(trimws(sent)), "identity")) {
    c(0, 1)
  } else {
    suppressWarnings(as.numeric(strsplit(sent, ",", fixed = TRUE)[[1L]]))
  }
  if (anyNA(given)) {
    return(paste0("--link ", sent, " is not numbers"))
  }
  names <- c("intercept", "slope", "depth_slope", "overdispersion")
  for (i in seq_along(given)) {
    recorded <- echo[[names[[i]]]]
    if (wire_null(recorded) || !isTRUE(all.equal(as.numeric(recorded), given[[i]], tolerance = 1e-12))) {
      return(paste0("recorded with another link, not --link ", sent))
    }
  }
  ""
}

# Why this recording cannot be the answer to these options, or "" when it can.
replay_mismatch <- function(data, options) {
  if (is.list(data[["link"]]) && is.character(options[["link"]])) {
    why <- replay_link_mismatch(data[["link"]], options[["link"]])
    if (nzchar(why)) {
      return(why)
    }
  }
  presence <- replay_presence_mismatch(data, options)
  if (nzchar(presence)) {
    return(presence)
  }
  null_echo <- ""
  for (name in names(options)) {
    value <- options[[name]]
    key <- if (name %in% names(REPLAY_ECHOED_AS)) REPLAY_ECHOED_AS[[name]] else gsub("-", "_", name, fixed = TRUE)
    if (name %in% c("limit", "offset", "out")) {
      next
    }
    if (isTRUE(value)) {
      flag <- data[[paste0(key, "_included")]]
      if (is.null(flag)) flag <- data[[key]]
      if (!is.null(flag) && replay_scalar(flag) && !isTRUE(flag)) {
        return(paste0("--", name, " given, but the recording has ", key, "=", replay_show(flag)))
      }
      next
    }
    # An input file option is echoed back as <name>_file (--responses -> responses_file): hold the
    # call's file to the recording's, by file name.
    echoed <- data[[paste0(key, "_file")]]
    if (!key %in% names(data) && is.character(echoed) && length(echoed) == 1L) {
      if (!identical(replay_base(echoed), replay_base(value))) {
        return(paste0("recorded from '", replay_base(echoed), "', not '", replay_base(value), "'"))
      }
      next
    }
    recorded <- data[[key]]
    # As in pyMzLib, a null echo of an option the call sends never fits: the recording was made
    # without it (sdrf isobaric-design echoes plex_file_name_pattern null when none was given).
    # Kept for last, so a rule below that names the difference better is the reason given.
    if (key %in% names(data) && !is.list(recorded) && length(recorded) == 1L && is.na(recorded)) {
      if (!nzchar(null_echo)) null_echo <- paste0("recorded with ", key, "=null, not ", value)
      next
    }
    if (is.null(recorded) || !replay_scalar(recorded)) {
      next
    }
    if (identical(key, "path") || endsWith(key, "_file")) {
      if (!identical(replay_base(recorded), replay_base(value))) {
        return(paste0("recorded from '", replay_base(recorded), "', not '", replay_base(value), "'"))
      }
    } else if (is.numeric(recorded) && !is.na(suppressWarnings(as.numeric(value)))) {
      # A number is compared as a number: 0.94999999999999996 on the wire is the recorded 0.95.
      if (!isTRUE(all.equal(as.numeric(value), as.numeric(recorded), tolerance = 1e-12))) {
        return(paste0("recorded with ", key, "=", replay_show(recorded), ", not ", value))
      }
    } else if (!identical(replay_show(recorded), value)) {
      return(paste0("recorded with ", key, "=", replay_show(recorded), ", not ", value))
    }
  }

  # BULK.md: a --paths-stdin call has its own envelope (files[], read_count, ...). Neither shape
  # answers for the other, or a one-path example would "fit" a bulk recording that has no path.
  # pride files also has a `files` list - of PRIDE files, which carry no `path` - so it is not bulk.
  # An empty list says nothing either way (PRIDE's answer for an unknown accession is files = []),
  # so it counts as bulk only when BULK.md's read_count says so.
  files <- data[["files"]]
  bulk_recording <- is.list(files) && (length(files) > 0L || "read_count" %in% names(data)) &&
    all(vapply(files, function(f) is.list(f) && !is.null(f[["path"]]), logical(1L))) &&
    is.null(data[["path"]])
  if (!is.null(options[["paths-stdin"]]) != bulk_recording) {
    return(if (bulk_recording) "a bulk (--paths-stdin) recording" else "a one-document recording")
  }

  # A recording that wrote a file answers only a call that asked for one, and the reverse: the
  # payload's `written` block is the evidence, so an `out =` example cannot print a recording that
  # wrote nothing.
  if ("written" %in% names(data) && !is.null(options[["out"]]) != is.list(data[["written"]])) {
    return(if (is.list(data[["written"]])) "a recording that wrote out=" else "a recording without out=")
  }

  # proteins read: a filtered read (--accessions-stdin) and an unfiltered one never stand in for
  # each other, or an example would print rows its filter did not select.
  if ("accession_filter_count" %in% names(data) &&
    !is.null(options[["accessions-stdin"]]) != replay_scalar(data[["accession_filter_count"]])) {
    return("a recording with the other accession filter (filtered vs unfiltered)")
  }

  # proteins read --sequences adds a sequence column: a recording with it and one without never
  # stand in for each other, or an example would print sequences it did not ask for (or lack them).
  if ("accession_filter_count" %in% names(data) &&
    !is.null(options[["sequences"]]) != ("sequence" %in% wire_strings(data[["column_names"]]))) {
    return("a recording with the other --sequences choice")
  }
  # proteins align: --no-residue-map leaves residue_map null, and only then.
  if ("residue_map" %in% names(data) && !is.null(options[["no-residue-map"]]) != !is.list(data[["residue_map"]])) {
    return("a recording with the other --no-residue-map choice")
  }
  # genes orthology: --pair-status FROM,TO fills pair_status for that direction, and only then; and
  # only a --restrict-gene-trees call says which gene sets its gene trees were restricted to.
  if ("pair_status" %in% names(data)) {
    status <- data[["pair_status"]]
    if (!is.null(options[["pair-status"]]) != is.list(status)) {
      return("a recording with the other --pair-status choice")
    }
    if (is.list(status)) {
      direction <- paste0(status[["from_species"]], ",", status[["to_species"]])
      if (!identical(options[["pair-status"]], direction)) {
        return(paste0("recorded for --pair-status ", direction))
      }
    }
  }
  trees <- data[["gene_trees"]]
  if (is.list(trees) && "restricted_to_gene_sets" %in% names(trees) && !"alignment" %in% names(data) &&
    !is.null(options[["restrict-gene-trees"]]) != is.list(trees[["restricted_to_gene_sets"]])) {
    return("a recording with the other --restrict-gene-trees choice")
  }
  # genes map-residues: --proteins fills proteins, and only then; and the alignment file is the
  # recording's (by file name), as --path is.
  if ("proteins_included" %in% names(data) && !is.null(options[["proteins"]]) != isTRUE(data[["proteins_included"]])) {
    return("a recording with the other --proteins choice")
  }
  alignment <- data[["alignment"]]
  if (is.list(alignment) && is.character(options[["alignment"]]) &&
    !identical(replay_base(options[["alignment"]]), as.character(alignment[["source_file_name"]]))) {
    return(paste0("recorded from '", alignment[["source_file_name"]], "', not '", replay_base(options[["alignment"]]), "'"))
  }

  # A filter the recording applied that the call did not ask for.
  ms_order <- data[["ms_order"]]
  if (!is.null(ms_order) && replay_scalar(ms_order) && is.null(options[["ms-order"]])) {
    return(paste0("recorded with ms_order=", replay_show(ms_order), ", but the call has no ms-order"))
  }
  # isobaric kits: one kit and every kit never stand in for each other, either way round.
  if ("kit" %in% names(data)) {
    kit <- data[["kit"]]
    if (replay_scalar(kit) && is.null(options[["kit"]])) {
      return(paste0("recorded for kit=", replay_show(kit), ", but the call asks for every kit"))
    }
    if (!replay_scalar(kit) && !is.null(options[["kit"]])) {
      return("recorded for every kit, but the call asks for one")
    }
  }

  total <- data[["record_count"]]
  if (is.null(total) || !replay_scalar(total)) total <- data[["row_count"]]
  returned <- data[["returned_count"]]
  if (!is.null(total) && replay_scalar(total) && !is.null(returned) && replay_scalar(returned) &&
    is.null(options[["out"]])) {
    offset <- if (is.null(options[["offset"]])) 0 else as.numeric(options[["offset"]])
    expected <- max(0, as.numeric(total) - offset)
    if (!is.null(options[["limit"]])) {
      expected <- min(expected, as.numeric(options[["limit"]]))
    }
    recorded_offset <- data[["offset"]]
    recorded_offset <- if (is.null(recorded_offset) || !replay_scalar(recorded_offset)) 0 else as.numeric(recorded_offset)
    if (recorded_offset != offset || as.numeric(returned) != expected) {
      return(paste0(
        "recorded window offset=", format(recorded_offset), " returned=", format(returned),
        " of ", format(total), ", not what limit/offset ask for"
      ))
    }
  }

  if (isTRUE(data[["peaks_included"]]) && is.null(options[["peaks"]])) {
    return("recorded with peaks, but the call did not ask for them")
  }
  null_echo
}

replay_json_string <- function(text) {
  text <- gsub("\\", "\\\\", text, fixed = TRUE)
  text <- gsub("\"", "\\\"", text, fixed = TRUE)
  text <- gsub("\n", "\\n", text, fixed = TRUE)
  text <- gsub("\r", "\\r", text, fixed = TRUE)
  text <- gsub("\t", "\\t", text, fixed = TRUE)
  paste0("\"", text, "\"")
}

replay_usage_error <- function(message) {
  paste0(
    "{\"ok\":false,\"data\":null,\"error\":{\"type\":\"usage\",\"message\":",
    replay_json_string(paste0("replay bridge: ", message)), "}}"
  )
}

replay_read_text <- function(path) {
  text <- rawToChar(readBin(path, what = "raw", n = file.info(path)$size))
  Encoding(text) <- "UTF-8"
  text
}

replay_describe_options <- function(options) {
  if (length(options) == 0L) {
    return("no options")
  }
  paste(vapply(names(options), function(n) {
    value <- options[[n]]
    if (isTRUE(value)) paste0("--", n) else paste0("--", n, " ", value)
  }, character(1L)), collapse = " ")
}

# The envelope, as JSON text, that the stand-in bridge prints for `argv`.
replay_answer <- function(argv, dir = replay_dir(), stdin = character(0)) {
  call <- replay_parse_argv(argv)
  candidates <- replay_read_table(dir)[[call$verb]]
  if (length(candidates) == 0L) {
    return(replay_usage_error(paste0("no recording is staged for '", call$verb, "'")))
  }

  fits <- character(0)
  fitting_text <- character(0)
  reasons <- character(0)
  for (fixture in candidates) {
    text <- replay_read_text(file.path(dir, fixture))
    parsed <- json_parse(text)
    enveloped <- is.list(parsed) && !is.null(parsed[["ok"]]) && "data" %in% names(parsed)
    data <- if (enveloped) parsed[["data"]] else parsed
    why <- if (is.list(data) && !is.null(names(data))) replay_mismatch(data, call$options) else ""
    if (!nzchar(why) && !is.null(REPLAY_RECORDED_WITH[[fixture]])) {
      why <- replay_recorded_with_mismatch(REPLAY_RECORDED_WITH[[fixture]], call$options)
    }
    if (!nzchar(why) && call$verb %in% REPLAY_STDIN_ECHO && is.list(data)) {
      why <- replay_stdin_mismatch(data, stdin, call$verb)
    }
    if (nzchar(why)) {
      reasons <- c(reasons, paste0(fixture, ": ", why))
    } else {
      fits <- c(fits, fixture)
      fitting_text <- c(fitting_text, if (enveloped) {
        text
      } else {
        paste0("{\"ok\":true,\"data\":", text, ",\"error\":null}")
      })
    }
  }

  if (length(fits) == 1L) {
    return(fitting_text[[1L]])
  }
  if (length(fits) == 0L) {
    return(replay_usage_error(paste0(
      "no recording of '", call$verb, "' fits ", replay_describe_options(call$options), ": ",
      paste(reasons, collapse = "; ")
    )))
  }
  replay_usage_error(paste0(
    "'", call$verb, "' ", replay_describe_options(call$options), " fits several recordings: ",
    paste(fits, collapse = ", ")
  ))
}

# The body of the stand-in executable: answer this process's arguments and exit as the bridge would.
replay_main <- function(dir = replay_dir()) {
  argv <- commandArgs(trailingOnly = TRUE)
  stdin <- if (replay_parse_argv(argv)$verb %in% REPLAY_STDIN_ECHO) {
    connection <- file("stdin")
    on.exit(close(connection), add = TRUE)
    readLines(connection, warn = FALSE)
  } else {
    character(0)
  }
  answer <- replay_answer(argv, dir, stdin)
  cat(answer, "\n", sep = "")
  quit(save = "no", status = if (startsWith(answer, "{\"ok\":false")) 2L else 0L)
}

# Write the stand-in executable for this platform into `where`, and return its path.
#
# A two-line launcher around Rscript: a batch file on Windows, a shell script elsewhere. It runs an
# R file that puts THIS installation of mzLibR first on the library path, so the stand-in answers
# with the same package the example is running - which matters under `R CMD check`, where the
# package under test is not installed in the usual library.
replay_write_launcher <- function(where = tempfile("mzlibr-replay-"), dir = replay_dir()) {
  dir.create(where, recursive = TRUE, showWarnings = FALSE)
  library_path <- dirname(system.file(package = "mzLibR"))
  main <- file.path(where, "replay-main.R")
  writeLines(c(
    paste0(".libPaths(c(", deparse(library_path), ", .libPaths()))"),
    paste0("mzLibR:::replay_main(", deparse(dir), ")")
  ), main)

  windows <- .Platform$OS.type == "windows"
  rscript <- file.path(R.home("bin"), if (windows) "Rscript.exe" else "Rscript")
  if (windows) {
    launcher <- file.path(where, "mzlib-bridge.cmd")
    writeLines(c(
      "@echo off",
      paste0("\"", normalizePath(rscript), "\" --vanilla \"", normalizePath(main), "\" %*"),
      "exit /b %ERRORLEVEL%"
    ), launcher)
  } else {
    launcher <- file.path(where, "mzlib-bridge")
    writeLines(c(
      "#!/bin/sh",
      paste0("exec ", shQuote(rscript), " --vanilla ", shQuote(main), " \"$@\"")
    ), launcher)
    Sys.chmod(launcher, "0755")
  }
  launcher
}

# Point mzLibR at the stand-in bridge; returns the previous `mzlibr.bridge` option for
# `replay_bridge_stop()`. Used by the examples, inside `\dontshow{}`.
replay_bridge_start <- function(dir = replay_dir()) {
  old <- options(mzlibr.bridge = replay_write_launcher(dir = dir))
  invisible(old)
}

replay_bridge_stop <- function(old) {
  launcher <- getOption("mzlibr.bridge")
  options(old)
  # Only ever the directory replay_write_launcher() made, never whatever the option pointed at.
  if (!is.null(launcher) && startsWith(basename(dirname(launcher)), "mzlibr-replay-")) {
    unlink(dirname(launcher), recursive = TRUE)
  }
  invisible(NULL)
}
