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
#     and no top-level `path`), and a one-path call only a one-document recording.
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

# Why this recording cannot be the answer to these options, or "" when it can.
replay_mismatch <- function(data, options) {
  for (name in names(options)) {
    value <- options[[name]]
    key <- gsub("-", "_", name, fixed = TRUE)
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
    recorded <- data[[key]]
    if (is.null(recorded) || !replay_scalar(recorded)) {
      next
    }
    if (identical(key, "path") || endsWith(key, "_file")) {
      if (!identical(replay_base(recorded), replay_base(value))) {
        return(paste0("recorded from '", replay_base(recorded), "', not '", replay_base(value), "'"))
      }
    } else if (!identical(replay_show(recorded), value)) {
      return(paste0("recorded with ", key, "=", replay_show(recorded), ", not ", value))
    }
  }

  # BULK.md: a --paths-stdin call has its own envelope (files[], read_count, ...). Neither shape
  # answers for the other, or a one-path example would "fit" a bulk recording that has no path.
  # pride files also has a `files` list - of PRIDE files, which carry no `path` - so it is not bulk.
  files <- data[["files"]]
  bulk_recording <- is.list(files) && length(files) > 0L &&
    all(vapply(files, function(f) is.list(f) && !is.null(f[["path"]]), logical(1L))) &&
    is.null(data[["path"]])
  if (!is.null(options[["paths-stdin"]]) != bulk_recording) {
    return(if (bulk_recording) "a bulk (--paths-stdin) recording" else "a one-document recording")
  }

  # A filter the recording applied that the call did not ask for.
  ms_order <- data[["ms_order"]]
  if (!is.null(ms_order) && replay_scalar(ms_order) && is.null(options[["ms-order"]])) {
    return(paste0("recorded with ms_order=", replay_show(ms_order), ", but the call has no ms-order"))
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
  ""
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
replay_answer <- function(argv, dir = replay_dir()) {
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
  answer <- replay_answer(commandArgs(trailingOnly = TRUE), dir)
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
