# A test harness, in base R, in about a hundred lines.
#
# `testthat` is an external package, and this one has no dependencies at all — not for users,
# not for developers. The consequence is better than the compromise: the suite runs on any R
# from 3.5 up, on a machine with nothing installed, with no toolchain and no network. That
# matters more here than in most packages, because the thing being tested is a binding whose
# whole promise is that it works where it lands.
#
# The only feature genuinely worth importing is `skip()`, and it is nine lines.

mz_results <- new.env(parent = emptyenv())
mz_results$pass <- 0L
mz_results$fail <- 0L
mz_results$skip <- 0L
mz_results$messages <- character(0)

# Reach into the package for the internal functions. Most of what needs testing is internal by
# design — a user has no business calling `json_parse()` — and a test suite that could only see
# the exported surface would be testing the least interesting third of the package.
mz <- asNamespace("mzLibR")

# Abandon the current test without failing it.
#
# Signalled rather than returned, so it works from inside a helper several frames down. Not an
# `error` subclass, so a `tryCatch(error = )` in the code under test cannot swallow it.
skip <- function(reason) {
  stop(structure(
    class = c("mz_skip", "condition"),
    list(message = reason, call = NULL)
  ))
}

# Locate a recorded fixture. The suite runs both from tests/ (via run-tests.R) and from the
# package root, so both are tried rather than assuming a working directory.
#
# Defined here rather than in the first test file that needed it: files are sourced in
# alphabetical order, so a helper living in test-pride.R is invisible to test-peptidoform.R.
fixture_path <- function(name) {
  for (candidate in c(file.path("fixtures", name), file.path("tests", "fixtures", name))) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  stop("fixture not found: ", name)
}

skip_if <- function(condition, reason) {
  if (isTRUE(condition)) skip(reason)
  invisible(NULL)
}

# A recorded payload's `data`, from tests/fixtures/. Used by the stats, diff and ptmqtl files.
recorded_data <- function(name) {
  path <- fixture_path(name)
  text <- rawToChar(readBin(path, what = "raw", n = file.info(path)$size))
  Encoding(text) <- "UTF-8"
  parsed <- mz$json_parse(text)
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}

# Two tables equal column by column, except that numbers may differ by `relative` - the last digits
# a platform's math library moves (1e-9 for closed-form fits, 1e-6 for iterative optimisers, as
# pyMzLib's live tests hold them). List columns are compared cell by cell. Returns "" or what differs.
records_disagree <- function(live, recorded, relative = 1e-9) {
  if (!identical(names(live), names(recorded))) {
    return(paste0("columns ", paste(names(live), collapse = ","), " vs ", paste(names(recorded), collapse = ",")))
  }
  if (nrow(live) != nrow(recorded)) {
    return(paste0(nrow(live), " rows vs ", nrow(recorded)))
  }
  close <- function(a, b) {
    if (is.null(a) || is.null(b)) {
      return(is.null(a) && is.null(b))
    }
    if (is.numeric(a) && is.numeric(b)) {
      known <- !is.na(b)
      return(length(a) == length(b) && identical(is.na(a), is.na(b)) &&
        all(abs(a[known] - b[known]) <= relative * pmax(abs(b[known]), 1e-300)))
    }
    identical(a, b)
  }
  for (name in names(recorded)) {
    a <- live[[name]]
    b <- recorded[[name]]
    if (is.list(b)) {
      for (i in seq_along(b)) {
        if (!close(a[[i]], b[[i]])) return(paste0(name, "[", i, "]"))
      }
    } else if (!close(a, b)) {
      return(name)
    }
  }
  ""
}

# The live tests skip on this and nothing else. See `conditions.R` for why a timeout must not
# be widened into it.
skip_if_no_bridge <- function() {
  skip_if(!mz$bridge_available(), "no bridge staged")
}

# Skip a LIVE test unless the staged bridge dispatches `verb`.
#
# The weekly live job runs against the newest published bridge, which can be older than the verbs
# this package already projects: a verb new in pyMzLib 0.2.0 skips against a 0.1.x bridge rather
# than failing, and runs as soon as a bridge that has it is published.
skip_unless_bridge_has <- function(verb) {
  bridge <- Sys.getenv("MZLIB_BRIDGE", "")
  skip_if(!nzchar(bridge), "no bridge staged (set MZLIB_BRIDGE)")
  old <- options(mzlibr.bridge = bridge)
  on.exit(options(old), add = TRUE)
  verbs <- tryCatch(mzlibr_bridge_version()$verbs, error = function(e) NA_character_)
  skip_if(!verb %in% verbs, paste0("the staged bridge does not dispatch '", verb, "'"))
  invisible(NULL)
}

# A recorded payload, unwrapped from its envelope if it has one.
recorded_data <- function(name) {
  path <- fixture_path(name)
  text <- rawToChar(readBin(path, what = "raw", n = file.info(path)$size))
  Encoding(text) <- "UTF-8"
  parsed <- mz$json_parse(text)
  if (is.list(parsed) && !is.null(parsed$ok) && "data" %in% names(parsed)) parsed$data else parsed
}

# Run `call()` against a stand-in for the bridge, through the package's own transport and verb
# check, and return what it sent and what it returned.
#
# `version` answers the verb check with pyMzLib's recorded `version` payload (every verb of the
# bridge this package pins), and every other call is answered by `answer`: a recorded fixture's
# name, or a function(args, stdin) returning the bridge's stdout. The bridge path is a fresh file,
# so the verb check's per-bridge cache cannot leak between tests.
recorded_bridge_call <- function(answer, call, version = "bridge_version.json") {
  path <- tempfile("mzlib-bridge-", fileext = if (.Platform$OS.type == "windows") ".exe" else "")
  writeLines("not a real executable", path)
  seen <- new.env(parent = emptyenv())
  seen$calls <- list()
  envelope <- function(name) {
    path <- fixture_path(name)
    text <- rawToChar(readBin(path, what = "raw", n = file.info(path)$size))
    if (startsWith(trimws(text), "{\"ok\"")) text else paste0("{\"ok\":true,\"data\":", text, ",\"error\":null}")
  }
  runner <- function(exe, args, stdin = NULL, timeout = NULL) {
    seen$calls[[length(seen$calls) + 1L]] <- list(args = args, stdin = stdin, timeout = timeout)
    stdout <- if (identical(args, "version")) {
      envelope(version)
    } else {
      seen$args <- args
      seen$stdin <- stdin
      seen$timeout <- timeout
      if (is.function(answer)) answer(args, stdin) else envelope(answer)
    }
    list(stdout = stdout, stderr = "", status = 0L, timed_out = FALSE)
  }
  original <- mz$bridge_run
  utils::assignInNamespace("bridge_run", runner, "mzLibR")
  old <- options(mzlibr.bridge = path)
  on.exit(
    {
      utils::assignInNamespace("bridge_run", original, "mzLibR")
      options(old)
      unlink(path)
    },
    add = TRUE
  )
  seen$result <- call()
  seen
}

# A bridge error envelope, as `recorded_bridge_call()` takes an answer.
bridge_error_answer <- function(type, message) {
  function(args, stdin) {
    paste0("{\"ok\":false,\"data\":null,\"error\":{\"type\":\"", type, "\",\"message\":\"", message, "\"}}")
  }
}

test_that <- function(description, code) {
  outcome <- tryCatch(
    {
      force(code)
      "pass"
    },
    mz_skip = function(condition) {
      mz_results$messages <- c(
        mz_results$messages,
        paste0("SKIP  ", description, " — ", conditionMessage(condition))
      )
      "skip"
    },
    error = function(condition) {
      mz_results$messages <- c(
        mz_results$messages,
        paste0("FAIL  ", description, "\n      ", conditionMessage(condition))
      )
      "fail"
    }
  )
  mz_results[[outcome]] <- mz_results[[outcome]] + 1L
  invisible(outcome)
}

mz_fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

expect_true <- function(value, info = "") {
  if (!isTRUE(value)) {
    mz_fail("expected TRUE, got ", paste(format(value), collapse = " "), if (nzchar(info)) paste0(" — ", info))
  }
  invisible(TRUE)
}

expect_false <- function(value, info = "") {
  if (!isFALSE(value)) {
    mz_fail("expected FALSE, got ", paste(format(value), collapse = " "), if (nzchar(info)) paste0(" — ", info))
  }
  invisible(TRUE)
}

expect_equal <- function(actual, expected, info = "") {
  if (!isTRUE(all.equal(actual, expected))) {
    mz_fail(
      "expected ", paste(deparse(expected), collapse = " "),
      ", got ", paste(deparse(actual), collapse = " "),
      if (nzchar(info)) paste0(" — ", info)
    )
  }
  invisible(TRUE)
}

expect_identical <- function(actual, expected, info = "") {
  if (!identical(actual, expected)) {
    mz_fail(
      "expected ", paste(deparse(expected), collapse = " "),
      ", got ", paste(deparse(actual), collapse = " "),
      if (nzchar(info)) paste0(" — ", info)
    )
  }
  invisible(TRUE)
}

# Assert that `code` fails, and optionally that it fails as a particular condition class with a
# particular substring in its message.
#
# The message check is not decoration. Most of what this package promises a stuck user is *in*
# the error text — the three ways to stage a bridge, the reason a `.gz` filter matched nothing —
# so an assertion that something merely failed would pass while the useful half rotted away.
expect_error <- function(code, class = NULL, contains = NULL) {
  condition <- tryCatch(
    {
      force(code)
      NULL
    },
    error = function(e) e
  )

  if (is.null(condition)) {
    mz_fail("expected an error, got none")
  }
  if (!is.null(class) && !inherits(condition, class)) {
    mz_fail(
      "expected condition class '", class, "', got: ",
      paste(class(condition), collapse = ", ")
    )
  }
  if (!is.null(contains)) {
    for (needle in contains) {
      if (!grepl(needle, conditionMessage(condition), fixed = TRUE)) {
        mz_fail(
          "expected the message to contain '", needle, "', got:\n      ",
          conditionMessage(condition)
        )
      }
    }
  }
  invisible(condition)
}

# Assert that `code` warns, and optionally that the warning says a particular thing.
#
# The message check matters as much as the fact of the warning: a warning nobody can act on is
# noise, and the ones here name an issue number and a remedy.
expect_warning <- function(code, contains = NULL) {
  warned <- NULL
  withCallingHandlers(
    force(code),
    warning = function(w) {
      warned <<- w
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(warned)) {
    mz_fail("expected a warning, got none")
  }
  for (needle in contains) {
    if (!grepl(needle, conditionMessage(warned), fixed = TRUE)) {
      mz_fail(
        "expected the warning to contain '", needle, "', got:\n      ",
        conditionMessage(warned)
      )
    }
  }
  invisible(warned)
}

# Assert that `code` does NOT warn.
expect_no_warning <- function(code) {
  warned <- NULL
  withCallingHandlers(
    force(code),
    warning = function(w) {
      warned <<- w
      invokeRestart("muffleWarning")
    }
  )
  if (!is.null(warned)) {
    mz_fail("expected no warning, got: ", conditionMessage(warned))
  }
  invisible(NULL)
}

# Print the tally and fail the R CMD check run if anything failed.
mz_report <- function() {
  for (line in mz_results$messages) {
    cat(line, "\n", sep = "")
  }
  cat(sprintf(
    "\n%d passed, %d failed, %d skipped\n",
    mz_results$pass, mz_results$fail, mz_results$skip
  ))
  if (mz_results$fail > 0L) {
    stop(sprintf("%d test(s) failed.", mz_results$fail), call. = FALSE)
  }
  invisible(NULL)
}
