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

skip_if <- function(condition, reason) {
  if (isTRUE(condition)) skip(reason)
  invisible(NULL)
}

# The live tests skip on this and nothing else. See `conditions.R` for why a timeout must not
# be widened into it.
skip_if_no_bridge <- function() {
  skip_if(!mz$bridge_available(), "no bridge staged")
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
