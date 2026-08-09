# Locating and invoking the mzLib bridge executable.
#
# This is the only file in mzLibR that knows the bridge exists. Everything above it sees
# ordinary R functions and data.frames. That boundary is deliberate: the transport - today, a
# self-contained .NET executable invoked once per call - could be replaced by a long-lived
# local server without any exported function changing.
#
# It is the direct counterpart of pyMzLib's `_bridge.py` and mzLibRust's `src/bridge.rs`, and
# for the same reason: the wire contract is language-neutral by design, so a third binding is a
# transport file and some constructors, not a third implementation of mzLib.
#
# No `processx`. Base R's `system2()` has carried a real `timeout` since R 3.5.0, takes stdin
# through `input=`, and writes stdout and stderr to separate files, which is every capability
# this transport needs. Both halves of that were verified against the real bridge on Windows
# before the dependency was dropped - see `tests/test-bridge.R`.

# Wire-format version this package understands. The bridge reports its own; a mismatch means
# the two halves were built from different sources.
MZLIB_PROTOCOL_VERSION <- 1L

# Ways to point mzLibR at a bridge, in the order they are consulted.
MZLIB_BRIDGE_OPTION <- "mzlibr.bridge"
MZLIB_BRIDGE_ENV_VAR <- "MZLIB_BRIDGE"

# The error type the bridge uses for availability failures. Must equal
# `Program.ServiceUnavailableType` on the C# side - the two halves agree by this string and
# nothing else.
MZLIB_SERVICE_UNAVAILABLE_TYPE <- "ServiceUnavailable"

# `system2()` reports a timed-out command with this exit status.
MZLIB_TIMEOUT_STATUS <- 124L

# ---------------------------------------------------------------- locating the payload

# The name of the bridge executable on this platform.
bridge_executable_name <- function() {
  if (.Platform$OS.type == "windows") "mzlib-bridge.exe" else "mzlib-bridge"
}

# The staging subdirectory for a named platform.
#
# These strings must equal the .NET runtime identifiers `publish-bridge.ps1` stages under, or
# nothing is ever found. Split out from the live lookup so the mapping can be tested without a
# machine of each kind, and so all three bindings can be checked against one table.
bridge_platform_tag <- function(sysname = Sys.info()[["sysname"]],
                                machine = Sys.info()[["machine"]]) {
  prefix <- switch(tolower(sysname),
    windows = "win",
    linux = "linux",
    darwin = ,
    macos = "osx",
    NULL
  )
  if (is.null(prefix)) {
    stop(mzlib_bridge_not_found(paste0("Unsupported platform: ", sysname, " ", machine)))
  }

  suffix <- switch(tolower(machine),
    "x86-64" = ,
    "x86_64" = ,
    amd64 = ,
    x64 = "x64",
    arm64 = ,
    aarch64 = "arm64",
    tolower(machine)
  )

  paste0(prefix, "-", suffix)
}

# Where a downloaded bridge is cached.
#
# `tools::R_user_dir()` does exactly this, but it arrived in R 4.0 and this package supports
# 3.5, so the same locations are resolved by hand. The layout matches `R_user_dir()`'s on
# purpose: a user on old R and a user on new R share one cache rather than downloading 140 MB
# twice.
mzlibr_cache_dir <- function() {
  explicit <- Sys.getenv("R_USER_CACHE_DIR", "")
  if (nzchar(explicit)) {
    return(file.path(explicit, "R", "mzLibR"))
  }

  if (.Platform$OS.type == "windows") {
    local <- Sys.getenv("LOCALAPPDATA", "")
    if (nzchar(local)) {
      return(file.path(local, "R", "cache", "R", "mzLibR"))
    }
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    return(path.expand("~/Library/Caches/org.R-project.R/R/mzLibR"))
  } else {
    xdg <- Sys.getenv("XDG_CACHE_HOME", "")
    if (nzchar(xdg)) {
      return(file.path(xdg, "R", "mzLibR"))
    }
  }

  path.expand("~/.cache/R/mzLibR")
}

# What to tell someone who has no bridge at all.
#
# Split out from the lookup so it can be tested for its content without the test depending on
# whether the machine running it happens to have a bridge staged - the mistake that made
# mzLibRust's equivalent test pass in CI and fail for every developer who had one.
#
# Every remedy is named. A bare "not found" leaves the reader with nowhere to go, and the error
# message is the only documentation a stuck user reliably reads.
bridge_missing_message <- function(cache_path) {
  paste0(
    "No mzLib bridge executable found.\n\n",
    "Three ways to fix it, cheapest first:\n",
    "  1. mzlibr_install_bridge() - downloads one into ", cache_path, "\n",
    "  2. Sys.setenv(", MZLIB_BRIDGE_ENV_VAR, " = \"/path/to/", bridge_executable_name(),
    "\") - for a bridge you already have,\n",
    "     for example the one pyMzLib stages under pkg/python/src/pymzlib/_dotnet/<rid>/.\n",
    "  3. options(", MZLIB_BRIDGE_OPTION, " = \"/path/to/", bridge_executable_name(),
    "\") - same thing, scoped to this session.\n\n",
    "Overriding the bridge is also how you relink a modified mzLib without rebuilding this ",
    "package, which LGPL section 4 requires mzLibR to allow."
  )
}

#' Path of the bridge executable mzLibR will use
#'
#' Resolution order: the `mzlibr.bridge` option, then the `MZLIB_BRIDGE` environment variable,
#' then a bridge downloaded by [mzlibr_install_bridge()].
#'
#' @return A single file path.
#' @export
mzlibr_bridge_path <- function() {
  from_option <- getOption(MZLIB_BRIDGE_OPTION, default = NULL)
  if (!is.null(from_option)) {
    if (!file.exists(from_option)) {
      stop(mzlib_bridge_not_found(paste0(
        "options(", MZLIB_BRIDGE_OPTION, ") points at '", from_option,
        "', which is not a file."
      )))
    }
    return(from_option)
  }

  from_env <- Sys.getenv(MZLIB_BRIDGE_ENV_VAR, "")
  if (nzchar(from_env)) {
    if (!file.exists(from_env)) {
      stop(mzlib_bridge_not_found(paste0(
        MZLIB_BRIDGE_ENV_VAR, " points at '", from_env, "', which is not a file."
      )))
    }
    return(from_env)
  }

  cached <- file.path(mzlibr_cache_dir(), bridge_platform_tag(), bridge_executable_name())
  if (file.exists(cached)) {
    return(cached)
  }

  stop(mzlib_bridge_not_found(bridge_missing_message(cached)))
}

# Is there a bridge to talk to? Used by the live tests to skip rather than fail.
bridge_available <- function() {
  !inherits(try(mzlibr_bridge_path(), silent = TRUE), "try-error")
}

# ---------------------------------------------------------------- running the process

# Quote one argument for the shell `system2()` will hand the command to.
#
# `system2()` quotes the *command* for you and the arguments not at all, so this is the
# caller's job - and only for the arguments, since quoting the command twice breaks it. It
# matters the moment a path contains a space, which on Windows is the normal case: a project
# directory under "My Documents" would otherwise arrive at the bridge split into pieces.
#
# The two shells want different quoting, and `shQuote()` knows both.
bridge_quote <- function(argument) {
  shQuote(argument, type = if (.Platform$OS.type == "windows") "cmd" else "sh")
}

# Read a file back as UTF-8 whatever the machine's locale is.
#
# The bridge always writes UTF-8. `readLines()` would decode it in the native encoding, which
# on an older Windows is not UTF-8, and a file name with an accent in it would arrive mangled -
# then fail to match anything, in a way that looks like a PRIDE problem rather than an ours
# problem.
bridge_read_utf8 <- function(path) {
  if (!file.exists(path)) {
    return("")
  }
  size <- file.info(path)$size
  if (is.na(size) || size <= 0) {
    return("")
  }
  bytes <- readBin(path, what = "raw", n = size)
  text <- rawToChar(bytes)
  Encoding(text) <- "UTF-8"
  text
}

# Run the bridge once and collect everything it produced.
#
# This is the seam pyMzLib gets by monkeypatching `subprocess.run` and mzLibRust gets from a
# `Runner` trait. In R it is simply a function argument, which is the cheapest of the three:
# the failure paths that matter - a process that dies silently, output that is not JSON, a
# timeout, a binary that will not launch - are the code most likely to be wrong and the least
# convenient to provoke with a real executable.
#
# Returns a list of `stdout`, `stderr`, `status` and `timed_out`.
bridge_run <- function(exe, args, stdin = NULL, timeout = NULL) {
  out_file <- tempfile("mzlibr-out-")
  err_file <- tempfile("mzlibr-err-")
  on.exit(unlink(c(out_file, err_file)), add = TRUE)

  # `system2(timeout=)` raises a warning as well as returning 124. The warning is muffled here
  # and turned into a typed condition below, so a timeout does not arrive as console noise
  # alongside whatever the caller does next.
  timed_out <- FALSE
  status <- withCallingHandlers(
    tryCatch(
      system2(
        exe,
        args = vapply(args, bridge_quote, character(1L), USE.NAMES = FALSE),
        stdout = out_file,
        stderr = err_file,
        stdin = "",
        input = stdin,
        timeout = if (is.null(timeout)) 0 else timeout
      ),
      error = function(e) {
        stop(mzlib_bridge_not_found(paste0(
          "Could not run the mzLib bridge at '", exe, "': ", conditionMessage(e)
        )))
      }
    ),
    warning = function(w) {
      if (grepl("timed out", conditionMessage(w), fixed = TRUE)) {
        timed_out <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  list(
    stdout = bridge_read_utf8(out_file),
    stderr = bridge_read_utf8(err_file),
    status = as.integer(status),
    timed_out = timed_out || identical(as.integer(status), MZLIB_TIMEOUT_STATUS)
  )
}

# ---------------------------------------------------------------- the envelope

# Turn a completed invocation into the `data` payload, or raise the right condition.
#
# Pure: it takes what `bridge_run()` returned and nothing else, so every failure path can be
# tested by handing it a list, with no subprocess and no mocking framework.
bridge_decode <- function(output, timeout = NULL) {
  if (isTRUE(output$timed_out)) {
    stop(mzlib_timeout(if (is.null(timeout)) NA_real_ else timeout))
  }

  text <- trimws(output$stdout)
  if (!nzchar(text)) {
    # A silent non-zero exit means the process died before it could report anything. stderr is
    # the only evidence left, so it goes in the message rather than being discarded.
    stderr_text <- trimws(output$stderr)
    stop(mzlib_protocol_error(paste0(
      "mzLib bridge exited with code ",
      if (length(output$status) == 1L && !is.na(output$status)) output$status else "unknown",
      " and no output. stderr: ",
      if (nzchar(stderr_text)) stderr_text else "(empty)"
    )))
  }

  envelope <- json_parse(text)
  if (!is.list(envelope) || is.null(names(envelope))) {
    stop(mzlib_protocol_error(paste0(
      "mzLib bridge returned JSON that is not an envelope: ", substr(text, 1L, 400L)
    )))
  }

  if (isTRUE(envelope$ok)) {
    return(envelope$data)
  }

  error <- envelope$error
  if (!is.list(error)) {
    error <- list()
  }
  error_type <- if (is.character(error$type)) error$type else "Unknown"
  message <- if (is.character(error$message)) {
    error$message
  } else {
    "mzLib reported a failure with no message."
  }

  if (identical(error_type, "usage")) {
    stop(mzlib_usage_error(message))
  }
  if (identical(error_type, MZLIB_SERVICE_UNAVAILABLE_TYPE)) {
    stop(mzlib_service_unavailable(message, error_type = error_type))
  }
  stop(mzlib_bridge_error(message, error_type = error_type))
}

# Reject a timeout that cannot mean anything, before spawning a process.
#
# A zero timeout looks exactly like a service that never answered, and a negative or infinite
# one is a typo the caller wants told about now rather than in ten minutes.
bridge_check_timeout <- function(timeout) {
  if (is.null(timeout)) {
    return(invisible(NULL))
  }
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) ||
    !is.finite(timeout) || timeout <= 0) {
    stop(mzlib_usage_error(
      "timeout must be a single finite number of seconds greater than zero; pass NULL to wait indefinitely."
    ))
  }
  invisible(NULL)
}

# Run one bridge command and return the decoded `data` payload.
#
# `args` is the verb and its options, already quoted. `stdin` carries a payload that would not
# fit on the command line - argv has a hard ceiling of roughly 32 KB and a real experiment's
# worth of file names goes straight past it. `timeout` of `NULL` waits indefinitely, which is
# the right default for a download that may run for an hour.
# Both siblings guard here against an argument containing an embedded null, which would
# terminate the string early somewhere inside the C runtime and hand the bridge a truncated
# accession. There is no such guard in this version because R cannot express the input: a
# character vector may not contain a nul, `"a\0b"` is a parse error, and `rawToChar()` refuses
# to build one. The check would be unreachable, and unreachable checks rot.
bridge_invoke <- function(args, stdin = NULL, timeout = NULL, runner = bridge_run) {
  bridge_check_timeout(timeout)

  exe <- mzlibr_bridge_path()
  bridge_decode(runner(exe, args, stdin = stdin, timeout = timeout), timeout = timeout)
}

# ---------------------------------------------------------------- version handshake

#' Version information reported by the mzLib bridge
#'
#' The whole transport story end to end - locate the executable, run it, parse an envelope,
#' agree on a wire format - in one call with no network and no arguments. It is the first thing
#' to make work and the first thing to check when something else is wrong.
#'
#' @param runner The function used to run the bridge. Present so the transport's failure paths
#'   can be tested; not something a caller normally sets.
#' @return A list with `bridge`, `protocol`, `runtime` and `mzlib`.
#'
#'   `mzlib` is which mzLib the bridge was built against, as `1.0.0+<commit>`, and is
#'   `NA_character_` when the bridge did not report one - either because it predates the field,
#'   or because its build recorded no source commit. This package installs a bridge rather than
#'   building one, so it is the only way to ask which mzLib is actually running.
#'
#'   It is deliberately **not** a compatibility check. `protocol` is that, and it is what this
#'   function verifies. `mzlib` is for reporting a run, filing a bug, or tying a result to the
#'   library that produced it.
#' @export
mzlibr_bridge_version <- function(runner = bridge_run) {
  data <- bridge_invoke("version", timeout = 60, runner = runner)

  reported <- data$protocol
  if (!is.numeric(reported) || length(reported) != 1L ||
    !identical(as.integer(reported), MZLIB_PROTOCOL_VERSION)) {
    stop(mzlib_protocol_error(paste0(
      "mzLib bridge speaks protocol ",
      if (is.numeric(reported)) format(reported) else "nothing",
      ", but this mzLibR expects ", MZLIB_PROTOCOL_VERSION,
      ". The R package and the bridge were built from different sources."
    )))
  }

  list(
    bridge = as.character(data$bridge),
    protocol = as.integer(reported),
    runtime = as.character(data$runtime),
    # An absent field projects to NA, never to "" or to a made-up placeholder. A bridge built
    # before this field existed simply does not send it, and the honest R answer to "which mzLib?"
    # in that case is "not known", which is what NA means.
    mzlib = wire_field(data, "mzlib", "character", NA_character_)
  )
}
