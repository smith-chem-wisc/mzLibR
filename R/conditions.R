# Everything mzLibR can fail with.
#
# R's condition system is class-based, so the hierarchy pyMzLib builds out of exception
# subclasses is spelled here as a class *vector*: a usage failure is
# `c("mzlib_usage_error", "mzlib_error", "error", "condition")`. A caller who wants to handle
# one kind names it; a caller who wants to handle all of them names `mzlib_error`; a caller who
# does not care catches `error` like any other R failure. Nothing has to be exported for that
# to work, because `tryCatch()` dispatches on the class vector directly:
#
#     tryCatch(
#       pride_list_files("PXD000001"),
#       mzlib_service_unavailable = function(e) message("EBI is having a bad morning"),
#       mzlib_error = function(e) stop(e)
#     )
#
# The classification itself is made in the **bridge**, not here, so every binding over this
# wire inherits the same judgements rather than each re-deriving them.

#' The conditions mzLibR signals, and how to handle each
#'
#' Every failure mzLibR raises is an R condition whose class vector names **what kind** of
#' failure it is, then `mzlib_error`, then `error`. Handle one kind by naming its class in
#' [tryCatch()]; handle all of them by naming `mzlib_error`; or treat them as any other R error.
#' Nothing needs to be imported: `tryCatch()` dispatches on the class vector itself.
#'
#' The classification is made in the **bridge**, not in this package, so pyMzLib, mzLibRust and
#' mzLibR agree on it. Each verb's help page lists, under **Errors**, which of these it raises and
#' when - generated from the bridge's verb spec, where the three kinds are called `usage`,
#' `service_unavailable` and `correctness`.
#'
#' @section The classes:
#'
#' `mzlib_usage_error` - the call was malformed: a missing or invalid argument, a file that is not
#' there, a file with no such view. Raised before any work happens, either by mzLibR itself or by
#' the bridge (its `usage` kind, exit status 2). Fix the call; retrying will not help. pyMzLib
#' raises `UsageError`, mzLibRust returns `MzLibError::Usage`.
#'
#' `mzlib_service_unavailable` - an external service is down, rate-limited, timing out or
#' unreachable: HTTP 408, 429 and 5xx from PRIDE, UniProt or EBI. **Retry later.** A 404 or 400 is
#' deliberately not this class, because a wrong URL or a malformed request is a defect, and
#' excusing it as an outage would hide it. The live tests skip on this class and nothing else.
#' pyMzLib raises `ServiceUnavailableError`.
#'
#' `mzlib_bridge_error` - mzLib itself reported a failure (the spec's `correctness` kind): a file
#' it could not parse, a vendor library that will not load on this platform. Its `error_type`
#' field carries the .NET exception type, e.g. `"HttpRequestException"`, so a handler can tell
#' failures apart without parsing prose. pyMzLib raises `BridgeError`.
#'
#' `mzlib_project_not_found` - PRIDE has no project with that accession. PRIDE answers an unknown
#' accession with an empty list rather than a 404, and an empty list is indistinguishable from a
#' project with nothing matching, so this is raised rather than returning `list()`.
#'
#' `mzlib_timeout` - the bridge did not finish within `timeout` seconds. Deliberately **not**
#' `mzlib_service_unavailable`: a timeout can mean a slow service, but equally a wedged bridge, a
#' corrupt executable or a timeout set too short, and mzLibR will not guess which. Its `seconds`
#' field is the limit that was exceeded.
#'
#' `mzlib_bridge_not_found` - no bridge executable could be located, or it could not be run. The
#' message names the three ways to provide one; see [mzlibr_install_bridge()] and
#' [mzlibr_bridge_path()].
#'
#' `mzlib_protocol_error` - the bridge answered with something this version cannot read: no output
#' from a process that died (stderr is quoted in the message), text that is not JSON, or a wire
#' protocol other than the one mzLibR speaks. Usually a bridge and a package from different
#' releases; [mzlibr_bridge_version()] says which.
#'
#' @section Fields:
#'
#' Every condition has `message` and a `NULL` `call` - the useful context is in the message, not
#' in a call stack pointing at the transport. `mzlib_bridge_error` and `mzlib_service_unavailable`
#' add `error_type`; `mzlib_timeout` adds `seconds`.
#'
#' @seealso [mzlibr_bridge_version()], [mzlibr_install_bridge()]
#' @name mzlib_error
#' @aliases mzlib_usage_error mzlib_service_unavailable mzlib_bridge_error mzlib_project_not_found
#'   mzlib_timeout mzlib_bridge_not_found mzlib_protocol_error
#' @examples
#' # A malformed call is refused before anything runs.
#' tryCatch(
#'   readers_read_spectra("run.mzML", limit = -1),
#'   mzlib_usage_error = function(e) conditionMessage(e)
#' )
#'
#' # Every class also carries mzlib_error, so one handler can catch them all.
#' old <- options(mzlibr.bridge = file.path(tempdir(), "no-such-bridge"))
#' failure <- tryCatch(mzlibr_bridge_version(), mzlib_error = function(e) e)
#' class(failure)
#' options(old)
NULL

# Build one of ours.
#
# `call. = NULL` deliberately: the useful context is in the message - the accession, the file
# name, the three ways to fix it - and an R call stack pointing at `bridge_invoke()` only adds
# noise for a user who never typed that name.
mzlib_condition <- function(class, message, ...) {
  structure(
    class = c(class, "mzlib_error", "error", "condition"),
    list(message = message, call = NULL, ...)
  )
}

# A call was malformed - a missing or invalid argument. Raised before any work happens.
#
# Either mzLibR rejected the input, or the bridge answered `{"type": "usage"}` and exit 2.
mzlib_usage_error <- function(message) {
  mzlib_condition("mzlib_usage_error", message)
}

# An external service is unavailable - down, rate-limited, timing out, or unreachable.
#
# Its own class because the difference between "the repository is having a bad morning" and
# "something is broken" is the difference between retrying later and filing a bug. HTTP 408,
# 429 and 5xx count; 404 and 400 do not, because a wrong URL or a malformed request is our
# problem and excusing it as an outage would hide a real defect.
#
# This is the class the live tests skip on. Nothing else may be widened into it.
mzlib_service_unavailable <- function(message, error_type = "ServiceUnavailable") {
  mzlib_condition("mzlib_service_unavailable", message, error_type = error_type)
}

# mzLib itself reported a failure.
#
# `error_type` carries the .NET exception type name, e.g. `HttpRequestException`, so a caller
# can tell a network failure from a bad accession without parsing prose.
mzlib_bridge_error <- function(message, error_type) {
  mzlib_condition("mzlib_bridge_error", message, error_type = error_type)
}

# The bridge process did not finish in time.
#
# Deliberately **not** `mzlib_service_unavailable`, and the distinction is the entire point. A
# subprocess timeout has several possible causes and only one of them is a slow service: the
# bridge may be wedged, the executable may be corrupt, antivirus may be holding it, or the
# caller may simply have passed a timeout that was too short. Reporting all of that as "the
# repository is down" is how a real bug gets skipped by every live suite and never seen again -
# which is exactly what pyMzLib did until it was found.
mzlib_timeout <- function(seconds) {
  mzlib_condition(
    "mzlib_timeout",
    paste0(
      "mzLib bridge did not finish within ", format(seconds), "s. This may mean the service ",
      "is slow, but it can equally mean the bridge is wedged or the timeout was too short - ",
      "mzLibR will not guess which."
    ),
    seconds = seconds
  )
}

# The bridge executable could not be located, or could not be run.
mzlib_bridge_not_found <- function(message) {
  mzlib_condition("mzlib_bridge_not_found", message)
}

# The bridge produced output this version cannot interpret.
#
# Empty output from a process that died, something that is not JSON, or a protocol version
# mzLibR does not speak.
mzlib_protocol_error <- function(message) {
  mzlib_condition("mzlib_protocol_error", message)
}

# No project with that accession exists, or it has no files matching.
#
# PRIDE answers an unknown accession with HTTP 200 and an empty list rather than a 404, so a
# naive binding hands back `list()` and the caller reports "0 files, done". An empty list is
# indistinguishable from "this project genuinely has nothing matching", so a typo'd accession
# becomes a script that succeeds and is wrong. A wrong answer that looks like a right answer is
# worse than an error, so this is an error.
mzlib_project_not_found <- function(message) {
  mzlib_condition("mzlib_project_not_found", message)
}
