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
