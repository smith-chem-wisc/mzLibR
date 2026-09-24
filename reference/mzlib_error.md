# The conditions mzLibR signals, and how to handle each

Every failure mzLibR raises is an R condition whose class vector names
**what kind** of failure it is, then `mzlib_error`, then `error`. Handle
one kind by naming its class in
[`tryCatch`](https://rdrr.io/r/base/conditions.html); handle all of them
by naming `mzlib_error`; or treat them as any other R error. Nothing
needs to be imported:
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) dispatches on the
class vector itself.

## Details

The classification is made in the **bridge**, not in this package, so
pyMzLib, mzLibRust and mzLibR agree on it. Each verb's help page lists,
under **Errors**, which of these it raises and when - generated from the
bridge's verb spec, where the three kinds are called `usage`,
`service_unavailable` and `correctness`.

## The classes

`mzlib_usage_error` - the call was malformed: a missing or invalid
argument, a file that is not there, a file with no such view. Raised
before any work happens, either by mzLibR itself or by the bridge (its
`usage` kind, exit status 2). Fix the call; retrying will not help.
pyMzLib raises `UsageError`, mzLibRust returns `MzLibError::Usage`.

`mzlib_service_unavailable` - an external service is down, rate-limited,
timing out or unreachable: HTTP 408, 429 and 5xx from PRIDE, UniProt or
EBI. **Retry later.** A 404 or 400 is deliberately not this class,
because a wrong URL or a malformed request is a defect, and excusing it
as an outage would hide it. The live tests skip on this class and
nothing else. pyMzLib raises `ServiceUnavailableError`.

`mzlib_bridge_error` - mzLib itself reported a failure (the spec's
`correctness` kind): a file it could not parse, a vendor library that
will not load on this platform. Its `error_type` field carries the .NET
exception type, e.g. `"HttpRequestException"`, so a handler can tell
failures apart without parsing prose. pyMzLib raises `BridgeError`.

`mzlib_project_not_found` - PRIDE has no project with that accession.
PRIDE answers an unknown accession with an empty list rather than a 404,
and an empty list is indistinguishable from a project with nothing
matching, so this is raised rather than returning
[`list()`](https://rdrr.io/r/base/list.html).

`mzlib_timeout` - the bridge did not finish within `timeout` seconds.
Deliberately **not** `mzlib_service_unavailable`: a timeout can mean a
slow service, but equally a wedged bridge, a corrupt executable or a
timeout set too short, and mzLibR will not guess which. Its `seconds`
field is the limit that was exceeded.

`mzlib_bridge_not_found` - no bridge executable could be located, or it
could not be run. The message names the three ways to provide one; see
[`mzlibr_install_bridge`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md)
and
[`mzlibr_bridge_path`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_path.md).

`mzlib_protocol_error` - the bridge answered with something this version
cannot read: no output from a process that died (stderr is quoted in the
message), text that is not JSON, or a wire protocol other than the one
mzLibR speaks. Usually a bridge and a package from different releases;
[`mzlibr_bridge_version`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_version.md)
says which.

## Fields

Every condition has `message` and a `NULL` `call` - the useful context
is in the message, not in a call stack pointing at the transport.
`mzlib_bridge_error` and `mzlib_service_unavailable` add `error_type`;
`mzlib_timeout` adds `seconds`.

## See also

[`mzlibr_bridge_version`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_version.md),
[`mzlibr_install_bridge`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md)

## Examples

``` r

# A malformed call is refused before anything runs.
tryCatch(
  readers_read_spectra("run.mzML", limit = -1),
  mzlib_usage_error = function(e) conditionMessage(e)
)
#> [1] "limit must be a positive whole number, or NULL for every record; got -1."

# Every class also carries mzlib_error, so one handler can catch them all.
old <- options(mzlibr.bridge = file.path(tempdir(), "no-such-bridge"))
failure <- tryCatch(mzlibr_bridge_version(), mzlib_error = function(e) e)
class(failure)
#> [1] "mzlib_bridge_not_found" "mzlib_error"            "error"                 
#> [4] "condition"             
options(old)
```
