# Version information reported by the mzLib bridge

Report the bridge's version, its wire protocol, its .NET runtime, the
mzLib it was built from, and every verb it dispatches.

The whole transport story end to end - locate the executable, run it,
parse an envelope, agree on a wire format - in one call with no network
and no arguments. It is the first thing to make work and the first thing
to check when something else is wrong.

## Usage

``` r
mzlibr_bridge_version(runner = bridge_run)
```

## Arguments

- runner:

  The function used to run the bridge. Present so the transport's
  failure paths can be tested; not something a caller normally sets.

## Value

A list with `bridge`, `protocol`, `runtime` and `mzlib`.

`mzlib` is which mzLib the bridge was built against, as
`1.0.0+<commit>`, and is `NA_character_` when the bridge did not report
one - either because it predates the field, or because its build
recorded no source commit. This package installs a bridge rather than
building one, so it is the only way to ask which mzLib is actually
running.

It is deliberately **not** a compatibility check. `protocol` is that,
and it is what this function verifies. `mzlib` is for reporting a run,
filing a bug, or tying a result to the library that produced it.

## Wraps

Wire verb `version`. Generated from the bridge's verb spec
`version.yaml` (bridge commit `5db922d4cfe1`) by `scripts/build-man.R`;
the spec owns these facts, and all three bindings render the same ones.

- [`AssemblyInformationalVersionAttribute`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/UsefulProteomicsDatabases.csproj)
  in `mzLib/UsefulProteomicsDatabases/UsefulProteomicsDatabases.csproj`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

This verb takes no parameters.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `bridge`:

  string; never `NA`. The bridge assembly's version. Not a compatibility
  number: use protocol.

- `protocol`:

  int; never `NA`. The wire-format version; a binding is compatible with
  a bridge by this.

- `runtime`:

  string; never `NA`. The bundled .NET runtime version.

- `mzlib`:

  string; `NA` when the build recorded no mzLib commit. Which mzLib the
  bridge links, as 1.0.0+\<commit\>, read from the linked assembly.

## On the wire but not projected yet

The bridge sends these, and this version of mzLibR does not return them
yet:

- `verbs`:

  arrives with the mzLib 1.0.592 port.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- verbs lists what the bridge ROUTES, not which options each verb takes:
  a bulk option (`--paths-stdin`) on an older bridge that lists the verb
  still fails with that bridge's usage error.

- Bridges built before verbs existed omit the key; a binding must treat
  its absence as 'lists nothing', which correctly refuses every verb
  added since.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.bridge_version`

- Rust (mzLibRust): `mzlib::bridge_version`

- R (mzLibR): `mzlibr_bridge_version`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The Unknown-command message and verbs share one generated list
  (Program.Verbs, from the DeriveVerbList MSBuild task); check_verbs.py
  still parses the switch itself. Both read the same pattern, so they
  cannot disagree unless the switch stops being one arm per line.

- pyMzLib returns version as a plain dict, typed as the TypedDict
  pymzlib.BridgeVersion so its keys are documented and linted; Rust and
  R spellings are the intended ones.

## Examples

``` r

info <- mzlibr_bridge_version()
info$protocol   # the compatibility contract
#> [1] 1
info$mzlib      # which mzLib produced the results, for a methods section
#> [1] "1.0.0+23c2490e10d3ccce71c941bca31610725826ba83"
```
