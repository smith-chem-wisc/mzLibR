# How the help pages are generated and checked

pyMzLib, mzLibRust and mzLibR project the same bridge, so they document the same verbs. When each
binding wrote that documentation by hand, it drifted: Python said mzLib reads 32 file types and
Rust said 31, this package's help pages were hand copies of Python docstrings, `name-parity.md`
had silently lost the SDRF module, and no binding ran any of its examples.

So the facts about each wire verb are written once, in the bridge, and every binding renders
them and checks its own prose against them. This page is how that works in mzLibR.

## The pieces

```text
bridge (private)                mzLibR (this repository)
design/verbs/*.yaml  ── sync ──▶ docs/specs/*.yaml + *.json + SOURCE      scripts/sync-specs.R
                                    │
                                    ├─▶ man/*.Rd, the fact sections         scripts/build-man.R
                                    ├─▶ docs/name-parity.md, and the lint   scripts/name-parity.R
                                    └─▶ inst/replay/TABLE                   scripts/stage-replay.R
tests/fixtures/*.json ──────────────▶ inst/replay/*.json  (the recordings the examples replay)
```

**A spec** (`<module>.<verb>.yaml`) holds the facts about one wire verb: its parameters with
types, defaults, **units** and ranges; its result fields and columns with **units** and **what
`NA` means**; its error kinds and when each happens; caveats; the mzLib code it wraps at the pin;
DOIs; its name in each binding; recorded example fixtures; and `since`. This package owns the
idiom - function names, argument names, prose, articles. Where the two disagree about a fact,
the spec wins, and if the spec is wrong it is fixed in the bridge.

**The vendored copy.** The bridge repository is private, so `docs/specs/` carries a copy and
`docs/specs/SOURCE` records the bridge commit. Never edit it. Base R has no YAML reader and this
package has no dependencies, so the sync script, which only a maintainer runs, writes each spec a
second time as JSON, which the package's own reader parses. It needs the `yaml` package; nothing
else does.

```sh
Rscript scripts/sync-specs.R --from ../bridge/design/verbs   # wherever your bridge checkout is
```

## Rendering: `@spec` in a `#'` block

A function that projects a wire verb names its spec:

```r
#' @spec readers.read-spectra
#' @export
readers_read_spectra <- function(...)
```

`scripts/build-man.R` then adds to that function's page, in the shared reference-page order all
three bindings use: the spec's one-sentence summary (first in the description), **Wraps**,
**Parameters: units, ranges and defaults**, **Returned fields** and **Columns** (every one with
its type, unit and `NA` meaning), **Errors** (as `mzlib_*` condition classes), **Caveats**,
**Performance**, **Same verb in other bindings**, **Since**, anything the spec lists as not yet
verified, and **References**. `@spec <verb> bulk` renders the many-files shape for a `_many`
function; `@spec <verb> selection` a second function over the same verb, such as
`pride_download_files()`.

Where mzLibR spells a param or field differently, returns it somewhere else, or does not return it
yet, that is declared with a reason in `R_DEVIATIONS` in `scripts/spec-facts.R`, and the page says
so. A field the bridge sends that this package does not project yet is listed under **On the wire
but not projected yet**, so no page claims a field its function does not return.

`build-man.R` reads the INSTALLED package, so always:

```sh
R CMD INSTALL .
Rscript scripts/build-man.R
```

## Checking: `scripts/name-parity.R`

Besides comparing names with pyMzLib and columns with the wire, it checks every help page against
its spec:

- every spec parameter is an argument, or a declared deviation;
- every `@param` for a parameter with a unit names that unit ("Scans to skip" passes, "How many to
  skip" fails);
- every returned field with a unit is named in `@return`, with its unit close by;
- every field the spec says the R object carries is really in it, checked by parsing the verb's
  recorded fixture with the package's own parser;
- every `R_DEVIATIONS` entry names something its spec has, so stale entries cannot pile up;
- every pyMzLib public callable is mapped to an R function or listed with a reason.

```sh
Rscript scripts/name-parity.R <pyMzLib>/pkg/python/src/pymzlib > docs/name-parity.md
```

pyMzLib is compared at the ref in `scripts/pymzlib-ref`: the release whose bridge this package
pins, or the pyMzLib commit a port is running ahead to.

## Examples that run

Every exported function has an `@examples` block, and `R CMD check` runs them all. They need no
bridge, no .NET and no network, because each starts, inside `\dontshow{}`, a stand-in bridge
(`R/replay.R`, a port of pyMzLib's `tests/replay_bridge.py`) that answers from a recording of the
real one:

```r
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' scans <- readers_read_spectra("sliced_ethcd.mzML", limit = 3)
#' scans$records
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
```

The stand-in is an executable invoked exactly as the bridge is, so an example exercises the real
argument handling, transport and parsing. It answers **only a call its recording fits**: the same
file name, the same echoed options, and a `limit`/`offset` that reproduces the recorded window;
a many-files call fits only a many-files recording. Anything else fails the example with every
recording and why it did not fit.

The recordings are pyMzLib's fixtures, byte for byte from `tests/fixtures/`, staged into
`inst/replay/` for the ones an example uses (the specs' `examples`, plus `REPLAY_EXTRA` in the
script). After adding a fixture or an example:

```sh
Rscript scripts/stage-replay.R
```

An example that needs the network goes in `\donttest{}` and only runs when a bridge is configured.
The pkgdown articles in `vignettes/articles/` replay the same way, in a hidden setup chunk.

## In CI

| where | check | fails when |
|---|---|---|
| `ci.yml` `docs` | `sync-specs.R --check` | a `.json` no longer matches its `.yaml` |
| `ci.yml` `docs` | `stage-replay.R --check` | a staged recording or the table is stale |
| `ci.yml` `docs` | `build-man.R`, then `git diff --exit-code man/` | a help page is stale |
| `ci.yml` `docs` | `name-parity.R`, then `git diff --exit-code docs/name-parity.md` | any finding above, or a stale page |
| `ci.yml` `check` | `R CMD check --as-cran` on three OSes and oldrel | an example errors |
| `pkgdown.yml` | `pkgdown::build_site_github_pages()` | an article or example errors |

## When you project a new verb

1. The spec lands in the bridge first. Sync it here.
2. Copy its recorded fixtures from pyMzLib byte for byte into `tests/fixtures/`.
3. Write the function, with `@spec` and an `@examples` block that replays the spec's example.
   Declare any deviation in `scripts/spec-facts.R`.
4. `R CMD INSTALL .`, then `stage-replay.R`, `build-man.R` and `name-parity.R`, until the last is
   clean. Commit what they write.
5. Add the verb to `NEWS.md`, and report the R spelling and `since` version back to the bridge.
