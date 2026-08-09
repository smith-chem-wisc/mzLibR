# mzLibR — status

Written 2026-07-24, at the end of the session that built it.

## Where it stands

**The crate exists, works, and is public:** <https://github.com/smith-chem-wisc/mzLibR>.

mzLib callable from R over the same language-neutral bridge as pyMzLib and mzLibRust, with **no
package dependencies at all** — not `jsonlite`, not `processx`, nothing but base R. **214 tests**
pass with a bridge staged; **199** pass offline with no bridge and no network, the rest skipping
rather than failing. `R CMD check --as-cran` on the built tarball is clean but for the unavoidable
"New submission" note.

Four capabilities, one more than the Rust and Python siblings currently expose:

| module | functions |
|---|---|
| `pride_*` | list, download, filtered download, locations, total size |
| `peptidoform_*` | digest + fragment, census, fragments-by-series, m/z, truncation |
| `flashlfq_*` | quantify, MBR peaks/rescued-count, peptide/protein counts |
| `readers_*` | **identify 29 file types, read the 3 quantifiable ones** — the capability the plan did not know the bridge had; neither sibling implements it |

## Why no dependencies — it was a correctness decision

- **`jsonlite` is gone.** Its auto-simplification and `null`→`NULL` (which *deletes* a list element)
  is the single biggest R-specific hazard: a protein-intensity vector would come back silently
  *short*, and FlashLFQ emits `null` for exactly the protein it could not resolve. `R/json.R` is a
  hand-written reader that never simplifies and reads `null` as `NA`. The trap cannot occur.
- **`processx` is gone.** `system2()` has had a real `timeout` since R 3.5.0 and delivers stdin at
  the same time — verified against the live bridge on Windows before the dependency was dropped.

R floor is **3.5.0**, set solely by `system2(timeout=)`.

## What the port produced beyond the crate

**A wrong number in a sibling.** PLAN §6.2 (and mzLibRust's own `design/bakeoff/DESIGN.md`
ground-truth table) said albumin gives **254/269** distinct sequences at `min_length = 1`. It is
**243/257** — verified by an exhaustive shell sweep, and mzLibRust's *own arm log already had 243*.
The ground-truth table disagreed with the arm that measured against it. Filed as
[mzLibRust#1](https://github.com/smith-chem-wisc/mzLibRust/issues/1).

**One documentation lesson, from the bake-off.** The modification census asserts "22 of 24 are
glycation" — a UniProt modification-*name* count the census only exposes as feature-*type*
(`glycosylation site`, all 24). Asserting a number the tool can't reproduce is the §8.3 trap. Fixed
here; filed for the siblings: [pyMzLib#11](https://github.com/smith-chem-wisc/pyMzLib/issues/11),
[mzLibRust#2](https://github.com/smith-chem-wisc/mzLibRust/issues/2).

**One capability gap.** The PRIDE manifest is knowingly incomplete but there's no way to fetch the
true FTP list/size. Belongs in the bridge:
[pyMzLib#12](https://github.com/smith-chem-wisc/pyMzLib/issues/12).

No **new** mzLib defect — the arms confirmed the known ones (#1109, #1106, #1112).

## The bake-off

Six blinded biologist arms, three tasks × two toolchains ([design](design/bakeoff/DESIGN.md),
[results](design/bakeoff/RESULTS.md), logs and answers alongside):

| | ecosystem | mzLibR |
|---|---|---|
| external lookups | 9 | **2** |
| dead ends | 5 | **1** |
| answers they would publish | ~2 of 9 | **9 of 9** |

The quant Arm A is the thesis: no installable R tool does FlashLFQ-style MBR, so the ecosystem
biologist hand-rolled an XIC engine, found their decoy control non-specific, and **refused to
report an MBR number** — then wished for "an installable R package that does FlashLFQ/IonQuant-style
quant end-to-end," not knowing Arm B existed.

## Two things done differently from the siblings, both deliberate

- **`max_threads` defaults to 1**, not pyMzLib's −1. Above one thread the roll-up drops MBR
  intensities nondeterministically (mzLib#1111); a binding that ships unreproducible-by-default is
  worse than one that diverges from its parent in a documented way. Setting anything else warns.
- **Docs and man pages are generated from the `#'` blocks** by `scripts/build-man.R` (no roxygen2
  dependency), and **name parity is machine-checked** by `scripts/name-parity.R` against pyMzLib and
  the wire — it exits non-zero on any undeclared difference. `docs/name-parity.md` currently: none.

## What is NOT done

1. **Not verified on macOS, Linux, or R 3.5.** Only R 4.6.1 on Windows ran this. `bridge_platform_tag`
   is tested as a table, not on a machine; the `chmod`/Gatekeeper paths in `mzlibr_install_bridge()`
   have never run on their target OS. **Test there before CRAN.** (Memory: `mzlibr-unverified-platforms`.)
2. **Not on CRAN / not submitted.** On hold. The installer's SHA-256 pins point at a pyMzLib
   *pre-release*; a real submission wants a stable bridge release first, plus (1). The exact version
   is deliberately not repeated here — it lives in `R/install-bridge.R` between the generated-pins
   markers, and `.github/workflows/bridge-watch.yml` rewrites it. A version number in prose is one
   nothing regenerates and nothing tests, so it is stale from the day it is written.
3. **The four upstream issues are filed, not resolved** — two on pyMzLib, two on mzLibRust.

## Immediate next steps

1. **Get a non-Windows machine** and run `R CMD check --as-cran` + the offline suite. This is the
   single biggest gap between "works" and "portable."
2. **When a stable bridge release exists**, update the four wheel SHA-256 pins in
   `R/install-bridge.R` (`MZLIB_BRIDGE_WHEELS`) and bump `MZLIB_BRIDGE_VERSION`.
3. **Watch the four sibling issues.** A fix upstream (e.g. mzLib PR #1114 removing the ETD `y`
   series) will require updating the tests here that assert current behaviour.
