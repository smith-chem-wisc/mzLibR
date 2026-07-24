# Findings

Defects and documentation gaps surfaced while building mzLibR.

Two standing rules govern what happens to anything on this page:

1. **A bug in mzLib is filed as an issue on mzLib.** A workaround in the bridge or in a binding is
   not the fix — mzLib bugs need to be fixed in mzLib, or MetaMorpheus and every other consumer
   stays broken and the next binding rediscovers it from scratch.
2. **A documentation lesson is back-ported to every binding.** The traps live in mzLib's behaviour,
   not in any one language, so wording that rescues an R user rescues a Python and a Rust one too.

mzLibR is the third binding, so most of mzLib's own defects were already found and filed by pyMzLib
and mzLibRust. This page records what was **new** to the R port: two documentation defects (one in a
*sibling's* docs), one capability gap, and the confirmations and language-specific decisions worth
keeping. The known mzLib defects the R arms re-confirmed are listed at the end.

---

## 1. A sibling's bake-off ground truth is wrong: 254/269 should be 243/257 — mzLibRust docs

**Status:** filed as [smith-chem-wisc/mzLibRust#1](https://github.com/smith-chem-wisc/mzLibRust/issues/1).
**Found by:** establishing mzLibR's own bake-off ground truth independently (from the shell, never
through the package — the discipline that exists precisely to catch this).

mzLibRust's `design/bakeoff/DESIGN.md` ground-truth table gives albumin (P02768) at `min_length = 1`
as **254** (`trypsin|P`) and **269** (`trypsin`). Neither reproduces. An exhaustive sweep — both
proteases × missed-cleavages 0–2 × min-length 1/7 × modifications on/off, counting distinct
sequences *and* peptidoforms — gives **243** and **257** at the default `missed_cleavages = 2`, and
254/269 appear at no setting at all.

**The tell:** mzLibRust's own `pep_armB_log.jsonl` already recorded 243. The ground-truth table
disagreed with the arm that measured against it, and the arm was right — the "numbers in
documentation rot" lesson sitting inside the ground-truth file itself.

**Handling here:** mzLibR's `design/bakeoff/DESIGN.md` and `ground-truth.md` carry the verified
values and the exact commands. This also corrected mzLibR's own PLAN.md, which had inherited the 254.

---

## 2. The modification census asserts "22 of 24 glycation" — a count it cannot show — pyMzLib / mzLibRust docs

**Status:** filed as [pyMzLib#11](https://github.com/smith-chem-wisc/pyMzLib/issues/11) and
[mzLibRust#2](https://github.com/smith-chem-wisc/mzLibRust/issues/2).
**Found by:** both peptidoform arms of the mzLibR bake-off, independently.

The census documentation states "22 of the 24 excluded albumin features are N-linked (Glc)
(glycation) lysine". But the census only surfaces the coarse UniProt **feature type**, which for all
24 is `glycosylation site` — the finer per-feature modification *name* (the source of the "22") is
not in the census output. A user sees `24 × glycosylation site` in `by_type`, reads "22 of 24
glycation", and cannot reconcile them.

**Impact:** this is the same class as finding 1 — asserting a number the tool cannot reproduce. It
originated in pyMzLib's docs and was inherited by both later bindings.

**Handling here:** mzLibR's `?census_explain` now frames the two levels explicitly — the census
reports the feature *type* (`glycosylation site`, all 24), and *within that* 22 carry the specific
UniProt modification name, which you confirm from the UniProt entry, not from the census. The "22" is
no longer presented as census output.

---

## 3. PRIDE: no way to fetch the true FTP file list or project size — the bridge

**Status:** filed as [pyMzLib#12](https://github.com/smith-chem-wisc/pyMzLib/issues/12) (the wire
capability lives in the bridge, so all three bindings gain it at once).
**Found by:** the mzLibR bake-off — it was Arm B's single external lookup.

The PRIDE manifest comes from the REST API, which is knowingly incomplete: PXD000001 lists **8**
files (0.51 GB) where the FTP tree holds **13** (~1.44 GB). The bindings document this candidly but
offer no function to fetch the complete list or the true size — so the one question the biologist
could not answer inside the tool was "how big is this project", and they left it for the FTP tree by
hand.

**Suggested direction:** a `pride ftp-files` verb, or at minimum a `files_hidden_from_rest`
attribute on the manifest, plus a companion that reports actual (compressed) transfer size.

---

## 4. R-specific decisions worth recording (not defects)

These are not bugs; they are places where R's semantics changed what the port needed to do, recorded
so the next reader does not re-derive them. All are in `docs/test-parity.md` too.

- **The embedded-null argument guard both siblings carry is unreachable in R.** A character vector
  may not contain a nul — `"a\0b"` is a parse error and `rawToChar()` refuses to build one. The check
  would be dead code, so it is a test asserting the *language* refuses the input, not a guard.
- **`jsonlite` was dropped for correctness, not just weight.** Its `null`→`NULL` deletes list
  elements, so a protein-intensity vector would come back silently *short* — and FlashLFQ emits
  `null` for exactly the protein it could not resolve. `R/json.R` reads `null` as `NA`.
- **`substring()` errors on zero-length positions** rather than returning `character(0)`, which made
  empty bridge output arrive as a bare `simpleError` instead of a typed protocol error. Guarded.

---

## Known mzLib defects the R arms re-confirmed (already filed)

The bake-off surfaced **no new mzLib defect**. The biologist arms independently walked into, and
correctly caught, defects pyMzLib and mzLibRust had already filed:

- **ETD/ECD emit a spurious `y` series**, ~⅓ of the fragment list — [mzLib#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109)
  (PR #1114 open). The mzLibR tests assert what mzLib *currently does*, so they fail when the fix lands.
- **`trypsin|P` applies the Keil rule and plain `trypsin` does not**, the reverse of MaxQuant/Mascot —
  [mzLib#1106](https://github.com/smith-chem-wisc/mzLib/issues/1106).
- **Glycosylation-site annotations are dropped on feature type with no report** —
  [mzLib#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112). The exclusion is correct; the
  silence is the defect. (Finding 2 above is a *documentation* refinement of how mzLibR explains it.)
