# mzLibR — plan

**mzLib callable from R.** A separate repo and package, the third sibling of pyMzLib and mzLibRust,
built on the **same language-neutral bridge**. Written 2026-07-23 from everything the pyMzLib *and*
mzLibRust builds taught us — so this is not a green-field design, it is the third pass at a proven
one, and the second one took a single session.

> **Pick-up for a fresh chat:** read this whole file. It is long because it is meant to be the only
> thing you need. Then read the three reference sets in §14. The single most important idea is §1;
> the single most valuable content is §6 (every trap, with numbers) and §8 (the documentation
> doctrine), because those are what the two previous builds actually paid for.

---

## 1. Why this is cheap: the bridge is already language-neutral (D6)

pyMzLib's `.NET` **bridge** is a self-contained executable speaking a **versioned JSON envelope over
stdin/stdout** that assumes nothing about its caller. It already does every genuinely hard part:

- the mzLib interop and the composition of mzLib's own methods (digestion, fragmentation, PRIDE,
  `MakeIdentifications` + `FlashLfqEngine`);
- the availability-vs-correctness **error classification** (timeouts/sockets/408/429/5xx →
  `ServiceUnavailable`, everything else by .NET exception type);
- keeping engine chatter off stdout so the envelope stays clean;
- carrying its own .NET runtime, so no .NET install is needed.

**mzLibR does none of that again.** It is a thin, idiomatic R package that spawns the bridge, writes
stdin, reads one JSON line, and parses it. The port is: transport file + typed constructors + a
condition hierarchy + tests + roxygen docs.

**Evidence this is true, not aspirational:** mzLibRust reached full test parity with pyMzLib in one
session — 137 tests against pyMzLib's 123 — and found seven upstream defects on the way. Nothing in
the wire had to change to accommodate a second language. Expect the same here.

## 2. The wire contract — exact

Protocol version **1**. Every invocation writes exactly one JSON object to stdout and nothing else;
diagnostics go to stderr. Exit 0 success, 1 handled failure, 2 usage.

```json
{"ok": true,  "data": { … }}
{"ok": false, "error": {"type": "…", "message": "…"}}
```

`error.type` is `"usage"`, `"ServiceUnavailable"`, or a .NET exception type name.

### The five verbs

| verb | arguments | stdin |
|---|---|---|
| `version` | — | — |
| `pride files` | `--accession` `--page-size` | — |
| `pride download` | `--accession` `--dest` `[--category]` `[--ext]` `[--no-overwrite]` `[--names-from-stdin]` | file names, one per line, when `--names-from-stdin` |
| `peptidoform fragments` | `--accession --protease --dissociation --terminus --missed-cleavages --min-length --max-length --max-mods --max-isoforms [--no-modifications]` | — |
| `quant flashlfq` | `--psms [--normalize] [--ppm] [--isotope-ppm] [--integrate] [--mbr] [--mbr-ppm] [--mbr-q] [--shared-peptides] [--bayesian] [--use-pep-q] [--threads] [--out]` | `path[⇥condition[⇥biorep[⇥techrep[⇥fraction]]]]`, one run per line |

**Large or variadic input goes on stdin, never argv** — argv has a ~32 KB ceiling and a real
experiment blows past it. `--max-length 0` means unbounded. Field names are the snake_case of their
mzLib names; **do not rename them on the R side either**.

## 3. Architecture map

| pyMzLib | mzLibRust | mzLibR |
|---|---|---|
| `_bridge.py` | `src/bridge.rs` | `R/bridge.R` — the ONLY file that knows the bridge exists |
| exception hierarchy | one `MzLibError` enum | **condition hierarchy**: `mzlib_error` → `mzlib_usage_error`, `mzlib_service_unavailable`, `mzlib_bridge_error`, `mzlib_timeout`, `mzlib_bridge_not_found`, `mzlib_protocol_error`, `mzlib_project_not_found` |
| `PYMZLIB_BRIDGE` | `MZLIB_BRIDGE` | `MZLIB_BRIDGE` **and** `options(mzlibr.bridge=)` |
| dataclasses | serde structs | S3 classes + **data.frames** (§4) |
| `pymzlib.pride.list_files` | `mzlib::pride::list_files` | `pride_list_files()` |

## 4. R-specific design — and the four R-shaped landmines

### 4.1 `jsonlite` will silently mangle the wire unless you stop it

**This is the single biggest R-specific hazard and it has no analogue in Python or Rust.**

- `jsonlite::fromJSON()` auto-simplifies. A JSON array of objects becomes a data.frame, an array of
  scalars becomes a vector, and **the shape depends on the data** — one file in a manifest gives a
  different structure from two. Parse with
  `fromJSON(txt, simplifyVector = FALSE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)` and
  build the data.frames yourself, deliberately.
- **JSON `null` becomes `NULL` inside an R list, and `NULL` elements vanish when a list is
  flattened.** This is exactly the protein-intensity trap (§6.4) with an extra R twist: the value
  does not become wrong, it *disappears*, and the vector silently comes back short. Map `null` →
  `NA_real_` at the boundary, explicitly, per field. Write a test with a fixture containing a null.
- Large integers: `file_size_bytes` can exceed 2^31. `jsonlite` may hand back a double; that is fine
  up to 2^53 but do not let it become `integer`.

### 4.2 R's `NA` is the *best* of the three languages for the None/0 lesson

A **peptide** intensity is `0` when missing; a **protein** intensity is "could not be resolved".
Python has to warn about `None` in prose; Rust makes it `Option<f64>`. R has `NA_real_`, which is
better than both, because **arithmetic propagates it**: `mean(x)` on a protein column returns `NA`
rather than a confidently wrong number, and `mean(x, na.rm = TRUE)` is an explicit choice the analyst
makes visibly.

```r
peptides$intensity   # numeric, 0 means "not measured here"
proteins$intensity   # numeric with NA_real_, NA means "FlashLFQ could not resolve one"
```

Say this in the docs *and* make it true in the types. Do not `na.rm = TRUE` anywhere inside the
package on the user's behalf.

### 4.3 Return data.frames, not lists of S3 objects

R users expect rectangles they can pipe into `dplyr`/`ggplot2`. Return **tidy data.frames** for
`files`, `peptides`, `peaks`, `proteins`, `fragments`, and a small S3 object with a `print()` method
for the result envelope that carries the scalars and the frames. One long frame beats a nested list
in every downstream use. Keep `tibble` optional, not a hard dependency.

### 4.4 Dependencies

`processx` (subprocess with real timeout support — `system2()` cannot time out portably) and
`jsonlite`. That is the floor. Everything else is a judgement call; nothing heavy. R's culture
tolerates dependencies better than Rust's, but every one is a CRAN-check liability.

## 5. Distribution — harder in R than in either sibling, and it is the real work

Python solved ~130 MB with a wheel that carries the payload. Rust solved it with `build.rs`
downloading a checksum-verified binary. **R has neither escape hatch:** CRAN's package size limit is
~5 MB, and **CRAN policy forbids writing outside `tempdir()` without explicit user consent** and
forbids downloading during install or `.onLoad`.

So the design is forced, and it is *different from the other two*:

1. **`mzlibr_install_bridge()`** — an exported, user-invoked function that downloads a
   checksum-verified bridge into `tools::R_user_dir("mzlibr", "cache")`. Interactive consent when
   `interactive()`, an explicit argument otherwise. Never called automatically.
2. **Resolution order at runtime:** `options(mzlibr.bridge=)` → `Sys.getenv("MZLIB_BRIDGE")` →
   the cache directory → an informative error naming all three remedies plus
   `mzlibr_install_bridge()`.
3. **Everything must work without a bridge.** `R CMD check`, all offline tests, all examples, all
   vignettes. Guard live code with `skip_on_cran()` and `skip_if(!bridge_available())`. This is not
   optional politeness — CRAN will run your checks on a machine that has no bridge and no network.

mzLibRust's `build.rs` and `scripts/stage-bridge.ps1` are the reference for the download-and-verify
logic; the *policy* wrapper around it is R-specific.

Two payload-shrink levers already found: mzLib **#1103** (TorchSharp/libtorch is ~238 MB, dragged in
transitively) and the mzML-only native-reader prune (~20 MB). Both matter more here than anywhere.

## 6. Every trap, with numbers — implement these correctly from line one

pyMzLib discovered these by shipping and running bake-offs; mzLibRust re-confirmed them and found
more. **You start knowing all of them. There is no excuse for a mzLibR user hitting any of these.**

### 6.1 PRIDE

| trap | the number |
|---|---|
| The API manifest is **incomplete**. PXD000001: API says **8** files, the FTP tree holds **13**. The five omitted include the two largest — a 450 MB mzML and a 472 MB mzXML | true project size **1.44 GB**, not the 0.514 GB the API implies |
| `file_size_bytes` is the **decompressed** size for some compressed files | the MGF reports 16,448,103 and downloads 5,984,662 — **2.75×** |
| A compressed file's extension is `.gz`, **not what it is compressed from**. `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz` | filtering on `.mgf` matches **zero files and exits successfully** |
| An unknown accession returns **HTTP 200 with `[]`**, not 404 | a typo becomes "0 files, done" unless you raise |
| `fileCategory == "PEAK"` matches **2** files, not 1 — the MGF *and* a 243 MB mzXML | 40× the intended download |
| `publicFileLocations` order is **not stable** — the mztab lists FTP first, the MGF lists **Aspera** first | indexing `[[1]]` yields an unfetchable `prd_ascp@…` URL. mzLib's `TryGetHttpsDownloadUrl` searches; do not re-implement |

Raise a condition when a filter matches nothing, and put the `.gz` explanation *in the error text*.

### 6.2 Peptidoforms — digestion

| trap | the number |
|---|---|
| mzLib's `trypsin\|P` **applies** the Keil rule; plain `trypsin` does not. This is the **reverse** of the MaxQuant/Mascot convention | albumin: **195** vs **202** peptides. mzLib#1106 |
| `min_length = 7` silently discards shorter peptides | albumin **195 → 254** at min_length 1; roughly a third of a histone digest |
| `peptides` are **peptidoforms**, not distinct sequences | albumin at 2 mods: **303** peptidoforms over **195** sequences. Provide both, name them differently |
| The isoform cap (default 1024) **truncates silently** | H3.1 at 4 mods loses ~30%. Surface `peptides_at_isoform_cap` and a `truncated` accessor |
| `--no-modifications` also discards **proteolysis products**, so the peptide list changes | albumin loses 2 signal-peptide peptides. pyMzLib#8 |

### 6.3 Peptidoforms — fragmentation

| trap | the number |
|---|---|
| **ETD/ECD emit `y` ions with no `b` ions.** No fragmentation mechanism produces that | ~**⅓** of every ETD fragment list. mzLib#1109, **PR #1114 open** — check whether it merged before documenting |
| The `zDot` series runs **1..L**, not 1..L−1. The extra one is `M − NH₂`, the N–Cα cleavage at residue 1 — **correct and deliberate**, not a bug | 1 per peptide; absent when the peptide starts with proline |
| **z• is suppressed N-terminal to proline; the complementary c ion is not** | albumin: **138** proline sites, 138 suppressed z•, **0** suppressed c. mzLib#1110 |
| Fragments expose `neutral_mass`, **not m/z**, and must not gain an `mz()` — a c or z ion carries only the fixed charges within its own span | per-fragment charge accounting does not exist |
| Peptide `mz()` must use the **proton** mass `1.00727646677`, not hydrogen `1.007825`, and must not double-count `fixed_charges` | 1.1 ppm at m/z 500; half a Thomson on a 2+ trimethylated peptide |

### 6.4 Peptidoforms — the modification census

Albumin P02768 annotates **38** modification-like features. mzLib loads **14**. The other 24 are
dropped **on feature type** (`ProteinXmlEntry.ParseFeatureEndElement` handles only
`modified residue` and `lipid moiety-binding region`), and **nothing reports it**.

**The exclusion is correct.** Do not "fix" it. 22 of the 24 are `N-linked (Glc) (glycation) lysine`,
which UniProt's `ptmlist.txt` *does* give a formula (`C6H10O5`) and mass (`162.052823`) — but
glycation is a labile, heterogeneous adduct that progresses to AGEs and dissociates in preference to
the backbone, so an exact mass plus a clean fragment ladder describes an unobservable species.
Furthermore **14 of the 22 are annotated `; in vitro`** and both `N-linked (GlcNAc...)` sites exist
only in the Redhill and Casebrook variants. mzLib reads none of those qualifiers, for any feature
type. mzLib#1112 — **the defect is the silence, not the exclusion.**

**Also: `annotated_modification_sites` ≠ `annotated_modifications_loaded`.** A histone carries K9me1,
K9me2, K9me3 and K9ac at one residue — four modifications, one site. Conflating them made H3.1 look
as though 93 annotations had been dropped when all had loaded.

### 6.5 FlashLFQ — the most valuable tranche and the most dangerous

Ground truth on mzLib's own K562 pair (`AllPSMs.psmtsv` +
`20100614_Velos1_TaGe_SA_K562_{3,4}.mzML`): 594 identifications → **354** peptides, **943** protein
groups, **647** peaks.

| trap | the number |
|---|---|
| **The peptide roll-up drops most MBR transfers.** Read `peaks`, not `peptides` | **140** true transfers; the peptide table shows **52**. A 63% under-count |
| **A whole run's transfers can vanish.** Per run, from peaks: run_3 **62**, run_4 **78**. From the roll-up: run_3 **0**, run_4 52 | MBR appears not to have worked at all in half the experiment |
| "Peptides quantified in both runs" | **257** from peaks; **169** from the peptide table |
| **Protein intensity can be unresolvable** (median-polish NaN → `null` → `NA_real_`) | **2** proteins. But **847** are `0` in both runs, mostly because `use_shared_peptides_for_protein_quant` defaults FALSE. "No usable number" is **849**; "could not be resolved" is **2** |
| `mbr_rescued_peptide_count` is *distinct sequences among MBR peaks* — the prose "rescued" definition differs | **140** vs **135** strict; the 5 have both an MBR peak and a zero-intensity MSMS peak in the same run |
| **`max_threads = -1` makes results non-reproducible.** The roll-up nondeterministically drops MBR intensities | a borderline protein was unresolvable in **5 of 6** runs. **Default to 1 in mzLibR, or warn loudly.** mzLib#1111 |
| The peptide roll-up does **not sum** multiple peaks in a run; it reports one | pivoting `peaks` yourself will not reproduce `QuantifiedPeptides.tsv` intensities |
| MBR needs a **complete, balanced design** and `mbr_q_value_threshold` as its FDR control | ~80% false transfers without it |
| mzML only | reject `.raw`/`.d` up front with a clear message |

### 6.6 Post-digestion modifications produce impossible peptides

mzLib digests first and applies modifications after, without rechecking the cleavage site. Albumin
yields `EFNAETFTFHADICTLSEK` ending in a **succinylated lysine at zero missed cleavages** — trypsin
cannot cleave after a succinyl-K, the charge is neutralised. **17 of 303** peptidoforms are
chemically impossible; far worse on histones. mzLib#1113.

## 7. Known upstream defects — state, as of 2026-07-23

| issue | what | status |
|---|---|---|
| [mzLib#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109) | ETD/ECD emit `y` with no `b` | **PR #1114 open** — 5334 tests pass |
| [mzLib#1110](https://github.com/smith-chem-wisc/mzLib/issues/1110) | z• proline-suppressed, c ions not | open |
| [mzLib#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111) | FlashLFQ roll-up nondeterministic | open |
| [mzLib#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112) | Glycosylation annotations dropped silently | open; **read the correction comment** |
| [mzLib#1113](https://github.com/smith-chem-wisc/mzLib/issues/1113) | Mods applied after digestion | open |
| [mzLib#1106](https://github.com/smith-chem-wisc/mzLib/issues/1106) | `trypsin` / `trypsin\|P` naming inversion | open |
| [mzLib#1103](https://github.com/smith-chem-wisc/mzLib/issues/1103) | TorchSharp bloats the payload ~238 MB | open |
| [mzLib#1108](https://github.com/smith-chem-wisc/mzLib/issues/1108) | Duplicate peptidoforms at chain/initiator-Met boundary | worked around **in the bridge** |
| [pyMzLib#7](https://github.com/smith-chem-wisc/pyMzLib/issues/7) | `intensity()` returned `None` against its invariant | **PR #9 open** |
| [pyMzLib#8](https://github.com/smith-chem-wisc/pyMzLib/issues/8) | `--no-modifications` drops proteolysis products | open |

**Check each before writing docs** — some may have merged, and a doc describing a fixed bug is its
own kind of wrong.

## 8. Documentation doctrine — the three lessons that cost the most to learn

### 8.1 Put the warning where the mistake is made, not where the concept is defined

Confirmed three times. The `.gz` trap was documented on `PrideFile::extension` — but the person about
to make the mistake is reading the `extensions` filter argument, whose examples were both
uncompressed. A biologist walked straight in. Conversely, the MBR warning sits on the
`intensities` field itself, and a biologist reported **257 instead of 169** and named that doc
comment as the reason.

In R this means: on the **argument** in `@param`, not only in `@details`; and in the **error
message**, which is the only documentation a stuck user reliably reads.

### 8.2 A disclosure that states a false reason is worse than no disclosure

`ModificationCensus$explain()` was wrong **twice in one day** — first claiming the excluded sites had
no defined composition (false for 22 of 24), then implying they should therefore be loaded (also
false, §6.4). Both times, the feature built specifically to prevent silent wrongness *was* the wrong
thing. If you cannot state the reason accurately, state the fact and say the reason is unavailable.

### 8.3 Numbers in documentation rot, and readers quote them

`"the two differ by 37 peptides out of about 200"` was stale by 5× — the real figure is 7. It had
been copied into a second binding before anyone measured it. **Every number in a doc string should be
reproducible by a test, or removed.**

### 8.4 Corollaries worth writing down

- `fragment_count()`-style conveniences invite wrong answers when the underlying data has categories
  that matter. Provide `fragments_by_series()` alongside, and say which to use.
- Say what a function *actually computes*, in code terms, when prose and implementation could
  diverge. `mbr_rescued_peptide_count` needed exactly this.
- Do not describe an equality as reassurance unless it is one. "These are equal only when X" was
  false on real data and hid the discrepancy it was supposed to reveal.
- Calibrate the reader: `NA` is the *rare* protein outcome (2), `0` is the common one (847). Docs
  that sell `NA` hard imply the opposite.

## 9. Test parity

pyMzLib has **123** tests, mzLibRust **137**. Target ≥123 for mzLibR, mapped one-for-one, with a
`docs/test-parity.md` recording every test that maps, every one eliminated by R's semantics, and
every one deliberately not ported **with the reason**. Nothing silently dropped.

- **Reuse the fixtures verbatim**: `pride_PXD000001_files.json`, `peptidoform_P02768_small.json`,
  `flashlfq_small.json`. Note `flashlfq_small.json` contains a `null` peptide intensity — that is the
  §4.1 landmine and must have a test.
- Split each module into **pure functions the tests can call directly**: `build_args()` for argument
  assembly and validation, `parse_*()` for wire → data.frame. Then the offline suite needs no
  subprocess and no mocking framework.
- For `bridge.R`, inject the runner (a function argument defaulting to `processx::run`) so the
  process-level failures are testable: silent death, non-JSON output, timeout, unlaunchable binary.
- Live tests **skip, never fail**, on `mzlib_service_unavailable`. `cargo test` has no skip verdict;
  **testthat does** (`skip()`), so mzLibR can do this properly where mzLibRust could only print.
- Add the **cross-binding golden test**: the same fixture parsed by pyMzLib, mzLibRust and mzLibR
  must agree field for field. Three bindings over one bridge is only a selling point if that holds.

## 10. The bake-off — run it, and run it correctly

This is what found five of the seven defects. **The methodology is not optional and the persona is
the experiment.**

### 10.1 The biologist phenotype (verbatim, inherited)

> A proteomics biologist: ~8 years at the bench and on an Orbitrap, comfortable with label-free
> quantification and **knows what match-between-runs is in principle**. **Not an experienced coder or
> R user** — leans on docs and `?help`, gets stuck on language mechanics, **does not read package
> source**.

mzLibRust's first attempt used "a researcher who writes some Rust". Those agents read library source,
reverse-engineered a dependency's internals, and wrote 375–598 lines of verification scaffolding
each — and **the comparison had to be thrown away**, because a careful expert self-corrects
regardless of documentation quality. It measured *can an expert reach the right answer* instead of
*do the docs carry a biologist through*.

**Run both, in this order: expert arms are excellent defect scanners; biologist arms are the only
valid comparison.** The expert run found six of the seven defects; the biologist run proved the docs.

### 10.2 Instrumentation (mandatory)

One JSON object per attempt, appended **as you go**, to `<task>_arm{A,B}_log.jsonl`:

```json
{"n": 1, "action": "what you tried", "outcome": "worked|deadend|external_lookup", "note": "short"}
```

`external_lookup` — leaving the tool for something it should have supplied — is the cleanest single
signal. Also required: **"would you put it in a figure?"**, dead-end count, and the decisions the arm
had to make alone.

### 10.3 Arms and blinding

Six agents: three tasks × two toolchains. **Arm A is never told mzLibR exists.** Arm B is forbidden
from reading package source, tests, `docs/findings.md` and `design/`.

Arm A's toolchain for R: search CRAN/Bioconductor honestly first. Expect to find `MSnbase`, `mzR`,
`Spectra`, `protViz`, `rpx` (which *does* do PRIDE, unlike Rust — so the PRIDE arm will be a real
contest here, and that is worth knowing). There is no R equivalent of FlashLFQ's MBR.

### 10.4 Ground truth must be independent

Derive it by driving the **bridge executable directly from the shell**, never through the package.
If the package is wrong, ground truth taken through it is wrong in the same direction and the
bake-off scores itself. Two of nine answers in the mzLibRust run turned out to be wrong in ground
truth's own source, and that only became visible because a second toolchain disagreed.

### 10.5 Agreement between arms is evidence, not proof

Both mzLibRust arms independently concluded the 22 glycation sites were searchable. Both were wrong.
They made the same reasonable inference from the same file. **A domain expert caught it, not the
tooling.** Budget for a human read of any surprising scientific conclusion.

### 10.6 Result to beat

| | ecosystem | mzLibRust |
|---|---|---|
| external lookups | 13 | **2** |
| dead ends | 8 | **1** |
| answers they would publish | 4 of 11 | **10 of 11** |

## 11. Standing rules

1. **An mzLib bug becomes an mzLib issue.** A workaround in a binding is not a fix — MetaMorpheus and
   every other consumer stays broken, and the next binding rediscovers it. Search existing issues,
   then file with the observation, the mechanism, a **quantified** impact, and a suggested fix.
   Reference the issue number at the workaround site, and write the test to assert what mzLib
   *currently does* so it fails when the fix lands. **Verify before filing** — one mzLibRust finding
   (the extra `z_L`) was investigated and correctly dropped as intended behaviour.
2. **A documentation lesson is back-ported to every binding.** The traps live in mzLib's behaviour,
   not in one language. Wording that rescues an R user rescues a Python one.

## 12. Licence

**LGPL-3.0-or-later**, matching mzLib — the package redistributes mzLib in compiled form. mzLibRust
initially declared MIT, which would have been a real violation. Ship `LICENSE`, `LICENSE.GPL-3.0` and
a `NOTICE` attributing mzLib, the .NET runtime (MIT), Intel `libmmd`, Newtonsoft.Json, and the
reference data (Unimod, PSI-MOD, UniProt). In `DESCRIPTION`: `License: LGPL-3` and a
`LICENSE.note`. State that the bridge override is a **licence affordance** under LGPL §4 — it is how
a user relinks a modified mzLib without rebuilding the package.

## 13. Milestones

- **M0** — `bridge.R` + `mzlibr_bridge_version()` round-trips against a real bridge and checks the
  protocol version. Proves transport and distribution end to end.
- **M1** — PRIDE tranche + tests + live canaries.
- **M2** — peptidoform tranche + tests + live canaries.
- **M3** — FlashLFQ tranche + tests, with `peaks`, `NA_real_` proteins, `mbr_rescued_peptide_count`,
  and `max_threads = 1` correct from the first line.
- **Then** — `mzlibr_install_bridge()`, `R CMD check --as-cran` clean on all three platforms, pkgdown
  site, the six-arm bake-off, cross-binding golden test, CRAN submission.

## 14. First steps

1. `git init` here; **build locally, create the public repo when green** — the first public commit
   should be a working package. Mirror the pyMzLib repo split.
2. Read, in order:
   - `E:\CodeReview\mzLibRust\` — **the closest model**: `src/bridge.rs`, `src/pride.rs`,
     `src/peptidoform.rs`, `src/flashlfq.rs`, `docs/findings.md`, `docs/test-parity.md`,
     `design/bakeoff/DESIGN.md`, `design/bakeoff/RESULTS.md`, `STATUS.md`, `build.rs`
   - `E:\CodeReview\pyMzLib\code\pyMzLib\pkg\bridge\` — `Program.cs`, `Peptidoform.cs`,
     `Quantification.cs` (the wire contract itself)
   - `E:\CodeReview\pyMzLib\code\pyMzLib\docs\contributing\conventions.md` — §1–§10
   - `E:\CodeReview\pyMzLib\design\bakeoff-flashlfq\DESIGN.md` — the original methodology
3. Stage a bridge:
   `E:\CodeReview\pyMzLib\code\pyMzLib\pkg\python\src\pymzlib\_dotnet\win-x64\mzlib-bridge.exe`
   (or `pkg/build/publish-bridge.ps1` for another RID).
4. M0, then M1–M3. Copy the three fixtures verbatim.
5. Run the bake-off — expert arms first, then biologist arms.

The user writes no R (as with Python and Rust). **Claude owns all R authorship, idiom, packaging and
CRAN compliance.** Explain trade-offs in terms the user knows — C#, mzLib, mass spectrometry — the
way `PYTHON_PRIMER.md` did for Python.

---

## Appendix: the one-paragraph version

The bridge does the hard part; you are writing a transport file, some data.frame constructors, a
condition hierarchy, and tests. Watch `jsonlite`'s auto-simplification and its `null`→`NULL`
disappearing act. Use `NA_real_` for protein intensities and `0` for peptides. Read `peaks`, never
the peptide roll-up, for MBR. Default `max_threads` to 1. Put every warning on the argument that
causes the mistake and inside the error message. Never state a reason you have not verified. File
mzLib bugs upstream. Run the bake-off with biologists, not programmers.
