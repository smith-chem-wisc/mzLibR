# mzLibR bake-off — results

Six biologist arms, three tasks × two toolchains, run 2026-07-24 against the ground truth in
`ground-truth.md`. Each arm was a blinded agent in the biologist persona (PLAN.md section 10.1):
Arm A never told mzLibR exists and given the CRAN/Bioconductor stack; Arm B given mzLibR and
forbidden to read its source, tests, docs or fixtures. Instrumentation per PLAN.md section 10.2.

**On the expert phase.** PLAN.md section 10.1 asks for expert defect-scanner arms before the
biologist comparison. That job was done here by the exhaustive independent ground-truth sweep,
which is what a careful expert does — and it found the defect an expert arm would have: mzLibRust's
bake-off ground-truth table quotes **254/269** for albumin where the truth is **243/257**, a number
its own arm log already had right (see `DESIGN.md`). So the six arms below are the valid
*comparison*, not a re-run of defect scanning.

## Headline

| | ecosystem (Arm A) | mzLibR (Arm B) | sibling's ecosystem | sibling's mzLibRust |
|---|---|---|---|---|
| external lookups | **9** | **2** | 13 | 2 |
| dead ends | **5** | **1** | 8 | 1 |
| answers they would put in a paper | **~2 of 9** | **9 of 9** | 4 of 11 | 10 of 11 |

mzLibR's Arm B lands on the sibling's mzLibRust almost exactly (2 external lookups, 1 dead end).
The ecosystem arm had *less* raw friction than the Rust ecosystem (9 vs 13) for one telling
reason — on the quant task it could not even get far enough to rack up lookups, because the tool
to attempt it does not installably exist.

## The answers, scored

Correct = matches ground truth and the arm can say why. Bold = the ground-truth number.

### Task 1 — PRIDE PXD000001

| | Arm A (rpx) | Arm B (mzLibR) | truth |
|---|---|---|---|
| file count | 11 — wrong, and no sizes at all | 8 (API) **and** knew the true 13 | 8 API / **13** real |
| MGF download | got it, but **not through rpx** — `pxget` said "No files to download"; pulled raw FTP | one call, **5,984,662** bytes, name-selected to dodge the 243 MB mzXML | **5,984,662** |
| typo accession | bare HTTP 404 warning, "easy to mistake for a network error" | clean typed error naming both causes | error, not `[]` |

Arm A's `rpx` exposes **no file sizes**, silently truncates its listing, and cannot see or fetch
the `generated/` MGF the task asked for — so the biologist left the tool for the raw FTP tree to
answer at all. Arm B answered every part inside the tool, and its help had pre-warned the 8-vs-13
gap, the 2.75× decompressed-size inflation, and the `category="PEAK"`-grabs-the-mzXML trap.

### Task 2 — albumin digestion and ETD fragmentation

| | Arm A (OrgMassSpecR) | Arm B (mzLibR) | truth |
|---|---|---|---|
| digest count | 83 at 0 missed cleavages, **no length filter** — "won't publish bare" | **195** distinct / 303 peptidoforms, all params stated | **195** (defaults) |
| modifications | 14 usable, but no name→mass map; needs a manual Unimod step | **38 / 14 / 24**, exclusion explained | **38 / 14 / 24** |
| ETD series | c and z (textbook) | c, zDot **and the spurious y** — caught the trap | c/zDot/**y** |

The digest question separates the arms cleanly. `OrgMassSpecR`'s `Digest()` applies no
length filter and no missed-cleavage model, so the biologist got a convention-dependent 83 they
refused to publish. Arm B got 195 *and* the surrounding sensitivity (38 at zero missed cleavages,
243 at min-length 1, 202 for plain trypsin) — and, crucially, **noticed the spurious ETD y series
that Arm A could not, because the ecosystem tool only ever returns the textbook c/z.** A trap you
cannot hit is not a trap you avoided; it is a question you were never allowed to ask.

### Task 3 — K562 quant, MBR on

| | Arm A (hand-rolled on RaMS) | Arm B (mzLibR) | truth |
|---|---|---|---|
| both runs | ~121–127, "not paper-ready" | **257** (trusted `peaks` over the roll-up's 169) | **257** |
| MBR rescued | **could not answer** — decoy m/z control non-specific; refused to invent | **140** (135 strict, difference explained) | **140** |
| proteins no intensity | 2 of 299 — wrong grouping, threshold-sensitive | **849** (847 zero + 2 NA), both readings given | 849 / **2** NA |

This is the arm with no tool. `mzR`/`MSnbase`/`xcms` need a compiler this machine lacks; the
biologist installed pure-R `RaMS` to read the mzML and hand-built an XIC engine. It reached
ID-level counts but on MBR it did the honest thing: a decoy-m/z control matched almost as often as
the true m/z, so it **refused to report a rescue count**. Its closing wish, verbatim: *"an
installable R package that does FlashLFQ/IonQuant-style quant end-to-end … which is exactly what
the R side is missing a binding to."* That is the product thesis, written by someone who did not
know Arm B existed.

## What did NOT happen, and why it matters

The sibling's sharpest result was an ecosystem arm **confidently wrong** — 284 peptides against
257, a plausible number reported without hedging. **No mzLibR Arm A answer was confidently wrong.**
Every ecosystem miss here was *hedged* — "won't publish bare," "not paper-ready," "could not
answer." That is not the ecosystem being safer; it is the ecosystem being so tool-starved on these
tasks that it could not manufacture a confident number to be wrong with. The danger the bindings
exist to prevent — a plausible wrong answer that reaches a figure — needs a tool capable of
producing the number in the first place. Arm B produced nine, and nine were right.

## Findings to act on (per the standing rules)

**No new mzLib defect.** The arms confirmed known ones — the ETD y series (#1109), the trypsin|P
naming inversion (#1106), the glycosylation exclusion (#1112) — and surfaced nothing new upstream.

**One documentation lesson, back-portable to every binding.** Both peptidoform arms flagged the
same thing: the live census labels all 24 excluded features `"glycosylation site"`, while the
`?census_explain` help text asserts "22 of 24 are N-linked (Glc) glycation lysine." The 22 is a
UniProt-name distinction the census output cannot show — the tool reports only the coarse
feature-*type* bucket. Asserting a finer number the tool cannot reproduce is the section 8.3 trap
in miniature. Fixed in mzLibR by framing the 22 as a UniProt-inspection detail, not a census
output; the same wording should be corrected in pyMzLib and mzLibRust.

**One capability gap, worth a binding-wide decision.** Arm B's single external lookup was real:
the tool knows the PRIDE manifest is incomplete (it says so) but offers no function to fetch the
true FTP file list or project size. A `pride_ftp_files()` / `pride_project_size()`, or even a
returned "files hidden from REST" attribute, would keep the user inside the tool. This is a wire
capability all three bindings lack, so it belongs upstream in the bridge if it is built at all.

**The sibling's 254/269.** Already recorded in `DESIGN.md`; per standing rule 2 it should be
corrected in mzLibRust's `design/bakeoff/DESIGN.md` too.
