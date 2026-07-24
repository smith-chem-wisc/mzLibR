# Label-free quant across the two K562 runs, MBR ON — Answers (Arm B, R / mzLibR)

Data: `AllPSMs.psmtsv` + K562 run 3 and run 4 mzML.
Call: `flashlfq_quantify(psms, spectra, match_between_runs = TRUE, max_threads = 1)`
with `spectra` = both runs, `condition = "K562"`, `biological_replicate = 1 / 2`.
Identifications: 594. Protein groups: 943.

---

## Q1. How many peptides are quantified in BOTH runs?

**257** (distinct sequences with a positive intensity in each run, from the `peaks` table).

- From `peaks` (intensity > 0 in both runs): **257**
- From `peptides` roll-up (intensity > 0 in both runs): **169**

The two disagree by 88. **I trust 257 (peaks).** The help for `flashlfq_quantify`
explicitly warns the `peptides` roll-up drops most MBR transfers and undercounts,
"badly and unevenly." The detection_type tallies show it directly: `peaks` carries
140 MBR peaks; the `peptides` table carries only 52 MBR rows (plus 173 `NotDetected`
and 31 `MSMSIdentifiedButNotQuantified`). 169 is the undercount; 257 is the real answer.

**CONFIDENCE: paper-ready = yes.** Reproduces the documented K562 figure (257 vs 169)
exactly, and the peaks-vs-roll-up gap is a documented, understood behavior, not noise.

---

## Q2. How many peptides are rescued by match-between-runs?

**140** — `flashlfq_mbr_rescued_peptide_count(res)` = distinct modified sequences among
the MBR peaks. Cross-checks: `flashlfq_mbr_peak_count` = 140; per-run MBR peaks
run3 = 62, run4 = 78 (= 140); no peptide was rescued in more than one run, so peaks
and distinct-peptides coincide here.

Caveat the docs flag and I verified: 140 is "distinct sequences carrying an MBR peak."
The stricter prose reading — "peptides quantified in a run *only* by MBR" — is **135**,
because 5 of the MBR peaks sit in a run where that same peptide also has a
(zero-intensity) MSMS peak. **The tool's reported number is 140**; 135 is the answer to
a subtly different question. Report 140 unless your Methods define rescue as
"MBR-only in that run," in which case say 135.

Do NOT read MBR off the `peptides` table (it shows 52 — a 63% undercount).

**CONFIDENCE: paper-ready = yes.** Helper is purpose-built; both readings reproduced
exactly (140 / 135). State which definition you use.

---

## Q3. How many protein groups end up with no usable intensity?

**849** of 943 protein groups have no usable intensity (0 or NA in both runs).
Only **94** have a positive intensity in at least one run.

Breakdown of the 849:
- **847** are `0` in both runs — measured, no signal (mostly because their only evidence
  is shared peptides, and `use_shared_peptides_for_protein_quant` defaults to FALSE).
- **2** are `NA` in both runs — FlashLFQ's median-polish could not resolve a number
  (produced NaN). (Seen as 4 NA protein *rows* = 2 groups × 2 runs.)

Two defensible readings, and the docs name both:
- "No usable number" (0 or NA in both runs) = **849** ← what I report for this question.
- "Could not be resolved" (NA only) = **2**.

**I trust 849** as the answer to "no usable intensity," because a 0-in-both-runs protein
gives you nothing to quantify with regardless of why. If the question specifically means
"failed / unresolved," that's 2. `0` means "not measured," `NA` means "unresolvable";
mzLibR never silently applies `na.rm`.

**CONFIDENCE: paper-ready = yes**, with the caveat that 847 of these are an artifact of a
*parameter default*, not biology: `use_shared_peptides_for_protein_quant = FALSE`. If your
design intends shared peptides to count, re-run with it TRUE before quoting 849 in a paper.

---

## Parameters / choices I made alone
- `match_between_runs = TRUE` (as asked).
- `spectra` design: both runs `condition = "K562"`, `biological_replicate = 1` and `2`.
  The help says MBR "needs a complete, balanced design" and to set `condition` and
  `biological_replicate` so FlashLFQ knows which runs are comparable. Two runs of the same
  sample = one condition, two bioreps. This is a minimal balanced design (2 runs); with only
  two runs MBR FDR control is thin — fine for this bake-off, thinner than I'd want for a paper.
- `max_threads = 1` (the default). The help says >1 is nondeterministic and can drop MBR
  intensities ~1 run in 6 (mzLib#1111), so single-thread for reproducibility.
- Left at defaults: `normalize = FALSE`, `ppm_tolerance = 10`, `isotope_ppm_tolerance = 5`,
  `integrate = FALSE`, `mbr_ppm_tolerance = 10`, `mbr_q_value_threshold = 0.05`,
  `use_shared_peptides_for_protein_quant = FALSE`, `bayesian_protein_quant = FALSE`,
  `use_pep_q_value = FALSE`.
- Q1 "quantified" = positive intensity (MSMS or MBR both count as quantified).

## What the tool warned me about / what surprised me
- The `peptides` roll-up and the `peaks` table disagree on nearly everything (both-runs
  257 vs 169; MBR 140 vs 52). The help warns about this loudly, but it would be very easy
  to `read only $peptides` and silently publish "MBR did nothing in run 3" (roll-up shows
  run3 MBR = 0). **Always work from `$peaks`.**
- `0` vs `NA` in the protein table mean different things (not-measured vs unresolvable) and
  arithmetic propagates NA on purpose. Surprising but sensible.
- 847/943 protein groups being 0 is driven entirely by one parameter default
  (shared-peptides off), not the data.
- Everything I computed reproduced the numbers baked into the help text for this exact
  K562 pair — reassuring, and it's why confidence is high.

## What I wished the docs said
- A one-liner recipe: "For 'peptides in both runs' pivot `$peaks` on `file_name` and
  intersect `sequence` where intensity > 0" — I had to derive that.
- Whether MBR-transferred peaks should count as "quantified" for Q1 (I decided yes).
- An `mbr_q_value_threshold` sanity note for a 2-run design — with only two runs, is the
  0.05 FDR even meaningful? The help says MBR "needs a complete, balanced design" but a
  2-run experiment is the smallest possible; I'd like guidance on how much to trust the 140.
- `census_excluded(res)` errored on this result with no message; the help lists it but I
  couldn't use it. Would like to know when it applies.
