# mzLibR bake-off — design and ground truth

Methodology follows pyMzLib's `design/bakeoff-flashlfq/DESIGN.md` and mzLibRust's
`design/bakeoff/DESIGN.md`. Give an independent agent a real task, one toolchain, and no hints,
then measure whether it got a **correct** answer — not whether it got *an* answer.

The product claim for these bindings is not "faster" or "smaller". It is **"you are less likely
to be confidently wrong."** So trap-avoidance is the primary metric.

## Persona (all biologist arms) — inherited verbatim

**Getting this wrong invalidates the comparison.** From PLAN.md section 10.1:

> A proteomics biologist: ~8 years at the bench and on an Orbitrap, comfortable with label-free
> quantification and **knows what match-between-runs is in principle**. **Not an experienced
> coder or R user** — leans on docs and `?help`, gets stuck on language mechanics, **does not
> read package source**.

Run **both** phenotypes, in this order (PLAN.md section 10.1):

1. **Expert arms first** — a researcher who writes some R, allowed to read anything. Excellent
   **defect scanners**; the mzLibRust expert run found six of seven defects. But not a valid
   *comparison*: a careful expert self-corrects regardless of documentation quality, so an
   expert run that finds "nobody was confidently wrong" is an artifact and must be discarded.
2. **Biologist arms** — the phenotype above. The only valid comparison, because the
   documentation has to do the work.

## Blinding

- **Arm A is never told mzLibR exists.** Its toolchain is CRAN/Bioconductor: `rpx` (which *does*
  do PRIDE, so the PRIDE arm is a real contest), `MSnbase`/`mzR`/`Spectra`, `OrgMassSpecR` or
  `cleaver` for digestion. There is **no R equivalent of FlashLFQ's match-between-runs**, so the
  quant Arm A is expected to reach "could not answer" — itself a finding (PLAN.md section 10.6).
- **Arm B uses mzLibR** but is forbidden from reading package source, `tests/`, `docs/` and
  `design/`. A real user has none of the author's answers.
- Each arm works alone, in its own scratch directory, blind to the other's work.

## Instrumentation (mandatory, all arms)

One JSON object per attempt, appended as you go, to `<task>_arm{A,B}_log.jsonl`:

```json
{"n": 1, "action": "what you tried", "outcome": "worked|deadend|external_lookup", "note": "short"}
```

`external_lookup` — leaving the tool for something it should have supplied — is the cleanest
single signal. Also required per arm: **"would you put it in a figure?"**, the dead-end count,
and the decisions the arm had to make alone.

## The three tasks

| Task | Arm B (mzLibR) | Arm A (ecosystem) |
|---|---|---|
| Label-free quant across the two K562 runs, MBR on | `flashlfq_*` | hand-rolled on `mzR`/`Spectra`; no MBR exists |
| Digest and fragment serum albumin (P02768) | `peptidoform_*` | `OrgMassSpecR` / `cleaver` |
| PRIDE PXD000001 manifest and a filtered download | `pride_*` | `rpx` |

The arms are **not** told what the traps are. They are asked for numbers, their confidence, and
their surprises. Scoring is against the ground truth below.

---

## Ground truth

Established independently of the R package, by driving the bridge executable directly from the
shell (PLAN.md section 10.4). If the package were wrong, ground truth taken *through* it would be
wrong in the same direction and the bake-off would score itself.

Every number below was reproduced on the staged bridge (`bridge 1.0.0.0`, protocol 1, .NET
8.0.27) on 2026-07-24. The commands are in `ground-truth.md` beside this file.

### Task 1 — quant, the two K562 runs, MBR on

`AllPSMs.psmtsv` + `20100614_Velos1_TaGe_SA_K562_{3,4}.mzML` from mzLib's own FlashLFQ test data.
**594 identifications → 354 peptides, 943 protein groups, 647 peaks.**

| Question | Correct | The obvious-but-wrong route | 
|---|---|---|
| Peptides quantified in **both** runs | **257** (peaks with intensity > 0) | 169 (the peptide roll-up), or 266 (peaks *including* zero-intensity detections) |
| Peptides rescued by MBR | **140** (from `peaks`) | 52 (from the roll-up's `detection_type`) — a **63%** under-count |
| Protein groups with no obtainable intensity | **2** are NA (median-polish could not resolve them) | 0 (if NA read as zero) or 849 (if NA and 0 conflated) |

**The headline trap, verified.** MBR transfers per run:

| | run_3 | run_4 |
|---|---|---|
| from `peaks` (correct) | **62** | **78** |
| from the peptide roll-up | **0** | 52 |

An entire run's 62 transfers are absent from the roll-up. A user who builds their matrix from the
peptide table sees a result in which MBR appears not to have worked at all in half the experiment.

**Two subtleties this port surfaced, both my own first-pass errors:**

- **"Quantified" requires positive intensity.** Counting *any* peak gives 266; requiring
  `intensity > 0` gives 257. There are exactly **40** zero-intensity peaks — a detection is not a
  quantification, and conflating the two inflates the count by 9.
- **NA is counted per group, not per cell.** 2 protein groups are NA, which is 4 NA cells across
  the two runs. Reporting the cell count (4) for the group count (2) doubles it.

`no-usable-number` is **849** (all values NA or 0); `zero-in-both` is **847**; unresolvable is
**2**. Three different statements about 849 proteins, and the docs must not sell NA as the common
case — it is the rare one.

### Task 2 — digest and fragment serum albumin (P02768)

Distinct base sequences, exhaustive over both proteases × missed-cleavages 0–2 × min-length 1/7:

| protease | min_length | missed cleavages | distinct sequences | peptidoforms |
|---|---|---|---|---|
| `trypsin\|P` | 7 | 2 (default) | **195** | 303 |
| `trypsin` | 7 | 2 | **202** | 310 |
| `trypsin\|P` | 1 | 2 | **243** | 366 |
| `trypsin` | 1 | 2 | **257** | 381 |

- **The protease trap.** mzLib's `trypsin|P` *applies* the Keil rule; plain `trypsin` does not —
  the **reverse** of MaxQuant/Mascot. 195 vs 202. (smith-chem-wisc/mzLib#1106)
- **The min-length trap.** The default of 7 hides **48** distinct sequences here (195 → 243), a
  fifth of the digest.
- **The census.** UniProt annotates **38** modification-like features; **14** are applied; 24 are
  excluded on feature type, correctly. (smith-chem-wisc/mzLib#1112)
- **The ETD trap.** mzLib emits `c`, `zDot` **and `y`** for ETD; the `y` ions are spurious, about
  a third of the list. (smith-chem-wisc/mzLib#1109)

> **A correction the ground-truth phase forced, worth its own paragraph.** mzLibRust's
> `design/bakeoff/DESIGN.md` ground-truth table quotes **254** and **269** for the two
> `min_length = 1` rows. Neither reproduces at any setting — not any missed-cleavage count, not
> peptidoforms, not modifications off. The correct values are **243** and **257**, and
> mzLibRust's own `pep_armB_log.jsonl` *already recorded 243*: its ground-truth table disagrees
> with its own arm, and the arm was right. This is exactly the "numbers in documentation rot"
> lesson (PLAN.md section 8.3), sitting inside a sibling binding. It should be corrected there
> too — the standing rule is that a documentation lesson back-ports to every binding.

### Task 3 — PRIDE project PXD000001

**8 files, 514,278,049 bytes (0.514 GB)** by the API — the FTP tree holds 13, and the real
project is 1.44 GB.

- **The extension trap.** The MGF is `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`; its extension is
  **`.gz`**, not `.mgf`. Filtering on `.mgf` matches **zero** files and exits successfully.
- **The unknown-accession trap.** PRIDE answers a typo'd accession with HTTP 200 and `[]`, not a
  404, so a naive binding reports "0 files, done".
- **The location-order trap.** `publicFileLocations` is not stably ordered — the MGF lists Aspera
  first, so `locations[[1]]` yields an unfetchable `prd_ascp@` address for some files.

## Scoring

Per arm, per question: **correct** (and can say why), **confidently wrong** (a plausible number,
wrong, unhedged — the worst outcome), **hedged wrong**, or **could not answer** (legitimate for
Arm A, itself a finding).

## Result to beat (mzLibRust)

| | ecosystem | mzLibRust |
|---|---|---|
| external lookups | 13 | **2** |
| dead ends | 8 | **1** |
| answers they would publish | 4 of 11 | **10 of 11** |
