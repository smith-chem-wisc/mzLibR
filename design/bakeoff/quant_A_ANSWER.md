# Label-free quant with MBR across K562 runs 3 & 4 — R arm (arm A)

**Tools I actually had:** R 4.6.1 with OrgMassSpecR + jsonlite. I found and installed
**RaMS** (pure-R mzML reader, CRAN binary) to read the raw data. I could NOT install
mzR / MSnbase / Spectra / xcms (Bioconductor, C/C++ code, no compiler on this machine)
or MSstats. So there was **no installable R tool that does label-free quant with MBR**.
RaMS reads spectra but does zero quantification; MSstats only does statistics on a
feature table produced by an upstream tool (MaxQuant/OpenMS/Skyline) and does no MBR.
Everything below the ID counts is a **crude XIC engine I hand-built** — not a validated tool.

Data note: the two mzML files are **RT slices** (73 MS1 scans each, ~79–83 min, the
*same* window in both runs), not full gradients. That happens to make RT alignment ~1:1.

---

## Q1 — How many peptides are quantified in BOTH runs?
**~121–127 peptides** (127 with any signal, 121 at my detection threshold, 111 at a stricter one).
This is the set of modification-aware peptides that were **identified in both runs** and have
real precursor signal in both (263/264 run3 and 215/215 run4 IDed peptides carry clear own-run
signal, so essentially every co-identified peptide is quantified in both).

- **Paper-ready? NO.** The number is stable, but "quantified" here = apex of per-scan-summed
  intensity in a fixed 10 ppm / ±0.5 min box, not a real integrated peak area with detected
  elution boundaries and isotope-envelope confirmation. It also excludes any peptide quantified
  in a second run *only* via MBR (see Q2 — I couldn't validate those).

## Q2 — How many peptides are rescued by match-between-runs?
**COULD NOT ANSWER.**
There are **225 candidates** (137 peptides IDed only in run3 + 88 only in run4). My naive MBR
finds *some* signal for 104/137 and 76/88 of them in the other run — **but a decoy-m/z control
(m/z shifted +0.5 Th) lands almost as often as the true m/z** (e.g. 180 true vs 157 decoy hits
with no threshold; 109 vs 74 at my working threshold — a ~60–90% false-match rate). In these
dense MS1 slices a 10 ppm + 0.5 min window nearly always contains *something*, so my method
cannot separate a genuinely transferred peak from background. Real MBR needs isotope-envelope
fitting, peak-shape scoring, RT alignment modeling, and an FDR estimate — none of which RaMS
provides and none of which I could responsibly hand-code. Any rescue count I quoted would be
fabricated.

- **Paper-ready? NO — no number at all.** This is the core thing R could not do for me.

## Q3 — How many protein groups end up with no usable intensity?
**2 protein groups** (`E9PCY5|Q02880`, `Q13547|Q5TEE2`) out of 299, at my working threshold
(1 group at no-threshold; the count jumps to 19 only when I set an aggressively high threshold).

- **Paper-ready? NO.** Two big caveats: (1) it rides on the same crude own-run quant as Q1;
  (2) a "protein group" here is just each unique pipe-delimited Protein Accession string from the
  PSM table — I did **no parsimony/protein inference**, so the group count itself (299) is not a
  real protein-grouping and would differ under proper inference.

---

## Every parameter I had to choose alone
| Choice | What I picked | Why / caveat |
|---|---|---|
| Confident-ID filter | `Decoy/Contaminant/Target=="T"` & `QValue<=0.01` | standard 1% FDR; all rows were target anyway |
| Peptide identity | `Full Sequence` (modification-aware) | a plain-sequence choice would merge charge/mod variants |
| Precursor m/z | (mono mass + z·1.007276)/z, mono mass & charge from best PSM (min QValue) | ignores that a peptide may elute at several charges |
| ppm tolerance | **10 ppm** | guessed from typical Orbitrap/FlashLFQ defaults |
| RT window | **±0.5 min**, centered on donor-run ID RT | only defensible because both slices span the same 79–83 min; assumes ~1:1 RT, no alignment model |
| Peak "intensity" | apex of per-scan **summed** intensity in the box | not an integrated area; no monoisotopic-envelope selection |
| Detection threshold ("quantified") | apex **> 50,000** | ~10th percentile of genuine IDed-peptide apices — arbitrary; Q1/Q3 move with it |
| "Rescued" (MBR) | signal at donor m/z + donor RT in the run where the peptide was not IDed | shown by decoy control to be non-specific here |
| Protein group | unique Protein Accession string | no parsimony inference |

## What I wished existed
An **installable R package that does FlashLFQ/IonQuant-style quant end-to-end** from a
`.psmtsv` + mzML: XIC extraction with isotope-envelope confirmation, proper peak integration,
RT alignment, and **FDR-controlled match-between-runs** — returning a peptide×run intensity
matrix and a protein rollup. RaMS gets raw data into R but stops there; MSstats starts after
the feature table already exists; the tools that bridge the gap (mzR/xcms) won't compile here.
The one tool on this machine that actually implements this is FlashLFQ itself (the mzLib engine
these test files ship with) — which is exactly what the R side is missing a binding to.

## Bottom line
- Q1 ≈ **121–127** (co-identified & co-detected; not a validated integrated-area quant).
- Q2 = **could not answer** (no specificity for MBR; decoy control fails).
- Q3 = **2** (crude, threshold-sensitive, no protein inference).
None are paper-ready. The ID-level counts are solid; everything requiring true peak
quantification or MBR is beyond what installable R tooling can deliver on this machine.
