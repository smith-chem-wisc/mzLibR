# HSA (UniProt P02768) — bake-off answers, arm A (R / OrgMassSpecR)

Tool: R 4.6.1, OrgMassSpecR `Digest()` and `FragmentPeptide()`. Sequence pulled fresh
from UniProt (rest.uniprot.org) because the package's built-in `example.sequence` is HSA
but the 585-aa mature chain, not the 609-aa canonical accession.

---

## Q1 — Tryptic digest: how many peptides?

**Answer: 83 fully-tryptic peptides** from the canonical 609-aa sequence at 0 missed
cleavages. Only **42** of those are ≥6 aa (37 are ≥7 aa).

Exact settings I used (all are judgement calls — see below):
- **Protease rule:** trypsin, cleave C-terminal to K and R **except when followed by P**
  (`enzyme = "trypsin"`, the OrgMassSpecR default; the alternative `"trypsin.strict"`
  cleaves after every K/R and would give a different count).
- **Minimum peptide length:** the tool applies **none** — it returns every fragment down
  to single residues. The "83" is unfiltered. If I impose a realistic detectability floor
  of ≥6 aa the count is **42**; at ≥7 aa it is **37**.
- **Missed cleavages:** reported at **0** (the default). For context, missed=1 and missed=2
  on the same sequence give more (missed=2 → 246 peptides). A real LFQ/MBR search would
  typically allow up to 2, so 83 is a floor, not what you'd actually search.
- **Which sequence:** canonical 609-aa P02768. The mature circulating chain (585 aa, signal
  peptide + propeptide removed) gives **79** at missed=0. The signal/propeptide peptides are
  not present in serum, so biologically the mature-chain number is the honest one; I report
  the canonical because that is literally "P02768".

**CONFIDENCE: In a paper? Not as a bare number.** The digest math is reliable, but "how many
peptides" is meaningless without stating rule + min length + missed cleavages + which
sequence form — and those choices swing the answer from 37 to 246. I would publish it only
as "N fully-tryptic peptides ≥7 aa, 0 missed cleavages, canonical sequence."

## Q2 — UniProt modifications: how many, and usable for a mass search?

UniProt P02768 has, in PTM-type feature categories:
- **14 "Modified residue"** — 9 phospho (pS/pT), 4 N6-succinyllysine, 1 N6-methyllysine.
- **24 "Glycosylation"** — 22 non-enzymatic glycation on Lys (most flagged "in vitro"),
  + 2 N-linked GlcNAc on Asn that exist **only in disease variants (Redhill, Casebrook)**,
  not in the canonical protein.
- (Also 17 disulfide bonds — structural cross-links, not a variable mass mod in a normal
  reduced/alkylated tryptic search. Not counted as searchable PTMs.)

**Usable for a mass search:**
- The **14 modified residues** are usable — each has defined chemistry, so a known delta
  (phospho +79.966, succinyl +100.016, methyl +14.016). These I'd add as variable mods.
- The glycations have a defined delta too (+162.053 hexose) but most are "in vitro"
  artefacts and heterogeneous — I would not put them in a routine search. The 2 real
  N-glycans are variant-only and don't apply to canonical HSA.
- **Big caveat:** UniProt gives a modification *name and site*, **not a machine-readable
  mass delta**. OrgMassSpecR does no name→delta mapping. To actually search any of these I
  had to supply the deltas myself (or go to Unimod). So "usable" = yes, but only after a
  manual mapping step the tools don't do for you.

**Answer I'd give: 14 well-defined modified residues that are directly usable as search
mass deltas; ~24 glycosylation annotations that are mostly in-vitro/variant and not useful
for a routine canonical search.**

**CONFIDENCE: In a paper? Yes for the counts and the phospho/succinyl/methyl usability
(straight off the UniProt record). No for treating the glycations as real in-vivo sites —
that needs primary evidence, not a UniProt "in vitro" tag.**

## Q3 — ETD fragmentation: which ion series?

**Answer: c ions and z ions** (N–Cα backbone cleavage), vs. the b/y ions you get from CID/HCD.
Confirmed in-tool: `FragmentPeptide(peptide, fragments = "cz")` returns the c- and z-ion
series (verified on the HSA tryptic peptide LVNEVTEFAK → [c1]1+, [c2]1+ … plus z ions).

**Caveat:** physically ETD yields the **z-radical (z•, ≈ z+1.00 Da)**; the tool labels its
product simply "z". Before matching real ETD spectra I'd confirm whether its "z" is classic
z or z•. c/a• and the complementary series can also appear but c and z are the diagnostic pair.

**CONFIDENCE: In a paper? Yes** — c/z for ETD is textbook and the tool agrees. The only thing
I'd double-check before a mass match is the z vs z• (z+1) definition.

---

## Judgement calls I made alone
1. **Sequence form:** reported canonical 609-aa but flagged that mature 585-aa (79 peptides)
   is the biologically-present form. The package's own HSA example is the mature one, unlabeled.
2. **Protease convention:** used trypsin *with* the no-cleavage-before-proline rule (default),
   not `trypsin.strict`.
3. **Minimum length:** the tool has none; I chose to report both the raw 83 and length-filtered
   42 (≥6) / 37 (≥7). "A peptide" here = any tryptic fragment; a detectable peptide = ≥6–7 aa.
4. **Missed cleavages:** anchored on 0 (default) but noted a real search uses up to 2 (→246).
5. **What counts as a "modification":** counted UniProt "Modified residue" + "Glycosylation"
   as PTMs; excluded disulfides (structural) and variant-only glycans from "usable".

## What I wished the docs / tool told me
- That `example.sequence` is HSA **as the mature chain** and which accession it maps to.
- A **minimum-length / detectability filter** option on `Digest()`, and a stated convention
  for what counts as a "peptide" — the headline count is dominated by this.
- Any bridge from **UniProt modification names to mass deltas** (a Unimod lookup). Right now
  the digest tool and the annotation live in two different worlds and I glue them by hand.
- Whether `FragmentPeptide`'s "z" ion is z or the ETD z-radical (z•).
