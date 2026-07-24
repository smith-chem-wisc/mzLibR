# Human serum albumin (UniProt P02768) — mzLibR digest & fragmentation

Tool: R package `mzLibR`, one call does everything:
`peptidoform_fragments("P02768", protease="trypsin|P", dissociation="ETD", missed_cleavages=2, min_length=7, ...)`
It fetches UniProt, applies annotated mods, digests, and fragments in a single result object.

---

## Q1 — Tryptic digest: how many peptides?

**195 distinct base sequences** (equivalently **303 peptidoforms**) at the tool defaults.

- **Protease rule:** `"trypsin|P"`. In mzLib naming the `|P` means the **Keil rule is applied** — it does NOT cleave before proline. This is the *reverse* of MaxQuant/Mascot, where `/P` means "do cleave before proline". Getting it wrong is quiet: plain `"trypsin"` gives **202** distinct sequences vs **195** for `"trypsin|P"` — I confirmed both numbers.
- **Minimum length:** 7 (default). This silently discards everything shorter. Dropping to `min_length=1` raises the count to **243 distinct** — about a fifth of the digest lives below the default floor.
- **Missed cleavages:** **They matter enormously.** Default is `missed_cleavages=2` → 195 distinct. With `missed_cleavages=0` (fully specific) the count collapses to **38 distinct sequences** ≥7 residues. So "how many tryptic peptides" is meaningless without stating the missed-cleavage allowance.

**"Peptides" is ambiguous and the gap is large, not a rounding error:**
- `nrow(digest$peptides)` = **303** = peptidoforms (one row per sequence-AND-modification-placement).
- `digest_distinct_base_sequences()` = **195** = distinct sequences.
I report distinct base sequences (195) as the headline "how many peptides," and flag 303 as the peptidoform count.

**CONFIDENCE: Paper-ready? YES**, provided every parameter is stated in the methods (protease `trypsin|P`/Keil, min_length 7, missed_cleavages 2) and I say whether I mean distinct sequences (195) or peptidoforms (303). The number is worthless without those qualifiers.

---

## Q2 — UniProt modifications: how many, and usable for a mass search?

**38 annotated; 14 applied (usable); 24 excluded.**

- The **14 applied** modifications (across 14 residue positions) are the ones usable for a mass search — mzLib loads only `modified residue` and `lipid moiety-binding region` annotations and gives each an exact monoisotopic mass.
- The **24 excluded** are dropped on feature *type* alone — the live tool labels all 24 "glycosylation site" (the help text specifies 22 of them are N-linked (Glc) glycation lysine). **This exclusion is correct** and should not be defeated: glycation/glycosylation are labile, heterogeneous adducts; assigning one exact mass + a clean fragment ladder would describe a species you cannot actually observe.

**Caveat that blocks blind use:** mzLib reads no qualifiers. Among the 14 *applied* sites, some UniProt annotations are marked "in vitro" and some exist only in disease variants (Redhill, Casebrook). The census gives a count only — it cannot tell those apart. So the 14 are usable as *masses*, but you must read the UniProt entry before trusting any specific site for a real sample.

**CONFIDENCE: Paper-ready? PARTLY.** The counts (14 applied / 24 excluded / 38 annotated) are solid and reproducible. But "usable for a mass search" needs a manual UniProt pass over the 14 applied sites to strip in-vitro / variant-only annotations before I'd put a site list in a paper. The count is publishable; a per-site usable list is not, yet.

---

## Q3 — Fragment a peptide with ETD: which ion series?

**Three series: `c`, `zDot`, AND `y`** — this surprised me.

Whole-digest tally (`digest_fragments_by_series`): `c`=5338, `y`=5338, `zDot`=5390.
Single peptide DLGEENFK (8-mer, unmodified): `c`=7, `y`=7, `zDot`=8.

**What the tool actually returns vs the textbook:**
- Textbook ETD gives **c and z• only**. mzLib additionally emits a **y series** (~one third of every ETD fragment list). No dissociation mechanism produces y without b, so these y ions are spurious — a known issue (mzLib #1109; fix proposed in PR #1114). **Do not trust a bare fragment total for ETD** — filter by `product_type` with `digest_fragments_by_series`.
- The `zDot` series runs `1..length` (8 for an 8-mer), one more than `c`/`y` which run `1..length-1`. That extra full-length z-dot (whole peptide minus NH2, the N-Cα cleavage at residue 1) is **correct and deliberate**, not a bug; it's absent when the peptide starts with proline.
- Across albumin, zDot counts come in *below* c overall because z-dot is suppressed N-terminal to proline (138 proline sites → 138 suppressed z-dot, 0 suppressed c).
- Fragments carry `neutral_mass` only — there is deliberately **no m/z** on fragment rows.

**CONFIDENCE: Paper-ready? YES for the mechanism, with a mandatory caveat.** I can state "ETD in mzLib returns c and z• ions" — but I must NOT quote a total ETD fragment count, because it silently folds in the spurious y series (~1/3) plus one extra z-dot per peptide. Any reported fragment count must be per-series and must exclude y.

---

## Judgement calls I made alone
1. Reported **195 distinct base sequences** as the primary "how many peptides," with 303 peptidoforms flagged alongside — because "peptides" is genuinely ambiguous here and the two differ by >50%.
2. Kept **all tool defaults** (Keil trypsin, missed=2, min_length=7) as the reference answer, then ran a small matrix (missed 0, min_length 1, plain trypsin) to show how sensitive each number is.
3. Treated the **24 glycosylation exclusions as correct** and did not attempt to force them into the search.
4. Called the **y series spurious** and refused to report any aggregate ETD fragment total.

## What I wished the docs said (they actually did — unusually well)
Honestly the help pages pre-empted almost every trap: the trypsin|P/Keil reversal, peptidoforms-vs-distinct, the min_length floor, the spurious ETD y series, the glycation exclusion, and the proton-vs-H mass in `peptide_mz`. The one thing I had to reconcile myself: the **live census labels all 24 excluded features "glycosylation site"** while the help text says "22 of 24 are N-linked (Glc) glycation lysine" — a small wording mismatch I'd want resolved before quoting the breakdown, though the count of 24 excluded / 14 applied is consistent between them.

## Surprises (stated explicitly)
- **ETD returned a y series.** I expected c and z• only. A third of the fragment list was y ions that no ETD mechanism should produce.
- **Missed cleavages swing the digest from 38 to 195** distinct sequences — a far bigger effect than I'd have guessed before running it.
- **zDot has one more ion than c** per peptide (the deliberate full-length N-Cα ion), so the series counts are not equal even before proline suppression.
