# Digest a protein and fragment its peptidoforms

Fetch a UniProt entry with its annotated modifications, digest it, and
return every peptidoform with its fragment ions, plus a census of which
annotations were applied.

Fetches a UniProt entry, applies its annotated modifications, digests it
with the named protease, and computes fragment ions for every resulting
peptidoform.

## Usage

``` r
peptidoform_fragments(accession, protease = "trypsin|P", dissociation = "ETD",
  modifications = TRUE, missed_cleavages = 2, min_length = 7, max_length = NULL,
  max_modifications = 2, max_isoforms = 1024, terminus = "Both", timeout = 300)
```

## Arguments

- accession:

  A UniProt accession, e.g. `"P02768"` (serum albumin).

- protease:

  The protease, in **mzLib's naming**.

  **Read this if you are coming from MaxQuant or Mascot.** mzLib's
  `"trypsin|P"` \*applies\* the Keil rule - it does not cleave before
  proline - while plain `"trypsin"` cleaves everywhere. That is the
  **reverse** of the MaxQuant and Mascot convention, where the `/P`
  suffix means "do cleave before proline". On albumin the two give
  **195** and **202** peptides, so the mistake is quiet and small enough
  to survive review (smith-chem-wisc/mzLib#1106).

- dissociation:

  The dissociation type, e.g. `"ETD"`, `"HCD"`, `"CID"`, `"ECD"`.

  \*\*`"ETD"` and `"ECD"` return the `c` and `zDot` series\*\* (radical
  N-Ca cleavage yields c/z-dot, not the b/y of vibrational activation).
  mzLib PR \#1114 removed the spurious `y` series ETD used to emit.

- modifications:

  Whether to apply UniProt's annotated modifications.

  `FALSE` gives the bare sequences - a clean control, since the digest's
  distinct backbones are unchanged and only the modified variants of
  them go away. The peptidoform count drops a long way; on albumin from
  **303** to **195**.

  This once carried a caveat saying `FALSE` also discarded **proteolysis
  products**, so that the peptide list itself changed and albumin lost
  two signal-peptide peptides (smith-chem-wisc/pyMzLib#8). That was true
  and is no longer: verified against the published bridge, both peptides
  are present either way and albumin gives 195 distinct base sequences
  with modifications on or off.

- missed_cleavages:

  Maximum missed cleavage sites per peptide.

- min_length:

  Shortest peptide to keep, in residues.

  The default of 7 **silently discards** everything shorter. Albumin
  goes from **195** distinct sequences at `min_length = 7` to **243** at
  `min_length = 1` - a fifth of the digest lives below the default. If
  you are looking for a short peptide and not finding it, look here
  first.

- max_length:

  Longest peptide to keep, in residues, or `NULL` for no limit.

- max_modifications:

  Maximum modifications considered per peptide.

- max_isoforms:

  Maximum modification isoforms (peptidoforms) generated per peptide.

  **This cap truncates silently.** A truncated result and a genuinely
  short one look identical from the outside - histone H3.1 at four
  modifications loses about **30%**. Check
  [`digest_truncated`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_truncated.md)
  before treating a peptide list as exhaustive.

- terminus:

  Which terminus to fragment: `"N"`, `"C"`, `"Both"` or `"None"`.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_digest`: a list with the scalars describing the run, a
`census` (see
[`census_explain`](https://smith-chem-wisc.github.io/mzLibR/reference/census_explain.md)),
and three data.frames that join on `peptide_index` - `peptides`,
`fragments` and `modifications`.

The scalars: `sequence_length`, the protein's length in residues;
`max_modifications`, the cap in modifications per peptidoform;
`max_modification_isoforms`, the cap in peptidoforms per peptide
position; and `peptides_at_isoform_cap`, the peptide positions that
reached it - non-zero means the list is truncated. The census counts
`sites` in residues (distinct positions carrying a modification),
`applied` in modifications and `annotated` in UniProt features.

`peptides` has `monoisotopic_mass`, the neutral mass in Da; `length` in
residues; `missed_cleavages` in cleavage sites; `modification_count` in
modifications; and `fixed_charges`, the formal charge the intact peptide
carries before protonation.

\*\*`peptides` holds peptidoforms, not distinct sequences. **One row per
sequence-and-modification-placement, so albumin at two modifications
is** 303 **rows over** 195\*\* distinct sequences. Both are legitimate
answers to "how many peptides" and they are not interchangeable; quoting
one for the other is a large error, not a rounding one. See
[`digest_distinct_base_sequences`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_distinct_base_sequences.md).

## Wraps

Wire verb `peptidoform fragments`. Generated from the bridge's verb spec
`peptidoform.fragments.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`ProteinDbLoader.LoadProteinXML`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

- [`Loaders.LoadUniprot`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/Loaders.cs)
  in `mzLib/UsefulProteomicsDatabases/Loaders.cs` at mzLib `23c2490e`

- [`Protein.Digest`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/Protein.cs)
  in `mzLib/Proteomics/Protein/Protein.cs` at mzLib `23c2490e`

- [`DigestionParams`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/ProteolyticDigestion/DigestionParams.cs)
  in `mzLib/Proteomics/ProteolyticDigestion/DigestionParams.cs` at mzLib
  `23c2490e`

- [`ProteaseDictionary`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/ProteolyticDigestion/ProteaseDictionary.cs)
  in `mzLib/Proteomics/ProteolyticDigestion/ProteaseDictionary.cs` at
  mzLib `23c2490e`

- [`PeptideWithSetModifications.Fragment`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/ProteolyticDigestion/PeptideWithSetModifications.cs)
  in
  `mzLib/Proteomics/ProteolyticDigestion/PeptideWithSetModifications.cs`
  at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `accession`:

  string; required. A UniProtKB accession, e.g. P02768; fetched from
  rest.uniprot.org. The bridge checks only that it is not blank (the
  bindings upper-case it and check its shape first). UniProt's 400 or
  404 becomes a usage error.

- `protease`:

  string; default `trypsin|P`; range
  `a key of mzLib's ProteaseDictionary`. mzLib protease name. The
  default "trypsin\|P" is the Keil rule (cleave after K/R EXCEPT before
  P); mzLib's plain "trypsin" cleaves before P too. The bracket/bar
  marks a PREVENTING residue, the reverse of MaxQuant's and Mascot's
  "Trypsin/P".

- `dissociation`:

  string; default `ETD`; range
  `a mzLib DissociationType name, case-insensitive: HCD, CID, ETD, EThcD, ...`.
  Which fragment series to generate: HCD/CID give b and y; ETD gives c
  and zDot.

- `terminus`:

  string; default `Both`; range
  `a mzLib FragmentationTerminus name, case-insensitive: Both, N, C`.
  Which fragment termini to generate; also passed to DigestionParams.

- `modifications` (wire `--no-modifications`):

  flag; default `FALSE`. Digest the bare sequence: UniProt's
  modifications are discarded but its signal-peptide and propeptide
  boundaries are kept (pyMzLib#8), so the backbone peptide list is
  unchanged and only the modified variants go away. A clean control.

- `missed_cleavages` (wire `--missed-cleavages`):

  int; in **cleavage sites**; default `2`; range `>= 0`. Maximum missed
  cleavage sites per peptide.

- `min_length` (wire `--min-length`):

  int; in **residues**; default `7`; range `>= 1`. Shortest peptide
  kept. The default silently drops shorter peptides (about a third of a
  histone digest); pass 1 for every peptide.

- `max_length` (wire `--max-length`):

  int; in **residues**; default `0`; range `>= 0`. Longest peptide kept;
  0 means unbounded.

- `max_modifications` (wire `--max-mods`):

  int; in **modifications**; default `2`; range `>= 0`. Maximum
  modifications per peptidoform (DigestionParams maxModsForPeptides).
  Isoforms grow combinatorially: histone H3.1 gives 49 bare tryptic
  peptides, 2,563 at 2 and 7,040 at 3.

- `max_isoforms` (wire `--max-isoforms`):

  int; in **peptidoforms**; default `1024`; range `>= 1`. Cap on
  modification isoforms per peptide position (DigestionParams
  maxModificationIsoforms). mzLib truncates SILENTLY at the cap;
  peptides_at_isoform_cap reports it.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `accession`:

  string; never `NA`. The accession UniProt returned.

- `name`:

  string; `NA` when the entry records no entry name. UniProt entry name,
  e.g. ALBU_HUMAN.

- `full_name`:

  string; `NA` when the entry records no recommended name. Recommended
  protein name, e.g. Albumin.

- `organism`:

  string; `NA` when the entry records no organism. Scientific name.

- `sequence_length`:

  int; in **residues**; never `NA`. Length of the full precursor
  sequence.

- `modifications_applied`:

  bool; never `NA`. False under no-modifications.

- `census$sites` (wire `annotated_modification_sites`):

  int; in **residues**; never `NA`. Distinct residue positions carrying
  at least one loaded modification. Not a modification count:
  K9me1/me2/me3/ac are four modifications at one site.

- `census$applied` (wire `annotated_modifications_loaded`):

  int; in **modifications**; never `NA`. Modifications mzLib loaded onto
  the protein (summed over sites). Counted from the annotated entry, so
  unchanged by no-modifications.

- `census$annotated` (wire `uniprot_annotated_features`):

  int; in **features**; never `NA`. Modification-like features in the
  XML: modified residue, glycosylation site, lipid moiety-binding region
  and cross-link.

- `census$unresolved` (wire `unresolved_modifications`):

  string\[\]; never `NA`. Modification names the entry annotates that
  could not be resolved against UniProt's ptmlist, sorted. These are
  dropped silently by mzLib otherwise.

- `census$by_type` (wire `uniprot_features_by_type`):

  object\[\]; never `NA`. {type, count, loaded} per feature type, sorted
  by type. loaded is true only for 'modified residue' and 'lipid
  moiety-binding region': mzLib drops the others on feature TYPE alone,
  before any mass lookup.

- `protease`:

  string; never `NA`. The protease used.

- `dissociation`:

  string; never `NA`. DissociationType name used.

- `terminus`:

  string; never `NA`. FragmentationTerminus name used.

- `max_modifications`:

  int; in **modifications**; never `NA`. max-mods, echoed.

- `max_modification_isoforms`:

  int; in **peptidoforms**; never `NA`. max-isoforms, echoed.

- `peptides_at_isoform_cap`:

  int; in **peptide positions**; never `NA`. Peptide positions (start,
  end) whose distinct peptidoforms reached max-isoforms. Non-zero means
  the list is TRUNCATED.

- `peptide_count`:

  int; in **peptidoforms**; never `NA`. Entries in peptides.

- `peptides`:

  object\[\]; never `NA`. One entry per distinct peptidoform (start,
  end, full sequence), in digestion order; fields under result.columns.

## Columns of `peptides`

- `base_sequence`:

  string; never `NA`. Bare residues.

- `full_sequence`:

  string; never `NA`. mzLib full sequence with modifications inline,
  e.g. ...K\[UniProt:N6-succinyllysine on K\].

- `monoisotopic_mass`:

  float; in **Da**; never `NA`. Neutral monoisotopic mass, modifications
  included.

- `length`:

  int; in **residues**; never `NA`. Peptide length.

- `one_based_start`:

  int; never `NA`. First residue's position in the protein, 1-based.

- `one_based_end`:

  int; never `NA`. Last residue's position in the protein, 1-based,
  inclusive.

- `missed_cleavages`:

  int; in **cleavage sites**; never `NA`. Missed cleavage sites spanned.

- `modification_count`:

  int; in **modifications**; never `NA`. Entries in modifications.

- `fixed_charges`:

  int; in **charge**; never `NA`. Formal charges the INTACT peptide
  carries before protonation (e.g. trimethyllysine's quaternary
  ammonium), recovered from each modification's formula-vs-recorded mass
  deficit in whole electrons. m/z at charge z is (monoisotopic_mass +
  (z - fixed_charges) \* 1.00727646677) / z. Does NOT apply to
  fragments.

- `modifications`:

  object\[\]; never `NA`. {one_based_residue, terminus, id, mass,
  formal_charge} per modification: one_based_residue is 1-based within
  the PEPTIDE and null for a terminal modification, which sets terminus
  'N' or 'C' instead (terminus is null for a residue modification); id
  is mzLib's IdWithMotif; mass is its monoisotopic mass in Da (null only
  if the modification defines none); formal_charge in elementary
  charges.

- `fragments`:

  object\[\]; never `NA`. {product_type, fragment_number, neutral_mass,
  neutral_loss, residue_position} per ion: product_type e.g. b, y, c,
  zDot; fragment_number the position in its series; neutral_mass the
  monoisotopic NEUTRAL mass in Da (not an m/z; no proton added);
  neutral_loss in Da (0 when none); residue_position 1-based within the
  peptide.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  accession missing or blank; UniProt answers 400 or 404 (no such
  entry); the entry holds no protein sequence; protease not in
  ProteaseDictionary; dissociation or terminus not an enum name; an
  integer option that is not an integer

- `mzlib_service_unavailable` (service_unavailable):

  UniProt unreachable by the bridge's own classification: a timeout
  (HttpClient 100 s), a socket or TLS failure, a request that never got
  a response, an HTTP 408, 429 or 5xx, or a body cut off in transit

- `mzlib_bridge_error` (correctness):

  any other UniProt status (e.g. 403); UniProt's ptmlist.txt or
  PSI-MOD.obo.xml missing from the bridge payload
  (FileNotFoundException: every protein would otherwise look unmodified
  or carry wrong charges); mzLib failing to parse the entry or to digest
  it

## Caveats

- UniProt annotates more than mzLib applies: only 'modified residue' and
  'lipid moiety-binding region' features load. Serum albumin: 14 of 38
  annotated features applied, the 24 glycosylation-site features dropped
  on type (mzLib#1112). Read uniprot_features_by_type and
  unresolved_modifications before trusting a modification count.

- max-isoforms truncates silently in mzLib (histone H3.1 at four
  modifications loses about 30% of peptidoforms at 1024). Check
  peptides_at_isoform_cap before treating the list as exhaustive.

- Duplicate peptidoforms that mzLib's Digest emits where a chain
  boundary meets the initiator-Met cleavage site (mzLib#1108, open) are
  collapsed by the bridge on (start, end, full sequence).

- Fragment masses are neutral. Converting one to m/z needs the fixed
  charge within that fragment's span, which the wire does not report;
  fixed_charges is whole-peptide only.

- ETD's zDot ions are suppressed N-terminal to proline while the
  complementary c ions are not, leaving about 4% of the c series
  unobservable (mzLib#1110).

- Needs the network: UniProt is queried on every call.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.peptidoform.fragments`

- Rust (mzLibRust): `mzlib::peptidoform::fragments_with` with
  `FragmentOptions`

- R (mzLibR): `peptidoform_fragments`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- peptides_at_isoform_cap is computed on the DE-DUPLICATED list,
  although FragmentsAsync's own comment says it must count the raw
  digest (a truncated locus deduped below the cap then reads as
  untruncated); the correctly computed peptidesAtCap local is never
  used. Likely a bridge bug; the doc above describes what the wire does.

- The example fixture is hand-trimmed (2 peptides of a max-mods 1
  digest), so its peptide_count is not what the live call returns; the
  envelope keys are the live ones.

- Integer options are not range-checked by the bridge (negative
  min-length, max-isoforms 0); the bindings check them. Enum.TryParse
  also accepts numeric strings, so `--dissociation` 3 or `--terminus` 7
  pass as undefined members. Behaviour of mzLib on such values is
  unverified.

- The wire names max-mods and no-modifications project as
  max_modifications and modifications (inverted) in every binding;
  pyMzLib's spec lint will need PYTHON_DEVIATIONS entries for them.

- modifications\[\].mass nullability: Modification.MonoisotopicMass is
  double?; a loaded UniProt modification has not been seen without one.

## References

- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, the source of the entry, its sequence and its annotated
  modifications (UniProt Consortium 2025, NAR)

## See also

[`digest_fragments_by_series`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_fragments_by_series.md),
[`digest_truncated`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_truncated.md),
[`peptide_mz`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md)

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
digest
#> <mzlibr_digest> P02768 ALBU_HUMAN (Homo sapiens)
#>   609 residues, trypsin|P, ETD, terminus Both
#>   min_length 7, 2 missed cleavages, max 1 modifications
#>   2 peptidoforms over 2 distinct sequences
#>   127 fragments: c=64, zDot=63
#>   ! 14 of 38 annotated modifications applied. See ?census_explain.
head(digest$peptides[, c("base_sequence", "monoisotopic_mass", "modification_count")])
#>                         base_sequence monoisotopic_mass modification_count
#> 1     ALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK          3562.853                  0
#> 2 RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK          4136.902                  1
digest_fragments_by_series(digest)
#>   product_type  n
#> 1            c 64
#> 2         zDot 63
```
