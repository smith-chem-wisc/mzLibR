# Digest a protein and fragment its peptidoforms

Fetch a UniProt entry with its annotated modifications, digest it, and
return every peptidoform with its fragment ions, plus a census of which
annotations were applied.

Fetches a UniProt entry, applies its annotated modifications, digests it
with the named protease, and computes fragment ions for every resulting
peptidoform.

## Usage

``` r
peptidoform_fragments(accession, protease = NULL, dissociation = NULL,
  modifications = TRUE, missed_cleavages = NULL, min_length = NULL, max_length = NULL,
  max_modifications = NULL, max_isoforms = NULL, terminus = NULL, timeout = 300)
```

## Arguments

- accession:

  A UniProt accession, e.g. `"P02768"` (serum albumin).

- protease:

  A protease name, as mzLib names it. `NULL` means mzLib's default
  protease.

  **Read this if you are coming from MaxQuant or Mascot.** The names and
  what each one cuts are mzLib's. mzLib's `"trypsin"` cuts after K and R
  even before proline, and `"trypsin|P"` does not cut before proline:
  the `P` names a residue that prevents the cut. That is the **reverse**
  of MaxQuant's and Mascot's `Trypsin/P`, which ignores the proline
  rule. On serum albumin the two give **202** and **195** distinct
  peptides, a small count hiding a large difference: which peptides you
  get changes wherever K or R precedes a proline
  (smith-chem-wisc/mzLib#1106).

- dissociation:

  `"HCD"` or `"CID"` (b and y ions), `"ETD"`, `"ECD"`, and the rest of
  mzLib's dissociation types. `NULL` means `"HCD"`: mzLib names no
  default here, so HCD is the bridge's choice until mzLib names one
  (smith-chem-wisc/mzLib#1488).

  `"ETD"` and `"ECD"` **return the c and z-dot series**, `c` and `zDot`
  (radical N-Ca cleavage yields c/z-dot, not the b/y of vibrational
  activation). mzLib PR \#1114 removed the spurious `y` series ETD used
  to emit.

- modifications:

  Whether to apply UniProt's annotated modifications.

  `FALSE` gives the bare sequences - a clean control, since the digest's
  distinct backbones are unchanged and only the modified variants of
  them go away. The peptidoform count drops a long way; on albumin from
  **310** to **202**.

  This once carried a caveat saying `FALSE` also discarded **proteolysis
  products**, so that the peptide list itself changed and albumin lost
  two signal-peptide peptides (smith-chem-wisc/pyMzLib#8). That was true
  and is no longer: both peptides are present either way, and albumin
  gives 202 distinct base sequences with modifications on or off.

- missed_cleavages:

  Maximum missed cleavage sites per peptide. `NULL` means mzLib's
  default.

- min_length:

  Shortest peptide to keep, in residues. `NULL` means mzLib's default.

  mzLib's default **silently discards** everything shorter. Albumin goes
  from **202** distinct sequences at the default to **257** at
  `min_length = 1` - a fifth of the digest lives below it - and histone
  H3.1 keeps 49 of the 96 peptides, by start and end, that it gives at
  `min_length = 1`. If you are looking for a short peptide and not
  finding it, look here first, and pass `min_length = 1` when you mean
  every peptide.

- max_length:

  Longest peptide to keep, in residues. `NULL` means mzLib's default,
  which is unbounded; `0` asks for no bound.

- max_modifications:

  Maximum modifications considered per peptide. `NULL` means mzLib's
  default. Modification isoforms are enumerated combinatorially: histone
  H3.1's 49 peptides give 2,419 peptidoforms at two modifications and
  6,580 at three.

- max_isoforms:

  Maximum modification isoforms (peptidoforms) generated per peptide
  position. `NULL` means mzLib's default.

  **This cap truncates silently.** A truncated result and a genuinely
  short one look identical from the outside - on histone H3.1 at four
  modifications the default cap keeps 8,644 of 12,808 peptidoforms,
  losing about a third. Check
  [`digest_truncated`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_truncated.md)
  before treating a peptide list as exhaustive.

- terminus:

  Which terminus to fragment: `"N"`, `"C"`, `"Both"` or `"None"`. `NULL`
  means mzLib's default.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**The defaults are mzLib's, not this package's.** Every digestion option
you leave as `NULL` is not sent, and mzLib's own default (its
`DigestionParams`) applies; the Parameters table below gives each value
at the pinned mzLib. mzLib names no default dissociation type, so the
bridge uses HCD until it does (smith-chem-wisc/mzLib#1488). The result
says which protease, dissociation type, terminus and caps were used.

## Value

An `mzlibr_digest`: a list with the scalars describing the run, a
`census` (see
[`census_explain`](https://smith-chem-wisc.github.io/mzLibR/reference/census_explain.md)),
and three data.frames that join on `peptide_index` - `peptides`,
`fragments` and `modifications`.

The scalars: `sequence_length`, the protein's length in residues;
`protease`, `dissociation` and `terminus`, each the one passed or the
default used; `max_modifications`, the cap in modifications per
peptidoform, and `max_modification_isoforms`, the cap in peptidoforms
per peptide position, each the value passed or mzLib's default;
`peptides_at_isoform_cap`, the peptide positions that reached it -
non-zero means the list is truncated; and `min_length` and `max_length`
in residues and `missed_cleavages` in cleavage sites, as passed, or `NA`
when left to mzLib, since the bridge does not report those three back.
The census counts `sites` in residues (distinct positions carrying a
modification), `applied` in modifications and `annotated` in UniProt
features.

`peptides` has `monoisotopic_mass`, the neutral mass in Da; `length` in
residues; `missed_cleavages` in cleavage sites; `modification_count` in
modifications; and `fixed_charges` (charge), **deprecated and always
0**. A permanently charged modification (trimethyllysine,
dimethylproline, phosphocholine-Ser and the other UniProt entries
PSI-MOD marks charged) counts at mzLib's charge-corrected mass, the
neutral convention Unimod and MetaMorpheus use: trimethyllysine is C3H6,
42.046951 Da, not UniProt's 43.054227. So no peptide carries a charge on
top of `monoisotopic_mass`. The `modifications` table's `mass` is
corrected the same way, and its `formal_charge` is likewise always 0 and
deprecated. Both columns stay until the bridge's next protocol version
removes them.

`peptides` **holds peptidoforms, not distinct sequences.** One row per
sequence-and-modification-placement, so albumin's default digest is
**310** rows over **202** distinct sequences. Both are legitimate
answers to "how many peptides" and they are not interchangeable; quoting
one for the other is a large error, not a rounding one. See
[`digest_distinct_base_sequences`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_distinct_base_sequences.md).

## Wraps

Wire verb `peptidoform fragments`. Generated from the bridge's verb spec
`peptidoform.fragments.yaml` (bridge commit `84c3b984b6b9`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`ProteinDbLoader.LoadProteinXML`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs)
  in `mzLib/UsefulProteomicsDatabases/ProteinDbLoader.cs` at mzLib
  `23c2490e`

- [`Mods.UniprotModifications`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/Modifications/Mods.cs)
  in `mzLib/Omics/Modifications/Mods.cs` at mzLib `d10017ce`

- [`Protein.Digest`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Proteomics/Protein/Protein.cs)
  in `mzLib/Proteomics/Protein/Protein.cs` at mzLib `23c2490e`

- [`DigestionParams`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Proteomics/ProteolyticDigestion/DigestionParams.cs)
  in `mzLib/Proteomics/ProteolyticDigestion/DigestionParams.cs` at mzLib
  `d10017ce`

- [`ProteaseDictionary`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Proteomics/ProteolyticDigestion/ProteaseDictionary.cs)
  in `mzLib/Proteomics/ProteolyticDigestion/ProteaseDictionary.cs` at
  mzLib `d10017ce`

- [`proteases.tsv (each protease name's cleavage motif)`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Proteomics/ProteolyticDigestion/proteases.tsv)
  in `mzLib/Proteomics/ProteolyticDigestion/proteases.tsv` at mzLib
  `d10017ce`

- [`PeptideWithSetModifications.Fragment`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Proteomics/ProteolyticDigestion/PeptideWithSetModifications.cs)
  in
  `mzLib/Proteomics/ProteolyticDigestion/PeptideWithSetModifications.cs`
  at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `accession`:

  string; required. A UniProtKB accession, e.g. P02768; fetched from
  rest.uniprot.org. The bridge checks only that it is not blank (the
  bindings upper-case it and check its shape first). UniProt's 400 or
  404 becomes a usage error.

- `protease`:

  string; default absent; range `a key of mzLib's ProteaseDictionary`.
  mzLib protease name. Absent, the bridge sends none and
  DigestionParams' default applies ("trypsin" at d10017ce). The names
  and what each one cuts are mzLib's, defined in proteases.tsv. At
  d10017ce "trypsin" is K\|,R\|, which cuts after K and R even before P;
  "trypsin\|P" is K\[P\]\|,R\[P\]\|, where \[P\] marks a residue that
  prevents the cut, so it does not cut before P. These are not
  MaxQuant's or Mascot's names: their "Trypsin/P" ignores the proline
  rule. On serum albumin (P02768) "trypsin" gives 202 distinct peptides
  and "trypsin\|P" gives 195 (other options absent).

- `dissociation`:

  string; default `HCD`; range
  `a mzLib DissociationType name, case-insensitive: HCD, CID, ETD, EThcD, ...`.
  Which fragment series to generate: HCD and CID give b and y; ETD gives
  c and zDot. Absent, HCD. mzLib's Fragment requires a type and names no
  default, so HCD is the bridge's own stop-gap until mzLib names one.
  mzLib#1488 asks for HCD, which MetaMorpheus's CommonParameters already
  defaults to.

- `terminus`:

  string; default absent; range
  `a mzLib FragmentationTerminus name, case-insensitive: Both, N, C`.
  Which fragment termini to generate. It is also passed to
  DigestionParams. Absent, DigestionParams' default applies (Both at
  d10017ce), to the digest and to the fragments.

- `modifications` (wire `--no-modifications`):

  flag; default `FALSE`. Digest the bare sequence: UniProt's
  modifications are discarded but its signal-peptide and propeptide
  boundaries are kept (pyMzLib#8), so the backbone peptide list is
  unchanged and only the modified variants go away. A clean control.

- `missed_cleavages` (wire `--missed-cleavages`):

  int; in **cleavage sites**; default absent; range `>= 0`. Maximum
  missed cleavage sites per peptide (DigestionParams
  maxMissedCleavages). Absent, DigestionParams' default applies (2 at
  d10017ce).

- `min_length` (wire `--min-length`):

  int; in **residues**; default absent; range `>= 1`. Shortest peptide
  kept (DigestionParams minPeptideLength). Absent, DigestionParams'
  default applies (7 at d10017ce). That default drops shorter peptides
  silently: histone H3.1 (P68431) keeps 49 of the 96 peptides it gives
  at 1 (distinct start and end, other options absent). Pass 1 for every
  peptide.

- `max_length` (wire `--max-length`):

  int; in **residues**; default absent; range `>= 0`. Longest peptide
  kept (DigestionParams maxPeptideLength); 0 means unbounded. Absent,
  DigestionParams' default applies (unbounded at d10017ce).

- `max_modifications` (wire `--max-mods`):

  int; in **modifications**; default absent; range `>= 0`. Maximum
  modifications per peptidoform (DigestionParams maxModsForPeptides).
  Absent, DigestionParams' default applies (2 at d10017ce). Isoforms
  grow combinatorially: histone H3.1's 49 peptides give 2,419
  peptidoforms at 2 and 6,580 at 3, where 2 peptide positions already
  reach the isoform cap (other options absent).

- `max_isoforms` (wire `--max-isoforms`):

  int; in **peptidoforms**; default absent; range `>= 1`. Cap on
  modification isoforms per peptide position (DigestionParams
  maxModificationIsoforms). Absent, DigestionParams' default applies
  (1024 at d10017ce). mzLib truncates SILENTLY at the cap;
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

  string; never `NA`. The protease used: the one given, or
  DigestionParams' default.

- `dissociation`:

  string; never `NA`. DissociationType name used: the one given, or HCD.

- `terminus`:

  string; never `NA`. FragmentationTerminus name used: the one given, or
  DigestionParams' default.

- `max_modifications`:

  int; in **modifications**; never `NA`. max-mods as used: the value
  given, or DigestionParams' default.

- `max_modification_isoforms`:

  int; in **peptidoforms**; never `NA`. max-isoforms as used: the value
  given, or DigestionParams' default.

- `peptides_at_isoform_cap`:

  int; in **peptide positions**; never `NA`. Peptide positions (start,
  end) whose peptidoforms in mzLib's raw digest, counted before the
  bridge collapses the mzLib#1108 duplicates, reached max-isoforms.
  Non-zero means the list is TRUNCATED.

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
  included. A permanently charged modification counts at mzLib's
  charge-corrected mass, the neutral convention Unimod and MetaMorpheus
  use (trimethyllysine C3H6, 42.046951, not UniProt's 43.054227), so m/z
  at charge z is (monoisotopic_mass + z \* 1.00727646677) / z.

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

  int; in **charge**; never `NA`. DEPRECATED: always 0, and removed at
  the next protocol bump. mzLib loads each permanently charged
  modification at its neutral-convention mass, so no peptide carries a
  charge on top of monoisotopic_mass and (monoisotopic_mass + (z -
  fixed_charges) \* 1.00727646677) / z is the ordinary formula. Until
  2026-10-27 the bridge inferred it from each modification's
  formula-vs-recorded mass gap, while monoisotopic_mass carried
  UniProt's charged mass.

- `modifications`:

  object\[\]; never `NA`. {one_based_residue, terminus, id, mass,
  formal_charge} per modification: one_based_residue is 1-based within
  the PEPTIDE and null for a terminal modification, which sets terminus
  'N' or 'C' instead (terminus is null for a residue modification); id
  is mzLib's IdWithMotif; mass is its monoisotopic mass in Da, never
  null for a loaded modification (mzLib keeps one only if it has a
  formula or a mass, and fills a missing mass from the formula),
  charge-corrected as for monoisotopic_mass; formal_charge is
  DEPRECATED, always 0 for the same reason as fixed_charges.

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

  any other UniProt status (e.g. 403); mzLib failing to load its
  built-in modification lists, to parse the entry or to digest it

## Caveats

- Since mzLib 1.0.595 (#1320), a protease with a preventing residue
  holds its rule at the ends of the sequence. trypsin\|P used to cut a
  protein that begins RP or KP after its first residue (RPAAAAK became R
  and PAAAAK); it no longer does, and the other \|P proteases
  (chymotrypsin\|P, elastase\|P, Lys-C\|P, subtilisin\|P) change the
  same way. mzLib's default trypsin has no preventing residue, so a
  default digest is unaffected, and so is any protein that begins with
  M.

- UniProt annotates more than mzLib applies: only 'modified residue' and
  'lipid moiety-binding region' features load. Serum albumin: 14 of 38
  annotated features applied, the 24 glycosylation-site features dropped
  on type (mzLib#1112). Read uniprot_features_by_type and
  unresolved_modifications before trusting a modification count.

- max-isoforms truncates silently in mzLib. Histone H3.1 at four
  modifications keeps 8,644 of its 12,808 peptidoforms at the default
  cap of 1024, losing about a third. Check peptides_at_isoform_cap
  before treating the list as exhaustive.

- Duplicate peptidoforms that mzLib's Digest emits where a chain
  boundary meets the initiator-Met cleavage site (mzLib#1108, open) are
  collapsed by the bridge on (start, end, full sequence).

- Fragment masses are neutral, charged modifications included at their
  corrected mass, so a fragment's m/z at charge z is (neutral_mass + z
  \* 1.00727646677) / z.

- Since 2026-10-27 a peptide carrying a permanently charged modification
  (the UniProt entries PSI-MOD marks charged, e.g. trimethyllysine,
  dimethylproline, phosphocholine-Ser) reports monoisotopic_mass, and
  its fragments neutral_mass, one proton lower per charge, with
  fixed_charges 0. Its peptide m/z is unchanged. N,N,N-trimethylglycine
  stays at UniProt's uncorrected mass, one hydrogen heavy, until the
  bridge's mzLib includes mzLib#1415, whose charge table is the first to
  list it.

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

- fixed_charges and formal_charge are always 0 and deprecated. Removing
  them breaks the wire (CONTRACT.md) and each binding's public mz():
  mzLibR's peptide_mz would return NA if the column vanished. Remove
  them, and simplify the three mz() functions, at the next protocol bump
  (G-bridge-charge-fields-removal).

- Integer options are not range-checked by the bridge (negative
  min-length, max-isoforms 0); the bindings check them, not alike
  (min_length 0 passes pyMzLib and mzLibRust, not mzLibR). Enum.TryParse
  also accepts numeric strings: `--dissociation` 3 passes as SID, and
  `--dissociation` 99 or `--terminus` 7 pass as undefined members.
  mzLib's behaviour on those values is unverified.

- The envelope echoes the protease, dissociation, terminus and both caps
  it used, but not missed-cleavages, min-length or max-length. Since
  2026-10-21 those can be mzLib's defaults, so a binding cannot show the
  values used (mzLibR prints NA as 'mzLib's default'). Echoing them,
  read from the DigestionParams the bridge built, would be additive and
  optional (CONTRACT.md).

- Unmeasured: mzLib caps isoforms per proteolytic peptide
  (ProteolyticPeptide.cs:62), so a position that two proteolysis
  products both produce (the mzLib#1108 duplicates) is capped twice but
  counted once. It could then read as at the cap when neither copy was
  cut.

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
#>   609 residues, trypsin, HCD, terminus Both
#>   min_length mzLib's default, missed cleavages mzLib's default, max 1 modifications
#>   2 peptidoforms over 2 distinct sequences
#>   128 fragments: b=64, y=64
#>   ! 14 of 38 annotated modifications applied. See ?census_explain.
head(digest$peptides[, c("base_sequence", "monoisotopic_mass", "modification_count")])
#>                         base_sequence monoisotopic_mass modification_count
#> 1     ALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK          3562.853                  0
#> 2 RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK          4136.902                  1
digest_fragments_by_series(digest)
#>   product_type  n
#> 1            b 64
#> 2            y 64
```
