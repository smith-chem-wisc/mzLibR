# Read accessions into entry, isoform or version, and applied variants, with mzLib

Read protein accessions back into their entry, isoform or version, and
applied sequence variants, with mzLib's accession grammars, never
repairing one that does not fit.

Answers \*\*which entry, isoform and variants does `P04637-2_S70N` or
`NP_000537.3_R72P` name?\*\* Wraps mzLib's
`VariantApplication.ParseAccession` (mzLib \#1382), and through it
`ProteinAccession.Parse`. It is the inverse of the names mzLib gives
variant proteoforms: `P04637_S70N` is entry `P04637` with `S70N`
applied, and `NP_000537.3_R72P` is RefSeq `NP_000537`, version 3, with
`R72P`. Use it to link a search result's accession to anything you keep
per entry - genes, orthologs, GO terms - without losing the isoform or
the variants. Cutting at the first `_` would turn `NP_000537` into `NP`.

## Usage

``` r
proteins_parse_accessions(accessions, timeout = NULL)
```

## Arguments

- accessions:

  A character vector, one accession per result row, in order, duplicates
  included. Sent on stdin, unchanged: not trimmed, not split. A blank or
  `NA` entry, or one with a line break, is refused, because it would
  shift every later row.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait.

## Details

**It parses and never repairs.** An accession outside UniProt's and
RefSeq's protein grammars comes back `"unrecognized"` and verbatim -
lower case, a stray space, a FASTA header, a `DECOY_` prefix, a
`|`-joined protein group, and `P12345_2`, mzLib's loader counter for a
\*different\* entry whose accession collided. mzLibR sends each
accession exactly as you give it.

## Value

An `mzlibr_parsed_accessions`. `record_count` counts the accessions
parsed, duplicates included. `namespace_counts` is a named vector,
namespace to the number of accessions in it, every namespace present
(`uniprot`, `refseq`, `unrecognized`), zeros included: a large
`unrecognized` count is worth reading, see `caveats`.
`with_applied_variants_count` counts the accessions that name applied
variants. `caveats` says what mzLib leaves unrecognized rather than
repairs.

`records` is a data.frame with one row per accession, in the order
given: `input`, exactly as given; `entry_accession`, the UniProt entry
or unversioned RefSeq accession (`P04637` for `P04637-2_S70N`), equal to
`input` when unrecognized - an **entry**, not a sequence; `isoform`, the
UniProt `-n` number, `NA` without one (always for RefSeq) or past
2,147,483,647; `version`, the RefSeq `.n` number, `NA` without one
(always for UniProt) or past 2,147,483,647; `namespace`;
`parent_accession`, the entry the variants were applied to, isoform or
version kept (`P04637-2` for `P04637-2_Q5*`), equal to `input` without
variants; `applied_variants`, the suffix as written (`S70N_A80T`), `NA`
without one; and `has_applied_variants`.

## Joining on it

`entry_accession` answers "same entry?", **not** "same sequence?":
`P04637` and `P04637-2` share it, but they are different sequences. Join
per-entry data on it, and sequence-level data on `input`. `namespace` is
spelled as
[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
spells it, so the two tables join.
`setNames(parsed$records$entry_accession, parsed$records$input)` maps
each accession to its entry, for grouping.

## Wraps

Wire verb `proteins parse-accession`. Generated from the bridge's verb
spec `proteins.parse-accession.yaml` (bridge commit `4e180efe8b59`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`VariantApplication.ParseAccession`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/BioPolymer/VariantApplication.cs)
  in `mzLib/Omics/BioPolymer/VariantApplication.cs` at mzLib `d10017ce`

- [`ProteoformAccession`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/BioPolymer/ProteoformAccession.cs)
  in `mzLib/Omics/BioPolymer/ProteoformAccession.cs` at mzLib `d10017ce`

- [`ProteinAccession.Parse`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/MzLibUtil/ProteinAccession.cs)
  in `mzLib/MzLibUtil/ProteinAccession.cs` at mzLib `d10017ce`

- [`AccessionNamespace`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/MzLibUtil/ProteinAccession.cs)
  in `mzLib/MzLibUtil/ProteinAccession.cs` at mzLib `d10017ce`

## Parameters: units, ranges and defaults

- `accessions`:

  string\[\]; required; range `>= 1 line`. Read from stdin (not an
  option): one accession per line, given to mzLib exactly as read, with
  no trimming, case folding or splitting. Blank lines are skipped. One
  row per line, in input order, duplicates included.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `record_count`:

  int; in **accessions**; never `NA`. Rows: the input accessions.

- `namespace_counts`:

  map\<string,int\>; in **accessions**; never `NA`. Namespace -\>
  accessions in it; all three keys (unrecognized, uniprot, refseq)
  always present, zeros included.

- `with_applied_variants_count`:

  int; in **accessions**; never `NA`. Accessions that name applied
  variants (has_applied_variants true).

- `column_names`:

  string\[\]; never `NA`. The table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per input, in input order.

- `caveats`:

  string\[\]; never `NA`. Seven static caveats.

## Columns of `records`

- `input`:

  string; never `NA`. The line exactly as read
  (ProteoformAccession.Verbatim).

- `entry_accession`:

  string; never `NA`. The UniProt entry or unversioned RefSeq accession
  (ProteinAccession.EntryAccession): P04637 for P04637-2_S70N, NP_000537
  for NP_000537.3_R72P. Equal to input when unrecognized. Names an
  entry, not a sequence.

- `isoform`:

  int; `NA` when no UniProt -n isoform suffix (always null for RefSeq
  and unrecognized), or a number too large for a 32-bit integer. The
  UniProt isoform number (ProteinAccession.Isoform).

- `version`:

  int; `NA` when no RefSeq .n version (always null for UniProt and
  unrecognized), or a number too large for a 32-bit integer. The RefSeq
  version (ProteinAccession.Version).

- `namespace`:

  string; never `NA`. uniprot, refseq or unrecognized: mzLib's
  AccessionNamespace, lower case as genes resolve and mzLib's
  GeneResolutionTsv write it.

- `parent_accession`:

  string; never `NA`. The accession of the entry the variants were
  applied to, as written, isoform or version kept
  (ProteoformAccession.Entry.Verbatim): P04637-2 for P04637-2_Q5\*.
  Equal to input when there are no applied variants.

- `applied_variants`:

  string; `NA` when the accession names no applied variant. The variant
  suffix as written, without the leading '\_'
  (ProteoformAccession.AppliedVariants): S70N, or S70N_A80T for several.

- `has_applied_variants`:

  bool; never `NA`. ProteoformAccession.HasAppliedVariants:
  applied_variants is not null.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  no accessions on stdin

- `mzlib_bridge_error` (correctness):

  never expected: ParseAccession does not throw

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Parses, never repairs (mzLib's rule): an accession outside the UniProt
  and RefSeq protein grammars is unrecognized and kept verbatim,
  entry_accession equal to input. That includes lower case, surrounding
  spaces, a FASTA header (sp\|P04637\|P53_HUMAN), a decoy, contaminant
  or entrapment prefix (DECOY_P04637), and ProteinDbLoader's collision
  counter (P04637_2, which names a different entry whose accession
  collided).

- entry_accession names the entry, not a sequence: P04637 and P04637-2
  share it but are different sequences, and the isoform UniProt displays
  for an entry can carry any number. It answers 'same entry', never
  'same sequence'.

- A protein-group accession joined with '\|' (P04637\|P02768) is
  unrecognized as a whole: mzLib does not split it. Split protein groups
  first.

- Only RefSeq protein accessions (NP\_, XP\_, YP\_, WP\_, AP\_) are
  recognized; a RefSeq transcript (NM\_) or an Ensembl id (ENSP...) is
  unrecognized.

- A variant suffix is read only after an entry mzLib recognizes, in the
  form mzLib's GetAccession writes: '\_', the original residues, a
  1-based position, the variant residues (S70N; a deletion A80; a stop
  gain Q5\*), repeated for several (S70N_A80T). applied_variants is that
  suffix as written: it is not checked against any sequence.

- An isoform or version number too large for a 32-bit integer
  (over 2147483647) reads as null while the accession still counts as
  UniProt or RefSeq: mzLib keeps the entry and drops the number.

- Each stdin line is one accession, given to mzLib exactly as read;
  blank lines are skipped. One row per line, in input order, duplicates
  included.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.parse_accessions`

- Rust (mzLibRust): `mzlib::proteins::parse_accessions`

- R (mzLibR): `proteins_parse_accessions`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Whether to split '\|'-joined protein-group accessions before parsing:
  mzLib does not, so the verb does not either. Splitting is a grouping
  rule, and it belongs in mzLib if anywhere.

- An isoform or version past int.MaxValue is dropped with no signal
  while the namespace still reads uniprot or refseq
  (ProteinAccession.NumberOrNull). Whether mzLib should leave such an
  accession unrecognized instead is an mzLib question, not raised
  upstream yet.

- Not yet ported to mzLibRust or mzLibR; the binding names above are
  proposals.

## References

- [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010):
  UniProt, whose published accession grammar (and -n isoform suffix)
  mzLib matches (UniProt Consortium 2025, NAR)

- [doi:10.1093/nar/gkv1189](https://doi.org/10.1093/nar/gkv1189):
  RefSeq, whose protein accession prefixes and .n versions mzLib matches
  (O'Leary et al. 2016, NAR)

## See also

[`proteins_resolve_genes`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md),
[`proteins_read`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)

## Examples

``` r

parsed <- proteins_parse_accessions(c("P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P",
  "P04637_S70N_A80T", "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"))
parsed
#> <mzlibr_parsed_accessions> 9 accession(s), 2 with applied variants
#>   unrecognized 4, uniprot 4, refseq 1
#>   ! unrecognized ones are kept verbatim, never repaired: see $caveats
parsed$records[, c("input", "namespace", "entry_accession", "isoform", "version",
                   "applied_variants")]
#>                 input    namespace     entry_accession isoform version
#> 1              P04637      uniprot              P04637      NA      NA
#> 2            P04637-2      uniprot              P04637       2      NA
#> 3            Q13409-2      uniprot              Q13409       2      NA
#> 4    NP_000537.3_R72P       refseq           NP_000537      NA       3
#> 5    P04637_S70N_A80T      uniprot              P04637      NA      NA
#> 6            P04637_2 unrecognized            P04637_2      NA      NA
#> 7        DECOY_P02768 unrecognized        DECOY_P02768      NA      NA
#> 8 sp|P04637|P53_HUMAN unrecognized sp|P04637|P53_HUMAN      NA      NA
#> 9       P04637|P02768 unrecognized       P04637|P02768      NA      NA
#>   applied_variants
#> 1             <NA>
#> 2             <NA>
#> 3             <NA>
#> 4             R72P
#> 5        S70N_A80T
#> 6             <NA>
#> 7             <NA>
#> 8             <NA>
#> 9             <NA>
parsed$namespace_counts
#> unrecognized      uniprot       refseq 
#>            4            4            1 
```
