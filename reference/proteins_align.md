# Align protein sequences with mzLib, and map every residue onto the other sequence

Align protein sequences pairwise with mzLib's global affine-gap aligner,
and map every residue of each sequence to its counterpart in the other,
or to none.

Answers **where is residue 120 of this protein in its isoform, its
ortholog, or last year's entry?** A position means something only in the
sequence it was counted in, and counting the same number along the other
sequence is wrong after the first insertion or deletion. This aligns
each pair with mzLib's `PairwiseAligner.Align` (mzLib \#1428) and maps
every residue to its counterpart in the other sequence, or to `NA` when
it has none.

## Usage

``` r
proteins_align(pairs, gap_open = NULL, gap_extend = NULL, end_gaps = NULL,
  max_cells = NULL, residue_map = TRUE, timeout = NULL)
```

## Arguments

- pairs:

  A list of `c(source, target)` character vectors of one-letter
  sequences, one result row each, in order:
  `list(c(canonical, isoform))` for one pair. `Map(c, sources, targets)`
  builds the list from two vectors. Sent to mzLib exactly as given (case
  is kept, and compared case-insensitively). Either sequence may be
  empty, not both. mzLib refuses `-` and `.`; it scores any other letter
  BLOSUM62 does not name (`U`, `O`, a digit, a space) as `X`, so strip
  modifications first.

- gap_open:

  Cost of opening a gap, in BLOSUM62's half-bits, charged once per gap:
  a gap of `k` residues costs `gap_open + k * gap_extend`. `NULL` (the
  default) sends nothing, and mzLib's default applies: 11 half-bits.
  mzLib refuses a negative value.

- gap_extend:

  Cost of each residue in a gap, the first included, in half-bits.
  `NULL` (the default): mzLib's, 1 half-bit. mzLib refuses less than 1.

- end_gaps:

  `"free"` - a gap at either end of either sequence costs nothing, so a
  signal peptide or a longer tail on one side is not penalised - or
  `"charged"`. `NULL` (the default): mzLib's, `"free"`.

- max_cells:

  The largest table mzLib will build,
  `(nchar(source) + 1) * (nchar(target) + 1)` cells of one byte each. A
  larger pair is refused rather than allowed to exhaust memory. `NULL`
  (the default): mzLib's, 200,000,000 cells (about 200 MB). Two titin
  isoforms need about 1.2 billion; mzLib accepts up to 2,147,483,591.

- residue_map:

  `TRUE` (the default) returns `residue_map`, one row per residue of
  both sequences. `FALSE` sends `--no-residue-map` and returns only the
  pairs table, for a large batch where you want the scores and counts.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait.

## Details

It is a **global** alignment with affine gaps (Gotoh's form of
Needleman-Wunsch), scored with **BLOSUM62**. All the pairs are aligned
in one bridge call, in order, so a few thousand isoform pairs cost one
process start.

**Every option you leave out is mzLib's to set.** mzLibR sends nothing
for it, and the result says what the aligner used (`gap_open`,
`gap_extend`, `end_gaps`, `max_cells`). mzLib's defaults in the pinned
build are gap open 11, gap extend 1, free end gaps and a
200,000,000-cell limit: BLASTP's gap costs for BLOSUM62, and EMBOSS
needle's free end gaps.

## Value

An `mzlibr_alignments`. `aligner_id` is mzLib's name for the method and
every parameter, e.g.
`"global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"`: store it with
any position you map, so the mapping says how it was made. `matrix` is
always `"BLOSUM62"`. `gap_open` and `gap_extend` are the gap costs used,
in BLOSUM62's half-bits; `end_gaps` the end-gap rule used; `max_cells`
the table limit used, in cells: yours, or mzLib's default when you gave
none. `record_count` counts the pairs aligned. `caveats` says what the
score and the mapping do and do not mean.

`records` is a data.frame with one row per pair, in input order:
`pair_index` (**1-based in R**, so it indexes `pairs`); `source_length`
and `target_length` in residues; `score`, the BLOSUM62 score in
half-bits (see the traps); `identical` and `substituted`, the paired
residues that are the same letter and that differ, in residues;
`source_residues_facing_gap` and `target_residues_facing_gap`, the
residues with no counterpart, in residues; and `aligned_source` and
`aligned_target`, each sequence with `-` wherever it faces a residue of
the other across a gap, the same length, to print one above the other.

`residue_map` is a data.frame with one row per residue of each sequence,
pair by pair, each pair's source residues first and then its target
residues: `pair_index` (1-based); `sequence`, `"source"` or `"target"`;
`position`, 1-based in that sequence, as in mzLib and UniProt;
`residue`, as given; `mapped_position`, the 1-based position of its
counterpart in the other sequence, \*\*`NA` when it faces a gap, never
the nearest residue\*\*; and `mapped_residue`, the counterpart's letter
(`NA` facing a gap), a different letter from `residue` on a
substitution. It is `NULL` with `residue_map = FALSE`.

## Reading one position

[`proteins_map_source_to_target`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md)
and
[`proteins_map_target_to_source`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_target_to_source.md)
look positions up for one pair;
[`proteins_alignments`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_alignments.md)
splits a result into one `mzlibr_alignment` per pair for them. They read
the residue map mzLib's `TryMapSourceToTarget` and
`TryMapTargetToSource` filled, and compute nothing.

## The traps

\*\*`score` is not an identity, a bit score or an E-value.\*\* It is a
raw sum of BLOSUM62 scores in half-bits, less the gaps, and it grows
with length. Compare scores only between alignments with the same
`aligner_id`. mzLibR computes no percent identity, because there are
several definitions: divide `identical` by the length you mean.

**It is a global alignment.** Every residue of both sequences is placed,
paired or facing a gap. With free end gaps a short sequence is placed on
the best-scoring stretch of a long one, even when it has no real match
there.

**One optimal alignment, chosen the same way every time.** Where two
alignments score the same, mzLib's traceback picks one, always the same
one. DYNC1I2 isoform 2B is a real case: canonical V76 and V82 are both
V, either could pair with the isoform's V76, and mzLib pairs V82.

**Always check the residue at a mapped position** before carrying a
modification across: a phosphosite on S is not a site on A.
`mapped_residue` sits beside `residue` for that.

## Wraps

Wire verb `proteins align`. Generated from the bridge's verb spec
`proteins.align.yaml` (bridge commit `60b1d25968a1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PairwiseAligner`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SequenceAlignment/PairwiseAligner.cs)
  in `mzLib/Omics/SequenceAlignment/PairwiseAligner.cs` at mzLib
  `d10017ce`

- [`PairwiseAligner.Align`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SequenceAlignment/PairwiseAligner.cs)
  in `mzLib/Omics/SequenceAlignment/PairwiseAligner.cs` at mzLib
  `d10017ce`

- [`PairwiseAligner.Id`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SequenceAlignment/PairwiseAligner.cs)
  in `mzLib/Omics/SequenceAlignment/PairwiseAligner.cs` at mzLib
  `d10017ce`

- [`PairwiseAlignment.TryMapSourceToTarget`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SequenceAlignment/PairwiseAlignment.cs)
  in `mzLib/Omics/SequenceAlignment/PairwiseAlignment.cs` at mzLib
  `d10017ce`

- [`PairwiseAlignment.TryMapTargetToSource`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SequenceAlignment/PairwiseAlignment.cs)
  in `mzLib/Omics/SequenceAlignment/PairwiseAlignment.cs` at mzLib
  `d10017ce`

- [`SubstitutionMatrix.Blosum62`](https://github.com/smith-chem-wisc/mzLib/blob/d10017ce/mzLib/Omics/SequenceAlignment/SubstitutionMatrix.cs)
  in `mzLib/Omics/SequenceAlignment/SubstitutionMatrix.cs` at mzLib
  `d10017ce`

## Parameters: units, ranges and defaults

- `pairs`:

  string\[\]; required; range `>= 1 pair; each line exactly one tab`.
  Read from stdin (not an option): one pair per line, the source
  sequence, a tab, the target sequence, as one-letter residue codes.
  Each sequence reaches mzLib exactly as given, case and all. Either
  sequence may be empty, but not both, because a blank line is skipped.
  One row per line, in input order, duplicates included. Sequences
  always travel on stdin, never argv: two titin isoforms are about
  70,000 characters, past argv's 32 KB ceiling on Windows.

- `gap_open` (wire `--gap-open`):

  int; in **half-bits (BLOSUM62 score units)**; default `11`; range
  `>= 0 (mzLib refuses a negative value)`. The cost of opening a gap,
  charged once per gap. A gap of length k costs gap-open plus k times
  gap-extend, BLAST's convention. Absent: nothing is sent and mzLib's
  PairwiseAligner default applies, which is 11 at d10017ce (BLASTP's
  default with BLOSUM62).

- `gap_extend` (wire `--gap-extend`):

  int; in **half-bits (BLOSUM62 score units)**; default `1`; range
  `>= 1 (mzLib refuses 0 or less)`. The cost of each residue in a gap,
  the first included. Absent: mzLib's default, 1 at d10017ce.

- `end_gaps` (wire `--end-gaps`):

  string; default `free`; range `free | charged (case-insensitive)`.
  free: a gap at either end of either sequence costs nothing, so one
  sequence may overhang the other (a signal peptide, a longer
  C-terminus) without penalty, as EMBOSS needle scores end gaps by
  default. charged: every gap costs. mzLib's
  PairwiseAligner.FreeEndGaps, true or false. Absent: mzLib's default,
  free (true) at d10017ce. The spelling is mzLib's own, from its aligner
  id.

- `max_cells` (wire `--max-cells`):

  int; in **cells**; default `2e+08`; range
  `1 to 2147483591 (Array.MaxLength; mzLib refuses anything outside)`.
  The largest dynamic-programming table mzLib will build, (source
  length + 1) x (target length + 1) cells of one byte each. A larger
  pair is refused rather than allowed to exhaust memory. Absent: mzLib's
  default, 200,000,000 cells (about 200 MB) at d10017ce.

- `residue_map` (wire `--no-residue-map`):

  flag; default `FALSE`. Leave out the residue map (residue_map is then
  null). The pairs table is unchanged. For a batch where only the scores
  and counts are wanted: the map has one row per residue of both
  sequences.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `aligner_id`:

  string; never `NA`. mzLib's PairwiseAligner.Id, naming the method and
  every parameter, e.g.
  global-affine;BLOSUM62;open=11;extend=1;end-gaps=free. Store it with
  any mapped position.

- `matrix`:

  string; never `NA`. The substitution matrix's name, always BLOSUM62:
  mzLib's default, the only one the wire offers.

- `gap_open`:

  int; in **half-bits (BLOSUM62 score units)**; never `NA`. The gap-open
  cost the aligner used: the option's value, or mzLib's default when it
  was absent.

- `gap_extend`:

  int; in **half-bits (BLOSUM62 score units)**; never `NA`. The
  gap-extension cost the aligner used.

- `end_gaps`:

  string; never `NA`. free or charged: the end-gap rule the aligner
  used.

- `max_cells`:

  int; in **cells**; never `NA`. The table limit the aligner used.

- `record_count`:

  int; in **pairs**; never `NA`. Rows: the input pairs.

- `column_names`:

  string\[\]; never `NA`. The pairs table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per pair, in input order.

- `residue_map`:

  object; `NA` when `--no-residue-map` was given. {row_count,
  column_names, columns}: one row per residue of each sequence of each
  pair, pair by pair, each pair's source residues (positions 1 to its
  length) and then its target residues (typed under
  result.tables.residue_map).

- `caveats`:

  string\[\]; never `NA`. Seven static caveats.

## Columns of `records`

- `pair_index`:

  int; never `NA`. 0-based position of the pair among the non-blank
  stdin lines.

- `source_length`:

  int; in **residues**; never `NA`. The source sequence's length.

- `target_length`:

  int; in **residues**; never `NA`. The target sequence's length.

- `score`:

  int; in **half-bits (BLOSUM62 score units)**; never `NA`.
  PairwiseAlignment.Score: substitution scores summed over the paired
  residues, less every charged gap. A raw score, not a bit score, an
  E-value or an identity; compare it only between alignments with the
  same aligner_id.

- `identical`:

  int; in **residues**; never `NA`. Paired residues that are the same
  letter, compared case-insensitively (PairwiseAlignment.Identical).

- `substituted`:

  int; in **residues**; never `NA`. Paired residues that differ
  (PairwiseAlignment.Substituted).

- `source_residues_facing_gap`:

  int; in **residues**; never `NA`. Source residues with no target
  counterpart: an extra or deleted stretch, or an overhang.

- `target_residues_facing_gap`:

  int; in **residues**; never `NA`. Target residues with no source
  counterpart.

- `aligned_source`:

  string; never `NA`. The source as given with '-' wherever it faces a
  target residue across a gap; as long as aligned_target.

- `aligned_target`:

  string; never `NA`. The target as given with '-' wherever it faces a
  source residue across a gap.

## Fields of `residue_map`

- `pair_index`:

  int; never `NA`. The pair the residue belongs to.

- `sequence`:

  string; never `NA`. source or target: which sequence position is in.

- `position`:

  int; never `NA`. 1-based position in that sequence, as everywhere in
  mzLib.

- `residue`:

  string; never `NA`. The residue at position, one letter, as given.

- `mapped_position`:

  int; `NA` when the residue faces a gap, so it has no counterpart
  (mzLib's TryMap returned false); never the nearest residue. 1-based
  position of the paired residue in the other sequence:
  TryMapSourceToTarget for a source row, TryMapTargetToSource for a
  target row.

- `mapped_residue`:

  string; `NA` when the residue faces a gap. The paired residue, one
  letter, as given. Different from residue on a substitution.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  no pairs on stdin; a line that is not exactly two tab-separated fields
  (the message names the pair)

- `mzlib_usage_error` (usage):

  `--gap-open`, `--gap-extend` or `--max-cells` not an integer, or given
  without a value; a value mzLib's PairwiseAligner constructor refuses
  (gap-open below 0, gap-extend below 1, max-cells below 1 or above
  2147483591), with mzLib's message and the option named; `--end-gaps`
  not free or charged, or given without a value

- `mzlib_usage_error` (usage):

  mzLib refused a pair: a sequence holding a gap character ('-' or '.'),
  or a table larger than max-cells. The message names the pair and
  quotes mzLib, which gives the cell count. The whole call fails on the
  first such pair, in input order.

- `mzlib_bridge_error` (correctness):

  never seen at d10017ce: mzLib's aligner throws nothing else for
  non-null input

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- score is a raw substitution-matrix score in BLOSUM62's half-bit units,
  summed over the alignment: it grows with length and is not a bit
  score, an E-value or a percent identity. Compare scores only between
  alignments with the same aligner_id.

- Global alignment: every residue of both sequences is in the alignment,
  paired or facing a gap. With end gaps free (mzLib's default) an
  overhang at either end costs nothing, so a short sequence is placed on
  the best-scoring stretch of a long one even when it has no real match
  there.

- A residue facing a gap has no counterpart: its mapped_position is
  null, never the nearest residue.

- One optimal alignment is returned. Where several score the same,
  mzLib's traceback prefers a residue pair, then a gap in the target,
  then a gap in the source, so the same input always gives the same
  alignment; a residue beside an indel may sit on the other side of the
  gap in another aligner's equally good answer. DYNC1I2 isoform 2B lacks
  six residues of the canonical sequence: V76 and V82 of the canonical
  are both V, either could pair with the isoform's V76, and mzLib pairs
  V82.

- mzLib refuses only the gap characters '-' and '.'. Any other character
  BLOSUM62 does not name (selenocysteine U, pyrrolysine O, J, a digit, a
  space, a bracket) is scored as X, so pass bare one-letter sequences.
  Letters are compared case-insensitively; aligned_source,
  aligned_target and the residue columns keep the case given.

- Time and memory grow with the product of the two lengths: mzLib builds
  a (source length + 1) x (target length + 1) table of one byte per cell
  and refuses a pair larger than max_cells. Two identical sequences free
  of X and \* are aligned residue for residue without the table, so
  max_cells does not limit them.

- Each stdin line is one pair, aligned in input order; blank lines are
  skipped. mzLib's aligner takes no thread count, so pairs are aligned
  one after another.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.align_many`

- Rust (mzLibRust): `mzlib::proteins::align_with` with `AlignOptions`

- R (mzLibR): `proteins_align`

## Since

Wire protocol 1; pyMzLib 0.6.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- A substitution matrix other than BLOSUM62 is not offered. mzLib can
  read one (SubstitutionMatrix.Parse(name, TextReader), NCBI's text
  format), but it needs a name for the aligner id, which the bridge
  would have to invent, and the units of every score would then be that
  matrix's. Whether to add `--matrix` FILE, and where its name comes
  from, is open.

- mzLib scores any character BLOSUM62 does not name as X, digits and
  spaces included, so a sequence with a stray modification bracket or a
  number aligns without complaint. Whether PairwiseAligner.Align should
  refuse non-letters, as it refuses '-' and '.', is an mzLib question,
  not raised upstream yet.

- Pairs are aligned one at a time: PairwiseAligner takes no degree of
  parallelism, so the wire has no `--threads` (PARALLELISM.md). A batch
  aligner with one (for proteome-wide ortholog or isoform maps) would be
  an mzLib addition.

- residue_map is not validated by tools/check_verbs.py, which checks
  only 'columns' and 'bulk'. The bridge test
  ADeletion_LeavesItsResiduesUnmappedAndMapsAcrossTheGap pins its column
  order.

- Not yet ported to mzLibRust or mzLibR; the binding names above are
  proposals.

## References

- [doi:10.1016/0022-2836(70)90057-4](https://doi.org/10.1016/0022-2836(70)90057-4):
  Global alignment by dynamic programming (Needleman and Wunsch 1970, J
  Mol Biol 48:443)

- [doi:10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9):
  Affine gap costs in the same time and memory order, the form mzLib
  implements (Gotoh 1982, J Mol Biol 162:705)

- [doi:10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915):
  BLOSUM62, the substitution matrix and its half-bit units (Henikoff and
  Henikoff 1992, PNAS 89:10915)

## See also

[`proteins_map_source_to_target`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md),
[`proteins_read`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)
with `sequences = TRUE` for the sequences of a database

## Examples

``` r

# mzLib's own deletion case: QRQ (source 9 to 11) is missing from the target.
aln <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")))
aln
#> <mzlibr_alignments> 1 pair(s), global-affine;BLOSUM62;open=11;extend=1;end-gaps=free
#>   pair 1: score 80 half-bits, 19 identical, 0 substituted, 3 + 0 facing a gap
#>   residue map: 41 rows; NA mapped_position faces a gap
writeLines(c(aln$records$aligned_source, aln$records$aligned_target))
#> MKTAYIAKQRQISFVKSHFSRQ
#> MKTAYIAK---ISFVKSHFSRQ
aln$records[, c("score", "identical", "source_residues_facing_gap")]
#>   score identical source_residues_facing_gap
#> 1    80        19                          3
proteins_map_source_to_target(aln, 8:12)
#> [1]  8 NA NA NA  9

# Two pairs in one call; the free overhang MLLL faces a gap.
alns <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"),
                            c("MLLLPEPTIDE", "PEPTIDE")))
alns$records[, c("pair_index", "score", "aligned_target")]
#>   pair_index score         aligned_target
#> 1          1    80 MKTAYIAK---ISFVKSHFSRQ
#> 2          2    39            ----PEPTIDE
proteins_map_source_to_target(proteins_alignments(alns)[[2]], 1:5)
#> [1] NA NA NA NA  1

# Charging the end gaps changes the overhang's score, and the aligner id says so:
charged <- proteins_align(list(c("MLLLPEPTIDE", "PEPTIDE")), end_gaps = "charged")
charged$records$score
#> [1] 24
charged$aligner_id
#> [1] "global-affine;BLOSUM62;open=11;extend=1;end-gaps=charged"
```
