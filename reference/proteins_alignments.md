# One alignment per pair, each with its residue map

Splits a
[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
result into one `mzlibr_alignment` per pair, in input order, for
[`proteins_map_source_to_target`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md)
and
[`proteins_map_target_to_source`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_target_to_source.md).
Nothing is computed: each carries its row of `records` and its pair's
slice of `residue_map`.

## Usage

``` r
proteins_alignments(alignments)
```

## Arguments

- alignments:

  A
  [`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
  result.

## Value

A list with one `mzlibr_alignment` per pair. Each holds `pair_index`
(1-based), `aligner_id`, `score` (half-bits), `identical`,
`substituted`, `source_residues_facing_gap` and
`target_residues_facing_gap` (residues), `source_length` and
`target_length` (residues), `aligned_source`, `aligned_target`,
`caveats`, and `source_to_target` and `target_to_source`: for each
1-based position `p` of one sequence, element `p` is the paired position
in the other, or `NA` when that residue faces a gap. Those two are
`NULL` when the call asked for no residue map.

## See also

[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)

## Examples

``` r

alns <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"),
                            c("MLLLPEPTIDE", "PEPTIDE")))
each <- proteins_alignments(alns)
vapply(each, function(a) proteins_map_target_to_source(a, 1), numeric(1))
#> [1] 1 5
```
