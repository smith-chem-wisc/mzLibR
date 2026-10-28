# The target position paired with each source position

Looks positions up in the residue map mzLib's `TryMapSourceToTarget`
filled. A residue facing a gap has no counterpart, and its answer is
`NA` - never the nearest residue.

## Usage

``` r
proteins_map_source_to_target(alignment, position)
```

## Arguments

- alignment:

  One pair's alignment:
  [`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
  of one pair, or an element of
  [`proteins_alignments`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_alignments.md).

- position:

  1-based source positions, 1 to the source's length in residues; a
  vector is fine.

## Value

The 1-based target positions, one per `position`, `NA` where the source
residue faces a gap.

## See also

[`proteins_map_target_to_source`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_target_to_source.md),
[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)

## Examples

``` r

aln <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")))
proteins_map_source_to_target(aln, c(8, 9, 10, 11, 12))
#> [1]  8 NA NA NA  9
```
