# The source position paired with each target position

The reverse of
[`proteins_map_source_to_target`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md),
from mzLib's `TryMapTargetToSource`.

## Usage

``` r
proteins_map_target_to_source(alignment, position)
```

## Arguments

- alignment:

  One pair's alignment:
  [`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
  of one pair, or an element of
  [`proteins_alignments`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_alignments.md).

- position:

  1-based target positions, 1 to the target's length in residues.

## Value

The 1-based source positions, one per `position`, `NA` where the target
residue faces a gap.

## See also

[`proteins_map_source_to_target`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md),
[`proteins_align`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)

## Examples

``` r

aln <- proteins_align(list(c("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")))
proteins_map_target_to_source(aln, 9)
#> [1] 12
```
