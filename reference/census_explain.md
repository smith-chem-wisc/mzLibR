# What UniProt annotated, and what mzLib could actually use

A plain-language account, because the alternative is a number arriving
with no indication that a rule was ever applied. For serum albumin, 14
modifications are applied out of 38 annotated.

## Usage

``` r
census_explain(digest)
```

## Arguments

- digest:

  An
  [`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  result, or its `census`.

## Details

**The exclusion is correct - do not try to defeat it.** mzLib loads only
`modified residue` and `lipid moiety-binding region` annotations. On
albumin the census reports the 24 excluded features under one feature
\*type\*, `glycosylation site` - that is the label you will see in
`census$by_type`, and it is the only granularity the census has. At
UniProt's finer modification-\*name\* level, 22 of those 24 are
specifically `N-linked (Glc) (glycation) lysine`; the census does not
surface that, so the "22" is something you confirm by reading the
UniProt entry, not a number this tool reports. Either way the exclusion
is right: glycation and glycosylation are labile, heterogeneous adducts,
so assigning one an exact mass and a clean fragment ladder would
describe a species you cannot observe.

**The defect is the silence, not the exclusion**
(smith-chem-wisc/mzLib#1112). Note also that mzLib reads no qualifiers
for any feature type. It cannot tell an annotation marked `; in vitro`,
or one that exists only in a disease variant (albumin's Redhill and
Casebrook), from any other - different grounds for exclusion needing
different judgements, which the census cannot make for you. Read the
UniProt entry before concluding anything about a specific site.

## Value

A single string.

## Examples

``` r

digest <- peptidoform_fragments("P02768", max_modifications = 1)
census_explain(digest)
#> [1] "14 of 38 annotated modifications were applied, across 14 residue positions. Excluded by type: 24 x glycosylation site - mzLib loads only 'modified residue' and 'lipid moiety-binding region' annotations, so these were dropped on feature type alone. The exclusion is usually right: a glycation or glycosylation annotation describes a labile, heterogeneous adduct, so assigning it one exact mass and a clean fragment ladder would invent a species you cannot observe. But the reason is not reported, and the qualifier is not read - some annotations are marked 'in vitro' and some exist only in disease variants, which are different grounds for exclusion needing different judgements from you. Read the annotations on the UniProt entry before concluding anything about a specific site; this census can only tell you the count (smith-chem-wisc/mzLib#1112)."
```
