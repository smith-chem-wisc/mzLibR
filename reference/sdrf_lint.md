# Find the concepts a set of SDRF documents annotated inconsistently

Find the concepts a set of SDRF documents wrote inconsistently, one row
per finding and spelling.

Validity is a property of one file; comparability is a property of the
relationship between files. Every document can pass
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)
and the pooled table still be unusable because one writes
`"Homo sapiens"` and another `"homo sapiens"`. This is mzLib's
`SdrfDriftLint`.

## Usage

``` r
sdrf_lint(documents, timeout = 60)
```

## Arguments

- documents:

  A character vector of paths. **Name it** to choose each document's
  label, as for
  [`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md):
  `c(cohort = "a.sdrf.tsv", partner = "b.sdrf.tsv")`. Name every element
  or none; mzLib's default label depends on where the files sit.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

Columns unique per row by construction are never compared -
`source name`, `assay name`, `comment[data file]`,
`comment[searched data file]`, `comment[file uri]` and
`comment[source document]` - so differing file names are not drift.
Reserved words and empty cells are skipped too, so documents that say
nothing lint clean: use
[`sdrf_assess`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)
for that.

## Value

An `mzlibr_sdrf_drift`: `document_count` documents linted, `paths`,
`labels`, `finding_count` distinct findings (not rows), `caveats`,
`column_names`, and `records`, a long data.frame with **one row per
finding x variant**: `finding_index` (1-based, most impactful first),
`kind`, `concept`, `column_name` (`NA` for a `ColumnNameVariant`, which
is about a name rather than a column's values), `variant_rank` (1 for
the majority spelling), `value` (this spelling, as written),
`occurrences` (how many documents used it), and `documents`, a list
column of the labels of those documents.

## Kinds of drift

`AccessionNameConflict` - one accession written with several names.
`NameAccessionConflict` - one name with several accessions, the serious
one. `MixedTermAndFreeText` - a column mixing controlled-vocabulary
terms and free text. `ColumnNameVariant` - column names differing only
by case or spacing. `ValueCaseVariant` - values differing only by case.
The majority spelling is not advice: it is only the most common.

## Wraps

Wire verb `sdrf lint`. Generated from the bridge's verb spec
`sdrf.lint.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SdrfDriftLint.Analyze`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfDriftLint.cs)
  in `mzLib/Readers/Sdrf/SdrfDriftLint.cs` at mzLib `23c2490e`

- [`SdrfCollection`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Sdrf/SdrfCollection.cs)
  in `mzLib/Readers/Sdrf/SdrfCollection.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `documents` (wire `--stdin`):

  string\[\]; required; range `>= 1 line`. Exactly as sdrf pool: one
  document per line, path then optionally a TAB and a label; label all
  or none. Labels are how the documents column names each file.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `document_count`:

  int; in **documents**; never `NA`. Documents linted.

- `paths`:

  string\[\]; never `NA`. The paths in stdin order.

- `labels`:

  string\[\]; never `NA`. Each document's label: supplied, or mzLib's
  containing-folder/file-stem default.

- `finding_count`:

  int; in **findings**; never `NA`. Distinct findings (not rows).

- `column_names`:

  string\[\]; never `NA`. The table's columns.

- `records` (wire `columns`):

  table; never `NA`. One row per finding x variant, findings in mzLib's
  impact order, variants majority first.

- `caveats`:

  string\[\]; never `NA`. Five static caveats, plus one when no labels
  were supplied.

## Columns of `records`

- `finding_index`:

  int; never `NA`. 0-based finding; rows of one finding share it.

- `kind`:

  string; never `NA`. AccessionNameConflict \| NameAccessionConflict \|
  MixedTermAndFreeText \| ColumnNameVariant \| ValueCaseVariant
  (SdrfDriftKind member names).

- `concept`:

  string; never `NA`. What was written inconsistently: an accession, a
  term name, a normalised free-text value, or a normalised column name
  (lower-cased, whitespace runs collapsed), by kind.

- `column_name`:

  string; `NA` when a ColumnNameVariant finding, which is about a
  column's NAME rather than its values. The column whose values drift.

- `variant_rank`:

  int; never `NA`. 0 = the majority spelling, then by frequency.
  Descriptive, not advice.

- `value`:

  string; never `NA`. This spelling as written; for
  MixedTermAndFreeText, "controlled vocabulary term" or "free text".

- `occurrences`:

  int; in **documents**; never `NA`. How widely this spelling was used.
  At 23c2490e every kind counts documents.

- `documents`:

  string\[\]; never `NA`. Labels of the documents that used this
  spelling, a list cell (labels are caller text, so no delimiter could
  be undone).

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  stdin empty; a blank path; a file not found; some lines labelled and
  some not

- `mzlib_bridge_error` (correctness):

  mzLib cannot read one of the documents

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- The majority spelling describes what the documents do; it is not
  advice (organism is free text in 963 curated documents and a CV term
  in 275).

- Never compared, being unique per row by construction: source name,
  assay name, comment\[data file\], comment\[searched data file\] (mzLib
  1.0.592, \#1335), comment\[file uri\], comment\[source document\].

- Reserved words and empty cells are skipped, so documents that say
  nothing lint clean; see sdrf assess.

- No `--threads`: SdrfDriftLint.Analyze is one call over the whole
  collection and does not parallelise.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.sdrf.lint`

- Rust (mzLibRust): `mzlib::sdrf::lint_labelled`

- R (mzLibR): `sdrf_lint`

## Since

Wire protocol 1; pyMzLib not yet shipped; mzLibRust not yet shipped;
mzLibR not yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Rust and R spellings are the intended ones, mirroring pool /
  pool_labelled and sdrf_pool; not yet ported.

- documents is the first list-valued cell in a columnar table on the
  wire. R will need a list column; confirm mzLibR's table projection
  accepts one before porting.

- since.pymzlib is null until released.

## References

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics and its controlled-vocabulary key=value cells (Dai et
  al., Nat Commun 12, 5854, 2021)

## See also

[`sdrf_pool`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md),
[`sdrf_validate`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)

## Examples

``` r

drift <- sdrf_lint(c(cohort = "sdrf_cohort.sdrf.tsv", partner = "sdrf_cohort_partner.sdrf.tsv"))
drift
#> <mzlibr_sdrf_drift> 2 documents: cohort, partner
#>   4 finding(s)
#>   ! AccessionNameConflict: Exploris 480 | Orbitrap Exploris 480
#>   ! MixedTermAndFreeText: controlled vocabulary term | free text
#>   ! ColumnNameVariant: Characteristics[sex] | characteristics[sex]
#>   ! ValueCaseVariant: Homo sapiens | homo sapiens
#>   5 caveat(s) - read x$caveats once
drift$records[, c("finding_index", "kind", "value", "occurrences")]
#>   finding_index                  kind                      value occurrences
#> 1             1 AccessionNameConflict               Exploris 480           1
#> 2             1 AccessionNameConflict      Orbitrap Exploris 480           1
#> 3             2  MixedTermAndFreeText controlled vocabulary term           1
#> 4             2  MixedTermAndFreeText                  free text           1
#> 5             3     ColumnNameVariant       Characteristics[sex]           1
#> 6             3     ColumnNameVariant       characteristics[sex]           1
#> 7             4      ValueCaseVariant               Homo sapiens           1
#> 8             4      ValueCaseVariant               homo sapiens           1
drift$records$documents[[1]]
#> [1] "partner"
```
