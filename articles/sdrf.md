# SDRF-Proteomics experimental designs

Every other reader in mzLibR answers *what did the search find*. An
SDRF-Proteomics file answers *what was searched*: which sample, organism
part, replicate and instrument settings. That is the half you need to
group results across experiments.

## Read one document

``` r

doc <- sdrf_read("PXD000070.sdrf.tsv")
doc
#> <mzlibr_sdrf> PXD000070.sdrf.tsv
#>   31 columns (some names repeat - use sdrf_all())
#>   6 rows in the document, 6 returned
#>   ! Cells are RAW STRINGS, never interpreted. The SDRF key=value grammar ("NT=Oxidation;AC=UNIMOD:35") is left exactly as written, because it cannot be told apart from a cell that merely contains '=' and ';' - comment[file uri] routinely carries pre-signed URLs that do.
#>   ! column_names is a LIST, not a set: names may repeat, and their order is part of the document. comment[modification parameters] legitimately appears many times in one file, so a name identifies a POSITION rather than a column.
#>   ! An SDRF reserved word - "not available", "not applicable" - is a real value meaning the experiment stated an absence. It is NOT the same as a column the document does not have, and the two must not be collapsed.
```

An SDRF document is **row-major**, not a data.frame, for two reasons
that are both real:

- **Column names repeat.** `comment[modification parameters]`
  legitimately appears once per modification, so a name identifies a
  *position*, not a column.
- **Rows can be ragged.** A row may carry fewer cells than the header
  has names, and mzLib keeps that so the file round-trips.

So `columns` is a character vector that may repeat, `rows` a list of
character vectors, and the accessors do the looking up:

``` r

sdrf_value(doc, "characteristics[organism part]")
#> [1] "human erythrocytes" "human erythrocytes" "human erythrocytes"
#> [4] "human erythrocytes" "human erythrocytes" "human erythrocytes"
sdrf_indexes_of(doc, "comment[modification parameters]")
#> [1] 20 21 22 23 24 25 26 27
sdrf_all(doc, "comment[modification parameters]")[[1]]
#> [1] "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"           
#> [2] "NT=Oxidation;AC=UNIMOD:35;TA=M,W,H;MT=Variable"         
#> [3] "NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=Variable"    
#> [4] "NT=Acetyl;AC=UNIMOD:1;TA=K;MT=Variable"                 
#> [5] "NT=Ammonia-loss;AC=UNIMOD:385;TA=Q,C;MT=Variable"       
#> [6] "NT=Glu->pyro-Glu;AC=UNIMOD:27;TA=E;MT=Variable"         
#> [7] "NT=Deamidated;AC=UNIMOD:7;PP=Protein N-term;MT=Variable"
#> [8] "NT=Phospho;AC=UNIMOD:21;TA=S,T,Y,A;MT=Variable"
```

`NA` from
[`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
means the document has no such column, or the row is too short to reach
it. The reserved words `"not available"` and `"not applicable"` are
different: they are values the experiment chose to write, and they come
back as themselves.

When a data.frame is what you want,
[`sdrf_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
gives one - lossy when a name repeats, which
[`sdrf_has_repeated_columns()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_has_repeated_columns.md)
tells you first:

``` r

sdrf_has_repeated_columns(doc)
#> [1] TRUE
head(sdrf_records(doc)[, c("source name", "characteristics[organism]")])
#>   source name characteristics[organism]
#> 1    Sample 1     plasmodium falciparum
#> 2    Sample 2     plasmodium falciparum
#> 3    Sample 3     plasmodium falciparum
#> 4    Sample 4     plasmodium falciparum
#> 5    Sample 5     plasmodium falciparum
#> 6    Sample 6     plasmodium falciparum
```

A ragged document, from mzLib’s own fixtures:

``` r

ragged <- sdrf_read("PXD059974.sdrf.tsv")
sdrf_ragged_row_count(ragged)
#> [1] 17
```

## Pool several into one analysis table

Name the documents: the names become each row’s provenance label.

``` r

pooled <- sdrf_pool(c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv"), limit = 4)
pooled
#> <mzlibr_pooled_sdrf> 2 documents: malaria, colon
#>   38 columns (some names repeat - use sdrf_all())
#>   24 rows in the document, 4 returned
#>   ! truncated - rows were left behind
#>   ! Every cell a source document did not have is filled with the reserved word "not available". Those words are therefore partly an artefact of pooling rather than something an experiment said, and a fill-rate computed over this table is not the fill rate of the originals.
#>   ! Provenance is in the 'comment[source document]' column. It is stamped only where a source document had none, so pooling an already-pooled table keeps the inner labels rather than overwriting them.
#>   ! This is an ANALYSIS table, not something to deposit. source name + assay name + comment[label] is unique within one document, but two experiments may both have a "Sample 1", so the pooled table will usually violate SDRF's uniqueness rule. Use the source-document column as part of any key.
#>   ! Columns are the UNION over all 2 documents, ordered by SDRF's own block structure, and a name that repeats is carried at the highest multiplicity any single document used, so no values are dropped.
sdrf_source_documents(pooled)
#> [1] "malaria" "malaria" "malaria" "malaria"
```

Columns are the union over the documents, and a cell a document did not
have is filled with `"not available"` - so a fill rate computed on the
pooled table is not the fill rate of the originals. The pooled table is
for analysis, not for deposit: two experiments may both have a “Sample
1”, so keep
[`sdrf_source_documents()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_source_documents.md)
in any key.

## Is the document well-formed?

[`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)
checks the specification’s structural rules - required columns, ragged
rows, the uniqueness of `source name` + `assay name` +
`comment[label]` - and locates each finding by line:

``` r

skeleton <- sdrf_validate("sdrf_skeleton.sdrf.tsv")
skeleton
#> <mzlibr_sdrf_validation> sdrf_skeleton.sdrf.tsv: NOT valid
#>   7 error(s), 2 warning(s) over 2 rows
#>   ! Error RequiredColumn x7
#>   ! Warning RecommendedColumn x2
#>   4 caveat(s) - read x$caveats once
skeleton$records[, c("severity", "rule", "line_number", "column_name")]
#>   severity              rule line_number
#> 1    Error    RequiredColumn          NA
#> 2    Error    RequiredColumn          NA
#> 3    Error    RequiredColumn          NA
#> 4    Error    RequiredColumn          NA
#> 5    Error    RequiredColumn          NA
#> 6    Error    RequiredColumn          NA
#> 7    Error    RequiredColumn          NA
#> 8  Warning RecommendedColumn          NA
#> 9  Warning RecommendedColumn          NA
#>                                   column_name
#> 1                             technology type
#> 2                         comment[instrument]
#> 3                              comment[label]
#> 4             comment[cleavage agent details]
#> 5                comment[technical replicate]
#> 6                comment[fraction identifier]
#> 7 comment[proteomics data acquisition method]
#> 8                  characteristics[cell type]
#> 9            comment[modification parameters]
```

Warnings never make a document invalid. A document-wide finding, such as
a missing column, has `NA` for its line.

## Does it describe its samples?

A document can be perfectly valid and still say nothing: reserved words
in every sample column, one replicate number, no factor.
[`sdrf_assess()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)
is the gate for that, and its verdict is how many of three checks pass:

``` r

assessment <- sdrf_assess("sdrf_cohort.sdrf.tsv")
assessment
#> <mzlibr_sdrf_assessment> sdrf_cohort.sdrf.tsv: Informative
#>   factor value varies: yes; sample described: yes; biological replicate varies: yes
#>   4 caveat(s) - read x$caveats once
assessment$records[, c("role", "column_name", "filled", "distinct_values", "fill_rate")]
#>                    role                           column_name filled
#> 1          factor_value                 factor value[disease]     12
#> 2 sample_characteristic                  characteristics[age]     10
#> 3 sample_characteristic            characteristics[cell type]      0
#> 4 sample_characteristic              characteristics[disease]     12
#> 5 sample_characteristic        characteristics[organism part]     12
#> 6 sample_characteristic                  characteristics[sex]     12
#> 7  biological_replicate characteristics[biological replicate]     12
#>   distinct_values fill_rate
#> 1               2 1.0000000
#> 2               5 0.8333333
#> 3               0 0.0000000
#> 4               2 1.0000000
#> 5               1 1.0000000
#> 6               2 1.0000000
#> 7               6 1.0000000
```

For a corpus,
[`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)
reads every document in one bridge call:

``` r

batch <- sdrf_assess_many(c("sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv"))
batch$files[, c("path", "verdict")]
#>                     path     verdict
#> 1   sdrf_cohort.sdrf.tsv Informative
#> 2 sdrf_skeleton.sdrf.tsv    Skeleton
#> 3     PXD000070.sdrf.tsv     Partial
```

## Each sample, merged over its rows, with its age

[`sdrf_samples()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)
gives one row per sample x column. Where a sample’s rows disagree about
a column, mzLib withholds it rather than let the first row win, and
`status` says `"conflicting"`. Ages are parsed into years on the
`characteristics[age]` rows:

``` r

samples <- sdrf_samples("sdrf_cohort.sdrf.tsv")
samples$records[samples$records$status == "conflicting", c("source_name", "column_name")]
#>    source_name              column_name
#> 54          S6 characteristics[disease]
ages <- samples$records[samples$records$column_name == "characteristics[age]", ]
ages[, c("source_name", "value", "age_years", "age_precision", "age_refusal")]
#>    source_name         value age_years age_precision   age_refusal
#> 6           S1           58Y      58.0         Exact          <NA>
#> 15          S2       40Y-85Y      62.5         Range          <NA>
#> 24          S3         >=90Y      90.0    LowerBound          <NA>
#> 33          S4            63        NA          <NA>       no_unit
#> 42          S5 Not available        NA          <NA> reserved_word
#> 50          S6       40Y-40Y      40.0         Exact          <NA>
```

A bare `63` is refused as `"no_unit"`: 63 years and 63 days are both
plausible in one study.
[`sdrf_parse_ages()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md)
applies the same rules to any cells you have:

``` r

sdrf_parse_ages(c("58Y", "40Y-85Y", ">=90Y", "6-8 weeks", "63"))$records[
  , c("cell", "years", "min_years", "max_years", "refusal")]
#>             cell      years  min_years  max_years       refusal
#> 1            58Y 58.0000000 58.0000000 58.0000000          <NA>
#> 2          30Y6M 30.5000000 30.5000000 30.5000000          <NA>
#> 3        40Y-85Y 62.5000000 40.0000000 85.0000000          <NA>
#> 4        40Y-40Y 40.0000000 40.0000000 40.0000000          <NA>
#> 5          >=90Y 90.0000000 90.0000000         NA          <NA>
#> 6            <1Y  1.0000000  0.0000000  1.0000000          <NA>
#> 7      6-8 weeks  0.1341547  0.1149897  0.1533196          <NA>
#> 8             63         NA         NA         NA       no_unit
#> 9  not available         NA         NA         NA reserved_word
#> 10                       NA         NA         NA         empty
#> 11   about forty         NA         NA         NA    unreadable
```

## Do several documents agree?

Every document can pass
[`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)
and a pooled table still be unusable because one writes `"Homo sapiens"`
and another `"homo sapiens"`.
[`sdrf_lint()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md)
finds those, one row per spelling:

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
```

## Instead of `readers_read_records()`

[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
recognises SDRF too, but joins each row into one semicolon-separated
string, and SDRF’s own `NT=...;AC=...` grammar puts semicolons inside
cells. Use
[`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md).
