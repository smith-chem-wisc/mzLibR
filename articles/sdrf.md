# SDRF-Proteomics experimental designs

Every other reader in mzLibR answers *what did the search find*. An
SDRF-Proteomics file answers *what was searched*: which sample, organism
part, replicate and instrument settings. That is the half you need to
group results across experiments.

| question | function | mzLib |
|----|----|----|
| what does this SDRF say, row by row? | [`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md) | `SdrfDocument.LoadResults` |
| several SDRFs as one analysis table, with provenance | [`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md) | `SdrfCollection.Merge` |
| is it well-formed? | [`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md), [`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md) | `SdrfValidator.Validate` |
| does it describe its samples, or is it a valid skeleton? | [`sdrf_assess()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md), [`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md) | `SdrfSampleInformativeness.Assess` |
| each sample’s characteristics, with ages in years | [`sdrf_samples()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md), [`sdrf_parse_ages()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md) | `SdrfSampleBlock.BySourceName`, `SdrfAge.TryParse` |
| do several SDRFs spell their terms alike? | [`sdrf_lint()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md) | `SdrfDriftLint.Analyze` |
| can it drive a label-free quantification, and with what design? | [`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md) | `SdrfLabelFreeDesign.Read` (#1363) |
| which rows did my search read, when it read `X-calib.mzML` and the SDRF names `X.raw`? | [`sdrf_restrict()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md) | `SdrfSearchScope.Restrict` |
| can it drive a TMT, iTRAQ or DiLeu quantification, channel by channel? | [`sdrf_isobaric_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md) | `SdrfIsobaricDesign.Read` |

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

`out` writes the whole pooled document, whatever `limit` is. **mzLib
adds `.sdrf.tsv` to a name that does not already end in it**, so
`out = "x.tsv"` writes `x.tsv.sdrf.tsv`; `written$path` names the file
it wrote. Name the file `.sdrf.tsv` and the two agree:

``` r

# Not run: writes a file, which the recordings this page replays did not.
pooled <- sdrf_pool(c(malaria = "PXD000070.sdrf.tsv", colon = "PXD026824.sdrf.tsv"),
                    out = "pooled.tsv")
pooled$written$path   # ".../pooled.tsv.sdrf.tsv"
```

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

A real deposit, PXD000070, is `Partial`: it describes its samples and
varies its replicates, but no factor value varies, so it cannot say what
was compared. A valid skeleton says nothing:

``` r

pxd70 <- sdrf_assess("PXD000070.sdrf.tsv")
str(pxd70[c("verdict", "factor_value_varies", "sample_is_described", "biological_replicate_varies")])
#> List of 4
#>  $ verdict                    : chr "Partial"
#>  $ factor_value_varies        : logi FALSE
#>  $ sample_is_described        : logi TRUE
#>  $ biological_replicate_varies: logi TRUE
empty <- sdrf_assess("sdrf_skeleton.sdrf.tsv")
empty$verdict
#> [1] "Skeleton"
empty$records[, c("column_name", "filled", "rows")]
#>                             column_name filled rows
#> 1              characteristics[disease]      0    2
#> 2        characteristics[organism part]      0    2
#> 3 characteristics[biological replicate]      2    2
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

## Turn an SDRF into a quantification design

FlashLFQ and MetaMorpheus need an experimental design: which runs are
the same condition, which are biological replicates of one another,
which are fractions of one sample. The SDRF already says all of that. A
hand-written `ExperimentalDesign.tsv` is how most people still give it,
and an invalid one is worse than none: MetaMorpheus skips quantification
with one warning and no error.

[`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
reads the design out of the SDRF, or refuses and says why. It calls
mzLib’s `SdrfLabelFreeDesign`, which runs every check MetaMorpheus’s own
design validator runs and reports *all* the failures at once:

| question | function | mzLib |
|----|----|----|
| the design, or every reason there is none | `sdrf_design(path, condition_columns)` | `SdrfLabelFreeDesign.Read` (#1363) |
| only the runs one search read | `sdrf_design(..., searched_files = ...)` | `SdrfLabelFreeDesignOptions.SearchedFiles` |
| it as `flashlfq_quantify(spectra = )` | [`sdrf_design_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md) | `SdrfLabelFreeDesign.Files` |
| it as `flashlfq_median_polish(design = )` | [`sdrf_design_run_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_run_design.md) | `SdrfLabelFreeDesign.Files` |
| MetaMorpheus’s `ExperimentalDesign.tsv` | `sdrf_design(..., out = "ExperimentalDesign.tsv")` | `WriteExperimentalDesignTsv` |

The examples use mzLib’s own design files: PXD067622, a TurboID
proximity-labelling study (two genotypes by four treatments, three
biological replicates each), and PXD049018 (two pulldowns, ten SDS-PAGE
bands each).

### A valid design

PXD067622 has two factor columns. Say which make up the condition; their
values are joined with `_`, in the order you give:

``` r

both <- c("factor value[genotype]", "factor value[treatment]")
d <- sdrf_design("PXD067622.sdrf.tsv", condition_columns = both)
d
#> <mzlibr_sdrf_design> PXD067622.sdrf.tsv: 24 runs in 8 condition(s)
#>   files keyed on: comment[data file]
#>   condition from: factor value[genotype] + factor value[treatment]
head(d$records[, c("file_name", "condition", "biological_replicate")], 3)
#>                               file_name                       condition
#> 1 20240830_HF_LC3_MAA_RK_12032_CA_DMSO4 SPRTN-TurboID CA_DMSO (vehicle)
#> 2 20240830_HF_LC3_MAA_RK_12032_CA_DMSO5 SPRTN-TurboID CA_DMSO (vehicle)
#> 3 20240830_HF_LC3_MAA_RK_12032_CA_DMSO6 SPRTN-TurboID CA_DMSO (vehicle)
#>   biological_replicate
#> 1                    0
#> 2                    1
#> 3                    2
```

`biological_replicate`, `technical_replicate` and `fraction` are
**0-based** - mzLib’s `SpectraFileInfo` coordinates, exactly what the
FlashLFQ functions take - although every other position in mzLibR is
1-based. `out` writes MetaMorpheus’s file, which starts at 1; mzLib adds
the one when it writes, and nowhere else. Never add it yourself.

Hand the design straight to FlashLFQ:

``` r

head(sdrf_design_spectra(d), 2)
#>                                        path                       condition
#> 1 20240830_HF_LC3_MAA_RK_12032_CA_DMSO4.raw SPRTN-TurboID CA_DMSO (vehicle)
#> 2 20240830_HF_LC3_MAA_RK_12032_CA_DMSO5.raw SPRTN-TurboID CA_DMSO (vehicle)
#>   biological_replicate technical_replicate fraction
#> 1                    0                   0        0
#> 2                    1                   0        0
```

The SDRF names files without a directory. When your spectra live
elsewhere, or the search read only some of them, pass `searched_files`:
each must be named **exactly** (case and extension) by one row, rows for
other files are dropped and listed in `notes`, and your paths become
`full_path`.

### A refusal lists every reason

PXD049018’s depositors did not say which pulldown was treated: its
`factor value[treatment]` is `not available` on every row. A condition
built from an unknown factor would pair samples nobody said were alike,
so mzLib refuses, once per row, naming each file:

``` r

refused <- sdrf_design("PXD049018.sdrf.tsv", condition_columns = both)
refused$is_valid
#> [1] FALSE
c(refusals = length(refused$refusals), files = refused$file_count)
#> refusals    files 
#>       20        0
writeLines(refused$refusals[[1]])
#> Line 2 (MSB67868ABand_01.raw): 'factor value[treatment]' is 'not available'. A condition cannot be built from an unknown factor; fill it in, or leave the column out of the declared condition columns to pool these rows.
# The hand-offs refuse a refused design, naming every reason:
inherits(try(sdrf_design_spectra(refused), silent = TRUE), "try-error")
#> [1] TRUE
```

**A refusal is a result, not an error.**
[`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
returns normally with `is_valid` `FALSE`; only the two hand-offs stop,
so a refused design cannot reach FlashLFQ by accident. Without
`condition_columns`, mzLib uses the document’s only factor column - and
PXD067622 has two, so leaving them undeclared is itself a refusal:
joining every factor could split a condition on a nuisance factor.

### Replicate numbers are ranked within each condition

A drafted SDRF often copies the study-wide run index out of the file
names. MetaMorpheus wants replicates numbered within each condition, so
mzLib ranks them back, and records every renumbering in `notes`, the
only place it is kept:

``` r

ranked <- sdrf_design("PXD067622_studywide.sdrf.tsv", condition_columns = both)
writeLines(ranked$notes[[2]])
#> Condition 'SPRTN-TurboID CA_formaldehyde 1 mM, 1 h': biological replicates renumbered 22 -> 1, 23 -> 2, 24 -> 3.
identical(ranked$records, d$records)
#> [1] TRUE
```

**Label-free only.** An isobaric SDRF needs a channel design, which
[`sdrf_isobaric_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md)
reads, below. For the channels themselves, see the isobaric kits
article.

## Which rows did my search read?

A search reads `X-calib.mzML`; the SDRF names the acquired `X.raw`. The
design readers match `searched_files` by **exact** name, so the search’s
own names match nothing.
[`sdrf_restrict()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)
joins them with mzLib’s `SdrfSearchScope.Restrict`: **by stem**,
ignoring case, after removing the derivative suffixes MetaMorpheus
writes (`-calib`, `-averaged`, in any order and number). Each kept row
is carried cell for cell, with the searched file’s name in a new
`comment[searched data file]` column; nothing is renumbered.

The example is mzLib’s TMT10 fixture from PXD008841 (40 rows, 4 files).
The search read two calibrated runs, the raw file of one of them as
well, a renamed file whose acquisition you know, and one the SDRF does
not describe. **Name an element by the acquired file it came from** when
you know it; that overrides the suffix rule:

``` r

pool1 <- "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
pool5 <- "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
searched <- c(paste0("D:\\search\\", pool1, "1-calib.mzML"), paste0(pool1, "2-calib.mzML"),
              paste0(pool1, "2.raw"), "sample_A.mzML", "not_in_sdrf.mzML")
names(searched) <- c("", "", "", paste0(pool5, "1.raw"), "")
scoped <- sdrf_restrict("PXD008841.subset.sdrf.tsv", searched, limit = 2)
c(input = scoped$input_row_count, kept = scoped$row_count, returned = scoped$returned_count)
#>    input     kept returned 
#>       40       30        2
sdrf_searched_data_files(scoped)
#> [1] "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
#> [2] "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
```

The result is an SDRF document, so
[`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
and the other accessors work on it. Four reports say what the
restriction could not place, and none of them is repaired:

``` r

scoped$searched_without_row   # searched, but the SDRF has no row for it
#> [1] "not_in_sdrf.mzML"
scoped$dropped_data_files     # acquired, but this search did not read it
#> [1] "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr02.raw"
writeLines(strwrap(scoped$ambiguous, 80))   # two searched files claimed one acquisition
#> HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02:
#> HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02-calib.mzML,
#> HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02.raw
scoped$shared_stem            # one searched file whose rows name two different acquired files
#> character(0)
```

`fr02` was searched twice, as `-calib.mzML` and as `.raw`: the first in
name order kept its rows, and the other has none and is not counted as
missing. `out` writes the **whole** restricted document, as for
[`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md):
name it `.sdrf.tsv`, because mzLib adds that to a name that does not end
in it, and `written$path` names the file written.

**Then read the design from the restricted document**, with the search’s
own paths as `searched_files`:
[`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
and
[`sdrf_isobaric_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md)
key on `comment[searched data file]` when a document has it. Give those
paths as this platform writes them.
[`sdrf_restrict()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)
reduces either separator on any platform, but the design readers split
only the running platform’s, so `D:\search\x.mzML` is not a file name on
Linux. [`file.path()`](https://rdrr.io/r/base/file.path.html) writes
`/`, which both platforms read.

## An isobaric design, channel by channel

[`sdrf_isobaric_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md)
is the TMT, TMTpro, iTRAQ and DiLeu counterpart of
[`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md).
It calls mzLib’s `SdrfIsobaricDesign.Read`, which runs MetaMorpheus’s
own checks on a TMT design - one plex, fraction and technical replicate
per file; a channel described once per plex; no sample twice at one
fraction in a plex - and refuses rather than repairs, reporting every
failure at once.

**It guesses nothing.** Not the kit, because one kit’s labels can be a
subset of another’s (TMT11 with 131C unannotated reads exactly like
TMT10), and not the plex, which you declare with exactly one of
`plex_column`, `plex_file_name_pattern` and `single_plex`. PXD008841
writes no plex column; the plex survives only in its file names, spelled
two ways (`TMTpool1`, `TMT_pool5`), so one .NET regular expression reads
both:

``` r

d <- sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10",
                          plex_file_name_pattern = "TMT_?pool(?<plex>\\d+)")
d
#> <mzlibr_sdrf_isobaric_design> PXD008841.subset.sdrf.tsv: 4 files, 2 plex(es), kit TMT10, 40 channel rows (1-based numbers)
#>   plex from: file-name pattern 'TMT_?pool(?<plex>\d+)' (declared)
#>   files keyed on: comment[data file]
d$plex_ids
#> 1 5 
#> 1 2
d$records[9:11, c("full_path", "plex_id", "channel", "sample_name", "fraction")]
#>                                                  full_path plex_id channel
#> 9  HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01.raw       1    130C
#> 10 HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01.raw       1    131N
#> 11 HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02.raw       1     126
#>    sample_name fraction
#> 9      OSL.41B        1
#> 10        pool        1
#> 11     OSL.53E        2
```

One row per file and **annotated** channel, files in SDRF order and each
file’s channels in reporter m/z order. **Every number is 1-based, as the
SDRF writes it** - `fraction`, `technical_replicate`,
`biological_replicate`, `plex_id` - where
[`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)’s
are 0-based. The split is mzLib’s (`SdrfIsobaricDesign` against
`SpectraFileInfo`), and mzLibR converts neither, so never mix the two
tables without adding or removing 1. `out` names a **folder**: mzLib
writes MetaMorpheus’s `TmtDesign.txt` inside it, with the same 1-based
numbers.

Leave the plex undeclared and mzLib refuses rather than guess one:

``` r

refused <- sdrf_isobaric_design("PXD008841.subset.sdrf.tsv", kit = "TMT10")
refused$is_valid
#> [1] FALSE
writeLines(strwrap(refused$refusals[[1]], 80))
#> No plex source was declared. Name the column that holds the plex, give a
#> file-name pattern whose group captures it, or name the single plex. The plex is
#> never inferred (MAP-09): 'comment[sample preparation batch]' is often not a
#> plex.
```

A declared column is no guarantee either. PXD061609’s
`comment[sample preparation batch]` holds three fractionation arms of
**one** TMTpro plex. Read as the plex, it gives three plexes. Declared
as one plex, each sample sits at the same fraction in three files, which
a plex cannot hold, and mzLib refuses, every case at once:

``` r

arms <- sdrf_isobaric_design("PXD061609.subset.sdrf.tsv", kit = "TMT18",
                             plex_column = "comment[sample preparation batch]")
arms$plex_ids
#> HighSalt  LowSalt    NoSAX 
#>        1        2        3
one <- sdrf_isobaric_design("PXD061609.subset.sdrf.tsv", kit = "TMT18", single_plex = "tube1")
one$is_valid
#> [1] FALSE
length(one$refusals)
#> [1] 36
writeLines(strwrap(one$refusals[[1]], 80))
#> Plex 'tube1' has sample 'brain_female' biorep 1 fraction 1 techrep 1 more than
#> once: 'ea12578.raw' 126, 'ea12603.raw' 126, 'ea12628.raw' 126. Files that share
#> a fraction number need different technical replicates, or different plexes.
```

Which reading is right is the depositor’s to say: mzLib reports, and you
decide. The kit is named as
[`isobaric_kits()`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
takes it, and an unknown name is the same usage error there.

## Instead of `readers_read_records()`

[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
recognises SDRF too, but joins each row into one semicolon-separated
string, and SDRF’s own `NT=...;AC=...` grammar puts semicolons inside
cells. Use
[`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md).

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3):
  SDRF-Proteomics sample columns (characteristics, factor value) (Dai et
  al., Nat Commun 12, 5854, 2021).
