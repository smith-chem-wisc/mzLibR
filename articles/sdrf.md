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

## Instead of `readers_read_records()`

[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
recognises SDRF too, but joins each row into one semicolon-separated
string, and SDRF’s own `NT=...;AC=...` grammar puts semicolons inside
cells. Use
[`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md).
