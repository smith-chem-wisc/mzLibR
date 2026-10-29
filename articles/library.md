# Spectral libraries

A spectral library says, for each **precursor** (a peptide or proteoform
at one charge), where a search should expect to find it: its m/z, its
retention time, its ion mobility, and the fragment ions its MS2 spectrum
should hold, with their relative intensities. MetaMorpheus and mzLib
write libraries in mzLib’s binary `.msl` format, predicted (Prosit,
through Koina) or built from a search. mzLibR reads them with mzLib’s
own reader, so a library reads in R exactly as MetaMorpheus reads it.

| question | function | mzLib |
|----|----|----|
| what precursors does this library hold? | [`library_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md) | `MslFileTypeHandler.Open`, `MslLibrary.GetEntry` |
| is this retention time iRT or minutes? | the `retention_time_unit` column | `MslLibraryEntry.RtIsCalibrated` |
| how many targets, decoys and entrapment precursors? | the `decoy_contaminant_target` column | `DecoyContaminantTargetLabel.For` |
| what are a precursor’s fragment ions? | `library_read(fragments = TRUE)` | `MslFragmentIon` |
| can I read a library too large for memory? | `library_read(index_only_threshold = )` | `MslLibrary.LoadIndexOnly` |
| several libraries into one table | [`library_read_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md) | the same reader, in one bridge process |
| where should this library precursor elute in my run, and back? | [`library_calibrate_rt()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_calibrate_rt.md) | `IrtCalibration.Fit`, `IrtCalibrationModel` (#1387) |

Two small libraries are used. `myPrositLib.msl` is mzLib’s own test
library `myPrositLib.msp` (five Prosit-predicted precursors) converted
to `.msl` by mzLib. `synthetic_library.msl` holds six precursors written
to cover every field: a target, its decoy, an entrapment target whose
retention time is in minutes, an entrapment decoy, a predicted precursor
with nothing measured, and a proteoform. Both are in pyMzLib’s
`pkg/python/tests/fixtures/library/`.

## Every precursor

``` r

lib <- library_read("synthetic_library.msl")
lib
#> <mzlibr_library_records> E:\CodeReview\wt-py-library-read\pkg\python\tests\fixtures\library\synthetic_library.msl: 6 precursors (4 target, 2 decoy), .msl format 3
#>   precursor m/z 406.242 to 952.493
#>   6 precursor(s) returned
#>   ! retention_time is in iRT and minutes: read retention_time_unit on each row
c(lib$record_count, lib$target_count, lib$decoy_count, lib$proteoform_count)
#> [1] 6 4 2 1
lib$records[, c("precursor_index", "full_sequence", "charge_state", "precursor_mz",
                "decoy_contaminant_target")]
#>   precursor_index     full_sequence charge_state precursor_mz
#> 1               1          PEPTIDEK            2     465.7400
#> 2               2          KEDITPEP            2     465.7401
#> 3               3          PEPTLDEK            2     458.7462
#> 4               4          KEDLTPEP            2     458.7463
#> 5               5       ELVISLIVESK            3     406.2420
#> 6               6 MAHVRSGKHQIPAEALQ            2     952.4930
#>   decoy_contaminant_target
#> 1                        T
#> 2                        D
#> 3                       ET
#> 4                       ED
#> 5                        T
#> 6                        T
```

One row per precursor, **in the library’s stored order**.
`precursor_index` is the position in the file, 1-based in R as every
position is (mzLib’s `MslLibrary.GetEntry` takes it from 0). It is not
m/z order: mzLib’s own `GetAllEntries` enumerates by ascending m/z, so
sort on `precursor_mz` if that is the order you want.

### What each column is

The names are mzLib’s `MslLibraryEntry` properties. Where mzLib stores a
value it documents as “not available”, the column holds `NA`:

| column | what it is | `NA` means |
|----|----|----|
| `precursor_index` | position in the library, 1-based | never `NA` |
| `full_sequence`, `base_sequence` | mzLib bracket notation, e.g. `PEPTM[Common Variable:Oxidation on M]IDE`; the bare residues | never `NA` |
| `precursor_mz` | precursor m/z, (M + zH)/z | never `NA` |
| `charge_state` | precursor charge | never `NA` |
| `retention_time`, `retention_time_unit` | the time, and whether it is `"iRT"` or `"minutes"` | never `NA` |
| `ion_mobility` | 1/K0, in V s/cm2 | mzLib stored 0.0, “not available” |
| `is_decoy`, `is_entrapment`, `decoy_contaminant_target` | target, decoy and entrapment flags, and mzLib’s label `T`, `D`, `ET` or `ED` | never `NA` |
| `protein_accession`, `protein_name`, `gene_name` | the source protein | mzLib stored an empty string: no protein |
| `is_proteotypic` | the library writer’s proteotypic flag | never `NA` |
| `q_value` | the library’s q-value for the precursor | mzLib stored NaN: no q-value, as in a predicted library |
| `elution_group_id` | shared by every precursor of one base sequence | never `NA` |
| `molecule_type` | `Peptide`, `Proteoform`, `Oligonucleotide` or `Glycopeptide` | never `NA` |
| `source_type` | `Predicted`, `Empirical`, or `EmpiricalRefined` (predicted intensities, measured time or mobility) | never `NA` |
| `dissociation_type` | e.g. `HCD`; `Unknown` when the library does not say | never `NA` |
| `nce` | normalized collision energy, a whole number | mzLib stored 0: unknown |

The empirical target shows the values filled in, and the predicted
precursor shows them “not available”:

``` r

lib$records[c(1, 5), c("ion_mobility", "q_value", "nce", "protein_accession", "source_type")]
#>   ion_mobility q_value nce protein_accession source_type
#> 1         0.95   0.001  28            P02768   Empirical
#> 5           NA      NA  NA              <NA>   Predicted
```

A whole column can have no source at all. `myPrositLib.msl` has no
protein table and no ion mobility, so those columns are named in
`absent_fields` and are `NA` in every row:

``` r

prosit <- library_read("myPrositLib.msl")
prosit$absent_fields
#> [1] "ion_mobility"      "protein_accession" "protein_name"     
#> [4] "gene_name"
prosit$records$nce
#> [1] NA NA NA NA NA
prosit$records$dissociation_type[[1]]
#> [1] "Unknown"
```

The `.msp` it was converted from records `Collision_energy=30.0`, yet
`nce` is `NA`. That is not a reading error. mzLib converts a text
library through `LibrarySpectrum`, which has no field for NCE,
dissociation type, q-value, ion mobility or protein, so
`MslLibraryEntry.FromLibrarySpectrum` writes mzLib’s “not available”
value for each. A converted library has none of them.

## Retention time: iRT or minutes, row by row

A library stores each precursor’s retention time either as an **indexed
retention time** (iRT, a unitless scale anchored on reference peptides)
or as **minutes** in one particular run. Which one is a flag on each
precursor (`MslLibraryEntry.RtIsCalibrated`), and one library can hold
both. So the unit is a column, not a property of the file:

``` r

lib$records[1:4, c("full_sequence", "retention_time", "retention_time_unit")]
#>   full_sequence retention_time retention_time_unit
#> 1      PEPTIDEK           35.2                 iRT
#> 2      KEDITPEP           35.2                 iRT
#> 3      PEPTLDEK           21.4             minutes
#> 4      KEDLTPEP           33.0                 iRT
```

**Never plot or compare `retention_time` without
`retention_time_unit`.** mzLibR converts nothing: iRT becomes minutes
only through a calibration against a run, which is not in the library
(see “Calibrate a run onto library iRT” below). Keep the rows of one
unit before you use the numbers:

``` r

irt <- lib$records[lib$records$retention_time_unit == "iRT", ]
c(irt = nrow(irt), all = nrow(lib$records))
#> irt all 
#>   5   6
```

iRT can be negative; Prosit predicts `-15.77` for one of the
`myPrositLib.msl` precursors. A `retention_time` of exactly 0 deserves
suspicion: mzLib’s conversion from a text library writes 0 when the
spectrum had no retention time, and documents no sentinel for it, so 0
is passed through as a value.

## Targets, decoys and entrapment

**The library is unfiltered.** Decoys and entrapment precursors are rows
like any other. `decoy_contaminant_target` is mzLib’s own label, written
by `DecoyContaminantTargetLabel.For`: `T` target, `D` decoy, `ET`
entrapment target, `ED` entrapment decoy. The `.msl` format has no
contaminant flag, so `C` never appears.

``` r

table(lib$records$decoy_contaminant_target)
#> 
#>  D ED ET  T 
#>  1  1  1  3
lib$records$full_sequence[!lib$records$is_decoy]
#> [1] "PEPTIDEK"          "PEPTLDEK"          "ELVISLIVESK"      
#> [4] "MAHVRSGKHQIPAEALQ"
```

Entrapment is a target for the search: it competes with the real
targets, so `target_count` includes the entrapment targets and
`decoy_count` the entrapment decoys. To count only the real targets,
keep `decoy_contaminant_target == "T"`.

## Fragment ions

With `fragments = TRUE` the table is **one row per fragment ion**, the
precursor’s columns repeated on each row, then the fragment’s own.
`limit` and `offset` still count precursors:

``` r

ions <- library_read("synthetic_library.msl", fragments = TRUE, limit = 1)
c(precursors = ions$returned_count, rows = ions$row_count)
#> precursors       rows 
#>          1          6
ions$records[, c("annotation", "fragment_mz", "fragment_intensity", "residue_position", "neutral_loss")]
#>   annotation fragment_mz fragment_intensity residue_position neutral_loss
#> 1         y1    147.1128               0.40                8      0.00000
#> 2         b2    227.1026               0.30                2      0.00000
#> 3         y2    276.1554               1.00                7      0.00000
#> 4   bIy[2-5]    326.1710               0.10                2      0.00000
#> 5 y3-18.0106    373.2082               0.20                6    -18.01056
#> 6 y4-43.0058    461.2242               0.15               NA    -43.00581
```

What the numbers are:

- **`fragment_intensity` is relative within its precursor, 0 to 1.**
  mzLib’s writer scales each precursor’s most intense fragment to 1, so
  intensities are comparable within a precursor and not across
  precursors.
- **`neutral_loss` is in Da, negative for a loss, and 0 means no loss.**
  It is 0, not `NA`, because 0 Da is the value of “no loss”, so
  `fragment_mz` and `neutral_loss` stay usable together.
- **`annotation` is mzLib’s library rendering**
  (`MslFragmentIon.Annotation`): the loss signed to four decimals and,
  for a charge above 1, the charge appended (`y4^2+`). It differs from
  `Product.Annotation` elsewhere in mzLib, so do not match the two
  strings against each other.
- **An internal ion** (`bIy[2-5]`) has `is_internal` `TRUE`, its
  C-terminal type in `secondary_product_type` and its end residue in
  `secondary_fragment_number`. Both are `NA` on terminal ions.
- **`residue_position` is `NA` for a fragment with a custom neutral
  loss**, like the `-43.0058` above. The `.msl` format reuses that field
  to index its table of custom losses, and mzLib reads it back as 0. The
  loss itself survives. Named losses (water, ammonia, phosphate) keep
  their position.
- **`fragment_count`** on each row is how many fragment ions the
  precursor has. A precursor with none is still one row, with
  `fragment_count` 0 and every fragment cell `NA`.

Fragments are stored, and returned, in m/z order: mzLib’s writer sorts
them. The table is a spectrum to plot as it is:

``` r

plot(ions$records$fragment_mz, ions$records$fragment_intensity, type = "h", lwd = 2,
     xlab = "m/z", ylab = "relative intensity", ylim = c(0, 1.1))
text(ions$records$fragment_mz, ions$records$fragment_intensity, ions$records$annotation,
     pos = 3, cex = 0.8)
```

![Stick plot of one precursor's six library fragment ions, relative
intensity against m/z, each labelled with its
annotation.](library_files/figure-html/unnamed-chunk-8-1.png)

A library has tens of fragment ions per precursor. For a whole library,
write the table to disk with `out` (a `.tsv`) and read that with
[`read.delim()`](https://rdrr.io/r/utils/read.table.html), rather than
returning it as one document.

## Large libraries: index-only

mzLib opens a library one of two ways. **Loaded whole**, every precursor
and every fragment is in memory. **Index-only**, the precursors are in
memory and each precursor’s fragment block is read from disk when it is
asked for. mzLib chooses by file size (`MslFileTypeHandler.Open`): a
library at or above `MslFileTypeHandler.DefaultIndexOnlyThresholdBytes`,
1 GiB in the mzLib this bridge is built from, opens index-only. mzLibR
sends nothing unless you pass `index_only_threshold`, a size in bytes;
`0` opens every library index-only:

``` r

small <- library_read("synthetic_library.msl", index_only_threshold = 0, limit = 2)
small$is_index_only
#> [1] TRUE
small$index_only_threshold
#> [1] 0
small$records$full_sequence
#> [1] "PEPTIDEK" "KEDITPEP"
```

**The rows are the same in both modes**; only memory and speed differ. A
**compressed** library (`is_compressed`) always loads whole: its
fragment offsets point into the decompressed data, so mzLib decompresses
it all and `is_index_only` stays `FALSE` whatever the threshold. And
index-only **with fragments** reads every fragment block from disk
twice, once to count the rows and once to write them; without fragments
it reads no fragment block at all, which is the fast way to list the
precursors of a very large library.

## Many libraries at once

[`library_read_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md)
reads a list of libraries into one long table, in **one** bridge
process. Its first two columns say where each row came from
(`source_index`, 1-based, and `source_path`), and `files` has one row
per library:

``` r

batch <- library_read_many(c("synthetic_library.msl", "myPrositLib.msl"))
batch
#> <mzlibr_library_batch> library read
#>   2 inputs: 2 read, 0 failed (on_error = "fail")
#>   11 rows in `records`, 11 records
batch$files$record_count
#> [1] 6 5
table(batch$records$source_index)
#> 
#> 1 2 
#> 6 5
batch$files$absent_fields
#> [[1]]
#> character(0)
#> 
#> [[2]]
#> [1] "ion_mobility"      "protein_accession" "protein_name"     
#> [4] "gene_name"
```

`threads` reads that many libraries at once. The answer is identical at
any value; the default is 1 because a library loaded whole is held in
memory while its rows are written. `on_error = "skip"` records an
unreadable library in its row of `files` and reads the rest. Never loop
[`library_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md)
over files: each call starts a bridge process.

## Calibrate a run onto library iRT

A library’s iRT is a run-independent scale (Escher et al. 2012): it says
that one peptide elutes after another, not at which minute. Each run has
its own column, gradient and dead volume, so the same iRT falls at a
different minute in every run. A **calibration** is the map between the
two for one run, fitted from **anchors**: precursors you identified
confidently in that run, each with the minute it was seen and its
library iRT. With it you can ask where any library precursor should
elute in that run, or put the run’s own times on the library’s scale.

| You want to | Call | mzLib does it with |
|----|----|----|
| Fit one run’s map from anchors | `library_calibrate_rt(anchors)` | `IrtCalibration.Fit` |
| Predict the minute a library precursor elutes at | a query with `irt` and `rt_minutes = NA` | `IrtCalibrationModel.ToRtMinutes` |
| Put an observed minute on the library’s iRT scale | a query with `rt_minutes` and `irt = NA` | `IrtCalibrationModel.ToIrt` |
| Know whether a time lies outside the anchors | the `is_extrapolated` column | `IrtCalibrationModel.IsExtrapolated` |
| See how well the run follows the library | `groups$residual_sd`, and `irt - fitted_irt` per anchor | `IrtCalibrationModel.ResidualSd`, `ToIrt` |
| Calibrate several runs in one call | one group per run | one `IrtCalibration.Fit` per group |
| Choose a straight line, or tune the curve | `kind`, `bandwidth`, `knots`, `robustness_iterations`, `minimum_anchors` | `IrtCalibrationOptions` |

**You choose the anchors.** mzLib asks for the observed apex time of
confidently identified precursors and is robust to some wrong ones, but
it applies no q-value cutoff, and neither does mzLibR. **Send the
queries in the same call as their anchors**: mzLib keeps the fitted
curve’s points private, so a calibration cannot be returned and used
later. Refitting the same anchors gives the same map, because mzLib
sorts them before it fits.

### One run, from a DIA-NN report

DIA-NN reports both numbers for every precursor it identifies: `RT`, the
minute it was seen, and `iRT`, its value in the spectral library. The
file used here is one HeLa run from a DIA-NN 2 report in mzLib’s own
test data, read by
[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
and cut to six columns; it ships in `extdata/calibration/`, with its
provenance. The precursors at a q-value of 0.01 or below are the
anchors; the rest are what we ask about.

``` r

hela <- read.delim(system.file("extdata", "calibration", "diann2_hela_rt.tsv", package = "mzLibR"),
                   stringsAsFactors = FALSE)
confident <- hela$q_value <= 0.01
others <- hela[!confident, ]
c(anchors = sum(confident), others = nrow(others))
#> anchors  others 
#>     418     114
anchors <- data.frame(group = hela$run[confident], rt_minutes = hela$retention_time[confident],
                      irt = hela$indexed_retention_time[confident])
queries <- rbind(
  data.frame(group = others$run, rt_minutes = NA, irt = others$indexed_retention_time),
  data.frame(group = others$run, rt_minutes = others$retention_time, irt = NA)
)
calibration <- library_calibrate_rt(anchors, queries)
calibration
#> <mzlibr_rt_calibration> lowess: 1 of 1 group(s) fitted from 418 anchors, 228 queries converted
#>   Fig2HeLa-0-5h_MHRM_R01_T0: lowess, 418 anchors over 15.99 to 40.9 min, residual SD 5 iRT
run <- calibration$groups
c(first = round(run$first_anchor_minutes, 1), last = round(run$last_anchor_minutes, 1))
#> first  last 
#>  16.0  40.9
```

Each anchor is a row with its `group` (the run), the minute it was seen
and its library iRT. A query gives one time and `NA` for the other, and
the `NA` is what the calibration fills in. With no options, mzLib fits
its default, a robust LOWESS curve (Cleveland 1979), so a gradient that
is not a straight line is followed.

The rows come back in the order you sent them: the anchors, then each
query. A `"to_minutes"` row predicts where a library precursor should
elute in this run. For the lower-confidence precursors, compare the
prediction with where DIA-NN saw them:

``` r

rows <- calibration$records
predicted <- rows$rt_minutes[rows$role == "to_minutes"]
round(median(abs(predicted - others$retention_time)), 2)
#> [1] 0.53
round(median(abs(others$predicted_retention_time - others$retention_time)), 2)
#> [1] 0.38
```

Half of the predictions are within about half a minute of the observed
time. DIA-NN’s own prediction for the same precursors, from its own
calibration of this run, is a little closer. A search uses this map to
decide where in the run to look for each library precursor, so the
prediction error sets how wide that retention-time window must be.

A `"to_irt"` row is the other direction: the run’s observed minute on
the library’s scale. That is how runs with different gradients are
compared, on one axis.

``` r

on_library_scale <- rows$irt[rows$role == "to_irt"]
round(median(abs(on_library_scale - others$indexed_retention_time)), 1)
#> [1] 3.4
```

### How good is the fit

`residual_sd` is mzLib’s robust spread of the anchors around the curve,
in **iRT units**, with the anchors that do not fit set aside. On the
anchor rows, `fitted_irt` is the curve’s iRT at the anchor’s minute, so
`irt - fitted_irt` is each anchor’s own residual. To read the spread in
minutes, divide by how many iRT units the run covers per minute:

``` r

anchor_rows <- rows[rows$role == "anchor", ]
per_minute <- diff(range(anchor_rows$fitted_irt)) / (run$last_anchor_minutes - run$first_anchor_minutes)
round(c(residual_sd = run$residual_sd, irt_per_minute = per_minute,
        minutes = run$residual_sd / per_minute), 2)
#>    residual_sd irt_per_minute        minutes 
#>           5.00           8.16           0.61
sum(abs(anchor_rows$irt - anchor_rows$fitted_irt) > 3 * run$residual_sd)
#> [1] 9
```

About 0.6 minutes of spread in a run where the anchors span 25 minutes.
Nine anchors sit more than three robust SDs from the curve: those are
the wrong identifications the robust fit set aside, and they did not
pull the curve toward them.

Beyond the first and last anchor the map is a straight line continued
from the curve’s ends, not a fit. `is_extrapolated` marks those rows,
and a caveat counts them:

``` r

sum(rows$is_extrapolated[rows$role != "anchor"])
#> [1] 5
calibration$caveats[length(calibration$caveats)]
#> [1] "5 queries fall outside the anchors' span, where the map is extrapolated (is_extrapolated)."
```

### Several runs, and when mzLib refuses one

Give every anchor and query its run as the group: each group is fitted
on its own anchors, and its queries use that fit. These are three
24-minute PBMC runs from a DIA-NN 1 report in mzLib’s test data. One
precursor, `AAAAAAAAGAFAGR` at charge 2, is held out of the anchors so
that each run can predict where it should elute:

``` r

pbmc <- read.delim(system.file("extdata", "calibration", "diann_pbmc_rt.tsv", package = "mzLibR"),
                   stringsAsFactors = FALSE)
held <- pbmc$precursor_id == "AAAAAAAAGAFAGR2"
anchors <- data.frame(group = pbmc$run[!held], rt_minutes = pbmc$retention_time[!held],
                      irt = pbmc$indexed_retention_time[!held])
queries <- data.frame(group = pbmc$run[held], irt = 11.886096)
refused <- library_calibrate_rt(anchors, queries, on_error = "skip")
data.frame(run = substr(refused$groups$group, 1, 8), anchors = refused$groups$anchor_count,
           error = refused$groups$error_message)
#>        run anchors                                          error
#> 1 20240316      13 Calibration needs at least 20 anchors; got 13.
#> 2 20240317      11 Calibration needs at least 20 anchors; got 11.
#> 3 20240315      10 Calibration needs at least 20 anchors; got 10.
```

mzLib refuses to fit a run from fewer than 20 anchors by default rather
than return a curve that a handful of points would decide. Without
`on_error = "skip"` the first refusal raises an `mzlib_usage_error`
naming the run. These runs have ten to thirteen anchors each. That is
too few for a LOWESS curve, which fits each point from a share of the
anchors, but enough for a straight line: say so, and lower mzLib’s
floor.

``` r

lines <- library_calibrate_rt(anchors, queries, kind = "linear", minimum_anchors = 10)
data.frame(run = substr(lines$groups$group, 1, 8), kind = lines$groups$kind,
           residual_sd = round(lines$groups$residual_sd, 1))
#>        run   kind residual_sd
#> 1 20240316 linear         4.1
#> 2 20240317 linear         3.6
#> 3 20240315 linear         0.7
predicted <- lines$records[lines$records$role == "to_minutes", ]
data.frame(run = substr(predicted$group, 1, 8), predicted = round(predicted$rt_minutes, 2),
           seen = round(pbmc$retention_time[held], 2))
#>        run predicted seen
#> 1 20240316      9.22 9.21
#> 2 20240317      9.45 9.23
#> 3 20240315      9.17 8.94
```

Each run’s line predicts the held-out precursor within a quarter of a
minute of where it was seen. With few anchors, prefer a line: a LOWESS
curve through ten points can pass through nearly all of them, and a
`residual_sd` near 0 then means an overfit, not a perfect run.

## Traps

- **Stored order, not m/z order.** Rows follow `precursor_index`. Sort
  on `precursor_mz` yourself.
- **The extra digits are the float.** `precursor_mz`, `retention_time`
  and `ion_mobility` are stored as 32-bit floats and widened by mzLib,
  so 465.74 reads back as `465.739990234375`. Round them for display;
  they carry no more precision than about seven significant digits.
- **`nce` is a whole number.** The format stores NCE x 10 and mzLib
  reads it back with integer division, so an NCE of 27.5 reads as 27.
- **Converted libraries are thin.** A library converted from `.msp` has
  no NCE, dissociation type, q-value, ion mobility or protein (see
  above).
- **A broken file is an error, not a short table.** mzLib checks the
  whole file (magic number, version, footer, precursor count and a
  CRC-32 checksum) before it returns a single precursor. A damaged
  library raises `mzlib_bridge_error` with `error_type`
  `"InvalidDataException"`, and a file that is not a library raises it
  with `"FormatException"`. A file whose name does not end in `.msl` is
  an `mzlib_usage_error` before anything is opened.

## What this does not do

- **Read a text library** (`.msp`, pDeep, ms2pip). mzLib can turn one
  into a library in memory (`MslConverter.FromMspFile`), but no wire
  verb projects it yet.
- **Write a library.** mzLib writes them (`MslLibrary.Save`); no wire
  verb projects that yet.
- **Query by m/z or retention-time window.** mzLib’s `QueryMzWindow` and
  `QueryWindow` are what a search engine calls. To find a precursor
  here, read the table and filter it.
- **Choose calibration anchors.**
  [`library_calibrate_rt()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_calibrate_rt.md)
  fits the anchors you send. A q-value cutoff, one anchor per precursor,
  apex or MS2-scan time: each is your choice, and mzLib makes none.
- **Keep a calibration for later.** mzLib keeps the fitted curve’s
  points private (`IrtCalibrationModel`), so nothing can carry it out of
  the call. Send the queries with their anchors, or send the same
  anchors again: the fit is deterministic.
- **Write calibrated times into a library.** mzLib’s
  `MslLibrary.WithCalibratedRetentionTimes` takes a slope and intercept
  and is not connected to `IrtCalibration`, and no wire verb projects
  it.

mzLib \#1383 added the entrapment flag and made `rt_is_calibrated` work,
\#1421 added `GetEntry` without fragments, and \#1420 made libraries
over 2 GB open.

## What to cite

The methods and resources behind the functions on this page, from their
specs:

- [doi:10.1002/pmic.201100463](https://doi.org/10.1002/pmic.201100463):
  the iRT concept: a normalized, run-independent retention time scale
  that a run is mapped onto by a calibration (Escher et al., Proteomics
  12:1111, 2012).
- [doi:10.1080/01621459.1979.10481038](https://doi.org/10.1080/01621459.1979.10481038):
  robust LOWESS with tricube local lines and bisquare reweighting at 6 x
  the median absolute residual, the curve IrtCalibrationKind.Lowess fits
  (Cleveland, J Am Stat Assoc 74:829, 1979).
