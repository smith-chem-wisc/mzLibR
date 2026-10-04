# Convert full sequences to another notation with mzLib, one row per input

Convert full sequences from one notation to another with mzLib's
sequence conversion service, one result per input: for example a
MetaMorpheus full sequence to Unimod accessions.

Hands each sequence to mzLib's `SequenceConversionService` (with
ProForma registered) and returns what mzLib made of it. The main use is
MetaMorpheus or mzLib full sequences to Unimod accessions:
`[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]SEQK`. Nothing is
converted in R or in the bridge; every output and every status is
mzLib's.

## Usage

``` r
peptidoform_convert(sequences, source = "mzLib", target = "Unimod", mode = "ReturnNull",
  threads = 1, timeout = 120)
```

## Arguments

- sequences:

  A character vector of full sequences, one result row each, in order,
  duplicates included, sent to mzLib exactly as given. An `NA`, a blank
  entry or one with a line break is refused, because it would shift
  every later row. Split an ambiguous MetaMorpheus full sequence
  (`|`-joined candidates) first: mzLib joins the candidates into one
  sequence, with only a warning.

- source:

  The notation the sequences are in: one of mzLib's registered source
  formats (the result's `source_formats`), matched case-insensitively.
  Sent as the wire's `--from`.

- target:

  The notation to write: one of mzLib's registered target formats (the
  result's `target_formats`), matched case-insensitively. Sent as the
  wire's `--to`.

- mode:

  mzLib's `SequenceConversionHandlingMode`, case-insensitive:
  `"ReturnNull"` (a sequence mzLib cannot convert is a `"failed"` row),
  `"RemoveIncompatibleElements"` or `"UsePrimarySequence"` (what the
  target cannot write is dropped and the row is
  `"converted_with_warnings"`), or `"ThrowException"` (the first such
  sequence, in input order, fails the whole call with an
  `mzlib_usage_error` naming it, and nothing is returned).

- threads:

  Sequences converted at once, or `-1` for every core. The rows are
  identical, in input order, at any value.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

**Choose Unimod, not ProForma, for UniProt-sourced modifications.**
mzLib's ProForma target does not resolve them and writes them back under
their mzLib name, with status `"converted"` (mzLib#1401). In ProForma
output, treat any bracket that is not a `UNIMOD:` term as unresolved.
The `pro_forma` column
[`readers_read_records`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
gives a `.psmtsv` comes from the same serializer and has the same gap.

## Value

An `mzlibr_sequence_conversions`:

\- `records`, a data.frame with one row per input, in input order:
`input`, the sequence as sent; `output`, mzLib's converted sequence,
`NA` when the row failed; `status`, `"converted"` (an output and nothing
recorded against it), `"converted_with_warnings"` (an output, but mzLib
recorded a warning, an error or an incompatible item) or `"failed"` (no
output); `failure_reason`, mzLib's `ConversionFailureReason` or `NA`
when it recorded none - which a failed row can be: the Unimod target
under `"ReturnNull"` names only the incompatible items; and three list
columns of character vectors, `character(0)` when empty:
`incompatible_items` (what the target could not write, as mzLib
describes it, `"Made Up:Not a modification on K @3(K)"`), `warnings` and
`errors`. - `source_format` and `target_format`, spelled as mzLib
registered them; `mode`, the handling mode used; `source_formats` and
`target_formats`, every notation mzLib has registered. - `record_count`
sequences, of which `converted_count`, `warned_count` and `failed_count`
have each status; `column_names`; and `caveats`, what a status does and
does not promise.

The failed and warned rows are
`records[records$status != "converted", ]`.

## Wraps

Wire verb `peptidoform convert`. Generated from the bridge's verb spec
`peptidoform.convert.yaml` (bridge commit `e76157831b15`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`SequenceConversionService.Default`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/SequenceConversion/SequenceConversionService.cs)
  in `mzLib/Omics/SequenceConversion/SequenceConversionService.cs` at
  mzLib `0a808fec`

- [`SequenceConversionService.GetConverter`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/SequenceConversion/SequenceConversionService.cs)
  in `mzLib/Omics/SequenceConversion/SequenceConversionService.cs` at
  mzLib `0a808fec`

- [`ISequenceConverter.Convert`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/SequenceConversion/SequenceConverter.cs)
  in `mzLib/Omics/SequenceConversion/SequenceConverter.cs` at mzLib
  `0a808fec`

- [`ProFormaSequenceConversion.RegisterWithDefault`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Readers/ProForma/SequenceConversion/ProFormaSequenceConversion.cs)
  in
  `mzLib/Readers/ProForma/SequenceConversion/ProFormaSequenceConversion.cs`
  at mzLib `0a808fec`

- [`SequenceConversionHandlingMode`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/SequenceConversion/Util/SequenceConversionHandlingMode.cs)
  in
  `mzLib/Omics/SequenceConversion/Util/SequenceConversionHandlingMode.cs`
  at mzLib `0a808fec`

- [`ConversionWarnings`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/Omics/SequenceConversion/Util/ConversionWarnings.cs)
  in `mzLib/Omics/SequenceConversion/Util/ConversionWarnings.cs` at
  mzLib `0a808fec`

## Parameters: units, ranges and defaults

- `source` (wire `--from`):

  string; default `mzLib`; range
  `a registered source format: MassShift, Modomics, ProForma, mzLib (at 0a808fec)`.
  The notation the sequences are in, matched case-insensitively against
  the service's AvailableSourceFormats. Any other name, or the option
  without a value, is a usage error that lists the registered names.

- `target` (wire `--to`):

  string; default `Unimod`; range
  `a registered target format: Chronologer, Essential, MassShift, ProForma, Unimod, mzLib (at 0a808fec)`.
  The notation to write, matched case-insensitively against the
  service's AvailableTargetFormats. Any other name, or the option
  without a value, is a usage error that lists the registered names.

- `mode`:

  string; default `ReturnNull`; range
  `ThrowException | ReturnNull | RemoveIncompatibleElements | UsePrimarySequence`.
  mzLib's SequenceConversionHandlingMode, by name, case-insensitive.
  ReturnNull: an input mzLib cannot convert is a failed row.
  RemoveIncompatibleElements and UsePrimarySequence: what the target
  cannot write is dropped and the row is converted_with_warnings.
  ThrowException: the first input mzLib cannot convert (in input order)
  fails the whole call as a usage error naming it. A number or any other
  name is a usage error.

- `threads`:

  int; default `1`; range `>= 1, or -1 for every core`. Sequences
  converted at once. Rows are identical, in input order, at any value.

- `sequences`:

  string\[\]; required; range `non-empty`. Read from stdin (not an
  option): one full sequence per line, given to mzLib exactly as read.
  Blank lines are skipped. One row per line, in order, duplicates
  included. No sequences is a usage error.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `source_format`:

  string; never `NA`. `--from`, spelled as mzLib registered it.

- `target_format`:

  string; never `NA`. `--to`, spelled as mzLib registered it.

- `mode`:

  string; never `NA`. The SequenceConversionHandlingMode used.

- `source_formats`:

  string\[\]; never `NA`. Every source format the service has
  registered, ordinal order.

- `target_formats`:

  string\[\]; never `NA`. Every target format the service has
  registered, ordinal order.

- `record_count`:

  int; in **sequences**; never `NA`. Rows: the input sequences.

- `converted_count`:

  int; in **sequences**; never `NA`. Rows with status converted.

- `warned_count`:

  int; in **sequences**; never `NA`. Rows with status
  converted_with_warnings.

- `failed_count`:

  int; in **sequences**; never `NA`. Rows with status failed.

- `column_names`:

  string\[\]; never `NA`. The table's columns in order.

- `records` (wire `columns`):

  table; never `NA`. Column name to one value per input, in input order.

- `caveats`:

  string\[\]; never `NA`. Six static caveats.

## Columns of `records`

- `input`:

  string; never `NA`. The line exactly as read from stdin.

- `output`:

  string; `NA` when mzLib returned no output for this input (status
  failed). mzLib's converted sequence.

- `status`:

  string; never `NA`. converted (an output and
  ConversionWarnings.IsClean) \| converted_with_warnings (an output, but
  mzLib recorded a warning, error or incompatible item) \| failed (no
  output).

- `failure_reason`:

  string; `NA` when mzLib recorded no ConversionFailureReason: always
  for converted rows, and for a row the Unimod target fails under
  ReturnNull. mzLib's ConversionFailureReason: InvalidSequence \|
  IncompatibleModifications \| UnsupportedDirection \| UnknownFormat.

- `incompatible_items`:

  string\[\]; never `NA`. ConversionWarnings.IncompatibleItems: what the
  target could not write, as mzLib describes it ('Made Up:Not a
  modification on K @3(K)'). Empty when none.

- `warnings`:

  string\[\]; never `NA`. ConversionWarnings.Warnings: mzLib's non-fatal
  messages for this input.

- `errors`:

  string\[\]; never `NA`. ConversionWarnings.Errors: mzLib's error
  messages for this input.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  no sequences on stdin; `--from` or `--to` not a registered format, or
  given without a value (the message lists the registered names);
  `--mode` not a SequenceConversionHandlingMode name; `--threads` 0 or
  below -1

- `mzlib_usage_error` (usage):

  `--mode` ThrowException and mzLib could not convert an input (the
  first in input order; the message names it and mzLib's reason)

- `mzlib_bridge_error` (correctness):

  mzLib threw something other than its SequenceConversionException
  (never seen at 0a808fec)

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- ProForma does not resolve UniProt-sourced modifications at 0a808fec:
  '\[UniProt:N-acetylserine on S\]SEQK' comes out as
  '\[UniProt:N-acetylserine on S\]-SEQK', status converted, not
  '\[UNIMOD:1\]-SEQK'. mzLib's ProForma serializer looks modifications
  up in the MetaMorpheus list only (ProFormaSequenceSerializer.cs:34,
  MzLibModificationLookup.cs:54); the Unimod target looks them up in
  every list mzLib loads, UniProt's included
  (GlobalModificationLookup.cs:23). Convert to Unimod when the sequences
  carry UniProt modifications. mzLib#1401.

- ProForma writes a modification it cannot resolve under its mzLib name
  and the row is still converted, so status converted does not mean
  every bracket in a ProForma output is a controlled-vocabulary term.
  Unimod either names a modification by its UNIMOD accession or reports
  it in incompatible_items.

- A modification the target cannot write: under ReturnNull the row fails
  (output null); under RemoveIncompatibleElements and UsePrimarySequence
  it is dropped from output, named in incompatible_items, and the row is
  converted_with_warnings; under ThrowException the first such input
  fails the whole call as a usage error.

- status is mzLib's verdict read off ConversionWarnings; failure_reason
  can be null on a failed row (the Unimod target under ReturnNull
  records only incompatible_items).

- An ambiguous MetaMorpheus full sequence ('\|'-joined candidates) is
  not refused: mzLib's parser skips each '\|' with a warning and joins
  the candidates into one sequence, status converted_with_warnings.
  Split ambiguous full sequences before converting them.

- Each stdin line is one input, given to mzLib exactly as read; blank
  lines are skipped; rows are in input order at any `--threads`.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.peptidoform.convert`

- Rust (mzLibRust): `mzlib::peptidoform::convert_with` with
  `ConvertOptions`

- R (mzLibR): `peptidoform_convert`

## Since

Wire protocol 1; pyMzLib 0.4.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The ProForma target's UniProt gap is an mzLib defect (the serializer's
  default lookup), mzLib#1401; close this when a pin carries its fix. A
  bridge-side substitute
  (ProFormaSequenceSerializer.WithLookup(GlobalModificationLookup))
  composes serializer internals and is deliberately not used.

- mzLib's MzLibSequenceParser joins '\|'-separated ambiguous candidates
  into one sequence with only a warning; whether it should refuse them
  (as SpectrumMatchFromTsv.ProFormaFromFullSequence does before calling
  it) is an mzLib question: mzLib#1405.

- AutoDetect (SequenceConversionService.ConvertAutoDetect) is not
  projected: `--from` is always explicit.

- Ported, unreleased: mzLibRust \#32 and mzLibR \#32 (drafts until
  pyMzLib 0.4.0 publishes the bridge). Set since.mzlibrust /
  since.mzlibr when they release.

## References

- [doi:10.1002/pmic.200300744](https://doi.org/10.1002/pmic.200300744):
  Unimod, the modification database whose accessions the Unimod target
  writes (Creasy and Cottrell 2004, Proteomics)

- [doi:10.1021/acs.jproteome.1c00771](https://doi.org/10.1021/acs.jproteome.1c00771):
  ProForma 2.0, the ProForma target's notation (LeDuc et al. 2022, J
  Proteome Res)

## See also

[`peptidoform_fragments`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md),
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)

## Examples

``` r

# UniProt and MetaMorpheus modifications to Unimod accessions, and one mzLib cannot map:
sequences <- c(
  "[UniProt:N-acetylserine on S]SEQK",
  "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
  "PEPM[Common Variable:Oxidation on M]K",
  "PEPK[Made Up:Not a modification on K]R"
)
result <- peptidoform_convert(sequences)
result
#> <mzlibr_sequence_conversions> 4 sequences, mzLib to Unimod (ReturnNull)
#>   3 converted, 0 with warnings, 1 failed
result$records[, c("output", "status")]
#>             output    status
#> 1   [UNIMOD:1]SEQK converted
#> 2 PEPK[UNIMOD:36]R converted
#> 3 PEPM[UNIMOD:35]K converted
#> 4             <NA>    failed
failed <- result$records[result$records$status != "converted", ]
failed$incompatible_items[[1]]
#> [1] "Made Up:Not a modification on K @3(K)"
```
