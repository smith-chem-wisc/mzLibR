# Identify a result file without parsing its contents

Say what kind of file a path is - the mzLib file type, its extension,
its reader and the cross-format views it offers - without parsing its
contents.

Cheap by design: mzLib resolves the type and stops, so identifying a
million-row file costs no more than identifying an empty one.

## Usage

``` r
readers_identify(path, timeout = 60)
```

## Arguments

- path:

  Path to a result file.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Details

It is not, however, \*pure\*. mzLib disambiguates a bare `.tsv` by
reading its first line, a `.mztab` by its first five, and a Bruker `.d`
by which analysis file the directory holds, so an unreadable file raises
rather than returning a guess.

## Value

An `mzlibr_file_info`: `path`, `file_type`, `extension`, `reader`,
`views` and `is_quantifiable`.

## Ask this before quantifying

`is_quantifiable` is the precondition for
[`flashlfq_quantify`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
— when `FALSE`, mzLib can still read the file, it simply has no uniform
view, and quantification would fail on it.

**But `TRUE` is not permission.** It reports what mzLib's \*interface\*
offers, not that the numbers are comparable. `MsFraggerPsm` is
quantifiable by interface and should not be quantified: among other
things its retention times are in seconds while MetaMorpheus's are in
minutes, and mzLib does not normalise them. See
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
and the `caveats` it returns.

## Wraps

Wire verb `readers identify`. Generated from the bridge's verb spec
`readers.identify.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`FileReader.ReadResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`SupportedFileTypeExtensions.ParseFileType`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Util/SupportedFileTypes.cs)
  in `mzLib/Readers/Util/SupportedFileTypes.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `path`:

  path; required. Any file (or Bruker .d directory). For many paths see
  the caveat on `--paths-stdin`.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `path`:

  string; never `NA`. Absolute path identified.

- `file_type`:

  string; never `NA`. The mzLib SupportedFileType.

- `extension`:

  string; `NA` when mzLib maps no extension to this type. The extension
  mzLib dispatched on.

- `reader`:

  string; never `NA`. The mzLib reader class that would parse it.

- `views`:

  string\[\]; never `NA`. The cross-format views it offers; often empty.

- `error`:

  object; `NA` when always, for a single path; in a paths-stdin entry,
  the path was identified. {kind, type, message} for a path that failed
  under on-error skip.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).
A condition that mentions `paths-stdin`, `threads` or `on-error` belongs
to the verb's many-files form.

- `mzlib_usage_error` (usage):

  path missing, blank or not found; mzLib does not recognise the file
  type

- `mzlib_usage_error` (usage):

  paths-stdin with path, offset or limit; paths-stdin with a value; no
  paths on stdin; a repeated path; threads 0 or below -1; on-error not
  fail or skip; on-error skip with a single path; out naming an input

- `mzlib_usage_error` (usage):

  with on-error fail, an input that is missing or not this view: the
  usage error of that input, prefixed 'Input \<i\> (\<path\>)'

- `mzlib_usage_error` (usage):

  out with paths-stdin: identify returns no table

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Cheap but not pure: mzLib reads a bare .tsv's first line, a .mztab's
  first five and a .d directory's listing to dispatch them.

- No software field: mzLib's Software property is Unspecified for
  everything its factory builds; file_type already names the tool.

- Many paths: `--paths-stdin` (newline-delimited paths on stdin) with
  `--threads` N (default 1, -1 = one per core; the answer is identical
  at any value) and `--on-error` fail\|skip answers {file_count,
  read_count, failed_count, on_error, files}, files\[i\] being path i's
  answer in input order with error set when it could not be identified
  under skip (readers_many_identify.json). There is no table and no
  `--out`.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.identify`; many files:
  `pymzlib.readers.identify_many`

- Rust (mzLibRust): `mzlib::readers::identify`; many files:
  `mzlib::readers::identify_many`

- R (mzLibR): `readers_identify`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The \_many spellings (py, rust, r) follow the 2026-09-23 cross-binding
  decision; Rust and R names here are intended, not yet ported, and
  since.pymzlib stays null until release.

- The `--paths-stdin` answer is not declared as result.bulk because
  check_verbs.py requires a record_count there (BULK.md section 2),
  which has no meaning for identify. Its shape is in the caveats and in
  readers_many_identify.json until the checker allows a table-less bulk
  verb.

## See also

[`readers_formats`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md),
[`readers_read_results`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)

## Examples

``` r

info <- readers_identify("PXD078927_msgf_1_1_0.mzid")
info
#> <mzlibr_file_info> C:\Users\trish\AppData\Local\Temp\claude\E--CodeReview-bridge\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\scratchpad\wt_c1\code\mzLib\mzLib\Test\DataFiles\PXD078927_msgf_1_1_0.mzid
#>   MzIdentML (.mzid), read by MzIdentMLResultFile
#>   views: spectral_match
info$views
#> [1] "spectral_match"
```
