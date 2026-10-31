# Identify many files in one bridge call

Say what kind of file a path is - the mzLib file type, its extension,
its reader and the cross-format views it offers - without parsing its
contents.

The many-files form of
[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md):
one bridge process resolves every path, `threads` at a time, and answers
in input order.

## Usage

``` r
readers_identify_many(paths, threads = 1, on_error = "fail", timeout = 60)
```

## Arguments

- paths:

  A character vector of paths. A Bruker `.d` directory is accepted.

- threads:

  Paths identified at once, in files. `1`, the default, or `-1` for one
  per core. The answer does not depend on it.

- on_error:

  `"fail"`, the default, stops at the first path that cannot be
  identified. `"skip"` records why in that path's row and identifies the
  rest.

- timeout:

  Seconds to allow, or `NULL` to wait indefinitely.

## Value

An `mzlibr_identify_batch`. `files` is a data.frame with one row per
path, in input order: `path`, `file_type`, `extension`, `reader`, a
`views` list column, `is_quantifiable`, and `error_kind`, `error_type`
and `error_message`, which are `NA` for a path that was identified.
`file_count`, `read_count` and `failed_count` count paths.

## Wraps

Wire verb `readers identify` with `--paths-stdin`. Generated from the
bridge's verb spec `readers.identify.yaml` (bridge commit
`84c3b984b6b9`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`FileReader.ReadResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`SupportedFileTypeExtensions.ParseFileType`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Util/SupportedFileTypes.cs)
  in `mzLib/Readers/Util/SupportedFileTypes.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `paths` (wire `--paths-stdin`):

  flag; default absent. Identify MANY paths: newline-delimited paths on
  stdin (blank lines ignored; a repeated path is a usage error).
  Exclusive with path; offset, limit and out are refused with it.

- `threads`:

  int; in **files**; default `1`; range `>= 1, or -1 for one per core`.
  Paths identified at once. The answer is identical at any value (files
  are reported in input order).

- `on_error` (wire `--on-error`):

  string; default `fail`; range `fail | skip`. fail: the first path that
  cannot be identified stops the batch with its usage or correctness
  error, the message starting 'Input \<i\> (\<path\>)'. skip: the
  failure goes in that path's files\[i\].error and the rest are
  identified. skip with a single path is a usage error.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `file_count`:

  int; in **files**; never `NA`. Paths given.

- `read_count`:

  int; in **files**; never `NA`. Paths identified.

- `failed_count`:

  int; in **files**; never `NA`. Paths not identified; non-zero only
  under on-error skip.

- `on_error`:

  string; never `NA`. fail or skip, as requested. `--threads` is
  deliberately not echoed: the envelope is identical at every thread
  count.

- `files`:

  object\[\]; never `NA`. One per input, in input order: the single-path
  answer (path, file_type, extension, reader, views, error). A path not
  identified has error set and file_type, extension and reader null,
  views empty.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

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

- `mzlib_bridge_error` (correctness):

  a .parquet that exists but cannot be opened, because another process
  holds it or the user may not read it: IOException or
  UnauthorizedAccessException from mzLib's detection (since mzLib
  \#1478), not the not-supported usage error

- `mzlib_service_unavailable` (service_unavailable):

  Never raised by this verb.

## Caveats

- Cheap but not pure: mzLib reads a bare .tsv's first line, a .txt's and
  a .mztab's first five, a .parquet's column schema and a .d directory's
  listing to dispatch them. A .parquet mzLib can open and read without
  the DIA-NN report columns (Precursor.Id, Stripped.Sequence, Run), or
  bytes that are not Parquet at all, is a usage error, 'Parquet file
  type not supported' (mzLib \#1384). One it cannot open is an I/O error
  instead (#1478).

- No software field: mzLib's Software property is Unspecified for
  everything its factory builds; file_type already names the tool.

- Many paths (result.bulk): there is no table and no `--out`, and no
  record_count, because nothing is read.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work. This bulk form pays it once for the whole list, which is the
reason it exists: never loop the one-path function over many files.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.readers.identify`; many files:
  `pymzlib.readers.identify_many`

- Rust (mzLibRust): `mzlib::readers::identify`; many files:
  `mzlib::readers::identify_many`

- R (mzLibR): `readers_identify`; many files: `readers_identify_many`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust 0.1.0; mzLibR 0.1.0.

## See also

[`readers_identify`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md),
[`readers_formats`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)

## Examples

``` r

found <- readers_identify_many(
  c("PXD078927_msgf_1_1_0.mzid", "no-such-run.mzML",
    "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"),
  on_error = "skip"
)
found$files[, c("file_type", "is_quantifiable", "error_message")]
#>                             file_type is_quantifiable
#> 1                           MzIdentML           FALSE
#> 2                                <NA>           FALSE
#> 3 MetaMorpheusQuantifiedProteinGroups           FALSE
#>                                                                                                                                                        error_message
#> 1                                                                                                                                                               <NA>
#> 2 File not found: 'C:\\Users\\trish\\AppData\\Local\\Temp\\claude\\E--CodeReview-bridge\\eebe6abc-2c7e-4422-9d4e-f0aaa03b3152\\scratchpad\\wt_c1\\no-such-run.mzML'.
#> 3                                                                                                                                                               <NA>
```
