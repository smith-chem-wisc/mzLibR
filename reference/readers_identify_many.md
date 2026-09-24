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

  Paths identified at once. `1`, the default, or `-1` for one per core.
  The answer does not depend on it.

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
`5db922d4cfe1`) by `scripts/build-man.R`; the spec owns these facts, and
all three bindings render the same ones.

- [`FileReader.ReadResultFile`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/FileReader.cs)
  in `mzLib/Readers/FileReader.cs` at mzLib `23c2490e`

- [`SupportedFileTypeExtensions.ParseFileType`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/Readers/Util/SupportedFileTypes.cs)
  in `mzLib/Readers/Util/SupportedFileTypes.cs` at mzLib `23c2490e`

## Parameters: units, ranges and defaults

This verb takes no parameters.

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
