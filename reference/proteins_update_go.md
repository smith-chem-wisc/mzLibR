# Download the current Gene Ontology release to a go.obo, on purpose

Download the current Gene Ontology release (go.obo) to a path, on
purpose, keeping any different file already there as a timestamped
backup.

Wraps mzLib's `Loaders.UpdateGeneOntology`. The whole go.obo (tens of
megabytes) is streamed from GO's PURL, which always serves the
**current** release. When a file is already at the path, it is kept
beside the new one as `go.obo.<yyyyMMdd-HHmmss-fff>` if the download
differs, and left alone if it is the same, so earlier runs stay
reproducible. A failed download leaves any existing file untouched.

## Usage

``` r
proteins_update_go(go_obo, timeout = NULL)
```

## Arguments

- go_obo:

  Where to write, e.g. `"go.obo"`. Its folder must exist.

- timeout:

  Seconds to allow, or `NULL` (the default) to wait; mzLib itself gives
  up after two minutes without data.

## Details

This is the only function in mzLibR that fetches a go.obo.
[`proteins_annotate_go`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)
never does, so the release a result was computed against is always a
file you chose to keep.

## Value

An `mzlibr_go_update`: `go_obo_file` written; `url` fetched from;
`existed_before`; `previous_sha256` of the file that was there, `NA`
when there was none; `changed`, whether the file on disk is now
different; `go`, the release now at the path (`source_file_name`,
`sha256`, `release`, `term_count`); and `caveats`, which name the backup
kept.

## Wraps

Wire verb `proteins update-go`. Generated from the bridge's verb spec
`proteins.update-go.yaml` (bridge commit `e76157831b15`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`Loaders.UpdateGeneOntology`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/Loaders.cs)
  in `mzLib/UsefulProteomicsDatabases/Loaders.cs` at mzLib `0a808fec`

- [`GeneOntologyGraph.Load`](https://github.com/smith-chem-wisc/mzLib/blob/0a808fec/mzLib/UsefulProteomicsDatabases/GeneOntology/GeneOntologyGraph.cs)
  in `mzLib/UsefulProteomicsDatabases/GeneOntology/GeneOntologyGraph.cs`
  at mzLib `0a808fec`

## Parameters: units, ranges and defaults

- `go_obo` (wire `--go-obo`):

  path; required. Where to write, e.g. go.obo. Its folder must exist.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `go_obo_file`:

  string; never `NA`. Absolute path written.

- `url`:

  string; never `NA`. Loaders.GeneOntologyUrl,
  https://purl.obolibrary.org/obo/go.obo: a moving PURL that always
  serves the current release.

- `existed_before`:

  bool; never `NA`. Whether a file was already at the path.

- `previous_sha256`:

  string; `NA` when no file was there. sha256 of the file that was
  there.

- `changed`:

  bool; never `NA`. Whether the file on disk now differs from what was
  there (true for a first download).

- `go`:

  object; never `NA`. {source_file_name, sha256, release (null when the
  file has no data-version), term_count}, read back with
  GeneOntologyGraph.Load, so a truncated download fails here rather than
  in a later annotate-go.

- `caveats`:

  string\[\]; never `NA`. The backup kept when a different file was
  replaced; a file with no data-version.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  the folder for go-obo does not exist

- `mzlib_service_unavailable` (service_unavailable):

  GO's PURL is unreachable, answers 408/429/5xx, or stalls for
  Loaders.GeneOntologyStallTimeout (2 min); any existing file is left
  untouched

- `mzlib_bridge_error` (correctness):

  a non-transient HTTP failure (e.g. 404); a downloaded file
  GeneOntologyGraph.Load refuses

## Caveats

- The whole file (~37 MB) is fetched every time; there is no conditional
  request. A file that is byte-identical is left in place.

- When the download differs, the old file is kept beside the new one as
  \<name\>.\<yyyyMMdd-HHmmss-fff\>, so earlier runs stay reproducible.

- This is the only verb that fetches go.obo; proteins annotate-go never
  downloads.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.proteins.update_go`

- Rust (mzLibRust): `mzlib::proteins::update_go`

- R (mzLibR): `proteins_update_go`

## Since

Wire protocol 1; pyMzLib 0.3.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- The backup's exact file name is not returned: mzLib's
  UpdateGeneOntology computes it internally and returns nothing.

## References

- [doi:10.1093/genetics/iyad031](https://doi.org/10.1093/genetics/iyad031):
  The Gene Ontology knowledgebase in 2023 (GO Consortium, Genetics 224,
  iyad031), the release fetched

## See also

[`proteins_annotate_go`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)

## Examples

``` r

update <- proteins_update_go("go.obo")
update
#> <mzlibr_go_update> go.obo changed
#>   release: releases/2026-07-26, 48340 terms
update$go[c("release", "term_count")]
#> $release
#> [1] "releases/2026-07-26"
#> 
#> $term_count
#> [1] 48340
#> 
```
