# Changelog

## mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb’s `since.mzlibr` in
the bridge’s verb specs, which every help page’s **Since** section
renders.

### Reading files

- [`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  and
  [`readers_identify()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
  say what a file is and which cross-format views it offers.
- [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  reads scans from mzML, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF
  and msalign: headers always, peaks on request.
- [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md),
  [`readers_read_features()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md)
  and
  [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  read the three cross-format views of search results;
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  reads any recognised file into that format’s own fields.
- [`readers_retention_time_in_minutes()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
  converts using each file’s declared unit, and refuses when there is
  none.

### SDRF-Proteomics

- [`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  reads one experimental-design file and
  [`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  pools several into one analysis table with provenance;
  [`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md),
  [`sdrf_all()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md),
  [`sdrf_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
  and the other `sdrf_` accessors read them.

### PRIDE, peptidoforms and quantification

- [`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
  [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md),
  [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  and
  [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  list and fetch PRIDE Archive projects.
- [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  digests a UniProt entry and fragments its peptidoforms.
- [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  runs FlashLFQ label-free quantification with match-between-runs.

### Documentation

- Every help page for a wire verb now carries that verb’s facts,
  generated from the bridge’s per-verb spec rather than retyped: what it
  wraps in mzLib, every parameter with its unit, range and default,
  every returned field and column with its unit and what `NA` means, the
  errors it raises as condition classes, caveats, the same verb in
  pyMzLib and mzLibRust, references and `since`. `scripts/name-parity.R`
  checks the hand-written `@param` and `@return` text against the same
  specs, and CI fails on any difference.
- Every exported function has an example, and `R CMD check` runs them
  all: a stand-in bridge answers each call from a recording of the real
  one, so the examples need no bridge, .NET or network.
- New
  [`?mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  documents the condition classes and how to handle each.
- A pkgdown site, with an article for each module.
