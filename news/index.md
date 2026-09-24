# Changelog

## mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb’s `since.mzlibr` in
the bridge’s verb specs, which every help page’s **Since** section
renders.

This release projects mzLib 1.0.592. The verbs marked *needs the pyMzLib
0.2.0 bridge* below are new in the bridge pyMzLib 0.2.0 publishes; with
an older bridge they are refused, naming the release they need.

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
  none - for a `_many` batch, each file’s rows by that file’s unit.
- Many files in one call:
  [`readers_identify_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify_many.md)
  and a `_many` form of every reader hand the whole list to one bridge
  process with `threads` and `on_error` stated, and return one long
  `records` table plus a `files` data.frame of per-file facts, one row
  per input (*needs the pyMzLib 0.2.0 bridge*).
- Every read reports the per-file block: `reader`, `rows_not_read`,
  `retention_time_unit`, `caveats`, and the three ways a column can be
  empty kept apart - `absent_fields` (the format has no source for it,
  so it is `NA` in every row), `failed_fields` and `excluded_fields`.
- [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  reports the run’s `source`: instrument model and its PSI-MS accession,
  serial number, and acquisition start time with whether it is UTC
  (mzLib
  [\#1349](https://github.com/smith-chem-wisc/mzLibR/issues/1349)).
- [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  reads mzIdentML (`.mzid`, `.mzid.gz`) with each match’s `q_value`,
  `rank` and `pass_threshold`, and with `scores = TRUE` one row per
  engine score; mzIdentML items mzLib does not represent are listed in
  `skipped` with the reason.
- New
  [`readers_read_protein_groups()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md),
  [`readers_read_quantified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
  and
  [`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
  read MetaMorpheus and FlashLFQ quantification tables long, one row per
  record per sample, and per-site PTM occupancy (mzLib
  [\#1347](https://github.com/smith-chem-wisc/mzLibR/issues/1347);
  *needs the pyMzLib 0.2.0 bridge*).
- [`mzlibr_bridge_version()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_version.md)
  reports `verbs`, every command the bridge dispatches.
- mzLib now recognises 36 file types (Pytheas, mzIdentML and its gzipped
  form, and the MetaMorpheus and FlashLFQ quantification tables are
  new), 17 of them with no cross-format view.
- `limit = 0` is accepted by every reader, as the bridge allows: the
  envelope alone.

### SDRF-Proteomics

- [`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  reads one experimental-design file and
  [`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  pools several into one analysis table with provenance;
  [`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md),
  [`sdrf_all()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md),
  [`sdrf_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
  and the other `sdrf_` accessors read them.
- [`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md),
  [`sdrf_lint()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md),
  [`sdrf_assess()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md),
  [`sdrf_samples()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)
  and
  [`sdrf_parse_ages()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md)
  ask mzLib whether an SDRF document is well-formed, whether several
  agree on their vocabulary, whether it describes its samples, and what
  its ages are in years;
  [`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md),
  [`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)
  and
  [`sdrf_samples_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md)
  read a corpus in one bridge call (mzLib
  [\#1207](https://github.com/smith-chem-wisc/mzLibR/issues/1207),
  [\#1325](https://github.com/smith-chem-wisc/mzLibR/issues/1325),
  [\#1326](https://github.com/smith-chem-wisc/mzLibR/issues/1326),
  [\#1333](https://github.com/smith-chem-wisc/mzLibR/issues/1333),
  [\#1335](https://github.com/smith-chem-wisc/mzLibR/issues/1335)).

### PRIDE, peptidoforms and quantification

- [`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
  [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md),
  [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  and
  [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  list and fetch PRIDE Archive projects.
- [`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md)
  finds PRIDE projects by keyword, one row per project, with dates as
  `Date`.
- [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  digests a UniProt entry and fragments its peptidoforms.
- [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  runs FlashLFQ label-free quantification with match-between-runs.
- [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
  re-runs FlashLFQ’s protein median polish from a
  `QuantifiedPeptides.tsv` under a new experimental design, without
  re-reading spectra.
- [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  no longer warns about mzLib#1111 when `max_threads` is not 1:
  mzLib#1155 fixed it inside the pinned bridge. The default stays 1 here
  (pyMzLib and mzLibRust use -1).

### Protein databases

- New
  [`proteins_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)
  reads UniProt XML or FASTA databases into one row per protein -
  organism, NCBI taxon, gene names, length, mass - with GO terms and
  Ensembl gene links on request;
  [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
  resolves proteins to stable Ensembl gene ids against a gene set you
  pin;
  [`proteins_classify_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md)
  calls each peptide `Unique`, `SharedWithinGene`, `SharedAcrossGenes`
  or `NotInDatabase`, with I and L one residue. One database or many, in
  one bridge call (mzLib
  [\#1336](https://github.com/smith-chem-wisc/mzLibR/issues/1336),
  [\#1338](https://github.com/smith-chem-wisc/mzLibR/issues/1338),
  [\#1348](https://github.com/smith-chem-wisc/mzLibR/issues/1348);
  *needs the pyMzLib 0.2.0 bridge*).

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
