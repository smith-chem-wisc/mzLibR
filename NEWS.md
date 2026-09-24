# mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb's `since.mzlibr` in the bridge's verb specs,
which every help page's **Since** section renders.

This release projects mzLib 1.0.592. The verbs marked *needs the pyMzLib 0.2.0 bridge* below are
new in the bridge pyMzLib 0.2.0 publishes; with an older bridge they are refused, naming the
release they need.

## Reading files

* `readers_formats()` and `readers_identify()` say what a file is and which cross-format views
  it offers.
* `readers_read_spectra()` reads scans from mzML, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF
  and msalign: headers always, peaks on request.
* `readers_read_results()`, `readers_read_features()` and `readers_read_matches()` read the three
  cross-format views of search results; `readers_read_records()` reads any recognised file into
  that format's own fields.
* `readers_retention_time_in_minutes()` converts using each file's declared unit, and refuses
  when there is none - for a `_many` batch, each file's rows by that file's unit.
* Many files in one call: `readers_identify_many()` and a `_many` form of every reader hand the
  whole list to one bridge process with `threads` and `on_error` stated, and return one long
  `records` table plus a `files` data.frame of per-file facts, one row per input
  (*needs the pyMzLib 0.2.0 bridge*).
* Every read reports the per-file block: `reader`, `rows_not_read`, `retention_time_unit`,
  `caveats`, and the three ways a column can be empty kept apart - `absent_fields` (the format has
  no source for it, so it is `NA` in every row), `failed_fields` and `excluded_fields`.
* `readers_read_spectra()` reports the run's `source`: instrument model and its PSI-MS accession,
  serial number, and acquisition start time with whether it is UTC (mzLib #1349).
* `readers_read_matches()` reads mzIdentML (`.mzid`, `.mzid.gz`) with each match's `q_value`,
  `rank` and `pass_threshold`, and with `scores = TRUE` one row per engine score; mzIdentML items
  mzLib does not represent are listed in `skipped` with the reason.
* New `readers_read_protein_groups()`, `readers_read_quantified_peptides()` and
  `readers_read_occupancy()` read MetaMorpheus and FlashLFQ quantification tables long, one row per
  record per sample, and per-site PTM occupancy (mzLib #1347; *needs the pyMzLib 0.2.0 bridge*).
* `mzlibr_bridge_version()` reports `verbs`, every command the bridge dispatches.
* mzLib now recognises 36 file types (Pytheas, mzIdentML and its gzipped form, and the
  MetaMorpheus and FlashLFQ quantification tables are new), 17 of them with no cross-format view.
* `limit = 0` is accepted by every reader, as the bridge allows: the envelope alone.

## SDRF-Proteomics

* `sdrf_read()` reads one experimental-design file and `sdrf_pool()` pools several into one
  analysis table with provenance; `sdrf_value()`, `sdrf_all()`, `sdrf_records()` and the other
  `sdrf_` accessors read them.
* `sdrf_validate()`, `sdrf_lint()`, `sdrf_assess()`, `sdrf_samples()` and `sdrf_parse_ages()` ask
  mzLib whether an SDRF document is well-formed, whether several agree on their vocabulary,
  whether it describes its samples, and what its ages are in years; `sdrf_validate_many()`,
  `sdrf_assess_many()` and `sdrf_samples_many()` read a corpus in one bridge call
  (mzLib #1207, #1325, #1326, #1333, #1335).

## PRIDE, peptidoforms and quantification

* `pride_list_files()`, `pride_list_ftp_files()`, `pride_download()` and
  `pride_download_files()` list and fetch PRIDE Archive projects.
* `pride_search()` finds PRIDE projects by keyword, one row per project, with dates as `Date`.
* `peptidoform_fragments()` digests a UniProt entry and fragments its peptidoforms.
* `flashlfq_quantify()` runs FlashLFQ label-free quantification with match-between-runs.
* `flashlfq_median_polish()` re-runs FlashLFQ's protein median polish from a
  `QuantifiedPeptides.tsv` under a new experimental design, without re-reading spectra.
* `flashlfq_quantify()` no longer warns about mzLib#1111 when `max_threads` is not 1: mzLib#1155
  fixed it inside the pinned bridge. The default stays 1 here (pyMzLib and mzLibRust use -1).

## Protein databases

* New `proteins_read()` reads UniProt XML or FASTA databases into one row per protein - organism,
  NCBI taxon, gene names, length, mass - with GO terms and Ensembl gene links on request;
  `proteins_resolve_genes()` resolves proteins to stable Ensembl gene ids against a gene set you
  pin; `proteins_classify_peptides()` calls each peptide `Unique`, `SharedWithinGene`,
  `SharedAcrossGenes` or `NotInDatabase`, with I and L one residue. One database or many, in one
  bridge call (mzLib #1336, #1338, #1348; *needs the pyMzLib 0.2.0 bridge*).

## Documentation

* Every help page for a wire verb now carries that verb's facts, generated from the bridge's
  per-verb spec rather than retyped: what it wraps in mzLib, every parameter with its unit, range
  and default, every returned field and column with its unit and what `NA` means, the errors it
  raises as condition classes, caveats, the same verb in pyMzLib and mzLibRust, references and
  `since`. `scripts/name-parity.R` checks the hand-written `@param` and `@return` text against the
  same specs, and CI fails on any difference.
* Every exported function has an example, and `R CMD check` runs them all: a stand-in bridge
  answers each call from a recording of the real one, so the examples need no bridge, .NET or
  network.
* New `?mzlib_error` documents the condition classes and how to handle each.
* A pkgdown site, with an article for each module.
