# Changelog

## mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb’s `since.mzlibr` in
the bridge’s verb specs, which every help page’s **Since** section
renders.

This release projects mzLib 1.0.593. The verbs marked *needs the pyMzLib
0.2.0 bridge*, *needs the pyMzLib 0.3.0 bridge* or *needs the pyMzLib
0.4.0 bridge* below are new in the bridge that release publishes; with
an older bridge they are refused, naming the release they need.
[`mzlibr_install_bridge()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md)
installs the bridge pyMzLib 0.4.0 publishes (mzLib 1.0.593, `0a808fec`),
verified against that release’s SHA256SUMS.

### mzLib 1.0.593

- [`readers_read_protein_groups()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md),
  [`readers_read_quantified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
  and
  [`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
  read an RNA search’s `AllQuantifiedTranscriptGroups.tsv` and
  `AllQuantifiedOligos.tsv` (mzLib
  [\#1388](https://github.com/smith-chem-wisc/mzLibR/issues/1388)).
  mzLib reads them with subclasses of the protein-group and peptide
  readers, so the columns keep their protein names: `protein_group_name`
  is the transcript group, `sequence` the oligonucleotide, and occupancy
  names RNA modifications.
- MetaMorpheus protein-group tables written without quantification now
  read (mzLib
  [\#1365](https://github.com/smith-chem-wisc/mzLibR/issues/1365)):
  `AllProteinGroups.tsv` and each file’s `<file>_ProteinGroups.tsv`
  failed with `Tsv file type not supported`, and now read as
  `MetaMorpheusQuantifiedProteinGroups`, with `intensity` in
  `absent_fields`.
- Four shipped modifications now write their Unimod accession in
  `pro_forma` (mzLib
  [\#1328](https://github.com/smith-chem-wisc/mzLibR/issues/1328)):
  `GG (Ubiquitination Site)` as `[UNIMOD:121]`, both `Myristoylation`
  entries as `[UNIMOD:45]` and `EQIGG` as `[UNIMOD:846]`, instead of by
  name. Masses are unchanged.
- A read fault on an existing `.mzid` is an `mzlib_bridge_error` naming
  the file (mzLib
  [\#1362](https://github.com/smith-chem-wisc/mzLibR/issues/1362)),
  where it was an `IOException`. A missing file is still an
  `mzlib_usage_error`.
- PRIDE download errors name the file and the host, never the URL (mzLib
  [\#1350](https://github.com/smith-chem-wisc/mzLibR/issues/1350)), so a
  reviewer token in a query string cannot reach a log. They are still
  `mzlib_service_unavailable`.
- The help pages no longer state how many file types mzLib recognises:
  [`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  asks mzLib, and a number written down goes stale at the next release.

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
- mzLib now recognises Pytheas, mzIdentML and its gzipped form, and the
  MetaMorpheus and FlashLFQ quantification tables.
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
- New
  [`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
  reads the label-free experimental design MetaMorpheus and FlashLFQ
  take out of an SDRF, with mzLib’s `SdrfLabelFreeDesign`, or every
  reason it refuses to - all at once, and as a result rather than an
  error.
  [`sdrf_design_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md)
  and
  [`sdrf_design_run_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_run_design.md)
  hand it to
  [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  and
  [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md);
  `out` writes MetaMorpheus’s `ExperimentalDesign.tsv`. Its replicate
  and fraction coordinates stay 0-based, as the quant functions take
  them (mzLib
  [\#1363](https://github.com/smith-chem-wisc/mzLibR/issues/1363);
  *needs the pyMzLib 0.3.0 bridge*).

### Isobaric kits

- New
  [`isobaric_kits()`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
  lists mzLib’s TMT, TMTpro, iTRAQ and DiLeu kits with every channel’s
  label, theoretical reporter-ion m/z and matching window, or one kit by
  MetaMorpheus’s name for it; `channel_index` is 1-based (mzLib
  [\#1375](https://github.com/smith-chem-wisc/mzLibR/issues/1375);
  *needs the pyMzLib 0.3.0 bridge*).

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
- New
  [`peptidoform_convert()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_convert.md)
  rewrites full sequences in another notation with mzLib’s
  `SequenceConversionService`. The main use is MetaMorpheus full
  sequences to Unimod accessions: `[UniProt:N-acetylserine on S]SEQK`
  becomes `[UNIMOD:1]SEQK`, and dimethyllysine becomes `UNIMOD:36`. One
  row per input, in order, in `records`, with mzLib’s own verdict
  (`converted`, `converted_with_warnings` or `failed`), the
  modifications it could not write, and its warnings. `source` and
  `target` take any notation mzLib has registered; every result lists
  them. `mode` is mzLib’s `SequenceConversionHandlingMode`. **The
  ProForma target does not resolve UniProt modifications** in this mzLib
  build (mzLib#1401) and writes them back by name; convert to Unimod for
  those. The `pro_forma` column of
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  has the same gap (pyMzLib
  [\#75](https://github.com/smith-chem-wisc/mzLibR/issues/75); *needs
  the pyMzLib 0.4.0 bridge*, which
  [`mzlibr_install_bridge()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md)
  now installs).
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
- [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
  with no design no longer hangs in terminal R. A call that sends the
  bridge nothing on stdin now gives it an empty stdin rather than R’s
  own, which a terminal never closes. Every verb goes through the same
  transport, so none can wait on it (pyMzLib
  [\#73](https://github.com/smith-chem-wisc/mzLibR/issues/73)).

### Differential abundance

- New
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  fits one linear model per feature of a feature-by-sample table and
  tests each named coefficient with limma’s empirical-Bayes moderated t,
  Benjamini-Hochberg adjusted: `lmFit()` then `eBayes(legacy = TRUE)`,
  computed by mzLib’s `LinearModel` and `EmpiricalBayes`. A feature that
  cannot be fitted is reported with its reason, never dropped, and
  `residual_df_differ` says when default limma would use a different
  prior estimator.
- New
  [`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
  Benjamini-Hochberg adjusts p-values from anywhere, keeping every
  position; `NA` is untested and not counted in m.
- New
  [`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
  pools one estimate per study into a DerSimonian-Laird random-effects
  estimate per feature, with direction agreement and leave-one-out
  sensitivity.
- The tests hold all three to limma, metafor and
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) to 1e-8
  relative - against their recorded output always, and against limma and
  metafor run in the test when they are installed (both are now in
  `Suggests`). The reference tables ship in `extdata/stats/` so the
  examples and the new article read real files (mzLib
  [\#1341](https://github.com/smith-chem-wisc/mzLibR/issues/1341),
  [\#1357](https://github.com/smith-chem-wisc/mzLibR/issues/1357);
  *needs the pyMzLib 0.3.0 bridge*).

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
- New
  [`proteins_annotate_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)
  annotates a stored MetaMorpheus protein-group table with Gene Ontology
  terms: one row per (group, term) that any member holds, directly or
  through an ancestor, naming the members that carry it, so consensus
  and direct-only views are filters on the rows. Every non-decoy group
  gets a row, a term-less one saying why. It reads the go.obo you name
  and never downloads one;
  [`proteins_update_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_update_go.md)
  is the one function that fetches the current release, keeping the
  previous file beside it (mzLib
  [\#1353](https://github.com/smith-chem-wisc/mzLibR/issues/1353),
  [\#1366](https://github.com/smith-chem-wisc/mzLibR/issues/1366);
  *needs the pyMzLib 0.3.0 bridge*).

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
- The articles are now vignettes, so `R CMD check` builds and runs every
  one of them, on every platform, against the replay bridge. Each opens
  with a question -\> function -\> mzLib table, ends with what to cite
  (rendered from the specs’ DOIs, each checked to resolve), and states
  no count of formats or verbs and no mzLib version in its prose;
  `scripts/docs-lint.R` holds all of that in CI.
- The articles teach on real data recorded from the real bridge: a
  FlashLFQ run with match-between-runs on mzLib’s K562 pair (and median
  polish reproducing its protein intensities for every group), the whole
  albumin digest, PXD000001’s live FTP listing, an RNA search’s
  transcript groups, and the MALAT1 dilution series for
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md).
- Corrected from those runs: the K562 peptide roll-up shows 21 of the
  140 MBR transfers, not 52; `use_pep_q_value` filters nothing, it
  changes the q-value FlashLFQ carries; PXD000001’s FTP tree now holds
  14 files.
- `knitr` and `rmarkdown` join `Suggests`, for the vignettes. The
  package still imports nothing but base R.
