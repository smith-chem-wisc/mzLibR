# mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb's `since.mzlibr` in the bridge's verb specs,
which every help page's **Since** section renders.

This release projects mzLib 1.0.595. The verbs marked *needs the pyMzLib 0.2.0 bridge*, *needs
the pyMzLib 0.3.0 bridge* or *needs the pyMzLib 0.4.0 bridge* below are new in the bridge that
release publishes; with an older bridge they are refused, naming the release they need. `mzlibr_install_bridge()` installs the bridge
pyMzLib 0.5.0 publishes (mzLib 1.0.595, `d10017ce`), verified against that release's SHA256SUMS.

## Defaults come from mzLib (the pyMzLib 0.5.0 bridge)

* **`peptidoform_fragments()` takes its defaults from mzLib, so a call that names no protease or
  dissociation type gives different peptides and ions.** Every digestion argument (`protease`,
  `dissociation`, `missed_cleavages`, `min_length`, `max_length`, `max_modifications`,
  `max_isoforms`, `terminus`) now defaults to `NULL`, mzLibR sends nothing for it, and mzLib's
  `DigestionParams` default applies. Two of the old defaults were mzLibR's own and change:
  - protease: mzLib's default is `"trypsin"`, which also cuts before proline. The old default was
    `"trypsin|P"`, which does not. Albumin's default digest goes from 303 to 310 peptidoforms, and
    from 195 to 202 distinct sequences.
  - dissociation: mzLib names no default, so the bridge uses HCD (b and y ions) until it does
    (smith-chem-wisc/mzLib#1488). The old default was ETD (c and z-dot ions).

  The other defaults (missed cleavages, the length window, the modification and isoform caps,
  terminus) have the same values as before, now set by mzLib. **To get the old answer, pass
  `protease = "trypsin|P", dissociation = "ETD"`.** A digest records `min_length`, `max_length`
  and `missed_cleavages` as `NA` when they were left to mzLib, since the bridge does not report
  those three back, and prints them as mzLib's default. A `max_length` below `min_length` is
  refused only when both are given. The help page's counts are re-measured at mzLib `d10017ce`.
* **`pride_list_files()` and `pride_search()` leave the page size to mzLib.** `page_size` now
  defaults to `NULL` and is sent only when given, so mzLib's own paging default applies (100
  today, as before). The results never depended on it. A given page size is checked as before.
* `?stats_fit` no longer says the bridge's spline default is 4. With `trend = TRUE` and no
  `spline_basis`, mzLibR has always sent none; the 0.5.0 bridge now asks mzLib for none either,
  so mzLib sizes the trend as limma does and a trend fit on fewer than 30 features gets limma's
  numbers. `prior$spline_basis_count` says what was used.
* The help pages are re-rendered from the bridge's specs, re-checked against mzLib `d10017ce`.
  PXD000001's REST manifest listed 9 of the 14 files in its FTP tree on 2026-10-21, and
  `?pride_list_files` says so rather than 8 of 13.

## mzLib 1.0.595

Built from mzLib 1.0.595 (tag commit `d10017ce`), which carries mzLib 1.0.594 as well.

* **DIA-NN 2.x's `report.parquet` reads** (mzLib #1384, file type `DiaNnReportParquet`). DIA-NN
  2.x writes its main report only as parquet. mzLib reads it into the same record as the 1.8/1.9
  `report.tsv`, so it is a fifth quantifiable type: `readers_read_results()`,
  `readers_read_records()` and `flashlfq_quantify()` take it. **`readers_formats()` lists 39
  types, up from 38, and 5 offer the `quantifiable` view, up from 4**; 19 still offer none. Like
  the TSV report it is recognised by its columns (`Precursor.Id`, `Stripped.Sequence`, `Run`), not
  its name, and any other `.parquet` is an `mzlib_usage_error`. **Its `retention_time_unit` is
  `"minutes"`**, as for the TSV report. `absent_fields` is always empty for it, because there is no
  text header to judge from, so the columns 2.x stopped writing (`File.Name`, `PG.Quantity`,
  `MS2.Scan` and others) are `NA` in `readers_read_records()` without being named there, and
  `ms2_scan_number` reads 0.
* **`peptidoform_convert(target = "ProForma")` names UniProt modifications by UNIMOD accession**
  (mzLib 1.0.594, #1403). `[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]-SEQK`, where it
  was written back by name; the caveat that said so is retired, so a result carries five caveats,
  not six. The `pro_forma` column `readers_read_records()` computes for a `.psmtsv` written
  without one gains the same accessions.
* The psmtsv caveat on `readers_read_results()` cites `SpectrumMatchFromTsv.cs:122` and `:197`;
  mzLib #1371 moved the lines it cited. The claim is unchanged.

## mzLib 1.0.593

* `readers_read_protein_groups()`, `readers_read_quantified_peptides()` and
  `readers_read_occupancy()` read an RNA search's `AllQuantifiedTranscriptGroups.tsv` and
  `AllQuantifiedOligos.tsv` (mzLib #1388). mzLib reads them with subclasses of the protein-group
  and peptide readers, so the columns keep their protein names: `protein_group_name` is the
  transcript group, `sequence` the oligonucleotide, and occupancy names RNA modifications.
* MetaMorpheus protein-group tables written without quantification now read (mzLib #1365):
  `AllProteinGroups.tsv` and each file's `<file>_ProteinGroups.tsv` failed with `Tsv file type not
  supported`, and now read as `MetaMorpheusQuantifiedProteinGroups`, with `intensity` in
  `absent_fields`.
* Four shipped modifications now write their Unimod accession in `pro_forma` (mzLib #1328):
  `GG (Ubiquitination Site)` as `[UNIMOD:121]`, both `Myristoylation` entries as `[UNIMOD:45]` and
  `EQIGG` as `[UNIMOD:846]`, instead of by name. Masses are unchanged.
* A read fault on an existing `.mzid` is an `mzlib_bridge_error` naming the file (mzLib #1362),
  where it was an `IOException`. A missing file is still an `mzlib_usage_error`.
* PRIDE download errors name the file and the host, never the URL (mzLib #1350), so a reviewer
  token in a query string cannot reach a log. They are still `mzlib_service_unavailable`.
* The help pages no longer state how many file types mzLib recognises: `readers_formats()` asks
  mzLib, and a number written down goes stale at the next release.

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
* mzLib now recognises Pytheas, mzIdentML and its gzipped form, and the MetaMorpheus and
  FlashLFQ quantification tables.
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
* New `sdrf_design()` reads the label-free experimental design MetaMorpheus and FlashLFQ take out
  of an SDRF, with mzLib's `SdrfLabelFreeDesign`, or every reason it refuses to - all at once,
  and as a result rather than an error. `sdrf_design_spectra()` and `sdrf_design_run_design()`
  hand it to `flashlfq_quantify()` and `flashlfq_median_polish()`; `out` writes MetaMorpheus's
  `ExperimentalDesign.tsv`. Its replicate and fraction coordinates stay 0-based, as the quant
  functions take them (mzLib #1363; *needs the pyMzLib 0.3.0 bridge*).

## Isobaric kits

* New `isobaric_kits()` lists mzLib's TMT, TMTpro, iTRAQ and DiLeu kits with every channel's
  label, theoretical reporter-ion m/z and matching window, or one kit by MetaMorpheus's name for
  it; `channel_index` is 1-based (mzLib #1375; *needs the pyMzLib 0.3.0 bridge*).

## PRIDE, peptidoforms and quantification

* `pride_list_files()`, `pride_list_ftp_files()`, `pride_download()` and
  `pride_download_files()` list and fetch PRIDE Archive projects.
* `pride_search()` finds PRIDE projects by keyword, one row per project, with dates as `Date`.
* `peptidoform_fragments()` digests a UniProt entry and fragments its peptidoforms.
* New `peptidoform_convert()` rewrites full sequences in another notation with mzLib's
  `SequenceConversionService`. The main use is MetaMorpheus full sequences to Unimod accessions:
  `[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]SEQK`, and dimethyllysine becomes
  `UNIMOD:36`. One row per input, in order, in `records`, with mzLib's own verdict (`converted`,
  `converted_with_warnings` or `failed`), the modifications it could not write, and its warnings.
  `source` and `target` take any notation mzLib has registered; every result lists them. `mode` is
  mzLib's `SequenceConversionHandlingMode`. The ProForma target names UniProt modifications by
  UNIMOD accession too, since mzLib #1403 (see mzLib 1.0.595 above; at mzLib 1.0.593 it wrote them
  back by name, mzLib#1401). (pyMzLib #75; *needs the pyMzLib 0.4.0 bridge*.)
* `flashlfq_quantify()` runs FlashLFQ label-free quantification with match-between-runs.
* `flashlfq_median_polish()` re-runs FlashLFQ's protein median polish from a
  `QuantifiedPeptides.tsv` under a new experimental design, without re-reading spectra.
* `flashlfq_quantify()` no longer warns about mzLib#1111 when `max_threads` is not 1: mzLib#1155
  fixed it inside the pinned bridge. The default stays 1 here (pyMzLib and mzLibRust use -1).
* `flashlfq_median_polish()` with no design no longer hangs in terminal R. A call that sends the
  bridge nothing on stdin now gives it an empty stdin rather than R's own, which a terminal never
  closes. Every verb goes through the same transport, so none can wait on it (pyMzLib #73).

## Differential abundance

* New `stats_fit()` fits one linear model per feature of a feature-by-sample table and tests each
  named coefficient with limma's empirical-Bayes moderated t, Benjamini-Hochberg adjusted:
  `lmFit()` then `eBayes(legacy = TRUE)`, computed by mzLib's `LinearModel` and `EmpiricalBayes`.
  A feature that cannot be fitted is reported with its reason, never dropped, and
  `residual_df_differ` says when default limma would use a different prior estimator.
* New `stats_adjust()` Benjamini-Hochberg adjusts p-values from anywhere, keeping every position;
  `NA` is untested and not counted in m.
* New `stats_meta()` pools one estimate per study into a DerSimonian-Laird random-effects estimate
  per feature, with direction agreement and leave-one-out sensitivity.
* The tests hold all three to limma, metafor and `stats::p.adjust()` to 1e-8 relative - against
  their recorded output always, and against limma and metafor run in the test when they are
  installed (both are now in `Suggests`). The reference tables ship in `extdata/stats/` so the
  examples and the new article read real files (mzLib #1341, #1357; *needs the pyMzLib 0.3.0
  bridge*).

## Protein databases

* New `proteins_read()` reads UniProt XML or FASTA databases into one row per protein - organism,
  NCBI taxon, gene names, length, mass - with GO terms and Ensembl gene links on request;
  `proteins_resolve_genes()` resolves proteins to stable Ensembl gene ids against a gene set you
  pin; `proteins_classify_peptides()` calls each peptide `Unique`, `SharedWithinGene`,
  `SharedAcrossGenes` or `NotInDatabase`, with I and L one residue. One database or many, in one
  bridge call (mzLib #1336, #1338, #1348; *needs the pyMzLib 0.2.0 bridge*).
* New `proteins_annotate_go()` annotates a stored MetaMorpheus protein-group table with Gene
  Ontology terms: one row per (group, term) that any member holds, directly or through an
  ancestor, naming the members that carry it, so consensus and direct-only views are filters on
  the rows. Every non-decoy group gets a row, a term-less one saying why. It reads the go.obo you
  name and never downloads one; `proteins_update_go()` is the one function that fetches the
  current release, keeping the previous file beside it (mzLib #1353, #1366; *needs the pyMzLib
  0.3.0 bridge*).

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
* The articles are now vignettes, so `R CMD check` builds and runs every one of them, on every
  platform, against the replay bridge. Each opens with a question -> function -> mzLib table,
  ends with what to cite (rendered from the specs' DOIs, each checked to resolve), and states no
  count of formats or verbs and no mzLib version in its prose; `scripts/docs-lint.R` holds all of
  that in CI.
* The articles teach on real data recorded from the real bridge: a FlashLFQ run with
  match-between-runs on mzLib's K562 pair (and median polish reproducing its protein intensities
  for every group), the whole albumin digest, PXD000001's live FTP listing, an RNA search's
  transcript groups, and the MALAT1 dilution series for `stats_fit()`.
* Corrected from those runs: the K562 peptide roll-up shows 21 of the 140 MBR transfers, not 52;
  `use_pep_q_value` filters nothing, it changes the q-value FlashLFQ carries; PXD000001's FTP
  tree now holds 14 files.
* `knitr` and `rmarkdown` join `Suggests`, for the vignettes. The package still imports nothing
  but base R.
