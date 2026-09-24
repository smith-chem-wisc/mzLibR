# Package index

## Reading files

What a file is, and its records: spectra, search results, features and
matches, or any recognised format in its own fields.

- [`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  : Every file type mzLib can recognise

- [`readers_identify()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
  : Identify a result file without parsing its contents

- [`readers_identify_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify_many.md)
  : Identify many files in one bridge call

- [`readers_read_features()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md)
  :

  Read deconvolved MS1 features, in the cross-format `ms1_features` view

- [`readers_read_features_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features_many.md)
  : Read MS1 features from many files, in one bridge call

- [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  :

  Read identifications, in the cross-format `spectral_match` view

- [`readers_read_matches_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches_many.md)
  : Read identifications from many files, in one bridge call

- [`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
  : Read PTM site occupancy from a MetaMorpheus protein-group table

- [`readers_read_occupancy_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy_many.md)
  : Read the PTM site occupancy of many protein-group tables, in one
  bridge call

- [`readers_read_protein_groups()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md)
  : Read a MetaMorpheus protein-group table, one row per group per
  sample group

- [`readers_read_protein_groups_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups_many.md)
  : Read many protein-group tables into one long table, in one bridge
  call

- [`readers_read_quantified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
  : Read a FlashLFQ peptide table, one row per peptide per sample

- [`readers_read_quantified_peptides_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides_many.md)
  : Read many peptide tables into one long table, in one bridge call

- [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  : Read any file mzLib recognises, into that format's own fields

- [`readers_read_records_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records_many.md)
  : Read many files, each into its format's own fields, in one bridge
  call

- [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
  : Read a result file into the uniform record view

- [`readers_read_results_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results_many.md)
  : Read many result files into the uniform record view, in one bridge
  call

- [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  : Read the scans of a spectra file: headers always, peaks on request

- [`readers_read_spectra_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra_many.md)
  : Read the scans of many spectra files, in one bridge call

- [`readers_retention_time_in_minutes()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
  : Retention times in minutes, whatever unit the format wrote

## SDRF-Proteomics experimental designs

Read one design file, or pool several into one analysis table.

- [`sdrf_all()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md)
  : Every cell under a column, one list element per row
- [`sdrf_assess()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)
  : Decide whether an SDRF file describes its samples, or is a valid
  skeleton
- [`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)
  : Assess many SDRF files in one bridge call
- [`sdrf_has_repeated_columns()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_has_repeated_columns.md)
  : Whether an SDRF document repeats a column name
- [`sdrf_index_of()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_index_of.md)
  : Where a column first sits in an SDRF document
- [`sdrf_indexes_of()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_indexes_of.md)
  : Every position a column name occupies in an SDRF document
- [`sdrf_lint()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md)
  : Find the concepts a set of SDRF documents annotated inconsistently
- [`sdrf_parse_ages()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md)
  : Read SDRF age cells into years, refusing anything that would need a
  guess
- [`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  : Merge several SDRF documents into one analysis table
- [`sdrf_ragged_row_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_ragged_row_count.md)
  : How many rows are shorter than the header
- [`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  : Read one SDRF-Proteomics experimental-design file
- [`sdrf_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
  : An SDRF document as a data.frame, when its shape allows one
- [`sdrf_samples()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)
  : Each sample's characteristics and factor values, merged over its
  rows, ages parsed
- [`sdrf_samples_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md)
  : The samples of many SDRF files, in one bridge call
- [`sdrf_source_documents()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_source_documents.md)
  : Which document each pooled row came from
- [`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)
  : Check one SDRF file against the specification's structural rules
- [`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md)
  : Validate many SDRF files in one bridge call
- [`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
  : The first cell under a column, one per row

## Protein databases

What an accession is, which gene it resolves to, and whether a peptide
is unique.

- [`proteins_classify_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md)
  : Classify peptides by how widely they are shared across protein
  databases, with I = L
- [`proteins_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)
  : Read protein databases: one row per protein, with GO terms and
  Ensembl genes on request
- [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
  : Resolve every protein to stable Ensembl gene ids, against a gene set
  you pin

## PRIDE Archive

List a project’s files and download the ones you choose.

- [`pride_approximate_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_approximate_total_size_bytes.md)
  : Sum the approximate sizes of some FTP files
- [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  : Download files from a PRIDE Archive project
- [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  : Download exactly the files you selected
- [`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  : List the files in a PRIDE Archive project
- [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
  : The complete file list of a PRIDE Archive project, from its FTP
  directory tree
- [`pride_locations()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_locations.md)
  : The published locations of each file, as controlled-vocabulary terms
- [`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md)
  : Find PRIDE Archive projects by keyword
- [`pride_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_total_size_bytes.md)
  : Total size of a set of PRIDE files

## Peptidoforms

Digest a UniProt entry and fragment its peptidoforms.

- [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  : Digest a protein and fragment its peptidoforms
- [`digest_distinct_base_sequences()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_distinct_base_sequences.md)
  : How many distinct base sequences a digest produced
- [`digest_fragments_by_series()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_fragments_by_series.md)
  : Fragment ions per product type
- [`digest_modified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_modified_peptides.md)
  : Only the peptidoforms carrying at least one modification
- [`digest_truncated()`](https://smith-chem-wisc.github.io/mzLibR/reference/digest_truncated.md)
  : Whether a digest hit the isoform cap, and is therefore incomplete
- [`census_excluded()`](https://smith-chem-wisc.github.io/mzLibR/reference/census_excluded.md)
  : Annotated features that could not be used
- [`census_explain()`](https://smith-chem-wisc.github.io/mzLibR/reference/census_explain.md)
  : What UniProt annotated, and what mzLib could actually use
- [`peptide_mz()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptide_mz.md)
  : The m/z of intact peptides at a given charge

## Quantification

Label-free quantification with FlashLFQ.

- [`flashlfq_mbr_peak_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_peak_count.md)
  : Total match-between-runs peaks across every run
- [`flashlfq_mbr_peaks()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_peaks.md)
  : Exactly the peaks transferred by match-between-runs
- [`flashlfq_mbr_rescued_peptide_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_mbr_rescued_peptide_count.md)
  : Distinct peptides rescued by match-between-runs
- [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
  : Roll a FlashLFQ peptide table up to protein intensities, under a new
  design
- [`flashlfq_peptide_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_peptide_count.md)
  : Number of quantified peptides
- [`flashlfq_protein_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_protein_count.md)
  : Number of quantified protein groups
- [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  : Quantify a search's peptides across mzML runs with FlashLFQ

## The bridge

Installing and locating the bridge every function runs, and asking it
what it is.

- [`mzlibr_bridge_path()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_path.md)
  : Path of the bridge executable mzLibR will use
- [`mzlibr_bridge_version()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_version.md)
  : Version information reported by the mzLib bridge
- [`mzlibr_install_bridge()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md)
  : Download a bridge executable into a local cache

## Errors

- [`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_usage_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_service_unavailable`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_bridge_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_project_not_found`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_timeout`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_bridge_not_found`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  [`mzlib_protocol_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  : The conditions mzLibR signals, and how to handle each

## Printing

- [`print(`*`<mzlibr_batch>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_batch.md)
  : Print a batch read from many files
- [`print(`*`<mzlibr_census>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_census.md)
  : Print a modification census
- [`print(`*`<mzlibr_digest>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_digest.md)
  : Print a digest
- [`print(`*`<mzlibr_feature_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_feature_records.md)
  : Print an MS1 feature table
- [`print(`*`<mzlibr_file_info>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_file_info.md)
  : Print a file identification
- [`print(`*`<mzlibr_gene_resolutions>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_gene_resolutions.md)
  : Print a gene resolution
- [`print(`*`<mzlibr_match_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_match_records.md)
  : Print a spectral-match table
- [`print(`*`<mzlibr_median_polish>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_median_polish.md)
  : Print a median-polish result
- [`print(`*`<mzlibr_native_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_native_records.md)
  : Print a native record table
- [`print(`*`<mzlibr_occupancy_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_occupancy_records.md)
  : Print a site-occupancy table
- [`print(`*`<mzlibr_peptide_classification>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_peptide_classification.md)
  : Print a peptide classification
- [`print(`*`<mzlibr_pooled_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_pooled_sdrf.md)
  : Print a pooled SDRF table
- [`print(`*`<mzlibr_protein_database>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_protein_database.md)
  : Print a protein database read
- [`print(`*`<mzlibr_protein_group_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_protein_group_records.md)
  : Print a protein-group table
- [`print(`*`<mzlibr_quant>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_quant.md)
  : Print a quantification result
- [`print(`*`<mzlibr_quantified_peptide_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_quantified_peptide_records.md)
  : Print a quantified-peptide table
- [`print(`*`<mzlibr_result_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_result_records.md)
  : Print a record view
- [`print(`*`<mzlibr_scan_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_scan_records.md)
  : Print a scan table
- [`print(`*`<mzlibr_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf.md)
  : Print an SDRF document
- [`print(`*`<mzlibr_sdrf_ages>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_ages.md)
  : Print parsed SDRF ages
- [`print(`*`<mzlibr_sdrf_assessment>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_assessment.md)
  : Print an SDRF informativeness assessment
- [`print(`*`<mzlibr_sdrf_drift>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_drift.md)
  : Print SDRF drift findings
- [`print(`*`<mzlibr_sdrf_samples>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_samples.md)
  : Print an SDRF document's samples
- [`print(`*`<mzlibr_sdrf_validation>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_validation.md)
  : Print an SDRF validation
