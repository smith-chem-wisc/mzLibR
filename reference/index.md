# Package index

## Reading files

What a file is, and its records: spectra, search results, features and
matches, or any recognised format in its own fields.

- [`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  : Every file type mzLib can recognise

- [`readers_identify()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
  : Identify a result file without parsing its contents

- [`readers_read_features()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md)
  :

  Read deconvolved MS1 features, in the cross-format `ms1_features` view

- [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  :

  Read identifications, in the cross-format `spectral_match` view

- [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  : Read any file mzLib recognises, into that format's own fields

- [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
  : Read a result file into the uniform record view

- [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  : Read the scans of a spectra file: headers always, peaks on request

- [`readers_retention_time_in_minutes()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
  : Retention times in minutes, whatever unit the format wrote

## SDRF-Proteomics experimental designs

Read one design file, or pool several into one analysis table.

- [`sdrf_all()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md)
  : Every cell under a column, one list element per row
- [`sdrf_has_repeated_columns()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_has_repeated_columns.md)
  : Whether an SDRF document repeats a column name
- [`sdrf_index_of()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_index_of.md)
  : Where a column first sits in an SDRF document
- [`sdrf_indexes_of()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_indexes_of.md)
  : Every position a column name occupies in an SDRF document
- [`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  : Merge several SDRF documents into one analysis table
- [`sdrf_ragged_row_count()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_ragged_row_count.md)
  : How many rows are shorter than the header
- [`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  : Read one SDRF-Proteomics experimental-design file
- [`sdrf_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
  : An SDRF document as a data.frame, when its shape allows one
- [`sdrf_source_documents()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_source_documents.md)
  : Which document each pooled row came from
- [`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
  : The first cell under a column, one per row

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

- [`print(`*`<mzlibr_census>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_census.md)
  : Print a modification census
- [`print(`*`<mzlibr_digest>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_digest.md)
  : Print a digest
- [`print(`*`<mzlibr_feature_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_feature_records.md)
  : Print an MS1 feature table
- [`print(`*`<mzlibr_file_info>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_file_info.md)
  : Print a file identification
- [`print(`*`<mzlibr_match_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_match_records.md)
  : Print a spectral-match table
- [`print(`*`<mzlibr_native_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_native_records.md)
  : Print a native record table
- [`print(`*`<mzlibr_pooled_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_pooled_sdrf.md)
  : Print a pooled SDRF table
- [`print(`*`<mzlibr_quant>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_quant.md)
  : Print a quantification result
- [`print(`*`<mzlibr_result_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_result_records.md)
  : Print a record view
- [`print(`*`<mzlibr_scan_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_scan_records.md)
  : Print a scan table
- [`print(`*`<mzlibr_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf.md)
  : Print an SDRF document
