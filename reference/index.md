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

## Spectral libraries

mzLib’s .msl spectral libraries, read with mzLib’s own reader: every
precursor, and its fragment ions on request; and a run’s retention times
calibrated onto a library’s iRT scale.

- [`library_calibrate_rt()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_calibrate_rt.md)
  : Fit each run's map from minutes onto a library's iRT scale, and
  convert times either way

- [`library_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md)
  :

  Read an mzLib spectral library (`.msl`): one row per precursor, or per
  fragment ion

- [`library_read_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md)
  : Read many spectral libraries into one long table, in one bridge call

## SDRF-Proteomics experimental designs

Read one design file or pool several, ask whether each is well-formed
and informative, and read the label-free or isobaric design a
quantification takes, restricted to one search’s files.

- [`sdrf_all()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md)
  : Every cell under a column, one list element per row

- [`sdrf_assess()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md)
  : Decide whether an SDRF file describes its samples, or is a valid
  skeleton

- [`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)
  : Assess many SDRF files in one bridge call

- [`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
  : Read a label-free experimental design out of an SDRF, or every
  reason it cannot be read

- [`sdrf_design_run_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_run_design.md)
  :

  A label-free design as
  [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
  takes its `design`

- [`sdrf_design_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md)
  :

  A label-free design as
  [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  takes its `spectra`

- [`sdrf_has_repeated_columns()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_has_repeated_columns.md)
  : Whether an SDRF document repeats a column name

- [`sdrf_index_of()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_index_of.md)
  : Where a column first sits in an SDRF document

- [`sdrf_indexes_of()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_indexes_of.md)
  : Every position a column name occupies in an SDRF document

- [`sdrf_isobaric_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md)
  : Read an isobaric (TMT, TMTpro, iTRAQ, DiLeu) design out of an SDRF,
  or every reason it cannot be read

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

- [`sdrf_restrict()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md)
  : Restrict an SDRF to the files one search read, and report what could
  not be placed

- [`sdrf_samples()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)
  : Each sample's characteristics and factor values, merged over its
  rows, ages parsed

- [`sdrf_samples_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md)
  : The samples of many SDRF files, in one bridge call

- [`sdrf_searched_data_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_searched_data_files.md)
  : The searched file each row of a restricted SDRF was joined to

- [`sdrf_source_documents()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_source_documents.md)
  : Which document each pooled row came from

- [`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md)
  : Check one SDRF file against the specification's structural rules

- [`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md)
  : Validate many SDRF files in one bridge call

- [`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md)
  : The first cell under a column, one per row

## Isobaric kits

TMT, TMTpro, iTRAQ and DiLeu channels, with mzLib’s reporter-ion m/z.

- [`isobaric_kits()`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
  : The isobaric kits mzLib can name, with every channel's label and
  reporter-ion m/z

## Protein databases

What an accession is, which gene it resolves to, whether a peptide is
unique, and what each protein group does in Gene Ontology terms; where a
residue is in another sequence, which entry an accession names, which
genes are orthologous across species, and which residue corresponds
across species in Compara’s gene-tree alignment.

- [`proteins_align()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md)
  : Align protein sequences with mzLib, and map every residue onto the
  other sequence
- [`proteins_alignments()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_alignments.md)
  : One alignment per pair, each with its residue map
- [`proteins_annotate_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)
  : Annotate MetaMorpheus protein groups with Gene Ontology terms,
  keeping every member
- [`proteins_build_orthology()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md)
  : Build an Ensembl Compara orthology snapshot for your species, and
  write it as Parquet
- [`proteins_classify_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md)
  : Classify peptides by how widely they are shared across protein
  databases, with I = L
- [`proteins_map_residues()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_residues.md)
  : Map residues across species through Ensembl Compara's gene-tree
  alignments
- [`proteins_map_source_to_target()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md)
  : The target position paired with each source position
- [`proteins_map_target_to_source()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_target_to_source.md)
  : The source position paired with each target position
- [`proteins_orthology_pair()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_orthology_pair.md)
  : One species pair's counts, in either order
- [`proteins_parse_accessions()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_parse_accessions.md)
  : Read accessions into entry, isoform or version, and applied
  variants, with mzLib
- [`proteins_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)
  : Read protein databases: one row per protein, with GO terms and
  Ensembl genes on request
- [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
  : Resolve every protein to stable Ensembl gene ids, against a gene set
  you pin
- [`proteins_update_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_update_go.md)
  : Download the current Gene Ontology release to a go.obo, on purpose

## PRIDE Archive

List a project’s files, download the ones you choose, and check a
download against the MD5 and exact size PRIDE recorded.

- [`pride_approximate_total_size_bytes()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_approximate_total_size_bytes.md)
  : Sum the approximate sizes of some FTP files
- [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  : Download files from a PRIDE Archive project
- [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  : Download exactly the files you selected
- [`pride_file_checksums()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_file_checksums.md)
  : The MD5 and exact size PRIDE recorded for each file deposited in a
  project
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

Digest a UniProt entry and fragment its peptidoforms, and convert full
sequences to Unimod or ProForma.

- [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  : Digest a protein and fragment its peptidoforms
- [`peptidoform_convert()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_convert.md)
  : Convert full sequences to another notation with mzLib, one row per
  input
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

## Differential abundance

limma’s moderated t-test with contrasts, Benjamini-Hochberg adjustment,
random-effects meta-analysis, combined p-values, Spearman correlation,
logistic and mixed models per feature (also on a latent scale, from
measured and counted samples), and entrapment FDP estimates, computed by
mzLib.

- [`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
  : Benjamini-Hochberg adjust p-values, keeping every position
- [`stats_combine_p()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_combine_p.md)
  : Combine several independent p-values per feature into one, by
  Fisher's or Stouffer's method
- [`stats_correlate()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_correlate.md)
  : Spearman-correlate every feature with a per-sample trait, optionally
  given covariates
- [`stats_entrapment_fdp()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_entrapment_fdp.md)
  : Estimate a search's false discovery proportion from its entrapment
  hits
- [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  : Fit one linear model per feature and test coefficients, and
  contrasts, with limma's moderated t
- [`stats_fit_latent_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_latent_mixed.md)
  : Fit a mixed model per feature on a latent scale, from samples
  measured on it, counted through it, or both
- [`stats_fit_link()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_link.md)
  : Fit the detection link: how a cell's latent value sets the chance
  that a trial succeeds
- [`stats_fit_logistic()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_logistic.md)
  : Fit one logistic regression per feature of a 0/1 table, such as a
  detection table
- [`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)
  : Fit one linear mixed model per feature, with a random intercept for
  each group of samples
- [`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
  : Pool one effect size per study into a random-effects estimate per
  feature

## Differential design and observations

A sample sheet to mzLib’s design matrix and contrasts per stratum, and a
peptide table to every peptide in every sample, with the reason for each
missing value.

- [`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md)
  : Build each stratum's design matrix and default contrasts from a
  sample sheet, with mzLib
- [`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)
  : Every peptide in every biological sample, as log2 intensity or NA
  with the reason, by mzLib
- [`factor_spec()`](https://smith-chem-wisc.github.io/mzLibR/reference/factor_spec.md)
  : A condition factor of a model, as mzLib's FactorSpec
- [`covariate_spec()`](https://smith-chem-wisc.github.io/mzLibR/reference/covariate_spec.md)
  : A numeric covariate of a model, as mzLib's CovariateSpec

## PTM site occupancy and PTM-QTL

Every site’s occupancy in every run, pairs of modifications, and a
trait’s effect on occupancy, computed by mzLib.

- [`ptmqtl_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md)
  : Compute every site's occupancy in every run, with mzLib's
  SiteOccupancyCalculator
- [`ptmqtl_pairs()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_pairs.md)
  : Find PTM-PTM pairs within one dataset, with mzLib's PtmPairEngine
- [`ptmqtl_trait_effects()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md)
  : Fit, per site, the effect of a per-run trait on the log odds of
  occupancy, with mzLib's SiteTraitEffects

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

- [`print(`*`<mzlibr_adjusted>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_adjusted.md)
  : Print Benjamini-Hochberg adjusted p-values
- [`print(`*`<mzlibr_alignment>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_alignment.md)
  : Print one pair's alignment
- [`print(`*`<mzlibr_alignments>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_alignments.md)
  : Print protein alignments
- [`print(`*`<mzlibr_batch>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_batch.md)
  : Print a batch read from many files
- [`print(`*`<mzlibr_census>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_census.md)
  : Print a modification census
- [`print(`*`<mzlibr_combined_p_values>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_combined_p_values.md)
  : Print combined p-values
- [`print(`*`<mzlibr_correlations>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_correlations.md)
  : Print Spearman correlations
- [`print(`*`<mzlibr_covariate_spec>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_covariate_spec.md)
  : Print a covariate of a model
- [`print(`*`<mzlibr_design_stratum>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_design_stratum.md)
  : Print one stratum's design
- [`print(`*`<mzlibr_detection_link>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_detection_link.md)
  : Print a detection link
- [`print(`*`<mzlibr_diff_design>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_diff_design.md)
  : Print a design
- [`print(`*`<mzlibr_digest>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_digest.md)
  : Print a digest
- [`print(`*`<mzlibr_entrapment_fdp_estimates>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_entrapment_fdp_estimates.md)
  : Print entrapment FDP estimates
- [`print(`*`<mzlibr_factor_spec>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_factor_spec.md)
  : Print a factor of a model
- [`print(`*`<mzlibr_feature_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_feature_records.md)
  : Print an MS1 feature table
- [`print(`*`<mzlibr_file_info>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_file_info.md)
  : Print a file identification
- [`print(`*`<mzlibr_gene_resolutions>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_gene_resolutions.md)
  : Print a gene resolution
- [`print(`*`<mzlibr_go_annotations>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_go_annotations.md)
  : Print a Gene Ontology annotation
- [`print(`*`<mzlibr_go_update>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_go_update.md)
  : Print a Gene Ontology update
- [`print(`*`<mzlibr_isobaric_kits>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_isobaric_kits.md)
  : Print isobaric kits
- [`print(`*`<mzlibr_latent_mixed_model_fit>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_latent_mixed_model_fit.md)
  : Print a latent-scale mixed-model fit
- [`print(`*`<mzlibr_library_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_library_records.md)
  : Print a spectral library read
- [`print(`*`<mzlibr_link_fit>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_link_fit.md)
  : Print a fitted detection link
- [`print(`*`<mzlibr_logistic_model_fit>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_logistic_model_fit.md)
  : Print a logistic fit
- [`print(`*`<mzlibr_match_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_match_records.md)
  : Print a spectral-match table
- [`print(`*`<mzlibr_median_polish>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_median_polish.md)
  : Print a median-polish result
- [`print(`*`<mzlibr_meta_analysis>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_meta_analysis.md)
  : Print a meta-analysis
- [`print(`*`<mzlibr_mixed_model_fit>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_mixed_model_fit.md)
  : Print a mixed-model fit
- [`print(`*`<mzlibr_moderated_contrasts>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_moderated_contrasts.md)
  : Print the contrasts of a moderated fit
- [`print(`*`<mzlibr_moderated_fit>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_moderated_fit.md)
  : Print a moderated fit
- [`print(`*`<mzlibr_native_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_native_records.md)
  : Print a native record table
- [`print(`*`<mzlibr_observations>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_observations.md)
  : Print an observation table
- [`print(`*`<mzlibr_occupancy_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_occupancy_records.md)
  : Print a site-occupancy table
- [`print(`*`<mzlibr_orthology_snapshot>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_orthology_snapshot.md)
  : Print an orthology snapshot
- [`print(`*`<mzlibr_parsed_accessions>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_parsed_accessions.md)
  : Print parsed accessions
- [`print(`*`<mzlibr_peptide_classification>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_peptide_classification.md)
  : Print a peptide classification
- [`print(`*`<mzlibr_pooled_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_pooled_sdrf.md)
  : Print a pooled SDRF table
- [`print(`*`<mzlibr_protein_database>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_protein_database.md)
  : Print a protein database read
- [`print(`*`<mzlibr_protein_group_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_protein_group_records.md)
  : Print a protein-group table
- [`print(`*`<mzlibr_ptm_pairs>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_ptm_pairs.md)
  : Print PTM pairs
- [`print(`*`<mzlibr_quant>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_quant.md)
  : Print a quantification result
- [`print(`*`<mzlibr_quantified_peptide_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_quantified_peptide_records.md)
  : Print a quantified-peptide table
- [`print(`*`<mzlibr_residue_mappings>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_residue_mappings.md)
  : Print residue mappings
- [`print(`*`<mzlibr_restricted_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_restricted_sdrf.md)
  : Print a restricted SDRF
- [`print(`*`<mzlibr_result_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_result_records.md)
  : Print a record view
- [`print(`*`<mzlibr_rt_calibration>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_rt_calibration.md)
  : Print a retention-time calibration
- [`print(`*`<mzlibr_scan_records>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_scan_records.md)
  : Print a scan table
- [`print(`*`<mzlibr_sdrf>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf.md)
  : Print an SDRF document
- [`print(`*`<mzlibr_sdrf_ages>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_ages.md)
  : Print parsed SDRF ages
- [`print(`*`<mzlibr_sdrf_assessment>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_assessment.md)
  : Print an SDRF informativeness assessment
- [`print(`*`<mzlibr_sdrf_design>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_design.md)
  : Print a label-free design
- [`print(`*`<mzlibr_sdrf_drift>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_drift.md)
  : Print SDRF drift findings
- [`print(`*`<mzlibr_sdrf_isobaric_design>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_isobaric_design.md)
  : Print an isobaric design
- [`print(`*`<mzlibr_sdrf_samples>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_samples.md)
  : Print an SDRF document's samples
- [`print(`*`<mzlibr_sdrf_validation>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sdrf_validation.md)
  : Print an SDRF validation
- [`print(`*`<mzlibr_sequence_conversions>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_sequence_conversions.md)
  : Print sequence conversions
- [`print(`*`<mzlibr_site_occupancy>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_site_occupancy.md)
  : Print site occupancy
- [`print(`*`<mzlibr_trait_effects>`*`)`](https://smith-chem-wisc.github.io/mzLibR/reference/print.mzlibr_trait_effects.md)
  : Print trait effects
