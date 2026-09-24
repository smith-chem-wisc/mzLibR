# Quantification: FlashLFQ. See scripts/spec-facts.R for what each list means.

R_DEVIATIONS[["quant flashlfq"]] <- list(
  "param.stdin" = as_r("spectra", "one stdin line per run, with its design columns"),
  "param.ppm" = as_r("ppm_tolerance", "mzLib's parameter name"),
  "param.isotope-ppm" = as_r("isotope_ppm_tolerance", "mzLib's parameter name"),
  "param.mbr" = as_r("match_between_runs", "mzLib's parameter name"),
  "param.mbr-ppm" = as_r("mbr_ppm_tolerance", "mzLib's parameter name"),
  "param.mbr-q" = as_r("mbr_q_value_threshold", "mzLib's parameter name"),
  "param.shared-peptides" = as_r("use_shared_peptides_for_protein_quant", "mzLib's parameter name"),
  "param.bayesian" = as_r("bayesian_protein_quant", "mzLib's parameter name"),
  "param.use-pep-q" = as_r("use_pep_q_value", "mzLib's parameter name"),
  "param.threads" = as_r("max_threads", "mzLib's parameter name; the default differs, see its argument"),
  "param.out" = as_r("output_directory", "mzLib's parameter name"),
  "field.peptide_count" = not_here("flashlfq_peptide_count() of the result"),
  "field.protein_count" = not_here("flashlfq_protein_count() of the result"),
  "field.peptides.intensities" = as_r("intensity", "unnested: one row per peptide per run, the run in file_name"),
  "field.peptides.detection_types" = as_r("detection_type", "unnested beside intensity"),
  "field.proteins.intensities" = as_r("intensity", "unnested: one row per protein group per sample, the sample in file_name")
)

PARENT_MAP[c(
  "flashlfq_quantify", "flashlfq_peptide_count", "flashlfq_protein_count",
  "flashlfq_mbr_peak_count", "flashlfq_mbr_peaks", "flashlfq_mbr_rescued_peptide_count"
)] <- c(
  "flashlfq.quantify", "FlashLfqResults.peptide_count", "FlashLfqResults.protein_count",
  "FlashLfqResults.mbr_peak_count", "FlashLfqResults.mbr_peaks",
  "FlashLfqResults.mbr_rescued_peptide_count"
)
PARENT_OMISSIONS[c(
  "Peptide.intensity", "Peptide.detection_type", "ProteinGroup.intensity", "Peak.is_mbr",
  "flashlfq.median_polish"
)] <- c(
  "a row of the long `peptides` frame",
  "a column of the long `peptides` frame",
  "a row of the long `proteins` frame",
  "a column comparison, `detection_type == \"MBR\"`",
  "parity debt: `quant median-polish` arrives with the mzLib 1.0.592 port"
)

FIELD_CHECKS[["quant flashlfq"]] <- list("flashlfq_small.json", function(d) mz$flashlfq_parse(d))
