# Isobaric labelling kits: isobaric kits. See scripts/spec-facts.R for what each list means.

R_TABLE["isobaric kits"] <- "records"
R_DEVIATIONS[["isobaric kits"]] <- list(
  "field.columns" = RECORDS,
  "field.channel_index" = as_r("channel_index", "1-based in R, like every position"),
  "field.kits" = as_r("kits", "a data.frame of kit and channel_count, one row per kit")
)

PARENT_MAP["isobaric_kits"] <- "isobaric.kits"
PARENT_OMISSIONS[c(
  "IsobaricKit.labels", "IsobaricKit.reporter_ion_mzs", "IsobaricKits.channels",
  "IsobaricKits.kits", "IsobaricKits.kit_named"
)] <- c(
  "`records$channel_label[records$kit == name]`",
  "`records$reporter_ion_mz[records$kit == name]`",
  "the channels are `records`, a data.frame",
  "`split(records, records$kit)`, and the per-kit summaries are `kits`",
  "`records[records$kit == name, ]`"
)

FIELD_CHECKS[["isobaric kits"]] <- list("isobaric_kits.json", function(d) mz$isobaric_parse_kits(d))
