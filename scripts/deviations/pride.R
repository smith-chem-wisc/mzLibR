# PRIDE Archive. See scripts/spec-facts.R for what each list means.

R_TABLE[c("pride files", "pride ftp-files", "pride search")] <- "(the returned data.frame itself)"

R_DEVIATIONS[["pride files"]] <- list(
  "field.accession" = as_r("project_accession", "stamped on every row, so a selection knows its project"),
  "field.file_count" = not_here("the result is the files data.frame; this is nrow() of it"),
  "field.total_size_bytes" = not_here("pride_total_size_bytes() of the result"),
  "field.files" = not_here("the result is this list itself, as a data.frame, one row per file")
)
R_DEVIATIONS[["pride ftp-files"]] <- list(
  "field.accession" = not_here("an FTP listing is not a download selection, so no row carries it"),
  "field.file_count" = not_here("the result is the files data.frame; this is nrow() of it"),
  "field.approximate_total_size_bytes" = not_here("pride_approximate_total_size_bytes() of the result"),
  "field.files" = not_here("the result is this list itself, as a data.frame, one row per file")
)
R_DEVIATIONS[["pride search"]] <- list(
  "field.keyword" = not_here("the result is the data.frame of hits; the keyword is the one you passed"),
  "field.result_count" = not_here("the result is the hits data.frame; this is nrow() of it"),
  "field.results" = not_here("the result is this list itself, as a data.frame, one row per project")
)
R_DEVIATIONS[["pride download"]] <- list(
  "param.dest" = as_r("destination", "spelled out"),
  "param.ext" = as_r("extensions", "spelled out; a character vector"),
  "param.no-overwrite" = as_r("overwrite", "stated positively: overwrite = FALSE sends --no-overwrite"),
  "param.names-from-stdin" = not_here("pride_download_files() sends a selection, and sets this itself"),
  "param.stdin" = not_here("pride_download_files() sends the selected file names"),
  "field.accession" = not_here("the result is the character vector of paths"),
  "field.destination_directory" = not_here("the result is the character vector of paths"),
  "field.downloaded_count" = not_here("length() of the result"),
  "field.paths" = not_here("the result is this character vector itself")
)
R_DEVIATIONS[["pride download (selection)"]] <- list(
  "param.accession" = not_here("taken from the rows' project_accession"),
  "param.category" = not_here("a selection is already filtered, by `[`"),
  "param.ext" = not_here("a selection is already filtered, by `[`"),
  "param.names-from-stdin" = not_here("always set: the selection travels on stdin"),
  "param.stdin" = as_r("files", "the file_name of each selected row, one per stdin line")
)

PARENT_MAP[c(
  "pride_list_files", "pride_list_ftp_files", "pride_download", "pride_download_files",
  "pride_total_size_bytes", "pride_approximate_total_size_bytes"
)] <- c(
  "pride.list_files", "pride.list_ftp_files", "pride.download", "pride.download_files",
  "pride.total_size_bytes", "pride.approximate_total_size_bytes"
)
PARENT_MAP["pride_search"] <- "pride.search"
PARENT_ADDITIONS["pride_locations"] <-
  "R-only: unnests the `locations` list column, which is pyMzLib's PrideFile.locations field"
PARENT_OMISSIONS[c(
  "PrideFile.size_mb", "PrideFile.extension", "PrideFile.downloadable", "PrideFile.as_dict",
  "PrideFtpFile.approximate_size_mb", "PrideFtpFile.extension", "PrideFtpFile.as_dict",
  "PrideProjectSearchResult.matched_fields", "PrideProjectSearchResult.as_dict"
)] <- c(
  "a column, `size_mb`, not a function",
  "a column, `extension`, not a function",
  "a column, `downloadable`, not a function",
  "meaningless in R: the data.frame is already the record",
  "a column, `approximate_size_mb`, not a function",
  "a column, `extension`, not a function",
  "meaningless in R: the data.frame is already the record",
  "a list column of pride_search(), `matched_fields`: the sorted names of `highlights`",
  "meaningless in R: the data.frame is already the record"
)

# pride ftp-files: the spec lists no example; this is the recording the offline tests use.
REPLAY_EXTRA[[length(REPLAY_EXTRA) + 1L]] <- c(verb = "pride ftp-files", fixture = "pride_ftp_PXD000001.json")

FIELD_CHECKS[["pride files"]] <- list("pride_PXD000001_files.json", function(d) mz$pride_parse_manifest(d, "PXD000001"))
FIELD_CHECKS[["pride ftp-files"]] <- list("pride_ftp_PXD000001.json", function(d) mz$pride_parse_ftp_files(d, "PXD000001"))
FIELD_CHECKS[["pride search"]] <- list("pride_search_plasmodium.json", function(d) mz$pride_parse_search(d))
