# PRIDE Archive access, backed by mzLib's `PrideArchiveClient`.
#
# The PRIDE Archive (https://www.ebi.ac.uk/pride/archive/) is EBI's public proteomics data
# repository. This module lists what is in a project and pulls files down, using the same
# paging, URL resolution and safe-download logic mzLib uses in C#.
#
# Everything here splits into pure functions the tests can call directly — `pride_build_*()`
# for argument assembly and `pride_parse_*()` for wire-to-data.frame — so the offline suite
# needs no subprocess and no mocking.

# ---------------------------------------------------------------- validation

# Validate and canonicalise an accession, failing loudly rather than returning nothing.
#
# Accessions are upper-cased because PRIDE's API is case-sensitive on the accession while its
# category matching is case-insensitive. Two rules pointing opposite ways is a trap, and this
# is the one that can be fixed without surprising anybody.
#
# The grammar is a short letter prefix and a run of digits, and it is **grammatical only**:
# "PXD0000019999" is well-formed and costs a live round trip before it fails, so no offline
# check will catch a transposed digit. That is deliberate — PXD accessions are not fixed-width
# forever, and rejecting a valid future accession would be worse than one wasted request.
pride_normalise_accession <- function(accession) {
  if (!is.character(accession) || length(accession) != 1L || is.na(accession)) {
    stop(mzlib_usage_error(
      "A PRIDE project accession is required, as a single string, e.g. 'PXD000001'."
    ))
  }

  candidate <- toupper(trimws(accession))
  if (!nzchar(candidate)) {
    stop(mzlib_usage_error(
      "A PRIDE project accession is required, e.g. 'PXD000001'."
    ))
  }

  if (!grepl("^[A-Z]{2,4}[0-9]{4,}$", candidate)) {
    stop(mzlib_usage_error(paste0(
      "'", accession, "' is not a valid repository accession. Expected a short letter prefix ",
      "followed by digits, e.g. 'PXD000001'."
    )))
  }

  candidate
}

# Reject a blank destination instead of quietly writing into the working directory.
#
# An empty path *is* the working directory, so a destination that came from an unset variable
# would spray a multi-gigabyte project across wherever R happens to be sitting.
pride_normalise_destination <- function(destination) {
  if (!is.character(destination) || length(destination) != 1L || is.na(destination) ||
    !nzchar(trimws(destination))) {
    stop(mzlib_usage_error(
      "A destination directory is required; got an empty path."
    ))
  }
  destination
}

# Accept a vector of extensions, refusing one that normalises to nothing.
#
# A caller who asked for a filter and whose list normalises to nothing would have had `--ext`
# omitted entirely, which the bridge reads as "no filter" and downloads the whole project.
# Asking for a filter and getting everything is never right.
pride_normalise_extensions <- function(extensions) {
  if (is.null(extensions) || length(extensions) == 0L) {
    return(character(0))
  }
  if (!is.character(extensions)) {
    stop(mzlib_usage_error("extensions must be a character vector, e.g. c('.raw', '.mzML')."))
  }
  if (any(grepl(",", extensions, fixed = TRUE))) {
    stop(mzlib_usage_error(paste0(
      "An extension may not contain a comma; got '",
      extensions[grepl(",", extensions, fixed = TRUE)][1L],
      "'. Pass separate vector elements."
    )))
  }

  kept <- trimws(extensions[!is.na(extensions)])
  kept <- kept[nzchar(kept)]
  if (length(kept) == 0L) {
    stop(mzlib_usage_error(paste0(
      "extensions was given but names no extensions; got ",
      paste(deparse(extensions), collapse = " "),
      ". Omit it to download every file type."
    )))
  }
  kept
}

# Refuse a value the bridge's parser would read as another option.
#
# The bridge treats `--a --b` as two flags, so a value beginning with '-' silently discards the
# option it belonged to — and can smuggle in a flag the caller never intended.
pride_reject_flag_like <- function(name, value) {
  if (startsWith(value, "-")) {
    stop(mzlib_usage_error(paste0(
      name, " may not begin with '-'; got '", value, "'. That would be read as another option."
    )))
  }
  value
}

# ---------------------------------------------------------------- derived values

# A file's lowercase extension including the dot, or "" if it has none.
#
# `"x.mgf.gz"` has the extension `".gz"`, which is the single most common way to end up with an
# empty result. See the `extensions` argument of [pride_download()], where the mistake is
# actually made.
pride_extension <- function(file_name) {
  base <- basename(file_name)
  at <- regexpr("\\.[^.]*$", base)
  extension <- ifelse(at > 1L, tolower(substring(base, at)), "")
  extension[is.na(file_name)] <- NA_character_
  extension
}

# Parse one ISO-8601 timestamp, treating anything unreadable as absent.
#
# A timestamp that will not parse is not worth failing a whole manifest over — the caller asked
# for a file list, not a date.
pride_parse_timestamp <- function(value) {
  absent <- .POSIXct(NA_real_, tz = "UTC")
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    return(absent)
  }

  text <- trimws(value)
  # `%z` wants "+0000". The bridge writes "+00:00" and sometimes "Z"; strptime accepts neither
  # on every platform, and the ones it silently fails on differ between Windows and Linux.
  text <- sub("Z$", "+0000", text)
  text <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", text)

  for (format in c("%Y-%m-%dT%H:%M:%OS%z", "%Y-%m-%dT%H:%M:%OS", "%Y-%m-%d")) {
    parsed <- as.POSIXct(text, format = format, tz = "UTC")
    if (!is.na(parsed)) {
      return(parsed)
    }
  }
  absent
}

# The vector form, kept as POSIXct rather than a list of them.
pride_parse_timestamps <- function(values) {
  seconds <- vapply(
    values,
    function(value) as.numeric(pride_parse_timestamp(value)),
    numeric(1L),
    USE.NAMES = FALSE
  )
  .POSIXct(seconds, tz = "UTC")
}

# Pull one scalar field out of a parsed file object, with an explicit type and an explicit
# answer for absent and for JSON `null`.
#
# This is where the reader's `NA` convention is turned into the `NA` of the right type. Doing it
# per field, in the module that knows what the field means, is the whole reason `json.R` refuses
# to guess: `checksum` absent means "the repository publishes none" and belongs as `""`, while
# `https_url` absent means "this file cannot be fetched over HTTPS" and must stay `NA` so that
# `downloadable` is honest.
pride_field <- function(entry, name, type, missing) {
  value <- entry[[name]]
  if (is.null(value) || (length(value) == 1L && is.logical(value) && is.na(value))) {
    return(missing)
  }
  switch(type,
    character = as.character(value),
    numeric = as.numeric(value),
    value
  )
}

# ---------------------------------------------------------------- parsing

# Turn one file's `locations` array into a small data.frame of controlled-vocabulary terms.
pride_parse_locations <- function(entry) {
  locations <- entry[["locations"]]
  if (!is.list(locations) || length(locations) == 0L) {
    return(data.frame(
      accession = character(0), name = character(0), value = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    accession = vapply(locations, pride_field, character(1L), "accession", "character", ""),
    name = vapply(locations, pride_field, character(1L), "name", "character", ""),
    value = vapply(locations, pride_field, character(1L), "value", "character", ""),
    stringsAsFactors = FALSE
  )
}

# Turn the `pride files` payload into a data.frame, refusing an empty manifest.
pride_parse_manifest <- function(data, accession) {
  entries <- if (is.list(data)) data[["files"]] else NULL
  if (!is.list(entries)) {
    entries <- list()
  }

  if (length(entries) == 0L) {
    # PRIDE answers an unknown accession with HTTP 200 and an empty list, not a 404. Handing
    # back a zero-row data.frame would let a typo produce "0 files, done" and a green exit.
    stop(mzlib_project_not_found(paste0(
      "PRIDE returned no files for '", accession, "'. Either the accession does not exist ",
      "(check for a typo) or the project is private. PRIDE does not distinguish the two, so ",
      "neither can mzLibR."
    )))
  }

  file_name <- vapply(entries, pride_field, character(1L), "file_name", "character", NA_character_)
  https_url <- vapply(entries, pride_field, character(1L), "https_url", "character", NA_character_)
  size <- vapply(entries, pride_field, numeric(1L), "file_size_bytes", "numeric", NA_real_)

  files <- data.frame(
    file_name = file_name,
    file_size_bytes = size,
    size_mb = size / 1e6,
    extension = pride_extension(file_name),
    category = vapply(entries, pride_field, character(1L), "category", "character", ""),
    category_accession = vapply(
      entries, pride_field, character(1L), "category_accession", "character", ""
    ),
    checksum = vapply(entries, pride_field, character(1L), "checksum", "character", ""),
    https_url = https_url,
    # `NA` and not `FALSE` for a missing URL would be defensible, but this column exists to be
    # used as a filter and `files[files$downloadable, ]` on an NA yields a row of NAs. A file
    # with no HTTPS location simply cannot be fetched, which is a fact, not an unknown.
    downloadable = !is.na(https_url),
    submission_date = pride_parse_timestamps(
      vapply(entries, pride_field, character(1L), "submission_date", "character", NA_character_)
    ),
    publication_date = pride_parse_timestamps(
      vapply(entries, pride_field, character(1L), "publication_date", "character", NA_character_)
    ),
    updated_date = pride_parse_timestamps(
      vapply(entries, pride_field, character(1L), "updated_date", "character", NA_character_)
    ),
    # Not on the wire — stamped on here so `pride_download_files()` can tell which project to
    # fetch from without the caller having to carry the accession alongside the frame.
    project_accession = accession,
    stringsAsFactors = FALSE
  )

  # A list column rather than an attribute. Attributes are dropped by `[.data.frame`, so a
  # biologist writing `files[files$size_mb < 5, ]` would silently lose the locations; a list
  # column subsets along with its rows like everything else.
  files$locations <- lapply(entries, pride_parse_locations)
  files
}

# Read the written paths out of a `pride download` payload.
pride_parse_paths <- function(data) {
  paths <- if (is.list(data)) data[["paths"]] else NULL
  if (!is.list(paths) || length(paths) == 0L) {
    return(character(0))
  }
  vapply(paths, function(path) as.character(path)[1L], character(1L), USE.NAMES = FALSE)
}

# Refuse to report success for a filter that matched nothing.
#
# A filter that matched nothing is nearly always a filter that does not mean what its author
# thought, and reporting success with an empty result lets a batch script carry on as though
# the work had been done.
pride_check_filter_matched <- function(written, accession, category, extensions) {
  filtered <- !is.null(category) || length(extensions) > 0L
  if (length(written) > 0L || !filtered) {
    return(invisible(NULL))
  }

  described <- character(0)
  if (!is.null(category)) {
    described <- c(described, paste0("category '", category, "'"))
  }
  if (length(extensions) > 0L) {
    described <- c(described, paste0("extensions ", paste0("'", extensions, "'", collapse = ", ")))
  }

  stop(mzlib_usage_error(paste0(
    "No file in ", accession, " matched ", paste(described, collapse = " and "),
    ". Use pride_list_files() to see what the project actually contains — note that ",
    "compressed files such as 'x.mgf.gz' have the extension '.gz', not '.mgf'."
  )))
}

# ---------------------------------------------------------------- argument assembly

pride_build_list_args <- function(accession, page_size) {
  canonical <- pride_normalise_accession(accession)

  if (!is.numeric(page_size) || length(page_size) != 1L || is.na(page_size) ||
    page_size != round(page_size)) {
    stop(mzlib_usage_error("page_size must be a single whole number."))
  }
  if (page_size < 1) {
    stop(mzlib_usage_error(paste0("page_size must be positive; got ", page_size, ".")))
  }
  if (page_size > 2147483647) {
    stop(mzlib_usage_error(paste0(
      "page_size is larger than the API allows; got ", format(page_size, scientific = FALSE), "."
    )))
  }

  c("pride", "files", "--accession", canonical,
    "--page-size", format(page_size, scientific = FALSE))
}

pride_build_download_args <- function(accession, destination, category, extensions, overwrite) {
  canonical <- pride_normalise_accession(accession)
  target <- pride_normalise_destination(destination)
  wanted <- pride_normalise_extensions(extensions)

  args <- c("pride", "download", "--accession", canonical, "--dest", target)

  if (!is.null(category)) {
    if (!is.character(category) || length(category) != 1L || is.na(category) ||
      !nzchar(trimws(category))) {
      stop(mzlib_usage_error(paste0(
        "category is empty. Omit it to download every category, rather than passing a blank ",
        "value — a filter that selects nothing must not silently select everything."
      )))
    }
    args <- c(args, "--category", pride_reject_flag_like("category", trimws(category)))
  }

  if (length(wanted) > 0L) {
    args <- c(args, "--ext", pride_reject_flag_like("extensions", paste(wanted, collapse = ",")))
  }

  if (!isTRUE(overwrite)) {
    args <- c(args, "--no-overwrite")
  }

  args
}

pride_build_download_files_args <- function(files, destination, overwrite) {
  target <- pride_normalise_destination(destination)

  if (!is.data.frame(files)) {
    stop(mzlib_usage_error(
      "files must be a data.frame from pride_list_files(), or a subset of one."
    ))
  }
  required <- c("file_name", "https_url", "downloadable", "project_accession")
  absent <- setdiff(required, names(files))
  if (length(absent) > 0L) {
    stop(mzlib_usage_error(paste0(
      "files is missing the column(s) ", paste(absent, collapse = ", "),
      ". Pass a data.frame from pride_list_files(), or a subset of one."
    )))
  }

  if (nrow(files) == 0L) {
    stop(mzlib_usage_error(paste0(
      "No files selected. An empty selection is almost always a filter that did not match what ",
      "you expected, so mzLibR refuses it rather than reporting success."
    )))
  }

  unreachable <- files$file_name[!files$downloadable]
  if (length(unreachable) > 0L) {
    stop(mzlib_usage_error(paste0(
      length(unreachable), " of ", nrow(files), " selected files have no HTTPS location and ",
      "cannot be downloaded (e.g. '", unreachable[1L], "'). Filter on `downloadable` first: ",
      "files[files$downloadable, ]."
    )))
  }

  accessions <- unique(files$project_accession[nzchar(files$project_accession) &
    !is.na(files$project_accession)])
  if (length(accessions) == 0L) {
    stop(mzlib_usage_error(paste0(
      "These rows carry no project accession, so mzLibR cannot tell which project to fetch ",
      "from. Obtain them from pride_list_files()."
    )))
  }
  if (length(accessions) > 1L) {
    stop(mzlib_usage_error(paste0(
      "All files must come from one project; got ", paste(accessions, collapse = ", "), "."
    )))
  }

  # The selection travels on stdin rather than argv: a few thousand names would blow the ~32 KB
  # command-line ceiling. The framing is newline-delimited, which is *almost* general — a POSIX
  # file name may legally contain a newline, so such a name would split in two and silently
  # select the wrong files. PRIDE has never published one, but "never seen it" is not a
  # contract, so it is refused explicitly rather than mis-parsed quietly.
  broken <- files$file_name[grepl("[\r\n]", files$file_name)]
  if (length(broken) > 0L) {
    stop(mzlib_usage_error(paste0(
      "Cannot select '", encodeString(broken[1L]), "': the file name contains a line break, ",
      "which the selection format cannot represent. Please open an issue — this is worth ",
      "fixing properly if a real repository ever publishes such a name."
    )))
  }

  args <- c(
    "pride", "download", "--accession", accessions[1L],
    "--dest", target, "--names-from-stdin"
  )
  if (!isTRUE(overwrite)) {
    args <- c(args, "--no-overwrite")
  }

  list(args = args, stdin = files$file_name)
}

# ---------------------------------------------------------------- the public surface

#' List the files in a PRIDE Archive project
#'
#' Paging is handled for you: however many pages the project spans, you get one data.frame.
#'
#' @section What this manifest is, and is not:
#'
#' **This is what PRIDE's REST API publishes, which is not always everything in the project.**
#' For PXD000001 the API returns **8** files while the FTP tree holds **13** — and the five it
#' omits include the two largest, a 450 MB `.mzML` and the matching 472 MB `.mzXML`, which are
#' exactly the modern open-format conversions most people want. The omission is PRIDE's, not
#' mzLib's.
#'
#' So a manifest that looks short may be short. If completeness matters — mirroring a project,
#' budgeting a download, proving you analysed everything — cross-check the FTP directory at
#' `https://ftp.pride.ebi.ac.uk/pride/data/archive/<year>/<month>/<accession>/`.
#'
#' @param accession A project accession, e.g. `"PXD000001"`. Case and surrounding whitespace are
#'   normalised. The check is grammatical only: a well-formed accession that does not exist
#'   costs one live request before it fails, because PXD accessions are not fixed-width forever
#'   and rejecting a valid future one would be worse.
#' @param page_size How many files to request per underlying API call. Affects only how the
#'   manifest is fetched, never what you get back.
#' @param timeout Seconds to allow for the whole fetch, or `NULL` to wait indefinitely.
#'
#' @return A data.frame with one row per file and the columns `file_name`, `file_size_bytes`,
#'   `size_mb`, `extension`, `category`, `category_accession`, `checksum`, `https_url`,
#'   `downloadable`, `submission_date`, `publication_date`, `updated_date`,
#'   `project_accession`, and a `locations` list column (see [pride_locations()]).
#'
#'   `file_size_bytes` is **the size PRIDE reports, which is not always what you will
#'   transfer**: for compressed files it is frequently the *decompressed* size. In PXD000001 the
#'   reported size of `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz` is 16,448,103 bytes and the
#'   actual download is 5,984,662 — **2.75x** smaller. The `.mztab.gz` behaves the same way; the
#'   `.xml.gz` does not. PRIDE's own metadata is inconsistent here, so neither mzLibR nor mzLib
#'   can correct it. Treat the sum as an upper bound on transfer, and see
#'   [pride_total_size_bytes()].
#'
#' @seealso [pride_download_files()], which is usually what you want next.
#' @export
pride_list_files <- function(accession, page_size = 100, timeout = 300) {
  args <- pride_build_list_args(accession, page_size)
  canonical <- pride_normalise_accession(accession)
  data <- bridge_invoke(args, timeout = timeout)
  pride_parse_manifest(data, canonical)
}

#' Download files from a PRIDE Archive project
#'
#' Files are streamed to a temporary name and moved into place only once complete, so an
#' interrupted download never leaves a truncated file behind.
#'
#' Prefer [pride_download_files()] when you can. `category` and `extensions` can only express
#' what they were built to express; "under 5 MB", "the three newest", or "everything except the
#' MGF" cannot be said in that vocabulary at all — and can all be said with `[`.
#'
#' @param accession A project accession, e.g. `"PXD000001"`.
#' @param dest Directory to write into. A blank path is refused rather than silently taken to
#'   mean the working directory.
#' @param category Keep only files of this category, e.g. `"RAW"`, `"PEAK"`, `"SEARCH"`,
#'   `"OTHER"`. `NULL` keeps all.
#'
#'   Categories are coarser than they look. In PXD000001 `"PEAK"` matches **2** files, not 1 —
#'   the 6 MB MGF *and* a 243 MB mzXML — so a download you expected to be small is 40x bigger.
#'   Check against [pride_list_files()] first.
#' @param extensions Keep only files with these extensions, e.g. `c(".raw", ".mzML")`. `NULL`
#'   keeps all. Combined with `category` as AND.
#'
#'   **A compressed file's extension is `.gz`, not what it is compressed from.** PXD000001's
#'   peak list is `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`, so `".mgf"` matches **nothing**,
#'   while `".gz"` over-matches to three unrelated files. To select one compressed type, combine
#'   `category = "PEAK"` with `extensions = ".gz"`, or skip these filters and pass the rows you
#'   want to [pride_download_files()].
#'
#'   A filter that matches nothing raises an error rather than reporting success, because an
#'   empty result here is nearly always a filter that does not mean what its author thought.
#' @param overwrite When `FALSE`, a file already present at the destination is left alone and
#'   not re-fetched — a cheap resume for a large project.
#' @param timeout Seconds to allow, or `NULL` to wait as long as it takes. `NULL` is the default
#'   because multi-gigabyte projects legitimately run for hours.
#'
#' @return A character vector of the paths where the files now are.
#' @export
pride_download <- function(accession, dest, category = NULL, extensions = NULL,
                           overwrite = TRUE, timeout = NULL) {
  args <- pride_build_download_args(accession, dest, category, extensions, overwrite)
  canonical <- pride_normalise_accession(accession)
  written <- pride_parse_paths(bridge_invoke(args, timeout = timeout))
  pride_check_filter_matched(
    written, canonical, category, pride_normalise_extensions(extensions)
  )
  written
}

#' Download exactly the files you selected
#'
#' The counterpart to [pride_list_files()], and usually the one you want: filter the data.frame
#' however you like, then pass the rows.
#'
#' ```r
#' files <- pride_list_files("PXD000001")
#' small <- files[files$size_mb < 5 & files$downloadable, ]
#' pride_download_files(small, "downloads")
#' ```
#'
#' @param files Rows of a [pride_list_files()] data.frame. They must all come from one project.
#'   Rows whose `downloadable` is `FALSE` are refused up front rather than failing halfway
#'   through a multi-gigabyte transfer.
#' @param dest Directory to write into.
#' @param overwrite When `FALSE`, files already present are left alone.
#' @param timeout Seconds to allow, or `NULL` to wait as long as it takes.
#'
#' @return A character vector of paths. These say **where each file is**, not what was
#'   transferred just now: with `overwrite = FALSE` a file already present is left alone and its
#'   path is still returned. Do not read `length()` of this as work done.
#' @export
pride_download_files <- function(files, dest, overwrite = TRUE, timeout = NULL) {
  prepared <- pride_build_download_files_args(files, dest, overwrite)
  pride_parse_paths(bridge_invoke(prepared$args, stdin = prepared$stdin, timeout = timeout))
}

#' The published locations of each file, as controlled-vocabulary terms
#'
#' One long data.frame rather than a nested list, so it pipes into anything.
#'
#' You rarely need this: `https_url` already carries the fetchable URL, resolved by mzLib's
#' `TryGetHttpsDownloadUrl`, which **searches** the locations rather than taking the first.
#' That distinction matters — **the order is not stable.** In PXD000001 the mztab lists FTP
#' first while the MGF lists **Aspera** first, so code that took `locations[[1]]` would get an
#' unfetchable `prd_ascp@fasp.ebi.ac.uk:...` address for some files and a working one for
#' others, in the same project. Do not re-implement the search.
#'
#' @param files A [pride_list_files()] data.frame.
#' @return A data.frame with `file_name`, `accession`, `name` and `value`, one row per location.
#' @export
pride_locations <- function(files) {
  if (!is.data.frame(files) || !"locations" %in% names(files)) {
    stop(mzlib_usage_error("files must be a data.frame from pride_list_files()."))
  }
  if (nrow(files) == 0L) {
    return(data.frame(
      file_name = character(0), accession = character(0),
      name = character(0), value = character(0), stringsAsFactors = FALSE
    ))
  }

  parts <- lapply(seq_len(nrow(files)), function(row) {
    locations <- files$locations[[row]]
    if (!is.data.frame(locations) || nrow(locations) == 0L) {
      return(NULL)
    }
    data.frame(
      file_name = files$file_name[row],
      accession = locations$accession,
      name = locations$name,
      value = locations$value,
      stringsAsFactors = FALSE
    )
  })

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(data.frame(
      file_name = character(0), accession = character(0),
      name = character(0), value = character(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, parts)
}

#' Total size of a set of PRIDE files
#'
#' @param files A [pride_list_files()] data.frame, or a subset.
#' @return A single number of bytes.
#'
#' @section Two reasons this is not the size of the project:
#'
#' **It over-reports compressed files.** PRIDE frequently gives the *decompressed* size — the
#' MGF in PXD000001 reports 16,448,103 bytes and downloads as 5,984,662, a factor of **2.75**.
#'
#' **It sums an incomplete manifest.** For PXD000001 this returns **0.51 GB**; the project on
#' disk is **1.44 GB**, because PRIDE's API omits five files including the two largest (see
#' [pride_list_files()]).
#'
#' The two errors run in opposite directions and do **not** cancel: compressed sizes are
#' inflated, whole files are missing entirely.
#' @export
pride_total_size_bytes <- function(files) {
  if (!is.data.frame(files) || !"file_size_bytes" %in% names(files)) {
    stop(mzlib_usage_error("files must be a data.frame from pride_list_files()."))
  }
  sum(files$file_size_bytes, na.rm = TRUE)
}
