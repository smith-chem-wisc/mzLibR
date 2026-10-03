# Isobaric labelling kits - TMT, TMTpro, iTRAQ, DiLeu - with every channel's reporter-ion m/z.
#
#   isobaric_kits()   mzLib's IsobaricMassTag for each kit it can name, or for one
#
# An isobaric experiment is read out of the low-m/z reporter ions: one ion per channel, a few
# millidaltons apart. Getting a channel's label or its m/z wrong mislabels every sample after it,
# so this package keeps no table of its own: every row is mzLib's.
#
# Ported from pyMzLib's `isobaric.py`. What changes is only the projection into R: the channels are
# a data.frame, `records`, with `channel_index` 1-based as every position in R is; the per-kit
# summaries are a data.frame, `kits`. pyMzLib's grouped views (IsobaricKits.kits, kit_named, ...)
# are subsets of `records` here.

# The pyMzLib release whose bridge first dispatches `isobaric kits`.
ISOBARIC_KITS_SINCE <- "0.3.0"

isobaric_parse_kits <- function(data) {
  records <- readers_parse_records_table(data)
  if (is.null(records)) {
    names <- wire_strings(data[["column_names"]])
    records <- as.data.frame(structure(rep(list(logical(0)), length(names)), names = names),
      stringsAsFactors = FALSE, optional = TRUE
    )
  }
  if ("channel_index" %in% names(records)) {
    records$channel_index <- as.numeric(records$channel_index) + 1
  }
  kit <- data[["kit"]]
  structure(
    list(
      kit = if (wire_null(kit)) NA_character_ else as.character(kit),
      kit_count = as.numeric(wire_field(data, "kit_count", "numeric", NA_real_)),
      record_count = as.numeric(wire_field(data, "record_count", "numeric", NA_real_)),
      absolute_tolerance = as.numeric(wire_field(data, "absolute_tolerance", "numeric", NA_real_)),
      kits = wire_objects(data[["kits"]], empty = c("kit", "channel_count")),
      column_names = wire_strings(data[["column_names"]]),
      records = records,
      caveats = wire_strings(data[["caveats"]])
    ),
    class = "mzlibr_isobaric_kits"
  )
}

#' The isobaric kits mzLib can name, with every channel's label and reporter-ion m/z
#'
#' Calls mzLib's `IsobaricMassTag.TryGetIsobaricMassTag` for each `IsobaricMassTagType`, or for the
#' one `kit` resolves to through `IsobaricMassTag.TryGetTagType`. No file, no network.
#'
#' **Nothing in the table is typed in.** mzLib holds only the channel labels. Each m/z is a
#' `DI HCD` diagnostic-ion line of the kit's `Multiplex Label` modification in mzLib's embedded
#' `TMT.txt`, plus one proton, sorted ascending and paired with the labels by position - the same
#' m/z MetaMorpheus uses to quantify a TMT search.
#'
#' @param kit One kit, matched by mzLib's whole-name rule: case-insensitive, optionally followed by
#'   `" on <motif>"` (`"TMT10"`, `"TMT6-plex"`, `"iTRAQ-4plex on K"`). Never a substring, so
#'   `"TMT10plex"` is refused rather than mistaken for a kit whose name it contains. `NULL`, the
#'   default, lists every kit.
#' @param timeout Seconds to allow, or `NULL` to wait indefinitely.
#'
#' @return An `mzlibr_isobaric_kits`: `kit`, the name you asked for exactly as given (`NA` when
#'   every kit was listed); `kit_count` kits listed; `record_count` channels listed over every kit;
#'   `absolute_tolerance`, the half-width in Da of each channel's matching window; `kits`, a
#'   data.frame of `kit` and `channel_count`, in mzLib's order; `caveats`; and `records`, a
#'   data.frame with one row per (kit, channel): `kit`, `channel_index` (1-based within the kit,
#'   ascending with the m/z), `channel_label` as MetaMorpheus spells it, `reporter_ion_mz` (m/z,
#'   charge 1, theoretical), and `mz_min` and `mz_max`, the matching window in m/z.
#'
#' @section What the m/z are not:
#'
#' `reporter_ion_mz` is theoretical, at charge 1, from mzLib's `TMT.txt`: not a calibrated or
#' observed value. mzLib reads a reporter intensity as the most intense peak between `mz_min` and
#' `mz_max`. `TMT16` is the lowest sixteen channels of TMTpro 18-plex, and `iTRAQ8` has no 120
#' channel (the phenylalanine immonium ion sits there); its eighth reagent is 121.
#'
#' @seealso [sdrf_design()], for a label-free design from an SDRF.
#' @examples
#' \dontshow{.mzlibr_example <- mzLibR:::replay_bridge_start()}
#' catalogue <- isobaric_kits()
#' catalogue$kits
#' catalogue$absolute_tolerance
#' catalogue$records$channel_label[catalogue$records$kit == "iTRAQ8"]
#'
#' # TMTpro 18-plex, whose top channel sits at 135.15:
#' tmtpro <- isobaric_kits("TMT18")
#' utils::tail(tmtpro$records[, c("channel_index", "channel_label", "reporter_ion_mz")], 3)
#' \dontshow{mzLibR:::replay_bridge_stop(.mzlibr_example)}
#' @spec isobaric.kits
#' @export
isobaric_kits <- function(kit = NULL, timeout = 60) {
  args <- c("isobaric", "kits")
  if (!is.null(kit)) {
    if (!is.character(kit) || length(kit) != 1L || is.na(kit) || !nzchar(trimws(kit))) {
      stop(mzlib_usage_error(paste0(
        "kit must be a kit name such as \"TMT10\" or \"iTRAQ-4plex\", or NULL; got ",
        paste(deparse(kit), collapse = " "), "."
      )))
    }
    args <- c(args, "--kit", kit)
  }
  bridge_require_verb("isobaric kits", ISOBARIC_KITS_SINCE)
  isobaric_parse_kits(bridge_invoke(args, timeout = timeout))
}

#' Print isobaric kits
#'
#' @param x An [isobaric_kits()] result.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mzlibr_isobaric_kits <- function(x, ...) {
  cat("<mzlibr_isobaric_kits> ", format(x$kit_count), " kit(s), ", format(x$record_count),
    " channels, matched within +/- ", format(x$absolute_tolerance), " Da\n",
    sep = ""
  )
  if (nrow(x$kits) > 0L) {
    cat("  ", paste0(x$kits$kit, " (", x$kits$channel_count, ")", collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
