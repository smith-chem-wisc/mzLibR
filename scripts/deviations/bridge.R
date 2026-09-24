# The transport: locating, installing and asking the bridge. See scripts/spec-facts.R for what
# each list means.

# Transport functions carry the package prefix, like mzlibr_bridge_path() and
# mzlibr_install_bridge(), so a user can find all three together.
R_NAMES["version"] <- "mzlibr_bridge_version"

R_DEVIATIONS[["version"]] <- list(
  "field.verbs" = pending()
)

PARENT_ADDITIONS["mzlibr_bridge_path"] <- "transport; pyMzLib's equivalent is private (_bridge)"
PARENT_ADDITIONS["mzlibr_bridge_version"] <-
  "transport; pyMzLib's equivalent is pymzlib.bridge_version, outside the modules compared here"
# Says why THIS package needs the function, not what the other bindings do instead. The earlier
# wording ("Rust downloads it from build.rs") was a claim about another repository, which nothing
# here can test and which goes false when that repository changes with nothing changing here -
# the mistake mzLibRust #16 was opened to stop making.
PARENT_ADDITIONS["mzlibr_install_bridge"] <- paste(
  "no pyMzLib counterpart: a wheel carries the payload, so Python never has to fetch one.",
  "CRAN forbids downloading at install time and writing outside tempdir() without consent,",
  "so here the download has to be a function the user calls."
)

FIELD_CHECKS[["version"]] <- list("bridge_version.json", function(d) {
  path <- file.path(fixtures, "bridge_version.json")
  text <- rawToChar(readBin(path, "raw", file.info(path)$size))
  old <- options(mzlibr.bridge = path)
  on.exit(options(old))
  mzlibr_bridge_version(runner = function(exe, args, stdin = NULL, timeout = NULL) {
    list(stdout = paste0("{\"ok\":true,\"data\":", text, "}"), stderr = "", status = 0L, timed_out = FALSE)
  })
})
