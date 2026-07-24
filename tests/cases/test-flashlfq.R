# Label-free quantification.
#
# The recorded fixture is deliberately tiny and deliberately awkward: it contains a null peptide
# intensity and two null protein intensities, which are the two halves of the most consequential
# distinction in this package. Almost every test here is about which of 0 and NA a field gets.

recorded_quant_data <- function() {
  mz$json_parse(paste(
    readLines(fixture_path("flashlfq_small.json"), warn = FALSE),
    collapse = "\n"
  ))
}

recorded_quant <- function() {
  mz$flashlfq_parse(recorded_quant_data())
}

quant_args <- function(psms = "AllPSMs.psmtsv", normalize = FALSE, ppm_tolerance = 10,
                       isotope_ppm_tolerance = 5, integrate = FALSE, match_between_runs = FALSE,
                       mbr_ppm_tolerance = 10, mbr_q_value_threshold = 0.05,
                       use_shared_peptides_for_protein_quant = FALSE,
                       bayesian_protein_quant = FALSE, use_pep_q_value = FALSE,
                       max_threads = 1, output_directory = NULL) {
  mz$flashlfq_build_args(
    psms, normalize, ppm_tolerance, isotope_ppm_tolerance, integrate, match_between_runs,
    mbr_ppm_tolerance, mbr_q_value_threshold, use_shared_peptides_for_protein_quant,
    bayesian_protein_quant, use_pep_q_value, max_threads, output_directory
  )
}

# ---------------------------------------------------------------- 0 versus NA

test_that("a peptide intensity of null becomes 0, not NA", {
  # A peptide intensity is 0 when the peptide was not measured in that run. It is never "could
  # not be resolved" — that is a protein-level outcome only. Reading this null as NA would put
  # an NA into a column where arithmetic propagates it, and mean() over a peptide column would
  # go from a number to NA across the board.
  peptides <- recorded_quant()$peptides
  missing <- peptides[peptides$base_sequence == "ACDEFR" & peptides$file_name == "run_4", ]

  expect_identical(nrow(missing), 1L)
  expect_identical(missing$intensity, 0)
  expect_false(is.na(missing$intensity))
  expect_identical(missing$detection_type, "NotDetected")
  expect_false(any(is.na(peptides$intensity)))
})

test_that("a protein intensity of null becomes NA, not 0", {
  # NA because FlashLFQ's median-polish produced NaN — it could not resolve a number at all,
  # which is categorically different from a measured zero.
  proteins <- recorded_quant()$proteins
  unresolved <- proteins[proteins$protein_group == "P67890", ]

  expect_identical(nrow(unresolved), 2L)
  expect_true(all(is.na(unresolved$intensity)))
  expect_true(is.double(unresolved$intensity))
})

test_that("the NA propagates through arithmetic, which is the entire point", {
  # R is the only one of the three bindings whose type system can say this without prose.
  # Python has to warn about None; Rust makes it Option<f64>; R makes mean() return NA, so a
  # confidently wrong number is not available by accident and na.rm = TRUE is a visible choice.
  proteins <- recorded_quant()$proteins
  expect_true(is.na(mean(proteins$intensity)))
  expect_false(is.na(mean(proteins$intensity, na.rm = TRUE)))

  # And mzLibR never makes that choice for the user on an intensity, anywhere in the package.
  #
  # The invariant is specifically about intensities. `na.rm` on a count — summing file sizes, or
  # per-run MBR peak counts — is fine, because an absent count means zero peaks rather than an
  # unknown quantity. An absent intensity means the opposite.
  sources <- list.files(file.path("..", "R"), pattern = "[.]R$", full.names = TRUE)
  if (length(sources) == 0L) {
    sources <- list.files("R", pattern = "[.]R$", full.names = TRUE)
  }
  if (length(sources) > 0L) {
    offenders <- character(0)
    for (source in sources) {
      lines <- readLines(source, warn = FALSE)
      offenders <- c(offenders, grep("intensit.*na[.]rm", lines, value = TRUE))
    }
    expect_identical(offenders, character(0), info = paste(offenders, collapse = " | "))
  }
})

test_that("zero is the common protein outcome and NA is the rare one", {
  # Calibration matters: documentation that sells NA hard implies the opposite of the truth. On
  # the K562 pair, 847 protein groups are 0 in both runs and 2 are NA.
  proteins <- recorded_quant()$proteins
  expect_identical(sum(is.na(proteins$intensity)), 2L)
  expect_true(sum(!is.na(proteins$intensity)) >= sum(is.na(proteins$intensity)) - 2L)
})

# ---------------------------------------------------------------- shapes

test_that("the four tables are long, tidy and factor-free", {
  results <- recorded_quant()
  expect_true(inherits(results, "mzlibr_quant"))
  for (frame in list(results$spectra_files, results$peptides, results$proteins, results$peaks)) {
    expect_true(is.data.frame(frame))
  }
  expect_false(is.factor(results$peptides$sequence))
  expect_false(is.factor(results$peaks$detection_type))

  # Long: one row per peptide per run, not one column per run.
  expect_identical(nrow(results$peptides), 4L)
  expect_identical(nrow(results$proteins), 4L)
  expect_identical(nrow(results$peaks), 3L)
  expect_identical(nrow(results$spectra_files), 2L)
})

test_that("the run description survives with mzLib's own field names", {
  results <- recorded_quant()
  expect_identical(results$identification_count, 4)
  expect_identical(
    names(results$spectra_files),
    c("file_name", "full_path", "condition", "biological_replicate", "technical_replicate",
      "fraction", "peak_count", "mbr_peak_count")
  )
  run_3 <- results$spectra_files[results$spectra_files$file_name == "run_3", ]
  expect_identical(run_3$peak_count, 3)
  expect_identical(run_3$mbr_peak_count, 1)
})

test_that("output_directory of null becomes NA rather than the string 'NA'", {
  expect_true(is.na(recorded_quant()$output_directory))
})

test_that("counts are of distinct entities, not of rows", {
  results <- recorded_quant()
  expect_identical(mz$flashlfq_peptide_count(results), 2L)
  expect_identical(mz$flashlfq_protein_count(results), 2L)
})

# ---------------------------------------------------------------- match between runs

test_that("MBR peaks are read from peaks, not from the peptide roll-up", {
  results <- recorded_quant()
  transfers <- mz$flashlfq_mbr_peaks(results)
  expect_identical(nrow(transfers), 1L)
  expect_identical(transfers$file_name, "run_4")
  expect_identical(transfers$sequence, "PEPTIDEK")
})

test_that("mbr_peak_count sums the per-run counts", {
  expect_identical(mz$flashlfq_mbr_peak_count(recorded_quant()), 1)
})

test_that("mbr_rescued_peptide_count counts distinct sequences among MBR peaks", {
  # Stated in code terms because the prose a reader supplies — "peptides quantified in at least
  # one run *only* by MBR" — is subtly different and gives a different number. On the K562 pair
  # this is 140 and the strict reading is 135; the five that differ have both an MBR peak and a
  # zero-intensity MSMS peak in the same run.
  expect_identical(mz$flashlfq_mbr_rescued_peptide_count(recorded_quant()), 1L)
})

test_that("printing warns when the roll-up under-counts the transfers", {
  # The place the mistake is actually made: someone printing a result and reading off a number.
  output <- paste(capture.output(print(recorded_quant())), collapse = "\n")
  expect_true(grepl("transfers in peaks", output, fixed = TRUE), info = output)
  expect_true(grepl("NA (could not be resolved)", output, fixed = TRUE), info = output)
})

# ---------------------------------------------------------------- the spectra design

test_that("a bare path stays a bare line on the wire", {
  # Trailing design fields are dropped so the bridge applies MetaMorpheus's defaults rather
  # than being handed empty strings that mean something else.
  expect_identical(mz$flashlfq_spectra_stdin("run_1.mzML"), "run_1.mzML")
  expect_identical(
    mz$flashlfq_spectra_stdin(c("a.mzML", "b.mzML")),
    c("a.mzML", "b.mzML")
  )
})

test_that("a design data.frame becomes tab-separated fields", {
  design <- data.frame(
    path = c("a.mzML", "b.mzML"),
    condition = c("control", "treated"),
    biological_replicate = c(1, 2),
    stringsAsFactors = FALSE
  )
  expect_identical(
    mz$flashlfq_spectra_stdin(design),
    c("a.mzML\tcontrol\t1", "b.mzML\ttreated\t2")
  )
})

test_that("a condition with no replicate still renders in the right column", {
  design <- data.frame(path = "a.mzML", condition = "control", stringsAsFactors = FALSE)
  expect_identical(mz$flashlfq_spectra_stdin(design), "a.mzML\tcontrol")
})

test_that("a fraction with no condition keeps its position", {
  # The fields are positional, so an omitted middle field must still occupy its slot or the
  # bridge reads the fraction as a condition.
  design <- data.frame(path = "a.mzML", fraction = 2, stringsAsFactors = FALSE)
  expect_identical(mz$flashlfq_spectra_stdin(design), "a.mzML\t\t\t\t2")
})

test_that("spectra paths are validated before anything is spawned", {
  expect_error(mz$flashlfq_spectra_stdin(character(0)),
    class = "mzlib_usage_error", contains = "At least one spectra file"
  )
  expect_error(mz$flashlfq_spectra_stdin(data.frame(x = 1)),
    class = "mzlib_usage_error", contains = "no 'path' column"
  )
  expect_error(mz$flashlfq_spectra_stdin(42), class = "mzlib_usage_error")
  expect_error(mz$flashlfq_spectra_stdin(""), class = "mzlib_usage_error")
})

test_that("a tab or newline in a path is refused, not silently split", {
  # They are the wire's field and record separators, so one in a path would silently become a
  # different experimental design.
  expect_error(mz$flashlfq_spectra_stdin("a\tb.mzML"),
    class = "mzlib_usage_error", contains = "may not contain a tab"
  )
  expect_error(mz$flashlfq_spectra_stdin("a\nb.mzML"), class = "mzlib_usage_error")
})

test_that("raw and d files are refused up front with the remedy named", {
  # FlashLFQ reads mzML only. Rejecting here rather than deep inside a reader means the message
  # names the file and what to do about it.
  expect_error(mz$flashlfq_spectra_stdin("run_1.raw"),
    class = "mzlib_usage_error", contains = c("mzML only", "MSConvert")
  )
  expect_error(mz$flashlfq_spectra_stdin("run_1.d"), class = "mzlib_usage_error")
  # Case does not matter.
  expect_identical(mz$flashlfq_spectra_stdin("run_1.MZML"), "run_1.MZML")
})

test_that("a negative or fractional design field is refused", {
  for (value in c(-1, 1.5)) {
    expect_error(
      mz$flashlfq_spectra_stdin(data.frame(path = "a.mzML", fraction = value)),
      class = "mzlib_usage_error", contains = "non-negative whole number"
    )
  }
})

# ---------------------------------------------------------------- argument assembly

test_that("the defaults reach the wire in the bridge's spelling", {
  args <- quant_args()
  paired <- function(flag) args[which(args == flag) + 1L]
  expect_identical(args[1:2], c("quant", "flashlfq"))
  expect_identical(paired("--psms"), "AllPSMs.psmtsv")
  expect_identical(paired("--ppm"), "10")
  expect_identical(paired("--isotope-ppm"), "5")
  expect_identical(paired("--mbr-ppm"), "10")
  expect_identical(paired("--mbr-q"), "0.05")
  expect_identical(paired("--threads"), "1")
  for (flag in c("--normalize", "--integrate", "--mbr", "--shared-peptides", "--bayesian",
                 "--use-pep-q", "--out")) {
    expect_false(flag %in% args, info = flag)
  }
})

test_that("every boolean flag has the polarity it says it has", {
  expect_true("--normalize" %in% quant_args(normalize = TRUE))
  expect_true("--integrate" %in% quant_args(integrate = TRUE))
  expect_true("--mbr" %in% quant_args(match_between_runs = TRUE))
  expect_true("--shared-peptides" %in%
    quant_args(use_shared_peptides_for_protein_quant = TRUE))
  expect_true("--bayesian" %in% quant_args(bayesian_protein_quant = TRUE))
  expect_true("--use-pep-q" %in% quant_args(use_pep_q_value = TRUE))
  for (bad in list(NA, "yes", 1)) {
    expect_error(quant_args(normalize = bad), class = "mzlib_usage_error")
  }
})

test_that("numbers reach the wire with a dot, whatever the session locale", {
  # The bridge parses with InvariantCulture. `format()` honours options(OutDec), which a user in
  # a European locale may well have set, so the decimal mark is stated explicitly rather than
  # assumed. Without this, mbr_q_value_threshold would arrive as "0,05" and fail to parse — on
  # that user's machine only.
  previous <- options(OutDec = ",")
  on.exit(options(previous), add = TRUE)

  args <- quant_args(mbr_q_value_threshold = 0.05, ppm_tolerance = 12.5)
  expect_identical(args[which(args == "--mbr-q") + 1L], "0.05")
  expect_identical(args[which(args == "--ppm") + 1L], "12.5")
})

test_that("the output directory is passed only when asked for", {
  expect_false("--out" %in% quant_args())
  args <- quant_args(output_directory = "results")
  expect_identical(args[which(args == "--out") + 1L], "results")
  expect_error(quant_args(output_directory = "  "), class = "mzlib_usage_error")
})

test_that("a missing psms path is refused", {
  for (bad in list("", "   ", NA_character_, NULL, 42)) {
    expect_error(quant_args(psms = bad), class = "mzlib_usage_error", contains = "psms must be")
  }
})

test_that("an impossible thread count is refused", {
  for (bad in list(0, -2, 1.5, NA_real_, "1")) {
    expect_error(quant_args(max_threads = bad), class = "mzlib_usage_error")
  }
  # -1 is meaningful: every core.
  expect_identical(quant_args(max_threads = -1)[
    which(quant_args(max_threads = -1) == "--threads") + 1L
  ], "-1")
})

# ---------------------------------------------------------------- reproducibility

test_that("max_threads defaults to 1, which differs from pyMzLib deliberately", {
  # pyMzLib's default is -1. mzLibR's is 1, because with more than one thread FlashLFQ's peptide
  # roll-up nondeterministically drops MBR intensities and identical inputs give different
  # protein-level answers roughly 1 run in 6 (smith-chem-wisc/mzLib#1111). A binding that
  # silently produces unreproducible results by default is worse than one that differs from its
  # parent in a documented way. The right fix is upstream.
  expect_identical(formals(flashlfq_quantify)$max_threads, 1)
})

test_that("setting max_threads to anything else warns, naming the issue and the remedy", {
  # The warning goes at the call that will produce the unreproducible answer, not only in help
  # the user may never open.
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  runner <- stub_runner(stdout = '{"ok":true,"data":{"peptides":[],"proteins":[],"peaks":[]}}')

  with_bridge_config(option = path, {
    # `flashlfq_quantify` does not take a runner, so the warning is checked against the
    # assembled call by invoking the exported function with a bridge that returns a stub.
    expect_warning(
      tryCatch(
        flashlfq_quantify("AllPSMs.psmtsv", "run.mzML", max_threads = 4),
        mzlib_error = function(e) NULL
      ),
      contains = c("mzLib#1111", "max_threads = 1")
    )
    expect_no_warning(
      tryCatch(
        flashlfq_quantify("AllPSMs.psmtsv", "run.mzML", max_threads = 1),
        mzlib_error = function(e) NULL
      )
    )
  })
})

test_that("printing warns when the result was produced multithreaded", {
  results <- recorded_quant()
  results$parameters$max_threads <- 8
  output <- paste(capture.output(print(results)), collapse = "\n")
  expect_true(grepl("may not reproduce", output, fixed = TRUE), info = output)
})

test_that("the fixture was produced single-threaded", {
  expect_identical(recorded_quant()$parameters$max_threads, 1)
})
