# The transport: finding the bridge, running it, and reading the envelope.
#
# Almost all of it runs with no bridge and no subprocess, because `bridge_decode()` is pure and
# the runner is a function argument. The handful of tests that do need a real executable skip
# when there is none, so the offline suite passes on a machine that has never seen .NET.

# Captured before anything below clears it, so the live tests can still find a staged bridge.
live_bridge <- Sys.getenv("MZLIB_BRIDGE", "")

# Run `code` with mzLibR's two override mechanisms under our control.
#
# Both are cleared first. A developer with MZLIB_BRIDGE set in their shell would otherwise get
# different results from CI, which is the exact mistake mzLibRust made: a test that passed on
# the build machine and failed for anyone who had staged a bridge.
with_bridge_config <- function(code, option = NULL, env = "") {
  old_option <- options(mzlibr.bridge = option)
  old_env <- Sys.getenv("MZLIB_BRIDGE", NA_character_)
  if (nzchar(env)) Sys.setenv(MZLIB_BRIDGE = env) else Sys.unsetenv("MZLIB_BRIDGE")
  on.exit(
    {
      options(old_option)
      if (is.na(old_env)) Sys.unsetenv("MZLIB_BRIDGE") else Sys.setenv(MZLIB_BRIDGE = old_env)
    },
    add = TRUE
  )
  force(code)
}

# A file that exists and is not an executable. Enough for path resolution, which only asks
# whether something is there.
fake_bridge_file <- function() {
  path <- tempfile("mzlib-bridge-", fileext = if (.Platform$OS.type == "windows") ".exe" else "")
  writeLines("not a real executable", path)
  path
}

# A runner that returns a canned result and records what it was asked to do.
stub_runner <- function(stdout = "", stderr = "", status = 0L, timed_out = FALSE) {
  seen <- new.env(parent = emptyenv())
  seen$called <- FALSE
  list(
    seen = seen,
    run = function(exe, args, stdin = NULL, timeout = NULL) {
      seen$called <- TRUE
      seen$exe <- exe
      seen$args <- args
      seen$stdin <- stdin
      seen$timeout <- timeout
      list(stdout = stdout, stderr = stderr, status = status, timed_out = timed_out)
    }
  )
}

# `bridge_decode()` takes what a runner returned, so a canned result is the whole input.
decoded <- function(stdout, stderr = "", status = 0L) {
  mz$bridge_decode(list(stdout = stdout, stderr = stderr, status = status, timed_out = FALSE))
}

# ---------------------------------------------------------------- platform mapping

test_that("platform tags match the .NET runtime identifiers on every supported machine", {
  # These strings must equal the RIDs `publish-bridge.ps1` stages under, or nothing is ever
  # found. The table is the contract between three bindings and one build script, so it is
  # checked here rather than discovered on a Mac.
  expectations <- list(
    c("Windows", "x86-64", "win-x64"),
    c("Windows", "x86_64", "win-x64"),
    c("Windows", "AMD64", "win-x64"),
    c("Windows", "arm64", "win-arm64"),
    c("Linux", "x86_64", "linux-x64"),
    c("Linux", "aarch64", "linux-arm64"),
    c("Linux", "arm64", "linux-arm64"),
    c("Darwin", "x86_64", "osx-x64"),
    c("Darwin", "arm64", "osx-arm64"),
    # Lower-case spellings, because `Sys.info()` is not consistent across platforms and
    # neither are the other two bindings.
    c("darwin", "aarch64", "osx-arm64"),
    c("linux", "X86_64", "linux-x64")
  )
  for (case in expectations) {
    expect_identical(
      mz$bridge_platform_tag(case[1], case[2]), case[3],
      info = paste(case[1], case[2])
    )
  }
})

test_that("an unsupported platform is reported, not guessed", {
  expect_error(
    mz$bridge_platform_tag("Plan9", "x86_64"),
    class = "mzlib_bridge_not_found", contains = "Unsupported platform"
  )
})

test_that("the executable name follows the platform", {
  expect_true(mz$bridge_executable_name() %in% c("mzlib-bridge", "mzlib-bridge.exe"))
})

test_that("the cache directory is laid out the way R_user_dir would have", {
  # `tools::R_user_dir()` is R >= 4.0 and this package supports 3.5, so the path is built by
  # hand — but it must land in the same place, or a user who upgrades R downloads 140 MB again.
  path <- mz$mzlibr_cache_dir()
  expect_true(is.character(path) && length(path) == 1L)
  expect_true(grepl("mzLibR", path, fixed = TRUE))
})

# ---------------------------------------------------------------- locating the bridge

test_that("the option wins over the environment variable", {
  from_option <- fake_bridge_file()
  from_env <- fake_bridge_file()
  on.exit(unlink(c(from_option, from_env)), add = TRUE)

  with_bridge_config(option = from_option, env = from_env, {
    expect_identical(mzlibr_bridge_path(), from_option)
  })
})

test_that("the environment variable is used when no option is set", {
  from_env <- fake_bridge_file()
  on.exit(unlink(from_env), add = TRUE)

  with_bridge_config(env = from_env, {
    expect_identical(mzlibr_bridge_path(), from_env)
  })
})

test_that("an override pointing at nothing says so rather than falling through", {
  # Falling through to the next candidate would be worse than failing: the user asked for a
  # specific bridge, and silently running a different one is how two halves built from
  # different sources end up talking to each other.
  with_bridge_config(option = "/definitely/not/here/mzlib-bridge", {
    expect_error(
      mzlibr_bridge_path(),
      class = "mzlib_bridge_not_found", contains = "not a file"
    )
  })
  with_bridge_config(env = "/definitely/not/here/mzlib-bridge", {
    expect_error(
      mzlibr_bridge_path(),
      class = "mzlib_bridge_not_found", contains = "not a file"
    )
  })
})

test_that("the missing-bridge message names every way out", {
  # Asserted against the message builder rather than by calling `mzlibr_bridge_path()` with the
  # environment stripped, so the test does not depend on whether this machine happens to have a
  # bridge cached. A test whose result depends on the developer's machine is worse than none.
  message <- mz$bridge_missing_message("/some/cache/win-x64/mzlib-bridge.exe")
  expect_true(grepl("mzlibr_install_bridge()", message, fixed = TRUE))
  expect_true(grepl("MZLIB_BRIDGE", message, fixed = TRUE))
  expect_true(grepl("mzlibr.bridge", message, fixed = TRUE))
  expect_true(grepl("/some/cache/win-x64/mzlib-bridge.exe", message, fixed = TRUE))
  # LGPL section 4: the override is how a user relinks a modified mzLib without rebuilding
  # this package, so the licence affordance is stated where it is exercised.
  expect_true(grepl("LGPL", message, fixed = TRUE))
})

test_that("bridge_available() answers without raising", {
  with_bridge_config(option = "/definitely/not/here/mzlib-bridge", {
    expect_false(mz$bridge_available())
  })
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  with_bridge_config(option = path, {
    expect_true(mz$bridge_available())
  })
})

# ---------------------------------------------------------------- reading the envelope

test_that("success returns only the data", {
  data <- decoded('{"ok":true,"data":{"a":1},"error":null}')
  expect_identical(data$a, 1)
})

test_that("a usage failure becomes a usage error", {
  expect_error(
    decoded('{"ok":false,"data":null,"error":{"type":"usage","message":"Missing --accession."}}'),
    class = "mzlib_usage_error", contains = "Missing --accession"
  )
})

test_that("an availability failure becomes its own class", {
  # This is the class the live suites skip on. Nothing else may be widened into it.
  expect_error(
    decoded('{"ok":false,"error":{"type":"ServiceUnavailable","message":"EBI is down"}}'),
    class = "mzlib_service_unavailable", contains = "EBI is down"
  )
})

test_that("any other failure keeps the .NET exception type", {
  # `error_type` is what lets a caller tell a network blip from a bad accession without
  # parsing prose.
  condition <- expect_error(
    decoded('{"ok":false,"error":{"type":"HttpRequestException","message":"503"}}'),
    class = "mzlib_bridge_error"
  )
  expect_identical(condition$error_type, "HttpRequestException")
})

test_that("a failure with no message still reads as something", {
  expect_error(decoded('{"ok":false}'), class = "mzlib_error", contains = "no message")
})

test_that("a silent death surfaces the exit code and stderr", {
  # A process that dies before writing anything leaves stderr as the only evidence there is.
  expect_error(
    decoded("", stderr = "Segmentation fault", status = 139L),
    class = "mzlib_protocol_error", contains = c("139", "Segmentation fault")
  )
})

test_that("a silent death with no stderr does not produce an empty message", {
  expect_error(
    decoded("", stderr = "", status = 1L),
    class = "mzlib_protocol_error", contains = "(empty)"
  )
})

test_that("output that is not JSON is quoted back, not swallowed", {
  expect_error(
    decoded("Unhandled exception. System.Whatever"),
    class = "mzlib_protocol_error", contains = "System.Whatever"
  )
})

test_that("JSON that is not an envelope is refused", {
  expect_error(decoded("[1,2,3]"), class = "mzlib_protocol_error", contains = "not an envelope")
})

# ---------------------------------------------------------------- timeouts

test_that("a timeout becomes a typed condition", {
  expect_error(
    mz$bridge_decode(
      list(stdout = "", stderr = "", status = 124L, timed_out = TRUE),
      timeout = 5
    ),
    class = "mzlib_timeout", contains = "5s"
  )
})

test_that("a timeout is NOT reported as a service outage", {
  # The regression this guards was live in pyMzLib: every subprocess timeout raised
  # ServiceUnavailable, which the live suites turn into a skip. A wedged bridge, a corrupt
  # binary, or a caller passing too small a timeout all reported "the repository is down" and
  # the suite went green — precisely the failure the skip convention exists to prevent.
  condition <- tryCatch(
    mz$bridge_decode(list(stdout = "", stderr = "", status = 124L, timed_out = TRUE), timeout = 1),
    error = function(e) e
  )
  expect_true(inherits(condition, "mzlib_timeout"))
  expect_false(inherits(condition, "mzlib_service_unavailable"))
})

test_that("a timeout that cannot mean anything is rejected before anything is spawned", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)

  for (bad in list(0, -1, NA_real_, Inf, "30", c(1, 2))) {
    runner <- stub_runner()
    with_bridge_config(option = path, {
      expect_error(
        mz$bridge_invoke("version", timeout = bad, runner = runner$run),
        class = "mzlib_usage_error"
      )
    })
    expect_false(runner$seen$called, info = "must not have run the bridge")
  }
})

test_that("NULL means wait indefinitely and is not an error", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  runner <- stub_runner(stdout = '{"ok":true,"data":{}}')
  with_bridge_config(option = path, {
    mz$bridge_invoke("version", timeout = NULL, runner = runner$run)
  })
  expect_true(runner$seen$called)
  expect_true(is.null(runner$seen$timeout))
})

# ---------------------------------------------------------------- what reaches the process

test_that("R itself makes an embedded-null argument unconstructible", {
  # pyMzLib and mzLibRust both check for a null inside an argument, because a truncated
  # accession reaching the bridge is worse than a rejected one. mzLibR has no such check, and
  # this test is why: R cannot represent the input at all. A string literal containing one is a
  # parse error, and every route that builds a string from bytes refuses it too. Porting the
  # guard would have added an unreachable branch — recorded in docs/test-parity.md as
  # eliminated by R's semantics rather than silently dropped.
  expect_error(rawToChar(as.raw(c(0x50L, 0x00L, 0x44L))), contains = "embedded nul")
})

test_that("a bridge that will not launch still fails as one of ours", {
  # A quarantined binary, a missing execute bit, a file that is not executable at all. Whatever
  # the platform makes of it, it must not escape as a bare system error - the promise is that
  # every failure from this package is an mzlib_error.
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)

  with_bridge_config(option = path, {
    condition <- tryCatch(mzlibr_bridge_version(), error = function(e) e)
    expect_true(inherits(condition, "mzlib_error"),
      info = paste(class(condition), collapse = ", ")
    )
  })
})

test_that("every condition mzLibR raises is catchable as one type", {
  # So a caller can write a single handler and be done, without enumerating the hierarchy.
  conditions <- list(
    mz$mzlib_usage_error("x"),
    mz$mzlib_service_unavailable("x"),
    mz$mzlib_bridge_error("x", "SomeException"),
    mz$mzlib_timeout(5),
    mz$mzlib_bridge_not_found("x"),
    mz$mzlib_protocol_error("x"),
    mz$mzlib_project_not_found("x")
  )
  for (condition in conditions) {
    expect_true(inherits(condition, "mzlib_error"),
      info = paste(class(condition), collapse = ", ")
    )
    expect_true(inherits(condition, "error"))
    expect_true(nzchar(conditionMessage(condition)))
  }
})

test_that("stdin reaches the runner unchanged", {
  # Large or variadic input goes on stdin and never on argv, because argv has a ceiling of
  # roughly 32 KB and a real experiment's worth of file names goes straight past it.
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  payload <- c("run_3.mzML\tA\t1\t1\t1", "run_4.mzML\tA\t2\t1\t1")
  runner <- stub_runner(stdout = '{"ok":true,"data":{}}')

  with_bridge_config(option = path, {
    mz$bridge_invoke(c("quant", "flashlfq"), stdin = payload, runner = runner$run)
  })
  expect_identical(runner$seen$stdin, payload)
})

test_that("arguments are quoted for the shell that will parse them", {
  # `system2()` quotes the command and not the arguments, so a path with a space would arrive
  # at the bridge in pieces. On Windows that is the normal case, not an edge case.
  quoted <- mz$bridge_quote("C:/Users/someone/My Documents/AllPSMs.psmtsv")
  expect_true(grepl(" ", quoted, fixed = TRUE))
  expect_true(nchar(quoted) > nchar("C:/Users/someone/My Documents/AllPSMs.psmtsv"))
})

# ---------------------------------------------------------------- version handshake

test_that("a matching protocol returns the version information", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  runner <- stub_runner(
    stdout = '{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27"},"error":null}'
  )
  with_bridge_config(option = path, {
    info <- mzlibr_bridge_version(runner = runner$run)
    expect_identical(info$bridge, "1.0.0.0")
    expect_identical(info$protocol, 1L)
    expect_identical(info$runtime, "8.0.27")
    # This payload carries no `mzlib`, which is what a bridge built before that field existed
    # sends. It must still read, and the missing value must project to NA rather than to "" or a
    # placeholder: an added wire field may not strand a caller on an older bridge.
    expect_identical(info$mzlib, NA_character_)
  })
})

test_that("the mzLib build is reported when the bridge sends it", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  runner <- stub_runner(
    stdout = paste0(
      '{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27",',
      '"mzlib":"1.0.0+f6b0f0d17f32383918ef895006aaecb71cdb9a7e"},"error":null}'
    )
  )
  with_bridge_config(option = path, {
    info <- mzlibr_bridge_version(runner = runner$run)
    expect_identical(info$mzlib, "1.0.0+f6b0f0d17f32383918ef895006aaecb71cdb9a7e")
    # It answers "which mzLib is this?", not "may I talk to it?" - protocol is that, and stays 1
    # precisely because adding a field is backward compatible.
    expect_identical(info$protocol, 1L)
  })
})

test_that("a mismatched protocol fails loudly", {
  # Halves built from different sources must not quietly produce wrong results.
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  runner <- stub_runner(
    stdout = '{"ok":true,"data":{"bridge":"9.9","protocol":2,"runtime":"8.0.27"}}'
  )
  with_bridge_config(option = path, {
    expect_error(
      mzlibr_bridge_version(runner = runner$run),
      class = "mzlib_protocol_error", contains = "different sources"
    )
  })
})

test_that("the version verb is what gets sent, with no stdin", {
  path <- fake_bridge_file()
  on.exit(unlink(path), add = TRUE)
  runner <- stub_runner(
    stdout = '{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27"}}'
  )
  with_bridge_config(option = path, {
    mzlibr_bridge_version(runner = runner$run)
  })
  expect_identical(runner$seen$args, "version")
  expect_true(is.null(runner$seen$stdin))
})

# ---------------------------------------------------------------- against a real bridge

test_that("LIVE: the bridge round-trips a version handshake", {
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")
  with_bridge_config(option = live_bridge, {
    info <- mzlibr_bridge_version()
    expect_identical(info$protocol, 1L)
    expect_true(nchar(info$bridge) > 0L)
    expect_true(nchar(info$runtime) > 0L)
  })
})

test_that("LIVE: an argument containing a space survives to the bridge intact", {
  # End-to-end proof that `bridge_quote()` is doing its job: the bridge echoes the path it was
  # given back in the error message, so if quoting were wrong the path would come back split.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")

  # The bridge validates the spectra files from stdin *before* it looks at --psms, so reaching
  # the argument under test means giving it a spectra file that exists. That the directory has
  # a space in it exercises both paths at once: stdin is tab-separated so spaces are safe
  # there, while argv has to survive the shell.
  spacey_dir <- file.path(tempdir(), "a directory with spaces")
  dir.create(spacey_dir, showWarnings = FALSE, recursive = TRUE)
  spectra <- file.path(spacey_dir, "run.mzML")
  file.create(spectra)
  on.exit(unlink(spacey_dir, recursive = TRUE), add = TRUE)
  psms <- file.path(spacey_dir, "AllPSMs.psmtsv")

  with_bridge_config(option = live_bridge, {
    condition <- tryCatch(
      mz$bridge_invoke(c("quant", "flashlfq", "--psms", psms), stdin = spectra, timeout = 120),
      error = function(e) e
    )
    expect_true(inherits(condition, "mzlib_error"))
    expect_true(
      grepl("a directory with spaces", conditionMessage(condition), fixed = TRUE),
      info = conditionMessage(condition)
    )
    # Split by the shell, the path would arrive as "a" and the message would name that
    # instead — so the full file name is what proves the quoting held.
    expect_true(
      grepl("AllPSMs.psmtsv", conditionMessage(condition), fixed = TRUE),
      info = conditionMessage(condition)
    )
  })
})

test_that("LIVE: stdin is delivered to the bridge", {
  # The bridge reports a different failure when stdin arrives empty than when it arrives with
  # content, which makes stdin delivery observable without a real experiment. This is the test
  # that replaced the `processx` dependency: base R's `system2(input =, timeout =)` really does
  # deliver both at once, on Windows included.
  skip_if(!nzchar(live_bridge), "no bridge staged (set MZLIB_BRIDGE)")

  with_bridge_config(option = live_bridge, {
    with_content <- tryCatch(
      mz$bridge_invoke(
        c("quant", "flashlfq", "--psms", file.path(tempdir(), "nope.psmtsv")),
        stdin = "no-such-run.mzML", timeout = 120
      ),
      error = function(e) conditionMessage(e)
    )
    without <- tryCatch(
      mz$bridge_invoke(
        c("quant", "flashlfq", "--psms", file.path(tempdir(), "nope.psmtsv")),
        stdin = character(0), timeout = 120
      ),
      error = function(e) conditionMessage(e)
    )

    expect_true(grepl("no-such-run.mzML", with_content, fixed = TRUE), info = with_content)
    expect_true(grepl("provided on stdin", without, fixed = TRUE), info = without)
  })
})
