# The JSON reader.
#
# These are the tests that justify not depending on `jsonlite`. Most of them assert a *shape*
# rather than a value, because the shape is what the alternative gets wrong.

# ---------------------------------------------------------------- never simplify

test_that("an object is a named list", {
  parsed <- mz$json_parse('{"accession":"PXD000001","count":8}')
  expect_true(is.list(parsed))
  expect_false(is.data.frame(parsed))
  expect_identical(names(parsed), c("accession", "count"))
  expect_identical(parsed$accession, "PXD000001")
  expect_identical(parsed$count, 8)
})

test_that("one array element has the same shape as two", {
  # The `jsonlite` hazard in one test. `fromJSON()` gives a data.frame for an array of objects
  # and a *list* for an array of one, so a manifest with a single file comes back a different
  # shape from a manifest with two, and every caller has to defend against both. Here the shape
  # is decided by the reader, not by the data.
  one <- mz$json_parse('{"files":[{"name":"a.mgf"}]}')
  two <- mz$json_parse('{"files":[{"name":"a.mgf"},{"name":"b.mzML"}]}')

  expect_identical(class(one$files), class(two$files))
  expect_true(is.list(one$files))
  expect_false(is.data.frame(one$files))
  expect_true(is.null(names(one$files)), "an array is unnamed")
  expect_identical(length(one$files), 1L)
  expect_identical(one$files[[1]]$name, "a.mgf")
  expect_identical(two$files[[2]]$name, "b.mzML")
})

test_that("an array of scalars stays a list, not a vector", {
  parsed <- mz$json_parse("[1,2,3]")
  expect_true(is.list(parsed))
  expect_identical(length(parsed), 3L)
  expect_identical(parsed[[2]], 2)
})

test_that("empty containers keep their kind", {
  empty_object <- mz$json_parse("{}")
  expect_true(is.list(empty_object))
  expect_identical(names(empty_object), character(0))

  empty_array <- mz$json_parse("[]")
  expect_true(is.list(empty_array))
  expect_true(is.null(names(empty_array)))
  expect_identical(length(empty_array), 0L)
})

# ---------------------------------------------------------------- null

test_that("null keeps its place in an array instead of vanishing", {
  # The one that matters most. Read as R `NULL`, the middle element would be *deleted* and the
  # list would come back length 2 — a protein-intensity vector silently one short, with no
  # error anywhere. FlashLFQ emits `null` for exactly the protein it could not resolve, so the
  # value that means "no answer" is the one that would disappear.
  parsed <- mz$json_parse('{"intensities":[1.5,null,3.5]}')

  expect_identical(length(parsed$intensities), 3L)
  expect_true(is.na(parsed$intensities[[2]]))

  flattened <- as.numeric(unlist(parsed$intensities))
  expect_identical(length(flattened), 3L)
  expect_identical(flattened[1], 1.5)
  expect_true(is.na(flattened[2]))
  expect_identical(flattened[3], 3.5)
})

test_that("a null member keeps its name", {
  parsed <- mz$json_parse('{"checksum":null,"size":10}')
  expect_identical(names(parsed), c("checksum", "size"))
  expect_true(is.na(parsed$checksum))
})

test_that("null becomes NA of the right type when a field asks for it", {
  parsed <- mz$json_parse('{"intensity":null}')
  expect_identical(as.numeric(parsed$intensity), NA_real_)
  expect_identical(as.character(parsed$checksum), character(0))
})

test_that("the real version envelope carries a null and reads back whole", {
  # Verbatim from the bridge on this machine. `"error":null` is in every success response, so
  # the null path is not an edge case, it is every single call.
  parsed <- mz$json_parse(
    '{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27"},"error":null}'
  )
  expect_identical(names(parsed), c("ok", "data", "error"))
  expect_identical(parsed$ok, TRUE)
  expect_identical(parsed$data$protocol, 1)
  expect_true(is.na(parsed$error))
})

# ---------------------------------------------------------------- numbers

test_that("every number is a double, including whole ones", {
  # `file_size_bytes` runs past 2^31 on real data — PXD000001 holds a 472 MB mzXML and totals
  # 1.44 GB — so a field read as R's 32-bit integer would overflow to NA. Doubles are exact to
  # 2^53, which nothing on this wire approaches.
  parsed <- mz$json_parse('{"small":8,"bytes":1544000000,"huge":9007199254740991}')
  expect_true(is.double(parsed$small))
  expect_false(is.integer(parsed$small))
  expect_identical(parsed$bytes, 1544000000)
  expect_identical(parsed$huge, 9007199254740991)
})

test_that("negatives, decimals and exponents all read", {
  parsed <- mz$json_parse('[-1,0.5,1e3,1.5E-3,-2.25e+2]')
  expect_identical(parsed[[1]], -1)
  expect_identical(parsed[[2]], 0.5)
  expect_identical(parsed[[3]], 1000)
  expect_identical(parsed[[4]], 0.0015)
  expect_identical(parsed[[5]], -225)
})

test_that("a monoisotopic mass survives to full precision", {
  # The proton mass matters to 1.1 ppm at m/z 500, so a reader that rounded would be a defect
  # nobody would notice until a figure was wrong.
  parsed <- mz$json_parse('{"proton":1.00727646677}')
  expect_identical(parsed$proton, 1.00727646677)
})

# ---------------------------------------------------------------- strings

test_that("plain strings pass through untouched", {
  parsed <- mz$json_parse('["PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz"]')
  expect_identical(parsed[[1]], "PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz")
})

test_that("the escapes JSON defines all decode", {
  parsed <- mz$json_parse('["a\\"b","c\\\\d","e\\/f","g\\nh","i\\tj","k\\rl","m\\bn","o\\fp"]')
  expect_identical(parsed[[1]], "a\"b")
  expect_identical(parsed[[2]], "c\\d")
  expect_identical(parsed[[3]], "e/f")
  expect_identical(parsed[[4]], "g\nh")
  expect_identical(parsed[[5]], "i\tj")
  expect_identical(parsed[[6]], "k\rl")
  expect_identical(parsed[[7]], "m\bn")
  expect_identical(parsed[[8]], "o\fp")
})

test_that("the bridge's apostrophe escape decodes", {
  # Not hypothetical: the bridge writes `'` for every file name it quotes in an error
  # message, so this is the hot path for anything that went wrong.
  parsed <- mz$json_parse('{"message":"Spectra file not found: \\u0027run3.mzML\\u0027."}')
  expect_identical(parsed$message, "Spectra file not found: 'run3.mzML'.")
})

test_that("a character outside the basic plane survives its surrogate pair", {
  parsed <- mz$json_parse('["\\ud83d\\ude00"]')
  expect_identical(parsed[[1]], intToUtf8(0x1F600L))
})

test_that("non-ASCII text is not mangled", {
  parsed <- mz$json_parse('{"organism":"Sch\\u00e9ma \\u00b5g"}')
  expect_identical(parsed$organism, "Schéma µg")
})

# ---------------------------------------------------------------- nesting

test_that("nested structures keep their shape all the way down", {
  parsed <- mz$json_parse(
    '{"a":{"b":[{"c":[1,null]},{"c":[]}]},"d":true}'
  )
  expect_identical(parsed$d, TRUE)
  expect_identical(length(parsed$a$b), 2L)
  expect_identical(length(parsed$a$b[[1]]$c), 2L)
  expect_true(is.na(parsed$a$b[[1]]$c[[2]]))
  expect_identical(length(parsed$a$b[[2]]$c), 0L)
})

test_that("more members than the initial buffer holds still all arrive", {
  # The reader grows its buffers by doubling; this crosses the boundary twice.
  wide <- paste0(
    "{", paste0('"k', 1:40, '":', 1:40, collapse = ","), "}"
  )
  parsed <- mz$json_parse(wide)
  expect_identical(length(parsed), 40L)
  expect_identical(parsed$k40, 40)
  expect_identical(names(parsed)[9], "k9")

  long <- paste0("[", paste0(1:40, collapse = ","), "]")
  expect_identical(length(mz$json_parse(long)), 40L)
})

# ---------------------------------------------------------------- refusing bad input

test_that("output that is not JSON is refused, not guessed at", {
  expect_error(
    mz$json_parse("Unhandled exception. System.Whatever"),
    class = "mzlib_protocol_error"
  )
})

test_that("a stray character outside a value is caught", {
  # `gregexpr()` silently skips what it cannot match, which would turn this into a plausible
  # looking structural error somewhere else entirely.
  expect_error(mz$json_parse('{"a" @ 1}'), class = "mzlib_protocol_error", contains = "@")
})

test_that("trailing content is reported", {
  expect_error(
    mz$json_parse('{"a":1} {"b":2}'),
    class = "mzlib_protocol_error", contains = "trailing"
  )
})

test_that("a truncated document says so", {
  expect_error(
    mz$json_parse('{"a":[1,2'),
    class = "mzlib_protocol_error", contains = "unexpected end of input"
  )
})

test_that("empty output is a protocol error, not an empty result", {
  expect_error(mz$json_parse(""), class = "mzlib_protocol_error")
  expect_error(mz$json_parse("   "), class = "mzlib_protocol_error")
})

test_that("an unknown escape is named rather than swallowed", {
  expect_error(
    mz$json_parse('"a\\qb"'),
    class = "mzlib_protocol_error", contains = "unknown escape"
  )
})

test_that("a malformed object is reported in terms of what was expected", {
  expect_error(mz$json_parse('{1:2}'), class = "mzlib_protocol_error", contains = "member name")
  expect_error(mz$json_parse('{"a" 1}'), class = "mzlib_protocol_error", contains = "':'")
  expect_error(mz$json_parse('{"a":1 "b":2}'), class = "mzlib_protocol_error", contains = "','")
})
