# A JSON reader, in base R, written for this wire contract and no other.
#
# mzLibR takes no package dependencies, so there is no `jsonlite` here. That began as a
# portability decision and turned into a correctness one, because the way `jsonlite` reads
# JSON is the single largest hazard an R binding of this bridge faces:
#
#   * `fromJSON()` auto-simplifies. An array of objects becomes a data.frame, an array of
#     scalars becomes a vector, and *which one you get depends on the data* — a manifest with
#     one file comes back a different shape from a manifest with two. Every caller then has to
#     defend against both.
#   * JSON `null` becomes R `NULL`, and a `NULL` assigned into a list **deletes the element**.
#     A protein-intensity vector does not come back wrong, it comes back *short*, and the
#     shortness is silent. That is the worst failure mode in this whole package: FlashLFQ emits
#     `null` for a protein intensity it could not resolve (see `flashlfq.R`), so the one value
#     that means "no answer" is exactly the one that would vanish.
#
# This reader does neither thing. It never simplifies — an object is always a named list, an
# array is always an unnamed list, whatever is inside them — and `null` becomes `NA`, which
# occupies its slot like any other value. Building the actual data.frames is then a deliberate,
# per-field act in the module that knows what the field means, which is the only place that
# knowledge exists.
#
# Everything here is internal. `json_parse()` is the only entry point.

# ---------------------------------------------------------------- why `null` is NA

# JSON `null` is read as logical `NA`.
#
# `NA` was chosen over `NULL` because `NULL` deletes, and over a bespoke sentinel because the
# rest of R already knows what `NA` means: it survives `[[`, it survives `length()`, it prints
# as `NA`, and `is.na()` finds it. It is unambiguous here because nothing else in JSON can
# produce a logical `NA` — `true` and `false` are the only other logicals and both are known.
#
# Callers coerce it to the `NA` of the right type at the boundary (`as.numeric(NA)` is
# `NA_real_`), which is the explicit, per-field decision this whole file exists to force.
JSON_NULL <- NA

# ---------------------------------------------------------------- tokenising

# One regex for every JSON token, tried in this order.
#
# Strings come first so that a `{` or a `,` inside one is never mistaken for structure. The
# string pattern is deliberately permissive about *what* follows a backslash — `\\.` rather
# than the strict `["\\/bfnrtu]` — because a malformed escape should be reported by the
# unescaper, which can say which escape and where, rather than silently failing to match here
# and leaving the tokeniser to resynchronise in the middle of a string.
JSON_TOKEN_PATTERN <- paste0(
  '"(?:[^"\\\\]|\\\\.)*"',
  "|-?(?:0|[1-9][0-9]*)(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?",
  "|true|false|null",
  "|[][{}:,]"
)

# Split JSON text into tokens, refusing anything that is not one.
#
# `gregexpr()` silently skips text it cannot match, which would turn `{"a" @ 1}` into a
# perfectly good-looking `{"a" 1}` and then into a confusing structural error. So the gaps
# between matches are checked: every one must be whitespace. The bridge writes compact JSON, so
# in practice almost every gap is empty and the check costs nothing.
json_tokenize <- function(txt) {
  if (!is.character(txt) || length(txt) != 1L || is.na(txt)) {
    stop("json_tokenize() wants a single non-NA string.", call. = FALSE)
  }

  found <- gregexpr(JSON_TOKEN_PATTERN, txt, perl = TRUE)[[1]]
  starts <- as.integer(found)
  if (length(starts) == 1L && starts[1L] == -1L) {
    starts <- integer(0)
    lengths <- integer(0)
  } else {
    lengths <- attr(found, "match.length")
  }

  # Gap i runs from just after token i-1 to just before token i, with the last gap running to
  # the end of the input. Empty gaps (end < start) are dropped before the whitespace test.
  gap_start <- c(1L, starts + lengths)
  gap_end <- c(starts - 1L, nchar(txt))
  wanted <- gap_end >= gap_start
  if (any(wanted)) {
    gaps <- substring(txt, gap_start[wanted], gap_end[wanted])
    stray <- gaps[!grepl("^[ \t\r\n]*$", gaps)]
    if (length(stray) > 0L) {
      stop_json(paste0(
        "not valid JSON: unexpected ", encodeString(substr(trimws(stray[1L]), 1L, 40L),
          quote = "'"
        ),
        " outside any value"
      ))
    }
  }

  # Not merely an optimisation: `substring()` raises "invalid substring arguments" on
  # zero-length positions rather than returning `character(0)`, so empty or all-whitespace
  # input has to return before reaching it. The gap check above has already run, so input that
  # is *only* junk — with no tokens at all — is still refused rather than read as empty.
  if (length(starts) == 0L) {
    return(character(0))
  }

  substring(txt, starts, starts + lengths - 1L)
}

# ---------------------------------------------------------------- errors

# Every failure in this file is one condition class, so a caller can catch the lot.
#
# `mzlib_protocol_error` and not a plain `stop()`: unreadable output from the bridge is a
# protocol failure in exactly the same sense as a version mismatch, and a user should not have
# to know that the reason we could not read the answer was the JSON parser rather than the
# envelope check.
stop_json <- function(message) {
  stop(mzlib_condition(
    class = "mzlib_protocol_error",
    message = paste0("mzLib bridge returned unreadable JSON: ", message)
  ))
}

# ---------------------------------------------------------------- parsing

# The parse cursor. An environment, not a list, so advancing it does not copy the token vector
# at every step — with a FlashLFQ payload that vector has hundreds of thousands of entries.
json_state <- function(tokens) {
  state <- new.env(parent = emptyenv())
  state$tokens <- tokens
  state$i <- 1L
  state$n <- length(tokens)
  state
}

# The next token, consumed. Running off the end is its own message because "unexpected end of
# input" and "unexpected token" send a reader to completely different places.
json_take <- function(state) {
  if (state$i > state$n) {
    stop_json("unexpected end of input")
  }
  token <- state$tokens[state$i]
  state$i <- state$i + 1L
  token
}

# The next token, not consumed; `NA_character_` at the end of input.
json_peek <- function(state) {
  if (state$i > state$n) NA_character_ else state$tokens[state$i]
}

# How deep a nesting this reader will follow.
#
# The wire format nests four deep at its worst, so this only ever fires on input that is not
# ours. It exists because R's own expression-depth limit would otherwise turn a pathological
# document into a "evaluation nested too deeply" error that says nothing about JSON.
JSON_MAX_DEPTH <- 200L

# Read a complete JSON document.
#
# Returns a named list for an object, an unnamed list for an array, a length-1 vector for a
# scalar, and `NA` for `null`. Never a data.frame, never a simplified vector, whatever the
# shape of the data.
json_parse <- function(txt) {
  state <- json_state(json_tokenize(txt))
  if (state$n == 0L) {
    stop_json("the bridge produced no output to parse")
  }

  value <- json_value(state, 1L)
  if (state$i <= state$n) {
    stop_json(paste0(
      "trailing content after the value: ",
      encodeString(substr(state$tokens[state$i], 1L, 40L), quote = "'")
    ))
  }
  value
}

# One value of any kind, dispatched on its first character — which is enough, because the
# tokeniser has already decided what each token is.
json_value <- function(state, depth) {
  if (depth > JSON_MAX_DEPTH) {
    stop_json(paste0("nested more than ", JSON_MAX_DEPTH, " deep"))
  }

  token <- json_peek(state)
  if (is.na(token)) {
    stop_json("unexpected end of input")
  }

  head <- substr(token, 1L, 1L)
  if (head == "{") {
    return(json_object(state, depth))
  }
  if (head == "[") {
    return(json_array(state, depth))
  }

  state$i <- state$i + 1L
  if (head == '"') {
    return(json_unescape(token))
  }
  if (token == "true") {
    return(TRUE)
  }
  if (token == "false") {
    return(FALSE)
  }
  if (token == "null") {
    return(JSON_NULL)
  }
  if (head == "-" || (head >= "0" && head <= "9")) {
    return(json_number(token))
  }

  stop_json(paste0("unexpected ", encodeString(token, quote = "'")))
}

# Every JSON number becomes a double, including ones that would fit in an integer.
#
# Deliberate. `file_size_bytes` for a PRIDE project runs past 2^31 — the mzXML in PXD000001 is
# 472 MB and a whole project is 1.44 GB — so anything read as R's 32-bit `integer` would
# overflow to `NA` on real data. A double is exact to 2^53, which no field on this wire
# approaches, and it means no caller has to know which fields are large.
json_number <- function(token) {
  value <- as.numeric(token)
  if (is.na(value)) {
    stop_json(paste0("unreadable number ", encodeString(token, quote = "'")))
  }
  value
}

# An object, always as a named list — never simplified, never a data.frame.
#
# `{}` returns a zero-length *named* list, so `names()` is `character(0)` rather than `NULL`
# and code that iterates names does not have to special-case the empty case.
json_object <- function(state, depth) {
  state$i <- state$i + 1L
  if (identical(json_peek(state), "}")) {
    state$i <- state$i + 1L
    return(structure(list(), names = character(0)))
  }

  # Grown by doubling. Appending with `c()` or `out[[length(out) + 1L]]` copies the whole
  # accumulated list every time, which is quadratic; a FlashLFQ payload has enough members for
  # that to be the difference between instant and a visible pause.
  capacity <- 8L
  values <- vector("list", capacity)
  keys <- character(capacity)
  count <- 0L

  repeat {
    key <- json_take(state)
    if (substr(key, 1L, 1L) != '"') {
      stop_json(paste0("expected a member name, got ", encodeString(key, quote = "'")))
    }
    separator <- json_take(state)
    if (separator != ":") {
      stop_json(paste0("expected ':' after a member name, got ", encodeString(separator, quote = "'")))
    }

    value <- json_value(state, depth + 1L)
    count <- count + 1L
    if (count > capacity) {
      capacity <- capacity * 2L
      length(values) <- capacity
      length(keys) <- capacity
    }
    # Safe only because `null` was read as `NA`: assigning `NULL` here would delete the member
    # instead of storing it, which is the whole failure this file exists to prevent.
    values[[count]] <- value
    keys[count] <- json_unescape(key)

    delimiter <- json_take(state)
    if (delimiter == "}") {
      break
    }
    if (delimiter != ",") {
      stop_json(paste0("expected ',' or '}' in an object, got ", encodeString(delimiter, quote = "'")))
    }
  }

  length(values) <- count
  names(values) <- keys[seq_len(count)]
  values
}

# An array, always as an unnamed list — one element per entry, whatever the entries are.
#
# An array of numbers stays a list of length-1 doubles rather than collapsing to a numeric
# vector. That is more verbose to consume and it is the point: the caller that knows the field
# does the collapsing, and it gets the same shape for one element as for a hundred.
json_array <- function(state, depth) {
  state$i <- state$i + 1L
  if (identical(json_peek(state), "]")) {
    state$i <- state$i + 1L
    return(list())
  }

  capacity <- 8L
  values <- vector("list", capacity)
  count <- 0L

  repeat {
    value <- json_value(state, depth + 1L)
    count <- count + 1L
    if (count > capacity) {
      capacity <- capacity * 2L
      length(values) <- capacity
    }
    values[[count]] <- value

    delimiter <- json_take(state)
    if (delimiter == "]") {
      break
    }
    if (delimiter != ",") {
      stop_json(paste0("expected ',' or ']' in an array, got ", encodeString(delimiter, quote = "'")))
    }
  }

  length(values) <- count
  values
}

# ---------------------------------------------------------------- strings

JSON_ESCAPES <- list(
  '"' = '"', "\\" = "\\", "/" = "/",
  "b" = "\b", "f" = "\f", "n" = "\n", "r" = "\r", "t" = "\t"
)

# Turn a quoted token into the string it denotes.
#
# The bridge escapes generously — it writes `'` for an apostrophe, which shows up in every
# error message that quotes a file name — so `\u` is a hot path, not an edge case.
json_unescape <- function(token) {
  body <- substr(token, 2L, nchar(token) - 1L)
  if (!grepl("\\", body, fixed = TRUE)) {
    return(body)
  }

  chars <- strsplit(body, "", fixed = TRUE)[[1L]]
  n <- length(chars)
  out <- character(n)
  written <- 0L
  i <- 1L

  while (i <= n) {
    if (chars[i] != "\\") {
      written <- written + 1L
      out[written] <- chars[i]
      i <- i + 1L
      next
    }

    if (i == n) {
      stop_json("a string ends with a lone backslash")
    }
    escape <- chars[i + 1L]

    if (escape == "u") {
      if (i + 5L > n) {
        stop_json("a \\u escape is missing its four hex digits")
      }
      point <- json_code_point(chars, i)
      i <- i + 6L
      # A character outside the basic plane is written as a surrogate pair, and the two halves
      # only mean anything together — `intToUtf8()` on either alone gives nonsense.
      if (point >= 0xD800L && point <= 0xDBFF) {
        if (i + 5L <= n && chars[i] == "\\" && chars[i + 1L] == "u") {
          low <- json_code_point(chars, i)
          if (low >= 0xDC00L && low <= 0xDFFF) {
            point <- 0x10000L + (point - 0xD800L) * 0x400L + (low - 0xDC00L)
            i <- i + 6L
          }
        }
      }
      written <- written + 1L
      out[written] <- intToUtf8(point)
      next
    }

    replacement <- JSON_ESCAPES[[escape]]
    if (is.null(replacement)) {
      stop_json(paste0("unknown escape ", encodeString(paste0("\\", escape), quote = "'")))
    }
    written <- written + 1L
    out[written] <- replacement
    i <- i + 2L
  }

  paste0(out[seq_len(written)], collapse = "")
}

# The code point of the `\uXXXX` beginning at `chars[at]`.
json_code_point <- function(chars, at) {
  hex <- paste0(chars[(at + 2L):(at + 5L)], collapse = "")
  point <- strtoi(hex, base = 16L)
  if (is.na(point)) {
    stop_json(paste0("bad \\u escape ", encodeString(paste0("\\u", hex), quote = "'")))
  }
  point
}
