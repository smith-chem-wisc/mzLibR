# Test parity with pyMzLib

pyMzLib is the parent binding. This records what mzLibR's suite covers against it: what maps,
what R's semantics made unnecessary, what is deliberately not ported, and — the part that
matters most — **what is genuinely missing**.

Nothing is silently dropped. If a parent test has no counterpart here, it is named below with a
reason.

| suite | tests |
|---|---|
| pyMzLib | **154** |
| mzLibRust | **137** |
| mzLibR | **214** with a bridge staged; **199** offline, the other 15 skipping |

The offline suite passes with no bridge, no .NET and no network. Live tests `skip()`, never fail,
on `mzlib_service_unavailable` — a verdict `cargo test` does not have, so mzLibR can do this
properly where mzLibRust could only print.

## `readers` — found missing, then implemented as M4

**The bridge exposes eight verbs. PLAN.md section 2 describes five.** The three it omits are
`readers formats`, `readers identify` and `readers read-results`, which pyMzLib implements in
`readers.py` with 31 tests. **mzLibRust still does not implement them** — so this was a gap the
whole downstream family shared, and mzLibR is now the only binding besides the parent to close
it.

It was not marginal. mzLib recognises 29 result-file types, and `readers_identify()` is what
tells you which one a path is — including whether it is *quantifiable*, the precondition for
`flashlfq_quantify()`. Exactly three of the 29 carry that view: MetaMorpheus `psmtsv` and
`osmtsv`, and `MsFraggerPsm`. Without it, a user pointing `flashlfq_quantify()` at a PSM file
had no supported way to ask what that file was first.

Two things the port established that are worth carrying back to the other bindings:

- **`is_quantifiable` being `TRUE` is not permission.** It reports what mzLib's *interface*
  offers. `MsFraggerPsm` has the view and should not be quantified — its retention times are in
  seconds while MetaMorpheus's are in minutes, and mzLib normalises neither.
- **`identify()` does not validate contents.** A plain text file containing
  `"this is not a proteomics result file"` is identified as `CruxResult`, because mzLib
  dispatches `.txt` there and `identify()` deliberately stops at resolving the type. The honest
  signal is `views`, which is empty. A live test asserts exactly this, so the behaviour is
  documented rather than discovered by a user.

**mzLibRust should get this module too.** Everything above is a property of mzLib and the wire,
not of R.

## Eliminated by R's semantics

These parent tests defend against something R cannot express. Porting them would have added
unreachable branches, and unreachable checks rot.

| pyMzLib test | why it is not ported |
|---|---|
| `test_null_bytes_are_rejected_rather_than_raising_from_subprocess` | R character vectors may not contain a nul: `"a\0b"` is a parse error and `rawToChar()` refuses to build one. The guard both siblings carry is unreachable here. A test asserts the *language* refuses it, so the elimination is recorded rather than assumed. |
| `test_a_bare_string_of_extensions_is_refused` | Python must distinguish `"raw"` from `["raw"]`, because a string is iterable and would silently become a list of characters. In R a length-1 character vector *is* the correct input, so there is nothing to refuse. |
| `test_as_dict_includes_the_computed_properties` | `PrideFile.as_dict()` has no meaning in R: the data.frame row already is the record, and the computed properties are columns. |
| `test_spectra_must_be_a_list_not_a_string` | Same reason: a character vector of paths is the idiomatic input, and one path is a vector of length one. mzLibR does still refuse a non-character, non-data.frame `spectra`. |
| `test_boolean_max_threads_raises` | Python's `bool` is an `int`, so `max_threads=True` silently means 1. R's `TRUE` is not numeric and `is.numeric()` rejects it, which the existing type check already covers. |

## Deliberately not ported

| pyMzLib test | why |
|---|---|
| `test_a_real_selection_can_be_downloaded_directly` | Downloads real files from EBI. Not run by default: a test suite that pulls megabytes on every invocation is a test suite people stop running. The argument assembly, the stdin framing and the empty-filter refusal are all covered offline. |
| `test_a_real_download_still_works_end_to_end` | As above. |
| — | (the 31 `test_readers.py` tests are now covered; see the M4 section above) |

## Present in mzLibR and not in pyMzLib

R-specific hazards, or checks the parent has no reason to make.

| mzLibR test | what it defends |
|---|---|
| `null keeps its place in an array instead of vanishing` | The largest R-specific hazard there is. `jsonlite` reads JSON `null` as `NULL`, and `NULL` assigned into a list **deletes** the element — a protein-intensity vector comes back *short*, silently. mzLibR's own reader makes it `NA`. |
| `one array element has the same shape as two` | `jsonlite` auto-simplifies, so a manifest with one file parses to a different shape from one with two. |
| `every number is a double, including whole ones` | `file_size_bytes` exceeds 2^31 on real projects; an R 32-bit integer would overflow to `NA`. |
| `the manifest is a plain data.frame with no factors` | `stringsAsFactors` flipped default in R 4.0. Relying on the default gives factors on old R and characters on new. |
| `locations survive subsetting the data.frame` | `[.data.frame` drops attributes, so `locations` had to be a list column or a filtered manifest would silently lose them. |
| `numbers reach the wire with a dot, whatever the session locale` | `format()` honours `options(OutDec)`. A European locale would send `"0,05"` as the MBR q-value and fail to parse — on that user's machine only. |
| `the NA propagates through arithmetic, which is the entire point` | R is the only one of the three whose types can state the 0-versus-NA distinction without prose. A second test greps the sources to prove `na.rm` is never applied to an intensity. |
| `the c and z-dot ladders close on the precursor mass` | `c_k + z_(L-k) = M + 1.00782503` to eight decimals. An off-by-one in fragment numbering could not survive it. Adapted from pyMzLib's live `test_fragment_series_close_on_the_precursor_mass`, but runs offline against the fixture. |
| `arguments are quoted for the shell that will parse them` | `system2()` quotes the command and not the arguments, which pyMzLib never faces because `subprocess` takes a list. |
| the `mzlibr_install_bridge()` refusals | R has neither Python's bundled wheel nor Rust's build script, so the download is a function the user calls — and every way it can refuse is new surface. |

## Where the numbers in the documentation are pinned

PLAN.md section 8.3: every number in a doc string should be reproducible by a test, or removed.
These are the tests that hold the documented figures.

| figure | test |
|---|---|
| PXD000001 publishes 8 files, 0.51 GB | `the recorded manifest holds the eight files PRIDE's API publishes` |
| a compressed file's extension is `.gz` | `a compressed file's extension is .gz, not what it is compressed from` |
| albumin: 195 sequences, 303 peptidoforms | `LIVE: albumin digests to the numbers the documentation quotes` |
| albumin at `min_length = 1`: **243** | `LIVE: min_length is what hides the short peptides` |
| plain `trypsin` gives 202 | `LIVE: the trypsin naming inversion changes the peptide count` |
| 14 of 38 modifications applied | `the census counts sites, applied and annotated separately` |
| ETD's spurious `y` series is about a third | `the spurious y series is about a third of an ETD fragment list` |
| z-dot suppressed at proline, `c` not | `z-dot ions are suppressed at proline while the complementary c ions are not` |

PLAN.md's own figure of **254** for albumin at `min_length = 1` was wrong; it is **243**. The
live test caught it, ground truth was taken by driving the bridge from the shell with no R in
the loop, and PLAN.md is corrected in place. It had never shipped — neither sibling quotes a
`min_length` figure at all.
