# mzLibR

**[mzLib](https://github.com/smith-chem-wisc/mzLib) for R — mass spectrometry and proteomics,
with no .NET to install and no package dependencies.**

The third sibling of [pyMzLib](https://github.com/smith-chem-wisc/pyMzLib) and
[mzLibRust](https://github.com/smith-chem-wisc/mzLibRust), over the same language-neutral
bridge, so all three agree field for field.

```r
files <- pride_list_files("PXD000001")
files[files$size_mb < 5 & files$downloadable, c("file_name", "size_mb", "category")]

digest <- peptidoform_fragments("P02768")           # serum albumin, tryptic, ETD
digest_fragments_by_series(digest)

quant <- flashlfq_quantify("AllPSMs.psmtsv", c("run_3.mzML", "run_4.mzML"),
                           match_between_runs = TRUE)
transfers <- flashlfq_mbr_peaks(quant)              # read peaks, never the peptide roll-up
```

## What it does

| | |
|---|---|
| `pride_*` | List and download from the PRIDE Archive, with mzLib's paging and URL resolution |
| `peptidoform_*` | Fetch a UniProt entry, apply its modifications, digest, and fragment |
| `flashlfq_*` | Label-free quantification with match-between-runs |
| `readers_*` | Identify **and read all 31** result-file types — `readers_read_records()` reads any of them into that format's own fields; `readers_read_results()`, `readers_read_features()`, `readers_read_matches()` and `readers_read_spectra()` project the four cross-format views |

## Installing

```r
install.packages("remotes")
remotes::install_github("smith-chem-wisc/mzLibR")
mzlibr_install_bridge()          # ~140 MB, asks first, verifies a SHA-256
```

The bridge is a self-contained executable carrying its own .NET runtime, so **no .NET
installation is required**. It cannot ship inside the package — CRAN's limit is about 5 MB — so
`mzlibr_install_bridge()` fetches it on your say-so and caches it. If you already have one,
point at it instead and skip the download:

```r
Sys.setenv(MZLIB_BRIDGE = "/path/to/mzlib-bridge")
```

That override is also how you relink a modified mzLib without rebuilding this package, which
LGPL section 4 requires mzLibR to permit.

## No dependencies, and why that made it better

mzLibR imports nothing but base R. `jsonlite` and `processx` would have been the obvious floor,
and dropping them was a correctness decision more than a portability one:

- **`jsonlite` auto-simplifies, and JSON `null` becomes R `NULL`.** A `NULL` assigned into a list
  *deletes* the element, so a protein-intensity vector would not come back wrong, it would come
  back **short**, silently — and FlashLFQ emits `null` for precisely the protein it could not
  resolve. mzLibR's own reader never simplifies and reads `null` as `NA`, which occupies its
  slot like any other value.
- **`processx` was for the subprocess timeout.** `system2()` has had one since R 3.5.0, and
  delivers stdin at the same time, which was verified against the real bridge before the
  dependency was dropped.

Requires **R >= 3.5.0**, set solely by `system2(timeout=)`.

## Things that will otherwise catch you

Every one of these is documented on the argument that causes it, and pinned by a test.

- **A compressed file's extension is `.gz`, not what it is compressed from.** PXD000001's peak
  list is `…pride.mgf.gz`, so filtering on `".mgf"` matches nothing.
- **PRIDE's API manifest is incomplete.** It publishes 8 files for PXD000001; the FTP tree holds
  13, including the two largest. And `file_size_bytes` is the *decompressed* size for some
  compressed files — 2.75x off for that MGF.
- **mzLib's `"trypsin|P"` applies the Keil rule and plain `"trypsin"` does not** — the reverse of
  the MaxQuant and Mascot convention. 195 peptides against 202 on albumin.
- **`min_length` defaults to 7 and silently discards everything shorter** — albumin goes from
  195 distinct sequences to 243 at `min_length = 1`.
- **ETD and ECD return a `y` series with no `b` ions**, about a third of the fragment list. No
  fragmentation mechanism produces that; use `digest_fragments_by_series()`.
- **Read `peaks`, never the peptide roll-up, for match-between-runs.** The roll-up drops most
  transfers — 140 true against 52 shown — and a whole run's can vanish.
- **A peptide intensity of `0` and a protein intensity of `NA` mean different things.** `0` is
  "not measured in this run"; `NA` is "FlashLFQ could not resolve a number". `NA` is the rare
  outcome (2 proteins) and `0` the common one (847).
- **`max_threads` defaults to 1 here, unlike pyMzLib.** Above one thread the roll-up
  nondeterministically drops MBR intensities, so identical inputs disagree roughly 1 run in 6
  ([mzLib#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111)).
- **`readers_identify()` does not validate contents**, and a `TRUE` from `is_quantifiable` is not
  permission — it reports what mzLib's interface offers, not that the numbers are comparable.

## Documentation

- [`docs/name-parity.md`](https://github.com/smith-chem-wisc/mzLibR/blob/main/docs/name-parity.md) — every function, parameter and column checked
  mechanically against pyMzLib and against the wire. Generated, not written.
- [`docs/test-parity.md`](https://github.com/smith-chem-wisc/mzLibR/blob/main/docs/test-parity.md) — what the suite covers against the parent's, what
  R's semantics made unnecessary, and what is deliberately not run.

`?pride_list_files`, `?peptidoform_fragments`, `?flashlfq_quantify` and `?readers_identify` carry
the detail.

## Upstream issues found by this port

An mzLib bug becomes an mzLib issue; a workaround in a binding leaves every other consumer
broken. See mzLib
[#1106](https://github.com/smith-chem-wisc/mzLib/issues/1106),
[#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109),
[#1110](https://github.com/smith-chem-wisc/mzLib/issues/1110),
[#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111),
[#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112),
[#1113](https://github.com/smith-chem-wisc/mzLib/issues/1113).

## Licence

LGPL-3.0-or-later, matching mzLib, which this package redistributes in compiled form. See
[`LICENSE.note`](https://github.com/smith-chem-wisc/mzLibR/blob/main/LICENSE.note).
