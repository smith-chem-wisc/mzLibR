# PRIDE PXD000001 via mzLibR — Arm B answers

Tool: `mzLibR` (R 4.6.1) + mzlib-bridge. Persona: bench proteomics biologist, not a coder.

## Q1 — How many files, and total size in GB?

**Two honest numbers, because the tool and the actual project disagree:**

| Source | Files | Total size |
|---|---|---|
| `pride_list_files()` manifest (what the tool returns) | **8** | **0.51 GB** (514,278,049 bytes, via `pride_total_size_bytes()`) |
| Actual project on the PRIDE FTP tree | **13** | **~1.44 GB** |

- The tool's manifest comes from PRIDE's REST API, which omits 5 files including the
  two largest — a ~450 MB `.mzML` and a ~473 MB `.mzXML` (the modern open-format
  conversions). This omission is PRIDE's, not mzLib's, and the `?pride_list_files`
  help says so outright and tells you to cross-check the FTP directory.
- I cross-checked `ftp.pride.ebi.ac.uk/.../2012/03/PXD000001/`: 11 files in the root
  + 2 in `generated/` (the `.mgf.gz` and `.mztab.gz`) = **13**, exactly matching the docs.
- The manifest total is *also* wrong in the other direction: `pride_total_size_bytes()`
  over-reports compressed files because PRIDE hands back the *decompressed* size (see Q2).

**CONFIDENCE:**
- "8 files / 0.51 GB from the tool" — **yes, paper-grade**, if you state it as "files exposed
  by the PRIDE REST API." Reproducible and exact.
- "13 files / 1.44 GB = the whole project" — **the 13 count: yes** (I verified it on FTP).
  **The 1.44 GB: I'd cite it as ~1.44 GB from the PRIDE docs/FTP, not as my own exact
  measurement** — the byte sizes I scraped from the FTP HTML page looked rounded
  (suspiciously clean multiples of 1024), so I would not put a to-the-byte project total
  in a table. If a paper needs the exact project size I would sum actual downloaded bytes.

## Q2 — Download JUST the peak-list MGF. Success? Actual bytes?

**Success.** File: `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`
**Actual bytes on disk: 5,984,662** (5.98 MB).

- Reported `file_size_bytes` was **16,448,103** — 2.75x larger — because PRIDE reports the
  decompressed size for this `.gz`. On-disk matches the docs' predicted 5,984,662 exactly.
- I selected the single MGF row by matching its name (`.mgf.gz`) and passed it to
  `pride_download_files()`. I deliberately did **not** use `category = "PEAK"`: the docs warn
  that "PEAK" also matches the 243 MB `.mzXML`, so a "just the peak list" download would have
  been ~40x bigger than intended.

**CONFIDENCE: yes, paper-grade.** Deterministic, verified with `file.info()$size`, and the
download streams to a temp name then moves into place, so a truncated file can't masquerade
as complete.

## Q3 — Sanity check "PXD0000019999" (likely typo)

The call **failed cleanly with a useful message**, not a crash:

> ERROR: PRIDE returned no files for 'PXD0000019999'. Either the accession does not exist
> (check for a typo) or the project is private. PRIDE does not distinguish the two, so
> neither can mzLibR.

- The accession is grammatically well-formed, so per the docs it costs exactly one live API
  request before failing (mzLibR intentionally does not hard-code accession width).
- Good behavior: it names the two possibilities (typo vs private) and is honest that PRIDE
  itself can't tell them apart.

**CONFIDENCE: yes** — clear, honest, non-crashing error handling.

## Decisions I made with no guidance
- **Interpreting Q1.** "How many files does the project contain" is ambiguous between "what the
  tool returns" and "what's actually there." I reported both rather than pick one, because the
  docs are explicit that they differ. I would not want a reader to think 8/0.51 GB is the whole project.
- **Selecting the MGF.** Chose a name match over `category = "PEAK"` to avoid dragging in the mzXML.
- **Trusting the docs over my FTP scrape** for the exact project byte total.

## What I wished the docs / tool said (or did)
- The `?pride_list_files` and `?pride_total_size_bytes` help are unusually candid — they told me
  in advance about the 8-vs-13 gap, the decompressed-size inflation, and the PEAK-matches-two trap.
  That saved me a lot of confusion. Credit where due.
- **But the tool makes me leave it to answer a basic question.** It knows the manifest is
  incomplete (it says so) yet gives me no function to fetch the true FTP file list or the true
  project size. A `pride_project_size()` / `pride_ftp_files()` that reads the FTP tree — or even
  a returned attribute like "N files hidden from REST" — would keep me inside the tool instead of
  hand-copying a URL from prose into a browser. That FTP cross-check was my one external lookup.
- I wish `pride_total_size_bytes()` had an option (or a companion) that returns the **actual
  transfer size** (compressed), since that's the number you actually budget a download against.
  Right now the only way to learn the real size is to download the file.
- Minor: `pride_download_files()` returns paths whether or not it transferred anything, and the
  help correctly warns not to read `length()` as work done — but there's no return field telling
  me bytes transferred this run. I verified size myself with `file.info()`.
