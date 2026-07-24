# PRIDE PXD000001 — Arm A (R / `rpx`) answers

Tool assigned: R + the `rpx` package (v2.20.0). I stayed inside `rpx`'s documented
interface (`?PXDataset`, `pxfiles`, `pxget`) and only left it when `rpx` did not
supply something I needed. Everywhere I left the package is flagged below and in the log.

---

## Q1 — How many files, and total size in GB?

**File count: ambiguous — the tools disagree.**
- `rpx` (my assigned tool, reading the ProteomeXchange record): **11 files**.
- PRIDE REST API v3 (had to leave `rpx` to get it): **8 files**.

They differ both in count *and* in which files. `rpx` lists three RAW/peak variants
(`...-20141210.mzML`, `...-20141210.mzXML`, `...01.mzXML`) plus a README and a second
mzTab; the PRIDE API lists a leaner 8-file set that instead includes a `.pride.mgf.gz`
peak list that `rpx` never shows. I could not reconcile this from the documentation.

**Total size: `rpx` cannot answer this at all.** Nothing in `rpx`'s documented
interface returns file sizes — `pxfiles(..., as.vector = FALSE)` gives
`ID | NAME | URI | TYPE | MAPPINGS | PX`, no bytes. I got sizes from the PRIDE API:

- Sum of `fileSizeBytes` over the 8 API files = **514,278,049 bytes ≈ 0.51 GB**
  (decimal, 10^9) / 0.48 GiB (2^30).
- **Caveat that matters:** for the gzipped files, the API's `fileSizeBytes` is the
  *uncompressed* size, not the downloadable size (proven below on the MGF: API says
  16,448,103 but the actual `.gz` is 5,984,662). So 0.51 GB is a mix of compressed-on-disk
  and uncompressed logical sizes — it is **not** a clean "download this many bytes" number.

**CONFIDENCE — would I put this in a paper?**
- 11 files (or 8): **No, not as a bare number.** I would only state it as
  "PXD000001, as listed by `rpx`/ProteomeXchange = 11 files; as listed by the PRIDE
  archive API = 8," because they genuinely disagree and I can't tell which is canonical.
- 0.51 GB total: **No.** It came from a second source `rpx` couldn't give me, and it
  conflates compressed and uncompressed sizes. I'd call it "roughly half a gigabyte" in
  prose, never a table value.

## Q2 — Download just the peak-list MGF; did it work; actual bytes?

**Succeeded — but NOT through `rpx`.** `rpx` never lists an MGF, so the documented
download call fails: `pxget(px, "...pride.mgf.gz")` returns **"No files to download."**
The MGF exists in PRIDE's `generated/` subfolder, which `rpx` doesn't expose. I
downloaded it directly from the FTP/HTTPS URL the PRIDE API gave me:

`.../PXD000001/generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`

- **Bytes downloaded: 5,984,662** (on disk). `curl -I` Content-Length = 5,984,662, and
  `gzip -t` passes → complete, valid gzip.
- Uncompresses to **16,448,103 bytes** of MGF text (`BEGIN IONS / TITLE=... / PEPMASS=...`),
  which equals the API's `fileSizeBytes` — that's why that number looked like a "truncation"
  at first. It wasn't; API size = uncompressed size.

**CONFIDENCE:** **High / Yes** that I have the correct, complete MGF peak list for this
project (valid gzip, real MGF content, size confirmed two ways). The only asterisk is that
this is a PRIDE-*generated* MGF, not an author-submitted peak list — worth stating if it
fed a result.

## Q3 — List files for "PXD0000019999" (likely a typo)

`PXDataset2("PXD0000019999")` in `rpx` returns a **warning and no file list**:
`cannot open URL '.../pride/ws/archive/v2/projects/PXD0000019999': HTTP status was '404 Not Found'`.
So a bad accession fails as a raw HTTP 404, not a friendly "invalid/unknown accession"
message. Good enough to know it doesn't exist; a novice could easily mistake it for a
network problem rather than a typo.

**CONFIDENCE:** High — it clearly does not resolve to a project.

---

## Decisions I had to make with NO guidance
1. **Which file count to trust** (rpx's 11 vs API's 8). I reported both rather than pick.
2. **`PXDataset` vs `PXDataset2`.** `?PXDataset` says v1 is deprecated, so I switched to
   `PXDataset2` — but the deprecated one is what every old tutorial shows.
3. **`pxfiles(n = 10)` truncates.** Default `n` is 10 and there are 11 files; I had to
   guess to raise `n` or I'd silently miss a file.
4. **GB vs GiB, and compressed vs uncompressed.** No steer on which "size" is meant. I
   reported decimal GB and flagged the compressed/uncompressed mismatch.
5. **How to get the MGF at all**, since `rpx` won't. I fell back to the raw FTP/HTTPS URL.
6. **Trusting a "success" return code.** `download.file` returned 0 while the byte count
   didn't match the metadata; I had to independently verify (HEAD + gzip test) that the
   file was actually whole.

## What I wish the documentation had told me
- **That `rpx` returns no file sizes at all**, and where to get them (so I'm not left
  guessing whether I missed an argument). This forced every size answer out to the raw API.
- **Why `rpx`'s file list differs from the PRIDE website/API** (11 vs 8), and that
  PRIDE-`generated/` files like the MGF exist but are invisible to `pxget`. The whole
  "download the MGF" task is impossible via the documented `rpx` path with no hint why.
- **That `fileSizeBytes` in the PRIDE API is the *uncompressed* size for `.gz` files.**
  This is a genuine trap — it made a perfectly complete download look truncated.
- **A clearer error for a bad accession** than a bare HTTP 404.
- That `pxfiles` **silently truncates at n = 10**.

## Files written
- `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz` (5,984,662 bytes, the MGF)
- `pride_armA_log.jsonl` (attempt-by-attempt instrumentation)
- `ANSWER.md` (this file)
