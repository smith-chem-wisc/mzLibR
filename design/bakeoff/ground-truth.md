# Ground truth — the exact commands

Every number in `DESIGN.md` is reproducible by driving the bridge executable directly, with no R
package in the loop. Reproduced 2026-07-24 on `bridge 1.0.0.0`, protocol 1, .NET 8.0.27.

Set the bridge and the FlashLFQ test data:

```sh
B="…/pymzlib/_dotnet/win-x64/mzlib-bridge.exe"
D="…/mzLib/Test/bin/Debug/net8.0-windows/FlashLFQ/TestData"
```

## Task 2 — albumin digestion (the exhaustive sweep)

```sh
for prot in "trypsin|P" "trypsin"; do
  for ml in 7 1; do
    for mc in 0 1 2; do
      "$B" peptidoform fragments --accession P02768 --protease "$prot" \
        --dissociation ETD --terminus Both --missed-cleavages "$mc" \
        --min-length "$ml" --max-length 0 --max-mods 2 --max-isoforms 1024 \
      | python -c 'import sys,json; d=json.load(sys.stdin)["data"]; \
          s=[p["base_sequence"] for p in d["peptides"]]; \
          print(prot,mc,ml,"distinct",len(set(s)),"peptidoforms",len(s))'
    done
  done
done
```

Result (distinct base sequences, default missed-cleavages 2):

| protease | min_length | distinct | peptidoforms |
|---|---|---|---|
| trypsin\|P | 7 | 195 | 303 |
| trypsin | 7 | 202 | 310 |
| trypsin\|P | 1 | 243 | 366 |
| trypsin | 1 | 257 | 381 |

**254 and 269 (mzLibRust's design table) appear at no setting.** The full sweep including
missed-cleavages 0 and 1, and modifications off, is in the session log; none of it yields 254 or
269.

## Task 1 — K562 quant, MBR on

```sh
printf '%s\tK562\t1\t1\t1\n%s\tK562\t2\t1\t1\n' \
  "$D/20100614_Velos1_TaGe_SA_K562_3.mzML" "$D/20100614_Velos1_TaGe_SA_K562_4.mzML" \
| "$B" quant flashlfq --psms "$D/AllPSMs.psmtsv" \
    --ppm 10 --isotope-ppm 5 --mbr --mbr-ppm 10 --mbr-q 0.05 --threads 1
```

From the `data` object:

- `identification_count` = 594; `peptides` distinct = 354; `proteins` distinct = 943;
  `peaks` = 647.
- **Both runs quantified**: sequences with an `intensity > 0` peak in each of the two files =
  **257**. Counting any peak (including the 40 with intensity 0) gives 266. From the peptide
  roll-up (`peptides` with `intensity > 0`) = 169.
- **MBR**: peaks with `detection_type == "MBR"` = **140**, all distinct sequences. Per file:
  run_3 = **62**, run_4 = **78**. From the roll-up's `detection_type`: 52 total, run_3 = **0**,
  run_4 = 52.
- **Proteins**: groups all-NA across both runs = **2** (4 NA cells / 2 runs). Groups zero in both
  = 847. Groups with no usable number (all NA or 0) = 849.

The MBR-per-run figures and the protein figures are **independent of the experimental design** —
one condition/two bioreps, two conditions, and bare paths all give the same 140 / 62 / 78 and
2 / 847 / 849.

## Task 3 — PRIDE PXD000001

```sh
"$B" pride files --accession PXD000001 --page-size 100
```

- 8 files, `total_size_bytes` = 514,278,049.
- The MGF is `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`; its extension is `.gz`.
- `"$B" pride files --accession PXD0000019999` returns `{"ok":true,"data":{"files":[]}}` — an
  empty list with a success exit, not a 404.
