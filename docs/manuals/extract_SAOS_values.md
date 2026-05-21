# `extract_SAOS_values.py` — SAOS Value Extractor (Strain + Frequency Sweeps)

> **STATUS:** `active`. Library module — imported by `extract_hmax_v2.py`.
> **EXTRACTION LAYER UNDER REVIEW** — see *Caveats*.

## Purpose

Parses Anton Paar / TA-style `.xls` workbooks from **amplitude (strain) sweeps**
and **frequency sweeps** and extracts the quantities that go into the manuscript
tables and into the `h_max` printability calculation:

**From the amplitude sweep:**
- `G'_LVR`, `G''_LVR`, `tan δ_LVR` — plateau values inside the linear
  viscoelastic region
- `γ_LVR` — strain at which G' drops to 90 % of the plateau (yield-strain proxy)

**From the frequency sweep:**
- `G'`, `G''`, `tan δ` at ω ≈ 1, 10, 100 rad/s

## Inputs

| Item | Path / format |
|---|---|
| Strain sweeps | `./Reologia/Strain/<sample> - Deformacao CP50 0,1mm.xls` |
| Frequency sweeps | `./Reologia/Frequency/<sample> - Freq CP50 0,1mm.xls` |

File paths are hard-coded in `main()` — edit the `strain_files` and
`freq_files` dicts (lines ~240–250) for a new sample set.

## Outputs

- **stdout only** when run standalone (`python extract_SAOS_values.py`)
- Two reusable functions that other scripts import:
  - `read_antpar_xls(filepath)` → `(headers, units, ndarray, sheet_idx)`
  - `extract_amplitude_sweep(filepath)` → dict with all LVR metrics

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/extract_SAOS_values.py
```

## Dependencies

`xlrd`, `numpy`. No pandas (talks directly to `xlrd` so it can detect header
rows by keyword scanning).

## Key design choice — header detection

Rather than assume fixed column indices, the parser **scans every row** of every
sheet looking for at least two SAOS keywords (`G'`, `G''`, `Omega`, `Strain`,
…). The first such row becomes the header, the row after it is treated as the
units row, everything below is numeric data.

This is **why this script tolerates moderate format drift** between rheometer
models — as long as the column labels contain recognisable substrings, it will
find them. New equipment that uses entirely different labels (e.g., German
abbreviations) will need the keyword list (line 39) extended.

## LVR detection (amplitude sweep)

1. Slide a 5-point window over G' vs strain.
2. Find the longest contiguous run where local CV(G') < 5 %.
3. Take that run as the LVR plateau; average G', G'', tan δ over it.
4. `γ_LVR` = strain at which G' first falls below 0.9·plateau (log-interpolated).

If no clean plateau is found, falls back to median ± 15 % over the low-strain
half of the curve.

## Caveats — UNDER REVIEW

This module is the **most format-tolerant** of all the extractors but still
assumes:

- Workbook is `.xls` (binary, `xlrd`-readable). `.xlsx` exports will fail —
  switch the engine.
- Column headers contain at least one of: `G'`, `Storage`, `G''`, `Loss`,
  `Omega`, `Frequency`, `Strain`, `Deforma`.
- Either an angular-frequency column (rad/s) or a frequency column (Hz) is
  present. Hz is auto-converted via ω = 2π·f.

If a new equipment export uses non-standard headers (e.g., `Modulus_real`
instead of `G'`), extend the `keywords` list at line 39 and the `find_col()`
calls in `extract_amplitude_sweep` / `extract_frequency_sweep`.
