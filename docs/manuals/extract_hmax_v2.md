# `extract_hmax_v2.py` — Maximum Scaffold Height & Layer Count

> **STATUS:** `active`. Supersedes `extract_hmax_v1.py` (deprecated, not shipped).
> Depends on `extract_SAOS_values.py` — keep both together.

## Purpose

Predicts how tall a printed scaffold can grow before its own weight collapses
it, under **three different yield criteria** plus a buckling-corrected version,
using only rheological inputs (no printing trial required).

## Method (one paragraph)

Self-weight stress at the base of a column of height `h` is `σ = ρ g h`. Solving
for `h_max` requires a yield-stress estimate `σ_max`. The script computes three:

| Criterion | `σ_max` formula | Physical meaning |
|---|---|---|
| (A) LVR endpoint | `G'_LVR · γ_LVR` | Stress at which the ink leaves linear viscoelasticity |
| (B) Practical | `G'(ω = 1 rad/s)` | Modulus at a deposition-relevant timescale |
| (C) Quasi-static | `G'(ω = 0.01 rad/s)` (extrapolated) | Modulus on long-dwell timescales (low ω power-law fit) |

Plus a buckling-corrected unsupported-span estimate (Smay 2002):
`h_max(L) = √(1.94 G' / ρg) · √L`.

Finally, layer counts `N_max = ⌊h_max / h_layer⌋` are reported for typical DIW
layer heights (0.3, 0.4, 0.5 mm).

## Inputs

Same files as `extract_SAOS_values.py`:

- `./Reologia/Strain/<sample> - Deformacao CP50 0,1mm.xls`
- `./Reologia/Frequency/<sample> - Freq CP50 0,1mm.xls`

Edit `FREQ_FILES` and `STRAIN_FILES` dicts (lines 58–67) for a new sample set.

## Outputs

Written to `./Analises/Python/Results/`:

| File | Contents |
|---|---|
| `SAOS_hmax_v2.txt` | Full report (stdout copy) |
| `printing_parameters_per_ink.csv` | Machine-readable per-ink table (all three criteria + buckling + N_max) |
| `Gprime_extrap_<ink>.png` | Low-ω power-law fit visualisation, one per ink |

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/extract_hmax_v2.py
```

## Dependencies

`numpy`, `matplotlib`, `xlrd`, plus the local module `extract_SAOS_values.py`.

## Constants you may want to change

| Constant | Default | Where |
|---|---|---|
| `RHO` (density) | 1000 kg/m³ | line 52 |
| `G_ACC` (gravity) | 9.81 m/s² | line 53 |
| `OMEGA_EXTRAP` | 0.01 rad/s | line 54 |
| `H_LAYER_TARGETS_MM` | [0.30, 0.40, 0.50] | line 55 |

If your ink is significantly denser or lighter than water, update `RHO` —
`h_max` scales inversely with ρ.

## Caveats — UNDER REVIEW

Inherits the data-ingestion caveats of `extract_SAOS_values.py`. The math
itself is equipment-agnostic.

The low-ω power-law extrapolation uses the lowest 6 measured points; if your
frequency sweep starts above ω = 0.1 rad/s, the extrapolation to ω = 0.01 rad/s
is large and the (C) criterion should be treated as an upper bound only.
