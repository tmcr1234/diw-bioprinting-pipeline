# `extract_hmax_v3.py` — Maximum Scaffold Height (CSV-only)

> **STATUS:** `active` (current). The final v3 in the Anton Paar migration
> chain. Closes out the last `xlrd` dependency in the pipeline.

## Purpose

Predicts how tall a printed scaffold can grow before its own weight
collapses it, under four yield criteria, using **only the CSV outputs of
the v2 SAOS extractors**. No raw `.xls` or `.csv` parsing happens here —
all the upstream parsing + fitting was done once by `extract_SAOS_amplitude_v2`
and `extract_SAOS_frequency_v2`.

## What changed vs v2

| Aspect | v2 | v3 |
|---|---|---|
| Reads | `.xls` amplitude + frequency sweeps via `xlrd` | `SAOS_amplitude_v2.csv` + `SAOS_frequency_v2.csv` |
| Re-fits the low-ω power law? | Yes (every run) | No — uses the fit the freq v2 reported |
| Hard dependency on matplotlib | Yes (top-level import) | No — lazy-imported only when plotting |
| Output filenames | `SAOS_hmax_v2.txt`, `printing_parameters_per_ink.csv` | `SAOS_hmax_v3.txt`, `printing_parameters_per_ink_v3.csv` |
| Physics | identical (4 criteria + buckling + N_max) | identical |

Numerical results match v2 to within float precision when the upstream
v2 CSVs were produced with the same `OMEGA_EXTRAP = 0.01 rad/s` setting.

## Yield criteria (same math as v2)

With ρ = 1000 kg/m³ and g = 9.81 m/s², the self-supported stacked height is

| Criterion | σ_max formula | Physical meaning |
|---|---|---|
| **(A) LVR endpoint** | `G'_LVR · γ_LVR` (γ fractional) | stress at which the ink leaves linear viscoelasticity |
| **(B) Practical** | `G'(ω = 1 rad/s)` (measured) | modulus at a deposition-relevant timescale |
| **(C) Quasi-static** | `G'(ω = 0.01 rad/s)` (extrapolated) | modulus on long-dwell timescales |
| **(D) Buckling (Smay 2002)** | `h_max(L) = √(1.94 G' / ρg) · √L` | unsupported-span limit; uses G'(ω=1) |

`h_max = σ_max / (ρ g)` for criteria A–C. For each criterion's `h_max` and
each layer height in `H_LAYER_TARGETS_MM`, the script also reports
`N_max = floor(h_max / h_layer)` — the predicted maximum layer count.

## Inputs

| Input | Required columns | Source |
|---|---|---|
| `SAOS_amplitude_v2.csv` | `sample`, `Gp_LVR_Pa`, `gamma_LVR_pct` | `extract_SAOS_amplitude_v2.py` |
| `SAOS_frequency_v2.csv` | `sample`, `Gp_omega1_Pa`, `lowomega_G0_Pa`, `lowomega_beta`, `Gp_extrap_Pa`, `omega_extrap_rad_s` | `extract_SAOS_frequency_v2.py` |

The two CSVs are joined on `sample`. Inks present in only one of them
are skipped with a clear warning in the report.

Edit `AMPLITUDE_CSV` and `FREQUENCY_CSV` at the top of the script if
your output folder differs from the default.

## Optional plotting

Set `FREQUENCY_RAW_FOLDER` to the folder containing the raw Rheocompass
frequency-sweep CSVs to enable per-ink low-ω extrapolation plots. The fit
overlaid on the plot uses the **exact (G0, β) reported by the v2 freq
CSV** — no re-fitting, no risk of plot/table drift.

If `matplotlib` is not installed in your environment, the script falls
back gracefully — it skips plotting and prints a one-line note. The core
h_max calculation has no matplotlib dependency at all (lazy import).

## Outputs

Written to `./Analises/Python/Results/`:

| File | Contents |
|---|---|
| `SAOS_hmax_v3.txt` | Per-sample report (inputs, yield-stress estimates, h_max per criterion, buckling at 3 spans, N_max layer-count grid) + summary table |
| `printing_parameters_per_ink_v3.csv` | Machine-readable per-ink table (all four h_max criteria, buckling at 3 spans, N_max for 3 layer heights × 3 criteria) |
| `Gprime_extrap_<sample>_v3.png` | Per-ink low-ω fit visualisation (optional — only if `FREQUENCY_RAW_FOLDER` is set AND matplotlib is installed) |

## How to run

```bash
cd "<PROJECT ROOT>"

# 1. Make sure the upstream v2 CSVs are up to date
python Export/01_Python/extract_SAOS_amplitude_v2.py
python Export/01_Python/extract_SAOS_frequency_v2.py

# 2. Run v3
python Export/01_Python/extract_hmax_v3.py
```

If you change the SAOS data, regenerate the upstream v2 CSVs first
(they're the single source of truth) — then re-run v3.

## Dependencies

`numpy`, `pandas`. `matplotlib` only if you want the plots.

## Caveats

- **The script enforces consistency with the v2 fits.** If you want
  different low-ω fit parameters (e.g. n_low = 10 instead of 6),
  change the setting in `extract_SAOS_frequency_v2.py`, regenerate the
  CSV, then re-run v3. Do **not** override the parameters here — that
  defeats the single-source-of-truth design.
- **Criterion (C) extrapolation quality** depends entirely on the
  upstream `lowomega_R2`. Inspect that column in the freq v2 CSV: R² < 0.95
  means the extrapolation is shaky and (C) should be treated as an
  upper bound only.
- **ρ = 1000 kg/m³ default**. If your ink is significantly denser or
  lighter than water, edit `RHO` at the top of the script. `h_max` scales
  inversely with density.
- **Geometry-agnostic.** This script does not know anything about needle
  or syringe dimensions — those are the MATLAB solver's concern. h_max
  comes from material properties alone.

## Hand-off

`printing_parameters_per_ink_v3.csv` is structured to drop into LaTeX
manuscript tables and to feed any future automated print-planning script
(e.g. choosing a needle gauge that won't exceed N_max for a given
shape's required height).

Per-criterion guidance for picking the right h_max in a print
specification:

| When you are printing... | Use criterion |
|---|---|
| Single fast print, < 1 min per scaffold | (B) G'(ω=1) — closest to deposition timescale |
| Multi-layer print, > 10 min total | (A) LVR or (C) extrapolated — long-dwell relevant |
| Lattice / scaffold with overhangs | (D) buckling, with L = your unsupported span |
| Conservative manuscript bound | (C) extrapolated — most pessimistic |

## Versioning note

This is **v3** — supersedes `extract_hmax_v2.py`. v2 will be moved to
`deprecated` in the registry once a few projects have used v3 end-to-end
and confirmed numerical parity. v2 stays in the Export folder for now as
the audit reference for any previously published numbers.
