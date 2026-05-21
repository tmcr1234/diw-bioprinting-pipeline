# `Fit_Muitos_Modelos_v4.py` — Rheological Model Fitting

> **STATUS:** `active`. Output of this script feeds directly into the MATLAB extrusion solver. Re-fit every new ink before simulating.
> **EXTRACTION LAYER UNDER REVIEW** — see *Caveats* below.

## Purpose

Fits **seven constitutive models** (Herschel–Bulkley, Power Law, Bingham, Casson,
Cross, Carreau, Ellis) to a steady-state flow curve and ranks them by AIC/BIC in
log space. The two parameters the rest of the pipeline needs are:

- **Power-Law:** `K`, `n` → fed into `bioprinting_algorithm_3.m` / `_v3.m`
- **Cross:** `eta0`, `etaInf`, `lambda`, `m` → fed into `bioprinting_algorithm_cross_v2.m` / `_v3.m`

## Inputs

| Item | Path / format | Required |
|---|---|---|
| Flow-curve workbook (one per sample) | `./Reologia/Viscosity/<sample> - Visco_Artur CP50 0,1mm.xls` | yes |
| Sheet name | `Ramp - 2` (final ramp) — change `SHEET_PREFIX` in script | yes |
| Sample list | edit `SAMPLES` (line 208) | yes |

**Column layout assumed inside the `.xls`** (after `skiprows=3`, 0-based):
`COL_STRESS=2`, `COL_SHEAR_RATE=3`, `COL_VISCOSITY=4`.

## Outputs

Written to `./Analises/Python/Results/`:

| File | Contents |
|---|---|
| `FitAll-Ramp2-v4.txt` | Human-readable per-sample fit table + AIC/BIC ranking |
| `FitAll-Ramp2-v4.csv` | Same data, machine-readable (for plotting & MATLAB hand-off) |

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/Fit_Muitos_Modelos_v4.py
```

The script uses relative paths (`./Reologia/Viscosity`, `./Analises/Python/Results`)
so it **must be executed from the project root**, not from `Export/01_Python/`.

## Dependencies

`numpy`, `pandas`, `scipy`, **`xlrd`** (required for legacy `.xls`; do NOT
substitute `openpyxl` — that handles `.xlsx` only).

## Why log-space fitting?

Rheological data spans several decades. Linear-space OLS over-weights high-stress
points; log-space minimises *relative* error and is equivalent to WLS with
`w = 1/y²`. AIC/BIC are also computed in log space so model rankings are
consistent with the fit objective.

## Caveats — rheology data ingestion is UNDER REVIEW

The script currently parses **`.xls` workbooks exported by the rheometer this lab
used historically** (TA Instruments format, sheet structure `Ramp - 1 / Ramp - 2`,
columns ordered stress | shear-rate | viscosity at fixed offsets).

When you start working with a different rheometer (e.g., a different Anton Paar
model, or any new equipment that exports a different workbook layout), the
**ingestion block (lines ~210–280) will need rewriting**, not the fitting math.

### Validation checklist before trusting a fit on new data

1. Open one of your `.xls` files manually. Confirm:
   - which sheet contains the flow curve (current default: `Ramp - 2`),
   - how many rows of header you need to skip (current default: `skiprows=3`),
   - the column indices for stress, shear rate, viscosity.
2. Update `FILE_PATH`, `FILE_SUFFIX`, `SHEET_PREFIX`, `COL_STRESS`,
   `COL_SHEAR_RATE`, `COL_VISCOSITY` at the top of the script.
3. Run on one sample and inspect the `.txt` report. The Power-Law `K` should be
   on the order of 10–1000 Pa·sⁿ for typical DIW gels and `n` ∈ (0.1, 0.6).
   Values far outside this band usually mean a column-mapping error.

If the new equipment exports `.xlsx` instead of `.xls`, the fix is two lines:
- `engine="xlrd"` → `engine="openpyxl"` in `pd.ExcelFile(...)`
- file glob `*.xls` → `*.xlsx`

## Hand-off to MATLAB

After running, copy the K, n (Power Law) and eta0, etaInf, lambda, m (Cross)
values from the `.txt` report into the `samples(...)` struct at the top of
`02_MATLAB/run_solver_v3.m`. Keep these two in sync at all times — the MATLAB
simulation results are only as good as the Python fit they ingest.
