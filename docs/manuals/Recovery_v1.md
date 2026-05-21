# `Recovery_v1.py` — Structural Recovery Calculation

> **STATUS:** `active`.
> **EXTRACTION LAYER UNDER REVIEW** — see *Caveats*.

## Purpose

Quantifies the **percentage of structural recovery** of an ink after a
high-shear disruption phase. Implements the classical three-step recovery
protocol:

1. **Ramp 1** — low shear (reference viscosity, η₁)
2. **Ramp 2** — high shear (disruption)
3. **Ramp 3** — low shear (recovered viscosity, η₃)

Recovery (%) = η₃ / η₁ × 100, with uncertainty propagated through the division.

## Inputs

| Item | Path / format | Required |
|---|---|---|
| Recovery workbooks | `./Reologia/Recovery/<sample> - Recovery_<freq>hz CP50 0,1mm.xls` | yes |
| Sample list | edit `samples` (line 8) | yes |
| Shear-rate set | edit `frequencies` (line 9) | yes |

Default sweep: 50, 100, 150, 200, 250, 300 s⁻¹ (the script labels them `hz` for
historical reasons — they are shear rates, not angular frequencies).

## Outputs

Written to `./Analises/Python/Results/`:

| File | Contents |
|---|---|
| `Dados_Recovery_Estruturado.csv` | One row per (sample, shear rate) with η₁, η₂, η₃, recovery%, error% |
| `Relatorio_Recovery.txt` | Human-readable summary + missing-file log |

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/Recovery_v1.py
```

Must be run from the project root (relative paths).

## Dependencies

`pandas`, `xlrd`.

## Hand-off

The single number this script produces that the **MATLAB solver consumes** is
`Rrec_pct` (recovery percentage at deposition shear). Pick the recovery value
closest to your expected wall shear rate inside the nozzle (typically the
150–300 s⁻¹ band) and paste it into `samples(...).Rrec_pct` in
`02_MATLAB/run_solver_v3.m`.

`k_flow = (1 + β_swell)² · √(Rrec/100)` — recovery enters the slicer flow
multiplier under the square root. A drop from 95 % → 50 % recovery cuts
`k_flow` by ~27 %.

## Caveats — UNDER REVIEW

Same as `Fit_Muitos_Modelos_v4.py`: the script assumes the historical `.xls`
layout with column index 4 holding viscosity and at least three `Ramp` sheets
present. On new rheometer exports:

1. Confirm column index for viscosity (script line 44: `df.iloc[:, 4]`).
2. Confirm sheet names begin with `Ramp` (`startswith("Ramp")`, line 39).
3. Confirm there are at least 3 ramps (the script needs the first and third).

If sheet names use a different convention (e.g., `Interval 1`, `Interval 2`),
update `sheet_name.startswith("Ramp")` accordingly.
