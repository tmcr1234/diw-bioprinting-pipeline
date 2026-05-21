# `run_solver_v3.m` — Unified PL + Cross Extrusion Solver Driver

> **STATUS:** `active`. **This is the primary MATLAB entry point** for the
> simulation stage. New students should start here.

## Purpose

Top-level driver that loops over every (ink × needle) combination, runs both
the **Power-Law** and **Cross** solvers, and writes per-sample slicer-lookup
CSVs plus a master summary table at a reference piston velocity.

## What it produces

Inside `output_v3/`:

| File | Contents |
|---|---|
| `slicer_lookup_<sample>_<needle>.csv` | Full sweep over piston velocities; PL & Cross side-by-side. Hand this to whoever is running the printer. |
| `plots_<sample>_<needle>.png` | 4-panel: ΔP, wall shear rate, head speed v_print, and k_flow |
| `master_summary_v3.csv` | One row per (sample × needle) at the reference piston velocity Vp = 0.01 mm/s |

## Inputs

**Edit the `samples(...)` struct at the top of the script** (lines 18–43):

| Field | From Python script | Notes |
|---|---|---|
| `K_PL`, `n_PL` | `Fit_Muitos_Modelos_v4.py` → Power-Law row | Required |
| `eta0`, `etaInf`, `lambda`, `m_Cross` | `Fit_Muitos_Modelos_v4.py` → Cross row | Required |
| `Rrec_pct` | `Recovery_v1.py` → recovery % at deposition shear (typ. 150 s⁻¹) | Required |

**Geometries** (lines 49–61): syringe radius `Rs`, needle inner radius `R_n`,
needle length `L_n`. Defaults are BD 10 mL + 21G/22G blunt needles. Edit if you
print with different hardware.

**Piston velocity sweep** (line 69): default `[0.003 … 0.04] mm/s`. This sweep
is calibrated to the lab DIW printer's mechanical envelope — head speeds of
~2–30 mm/s. Values above 0.05 mm/s are mechanically unreachable on the current
rig.

## How to run

```matlab
cd '<PROJECT ROOT>'
addpath('Export/02_MATLAB')
run_solver_v3
```

## Output interpretation

The two columns the operator actually needs from `master_summary_v3.csv`:

- **`Vp_mm_s`** — the *input* sent to the printer (piston / extruder motor speed)
- **`PL_v_print_mm_s`** or **`Cross_v_print_mm_s`** — the resulting printhead
  travel speed needed to deposit a single-strand line of width
  `w_line ≈ 2·R·(1+β_swell)` at the chosen layer height

`k_flow` is the **slicer flow factor**: this is the predicted deposition
efficiency. The slicer's **Extrusion Multiplier should be set to `1/k_flow`**
(see `kflow_inversion_convention.md` memory).

## Dependencies

MATLAB R2020a+ (uses `inputParser`, `fzero`, `trapz`, `writetable`,
`exportgraphics`). Calls `bioprinting_algorithm_v3.m` — both files must be on
the path.

## Caveats

- Solver assumes fully developed laminar Newtonian-like profile inside a
  **straight cylindrical needle**. For tapered nozzles use `run_solver_Cross_v2.m`
  + `bioprinting_algorithm_cross_v2.m` or `bioprinting_algorithm_conical.m`.
- No wall slip, no yield stress in the solver itself (yield is handled
  separately via the `h_max` calculation in `extract_hmax_v2.py`).
- `β_swell` is a heuristic die-swell proxy (`0.30·(1 − n)`). It is intentionally
  conservative and not validated against velocimetry for this ink system.
