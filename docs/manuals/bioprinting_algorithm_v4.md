# `bioprinting_algorithm_v4.m` — Unified Superset Extrusion Solver

> **STATUS:** `active` — the single solver function for straight-needle
> extrusion. It supersedes three older functions, whose physics are folded
> into it:
> - `bioprinting_algorithm_3.m` (Power-Law, straight cylinder)
> - `bioprinting_algorithm_cross_v2.m` (Cross, straight-needle degenerate case)
> - `bioprinting_algorithm_v3.m` (slicer-parameter generator)
>
> The three legacy functions are retained under `02_MATLAB/archive/scripts/`
> for audit reproducibility and for `validate_v4.m`; do not call them for new
> work.

Called by [`run_solver_v4.m`](run_solver_v4.md) — you normally do not invoke
this function directly. It is documented here so the physics and the output
struct are discoverable.

## What it computes

For each piston velocity in `Vp_vec`, and for **both** the Power-Law and the
Cross constitutive models, it returns the complete syringe + needle analysis:

- volumetric flow rate, mean / max velocities (syringe + needle);
- full pressure decomposition: syringe drop, needle drop, hydrostatic head,
  total system drop, and absolute exit / needle-inlet / syringe-inlet
  pressures;
- wall shear rate and stress (syringe + needle; PL also reports the
  Rabinowitsch-true wall shear rate);
- generalised Reynolds numbers, a critical-Re estimate, and a laminar/turbulent
  verdict;
- radial profiles `r`, `u(r)`, `γ̇(r)`, `τ(r)` (syringe + needle);
- FFF-slicer-equivalent parameters: `v_print`, `w_line`, `h_layer`,
  `beta_swell`, `k_flow`.

## Signature

```matlab
results = bioprinting_algorithm_v4(sample, geom, Vp_vec, ...)
```

### `sample` struct

| Field | Meaning |
|---|---|
| `name` | Sample identifier |
| `K_PL` | Power-Law consistency (Pa·sⁿ) |
| `n_PL` | Power-Law flow index (–) |
| `eta0` | Cross zero-shear viscosity (Pa·s) |
| `etaInf` | Cross infinite-shear viscosity (Pa·s) |
| `lambda` | Cross time constant (s) |
| `m_Cross` | Cross exponent (–) |
| `rho` | Density (kg/m³) |
| `Rrec_pct` | Structural recovery (%) at the deposition shear rate |

### `geom` struct

| Field | Meaning |
|---|---|
| `Rs` | Syringe inner radius (m) |
| `R_n` | Needle inner radius (m) |
| `L_n` | Needle length (m) |
| `Ls` | Syringe length (m) |
| `label` | e.g. `'21G'`, `'22G'` |
| `h_factor` | Layer height as a fraction of needle ID (default 0.7) |

### `Vp_vec`

Vector of piston velocities **in metres per second** (`run_solver_v4` passes
`cfg.Vp_mm_s * 1e-3`).

### Name–value options

| Option | Default | Purpose |
|---|---|---|
| `'OutputFolder'` | `'output_v4'` | Output directory |
| `'Orientation'` | `'horizontal'` | `horizontal` / `upward` / `downward` |
| `'IncludeHydrostatic'` | `true` | Toggle the ρgh term |
| `'SaveData'` | `true` | Per-Vp combined `_data.txt` |
| `'SaveFigures'` | `true` | Per-Vp + summary PNGs |
| `'SaveCSV'` | `true` | Across-Vp slicer CSV |
| `'PlotResults'` | `true` | Master switch for any figure |
| `'NumPoints'` | `200` | Radial discretisation (≥ 20) |
| `'f_slip'` | `1.0` | Syringe-piston slip factor |

## Output struct

`results` carries `.PL` and `.Cross` sub-structs — each field is a column
vector over `Vp_vec` (`Q`, `u_avg_needle`, `u_max_needle`, `tau_w_needle`,
`gamma_w_needle`, `dP_needle`, `dP_syr`, `dP_total`, `Re_needle`, …) — plus:

- `.profiles{iVp}` — the per-Vp radial data for both models;
- `.slicer` — `PL` and `CR` slicer parameter sets (`v_print_mm_s`, `w_line_mm`,
  `k_flow`, …);
- `.sample`, `.geom`, `.Vp_mm_s`.

## Flow assumptions

Fully developed laminar flow; no entrance/exit losses; no wall slip; no yield
stress in the solver (yield handled separately by `extract_hmax_v3.py`).
Geometry is a **straight cylindrical needle** (21G/22G blunt hardware). For a
conical tip, use `bioprinting_algorithm_conical.m`.

## Validation

`validate_v4.m` confirms v4 reproduces the legacy solvers: Power-Law fields
match `bioprinting_algorithm_3` exactly; Cross fields match
`bioprinting_algorithm_cross_v2` (straight-needle case) within solver
tolerance; the slicer fields match the v3 formulas. See
[`validate_v4.md`](validate_v4.md).
