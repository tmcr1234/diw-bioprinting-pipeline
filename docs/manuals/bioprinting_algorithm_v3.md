# `bioprinting_algorithm_v3.m` — Unified Extrusion Solver (Library Function)

> **STATUS:** `active`. Called by `run_solver_v3.m`. Not meant to be run
> standalone — it's a function, not a script.

## Purpose

Given an ink (rheology + recovery) and a straight-cylindrical needle geometry,
computes the full set of extrusion observables under **both** the Power-Law and
Cross models, plus the four FFF-slicer-equivalent parameters needed to drive a
DIW printer.

## Signature

```matlab
results = bioprinting_algorithm_v3(sample, geom, Vp_vec, varargin)
```

### Inputs

`sample` — struct with fields:
- `name` — sample label
- `K_PL`, `n_PL` — Power-Law parameters (Pa·sⁿ, dimensionless)
- `eta0`, `etaInf`, `lambda`, `m_Cross` — Cross parameters (Pa·s, Pa·s, s, dimensionless)
- `Rrec_pct` — structural recovery % at deposition shear

`geom` — struct with fields:
- `Rs` — syringe inner radius (m)
- `R_n`, `L_n` — needle inner radius, length (m)
- `label` — e.g. `'21G'`
- `h_factor` — layer height as fraction of needle ID (default 0.7)

`Vp_vec` — vector of piston velocities (m/s)

### Optional name-value pairs

| Name | Default | Purpose |
|---|---|---|
| `OutputFolder` | `'output_v3'` | Where to write CSV + PNG |
| `SaveCSV` | `true` | Disable for batch sweeps |
| `PlotResults` | `true` | Disable for headless runs |
| `f_slip` | `1.0` | Piston-syringe slip factor |
| `NumPoints` | `400` | Cross integration grid |

### Returns

`results` struct with:

```
.sample, .geom, .Vp_mm_s, .Q_mm3_s
.PL    — sub-struct: tau_w_Pa, gamma_w_invs, dP_kPa, v_mean_mm_s,
                     v_print_mm_s, w_line_mm, h_layer_mm,
                     beta_swell, k_flow
.Cross — same sub-struct, computed with Cross model via fzero on the
         Weissenberg-Rabinowitsch integral
```

## Side effects

If `SaveCSV` is true → writes `slicer_lookup_<name>_<needle>.csv`.
If `PlotResults` is true → writes `plots_<name>_<needle>.png`.

## Math summary

**Power-Law analytical:**
`τ_w = K · ((3n+1)Q / (πR³n))ⁿ`
`γ̇_w = (τ_w / K)^(1/n)`
`ΔP = 2 L τ_w / R`

**Cross numerical:**
`Q(τ_w) = (πR³/τ_w³) ∫₀^{τ_w} τ² γ̇(τ) dτ`, inverted with `fzero`. The
inner inversion `η(γ̇)·γ̇ = τ` is itself a `fzero` call. Total cost per Vp is
`NumPoints` viscosity inversions.

**Slicer parameters:**
- `β_swell = max(0, 0.30·(1 − flow_index))` — die-swell proxy
- `w_line = 2R·(1 + β_swell)` — single-strand line width
- `v_print = Q / (w_line · h_layer)` — required head travel speed
- `k_flow = (1+β)² · f_slip · √(R_rec/100)` — deposition efficiency

## Caveats

- Geometry: straight cylinder only. Use `_cross_v2.m` for tapered nozzles.
- Cross solver brackets τ_w in `[1e-3, 1e8]` Pa. If your ink is so viscous that
  the required τ_w exceeds 1e8 Pa, you'll get `NaN` — bump `tau_max` (line 107).
- Plots are saved as PNG (200 dpi). For publication-quality vector output,
  switch `exportgraphics` to `.pdf` (line 238).
