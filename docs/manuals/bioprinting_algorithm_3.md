# `bioprinting_algorithm_3.m` — Power-Law Solver (Detailed, Standalone)

> **STATUS:** `active` (legacy). Largely superseded by `bioprinting_algorithm_v3.m`,
> but kept because it produces richer per-run diagnostics (radial velocity
> profiles, orientation-dependent hydrostatic contributions, detailed TXT
> exports).

## Purpose

Steady-state, fully-developed Power-Law flow through a syringe + straight
needle, with orientation-aware hydrostatic correction.

## Signature

```matlab
results = bioprinting_algorithm_3(Rs, Rn, Ln, Ls, Vp, K, n, rho, varargin)
```

### Required positional arguments

| Arg | Meaning | Units |
|---|---|---|
| `Rs` | syringe inner radius | m |
| `Rn` | needle inner radius | m |
| `Ln` | needle length | m |
| `Ls` | syringe length | m |
| `Vp` | piston velocity (scalar) | m/s |
| `K`  | Power-Law consistency | Pa·sⁿ |
| `n`  | Power-Law index | dimensionless |
| `rho`| fluid density | kg/m³ |

### Optional name-value pairs

- `'PlotResults'` (default true)
- `'NumPoints'` (radial discretisation, default 200)
- `'Name'` (sample label)
- `'SaveData'`, `'SaveFigures'` (default true)
- `'Orientation'`: `'horizontal'` | `'upward'` | `'downward'` (default
  `'horizontal'`). Controls sign of hydrostatic head.
- `'OutputFolder'`
- `'IncludeHydrostatic'` (default true)

## What it returns

Struct with: `tau_w`, `gamma_w`, `dP_needle`, `dP_syringe`, `dP_total`,
radial profiles `r`, `v(r)`, `gdot(r)`, plus exported figures and a TXT report
if enabled.

## When to use this instead of `v3`

- You need a **radial velocity profile** for visualisation.
- You want **orientation-dependent** hydrostatic pressure (upward vs downward
  printing).
- You want a verbose per-sample TXT report rather than a CSV row.

## Assumptions (apply to all PL solvers)

- Incompressible, steady, laminar, fully-developed flow
- No wall slip, no entrance/contraction losses
- No yield stress (PL only — for HB use a separate solver)
- No viscoelasticity or thixotropy
- Exit pressure = 1 atm (101 325 Pa)

Driver: `run_solver_improved.m` (sample loop wrapper).
