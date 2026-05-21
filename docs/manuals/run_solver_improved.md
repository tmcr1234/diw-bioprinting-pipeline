# `run_solver_improved.m` — Power-Law Solver Driver (per-sample)

> **STATUS:** `active` (legacy). For new work, prefer `run_solver_v3.m` —
> it produces the slicer-lookup CSVs the SOP refers to.

## Purpose

Loops over a list of samples and calls `bioprinting_algorithm_3.m` for each.
Supports both **unsheared** (`Ramp 1`-derived parameters) and **pre-sheared**
(`Visco_Artur`-derived parameters) datasets in a single run.

## Inputs

Edit two struct arrays at the top of the script:

```matlab
samples         = struct('Name',{...}, 'K',{...}, 'n',{...}, 'Vp', 0.025, 'rho',{...});
samples_sheared = struct('Name',{...}, 'K',{...}, 'n',{...}, 'Vp', 0.025, 'rho',{...});
```

Plus shared geometry (lines 24–28): `Rs`, `Rn`, `Ln`, `Ls`.

## Outputs

Per-sample TXT reports, PNG figures, and radial-profile data, dumped to
`output_folder` (default `Vp25-Sheared/` in the script's cwd).

## How to run

```matlab
cd '<PROJECT ROOT>'
addpath('Export/02_MATLAB')
run_solver_improved
```

## When to use this instead of `v3`

- You only care about **Power-Law** results (not Cross).
- You want **radial velocity profiles** per sample.
- You want to compare unsheared vs pre-sheared (`Artur`) parameter sets
  side-by-side.

For everyday slicer-parameter lookup tables, use `run_solver_v3.m`.
