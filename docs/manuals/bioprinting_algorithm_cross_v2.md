# `bioprinting_algorithm_cross_v2.m` — Cross Model Solver, Tapered Nozzle

> **STATUS:** `active`. Use this when your nozzle is **tapered** (different
> inlet and outlet radii). For straight cylinders, prefer `_v3` (faster,
> handles PL + Cross in a single call).

## Purpose

Steady-state, fully-developed flow of a Cross-model fluid through a syringe
followed by a **conically tapered nozzle** (`Rn_in` → `Rn_out`). Solves the
radial momentum balance numerically along the nozzle axis.

## Signature

```matlab
results = bioprinting_algorithm_cross_v2(Rs, Rn_in, Rn_out, Ln, Ls, Vp, ...
                                          eta0, etainf, lambda, m, rho, ...)
```

### Required arguments

| Arg | Meaning | Units |
|---|---|---|
| `Rs` | syringe inner radius | m |
| `Rn_in` | nozzle **inlet** radius | m |
| `Rn_out` | nozzle **outlet** radius | m |
| `Ln` | nozzle (taper) length | m |
| `Ls` | syringe length | m |
| `Vp` | piston velocity | m/s |
| `eta0` | zero-shear viscosity | Pa·s |
| `etainf` | infinite-shear viscosity | Pa·s |
| `lambda` | Cross time constant | s |
| `m` | Cross exponent | dimensionless |
| `rho` | density | kg/m³ |

### Optional name-value pairs

Same set as `_3.m`: `'PlotResults'`, `'NumPoints'`, `'Name'`, `'SaveData'`,
`'SaveFigures'`, `'Orientation'`, `'OutputFolder'`, `'IncludeHydrostatic'`.

## Method

At each axial slice `z`, the local radius `R(z)` is interpolated linearly
between `Rn_in` and `Rn_out`. The wall shear stress `τ_w(z)` is found by
`fzero` on the Weissenberg–Rabinowitsch integral with the local Cross viscosity
law. The total pressure drop is `trapz` of `dP/dz` along the nozzle.

## Driver

`run_solver_Cross_v2.m` — sample loop, default Vp sweep, exports plots and
TXT reports per sample.

## When to use the cylindrical `_v3` instead

If `Rn_in == Rn_out` (straight cylinder), use `bioprinting_algorithm_v3` — it
has an analytical PL branch and a faster Cross branch (single τ_w inversion
instead of one per axial slice).
