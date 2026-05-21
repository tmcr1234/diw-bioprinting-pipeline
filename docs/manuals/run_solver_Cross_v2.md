# `run_solver_Cross_v2.m` — Cross Model Driver (tapered nozzle)

> **STATUS:** `active`. For straight cylinders, prefer `run_solver_v3.m`.

## Purpose

Sample-loop driver for the Cross-model tapered-nozzle solver
`bioprinting_algorithm_cross_v2.m`. Single-model, single-geometry sweep over a
list of inks.

## Inputs

Edit the sample struct at the top:

```matlab
samples = struct('Name',  {...}, ...
                 'eta0',  {...}, 'etainf', {...}, ...
                 'lambda',{...}, 'm',      {...}, ...
                 'Vp',    0.025, 'rho',    {...});
```

Plus tapered-nozzle geometry: `Rs`, `Rn_in`, `Rn_out`, `Ln`, `Ls`.

## Outputs

`Cross-Vp<XX>/` folder containing per-sample TXT + PNG.

## When to use

Tapered-nozzle setups, when you want only Cross-model results (no PL
comparison). If you want both PL and Cross under a straight cylinder, use
`run_solver_v3.m`.
