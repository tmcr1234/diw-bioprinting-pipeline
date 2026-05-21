# `bioprinting_algorithm_conical.m` — Power-Law Solver, Conical Tip

> **STATUS:** `active`. PL counterpart to `bioprinting_algorithm_cross_v2.m`.

## Purpose

Steady-state Power-Law flow through a syringe + a **conical transition tip**.
The transition length is hard-wired to equal the needle length (`Lt = Ln`) so
the geometry has a single length parameter.

## Signature

```matlab
results = bioprinting_algorithm_conical(Rs, Rn, Ln, Ls, Vp, K, n, rho, varargin)
```

### Required arguments

| Arg | Meaning | Units |
|---|---|---|
| `Rs` | syringe inner radius | m |
| `Rn` | needle (exit) radius | m |
| `Ln` | needle length AND transition length | m |
| `Ls` | syringe length | m |
| `Vp` | piston velocity | m/s |
| `K`, `n` | Power-Law parameters | Pa·sⁿ, dimensionless |
| `rho` | density | kg/m³ |

### Optional name-value pairs

Identical set to `bioprinting_algorithm_3.m`.

## When to use

If you print with a **conical tapered tip** (Nordson Optimum style) rather than
a blunt cylindrical needle. The conical contraction adds an extra pressure
drop the straight-cylinder solver misses.

## Caveats

- Assumes the transition is linear (conical, not exponential or hyperbolic).
- Hard-wired `Lt = Ln`. If your transition length differs from your needle
  length, edit the source (single line near the top of the algorithm block).
