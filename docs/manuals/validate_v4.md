# `validate_v4.m` — Regression Check for the v4 Superset Solver

> **STATUS:** `active` — a maintainer / reviewer regression test. It is **not**
> part of the student analysis workflow; you do not need to run it to use the
> pipeline.

## Purpose

Confirms that the unified superset solver
[`bioprinting_algorithm_v4.m`](bioprinting_algorithm_v4.md) reproduces the
results of the three solvers it replaced. Anyone who wants to trust v4 — or
verify it after editing — can run this and read the relative-error report.

## What it does

Using a single fixed case (C15 / 21G at `Vp = 0.01 mm/s`, horizontal so the
hydrostatic term drops out and the legacy totals are directly comparable) it:

1. runs the case through `bioprinting_algorithm_v4`;
2. runs the **Power-Law** legacy solver `bioprinting_algorithm_3`;
3. runs the **Cross** legacy solver `bioprinting_algorithm_cross_v2`
   (`Rn_in = Rn_out`, i.e. degenerate straight needle);
4. prints a field-by-field comparison with a relative-error flag:
   `OK` (< 1e-3), `ok~` (< 2e-2), or `FAIL`;
5. checks the slicer parameters against the inline v3 formulas
   (`v_print`, `k_flow`, `w_line`);
6. confirms the combined `_data.txt` contains both the Power-Law and Cross
   sections plus the needle and syringe radial profiles.

Expected result: Power-Law fields match **exactly** (`OK`); Cross fields match
within solver tolerance (`OK`/`ok~`).

## Dependencies

```
02_MATLAB/
├── bioprinting_algorithm_v4.m
└── archive/scripts/
    ├── bioprinting_algorithm_3.m         ← legacy Power-Law solver
    └── bioprinting_algorithm_cross_v2.m  ← legacy Cross solver
```

`validate_v4.m` adds `archive/scripts/` to the path itself. This is the only
reason those two deprecated solvers are shipped in this export — they are the
ground truth the regression test compares against. Do not call them directly
for new analyses (see `SCRIPT_REGISTRY.md`).

## How to run

```matlab
cd '<PROJECT ROOT>'
addpath('Export/02_MATLAB')
validate_v4
```

Output is printed to the command window; nothing persistent is written (it
uses a temporary output folder).
