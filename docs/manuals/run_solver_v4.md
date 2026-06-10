# `run_solver_v4.m` — One-Shot Master Driver

> **STATUS:** `active` — the **sole** MATLAB extrusion driver. It replaces
> `run_solver_v3.m`, `run_solver_improved.m`, and `run_solver_Cross_v2.m`,
> all of which are now deprecated (see `SCRIPT_REGISTRY.md`).
>
> **2026-06-01 rewrite:** the driver no longer orchestrates three separate
> solvers. It makes a single call to the unified superset solver
> [`bioprinting_algorithm_v4.m`](bioprinting_algorithm_v4.md) per
> `(ramp × ink × needle)`, which computes the Power-Law **and** Cross results,
> the radial profiles, and the slicer parameters together.

## Purpose

One execution produces, in a single `output_v4/` tree, every output that
previously required separate runs of three drivers: per-operating-point radial
profiles (Power-Law + Cross), the slicer-lookup CSV, an across-`Vp` trend
summary, and a long-format master summary spanning all inks, needles, ramps,
velocities, and both constitutive models.

## What it produces

Inside `output_v4/`:

| Path | Contents |
|---|---|
| `<Ramp>/Vp<v>mmps/<ink>_<needle>_data.txt` | Combined human-readable dump — Power-Law section + Cross section (legacy `_data.txt` format), including the radial profiles for syringe and needle |
| `<Ramp>/Vp<v>mmps/<ink>_<needle>_profiles.png` | Radial profiles `u(r)`, `γ̇(r)`, `τ(r)` (syringe + needle) for both models |
| `<Ramp>/Vp<v>mmps/<ink>_<needle>_system.png` | System overview — pressure decomposition, generalised Reynolds, viscosity across the needle exit |
| `<Ramp>/slicer_lookup_<ink>_<needle>.csv` | Across-`Vp` slicer table, PL & Cross side-by-side |
| `<Ramp>/<ink>_<needle>_summary.png` | Across-`Vp` trends (ΔP, wall shear rate, head speed, k_flow) |
| `master_summary_v4.csv` | **Long format** — one row per `{ramp, ink, needle, Vp, model}` |

Master summary columns: `ramp, ink, needle, Vp_mm_s, model, Q_mm3_s,
u_avg_needle_mm_s, u_max_needle_mm_s, tau_wall_needle_Pa, gamma_w_needle_invs,
dP_total_kPa, Re_needle, v_print_mm_s, k_flow`.

Two rows are written per operating point — one `PowerLaw`, one `Cross`. The
`u_max_needle_mm_s` column is the analytical centreline for the PL row; the
Cross row reports its profile maximum. Treat `gamma_w_needle_invs` and
`Re_needle` as "wall metrics with model-specific definitions" (PL applies the
Rabinowitsch correction; Cross integrates the discrete profile directly).
`Q_mm3_s`, `tau_wall_needle_Pa`, and `dP_total_kPa` are directly comparable
across the two rows.

### Folder-name convention

Each `Vp` value becomes a folder tag `Vp<value>mmps`, with any decimal point
replaced by the letter `p`:

| `cfg.Vp_mm_s(k)` | folder tag |
|---|---|
| `0.01` | `Vp0p01mmps` |
| `0.003` | `Vp0p003mmps` |
| `10` | `Vp10mmps` |

**The number is the value of `Vp` in mm/s, not a truncated form.**
`Vp0p01mmps` is Vp = 0.01 mm/s, not 10 mm/s.

## Inputs

**Edit the `inks(...)` struct at the top of the script** (lines ~62–78). Each
ink carries:

| Field | From | Notes |
|---|---|---|
| `name` | — | Canonical sample name (e.g. `C15`, `C15 Gira 5.5`, `Bozzano's Hair Gel`) |
| `rho` | — | Density (kg/m³) — used for the generalised Reynolds number |
| `Rrec_pct` | `extract_recovery_v2.py` (or `Recovery_v1.py`) at the deposition wall shear rate (~200 s⁻¹) | Used by the slicer layer for `k_flow` |
| `Ramp1` | `Fit_Muitos_Modelos_v5.py` (or v4) — **no-pre-shear** fit | Struct of `K_PL, n_PL, eta0, etaInf, lambda, m_Cross` |
| `Ramp2` | Same script — **post-pre-shear** fit (200 s⁻¹ × 300 s) | Same fields; thixotropy/rheopexy comparison vs Ramp 1. **Optional** |

### Inks with only one ramp

Omit (or set to `[]`) the missing field. The driver's `has_ramp(ink, ramp)`
guard skips any `{ink, ramp}` combination that has no data — no per-Vp
folder, no slicer CSV, no master-summary rows for it. If no ink declares a
ramp, that `Ramp1/`/`Ramp2/` subtree is never created.

## Configuration block (lines ~24–49)

| Knob | Default | Purpose |
|---|---|---|
| `cfg.output_root` | `'output_v4'` | Where everything is written |
| `cfg.orientation` | `'downward'` | `horizontal` / `upward` / `downward` (sign of the ρgh term) |
| `cfg.include_hydrostat` | `true` | Toggle the ρgh contribution to ΔP |
| `cfg.save_data` | `true` | Toggle the per-Vp combined `_data.txt` |
| `cfg.save_figures` | `true` | Toggle per-Vp + summary PNGs |
| `cfg.save_csv` | `true` | Toggle the across-Vp slicer CSV |
| `cfg.num_points` | `200` | Radial discretisation |
| `cfg.ramps_to_run` | `{'Ramp1','Ramp2'}` | Run one or both ramps |
| `cfg.show_figures` | `false` | Headless figure mode (restored on exit via `onCleanup`, even on error / Ctrl-C). Set `true` only for interactive debugging |
| `cfg.legacy_Vp_mode` | `false` | If `true`, uses the old `[3..25] mm/s` grid — that grid sits on a mechanically unreachable operating point for the CMC/NE rig and is provided **only** to audit-reproduce pre-2026-05 outputs |

Default `Vp_mm_s` sweep when `legacy_Vp_mode` is `false`:
`[0.003 0.005 0.007 0.01 0.015 0.02 0.03 0.04]` mm/s, which maps to head
speeds of ~2–30 mm/s via the `Vp·(Rs/Rn)²` amplification.

**Geometries** (lines ~80–85): `geom_21G`, `geom_22G` carry `Rs`, `R_n`,
`L_n`, `Ls`, `label`, `h_factor` (layer height as a fraction of needle ID,
default 0.7). Defaults: BD 10 mL syringe + 21G/22G blunt needles. Edit if your
hardware differs.

## How to run

```matlab
cd '<PROJECT ROOT>'
addpath('Export/02_MATLAB')
run_solver_v4
```

## Runtime estimate

Defaults: 2 ramps × 3 inks × 2 needles = **12 superset calls**, each sweeping
8 `Vp` values for both PL and Cross. Cross is the bottleneck (`fzero` over the
Weissenberg–Rabinowitsch integral). Expect ~3–6 minutes on a typical laptop.

Quick smoke-test: set `cfg.ramps_to_run = {'Ramp1'}` and shorten
`cfg.Vp_mm_s` to `[0.01]` — drops the work to ~6 calls (<30 s).

## Dependencies

MATLAB R2020a+. Calls a **single** solver:

- [`bioprinting_algorithm_v4.m`](bioprinting_algorithm_v4.md) — the unified
  superset solver (Power-Law + Cross + radial profiles + slicer layer).

Both must be on the path (handled automatically by the `addpath` at the top
of the driver). The deprecated solvers under `archive/scripts/` are **not**
needed to run `run_solver_v4` — they are retained only so `validate_v4.m` can
cross-check v4 against the legacy physics.

## Caveats

- Same flow assumptions as the underlying solver: fully developed laminar
  flow, no entrance/exit losses, no wall slip, no yield stress in the solver
  (yield handled separately by `extract_hmax_v3.py`).
- Geometry is **straight cylindrical needle only**. For a conical tip use
  [`bioprinting_algorithm_conical.m`](bioprinting_algorithm_conical.md)
  directly — it is outside v4's scope.
- `cfg.legacy_Vp_mode = true` reproduces the pre-2026-05 Vp grid for audit
  only. Do **not** report results at those Vps as if achievable on the current
  rig.

## Verifying the solver (optional)

`validate_v4.m` re-runs the C15 / 21G case through `bioprinting_algorithm_v4`
and the two legacy solvers (`bioprinting_algorithm_3`,
`bioprinting_algorithm_cross_v2` in `archive/scripts/`) and prints a
field-by-field relative-error report. PL matches the legacy solver exactly;
Cross matches within solver tolerance. See
[`validate_v4.md`](validate_v4.md).
