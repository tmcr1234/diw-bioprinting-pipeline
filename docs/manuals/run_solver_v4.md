# `run_solver_v4.m` — One-Shot Master Driver

> **STATUS:** `draft` (smoke-test pending). Once validated against `run_solver_v3.m`
> on a known-good ink set, this will be promoted to `active` and become the
> primary MATLAB entry point.
>
> **v4.1 iteration (2026-05-21):** added per-ink ramp optionality and the
> per-subset `all_samples_velocity_<orientation>.png` overlay. See the
> "Per-subset cross-ink overlay" and "Inks with only one ramp" sections below.
>
> **v4.2 iteration (2026-05-21):** headless-figures mode by default (no GUI
> windows pop up during the run; restored on exit via `onCleanup` even on
> Ctrl-C / error) and simplified `Vp` folder-name convention
> (`Vp10mmps` / `Vp0p01mmps` instead of `Vp_10_mmps` / `Vp_0p01_mmps`). See
> the "Configuration" and "Folder-name convention" sections.

## Purpose

Single execution that reproduces, in one folder tree, every output that
previously required manual reruns of three separate drivers:

- the legacy per-sample **radial-profile** workflow from `run_solver_improved.m`
  (Power-Law) and `run_solver_Cross_v2.m` (Cross) — 4 PNGs + `_data.txt` per
  sample, per piston velocity, and
- the legacy **cross-ink overlay** (`all_samples_velocity_<orientation>.png`)
  emitted per `{ramp, Vp, needle, model}` subset, and
- the **slicer-lookup** workflow from `run_solver_v3.m` (PL + Cross CSV +
  4-panel summary PNG), and
- a new **long-format master summary** with one row per
  `{ramp × ink × needle × Vp × model}`.

## What it produces

Inside `output_v4/`:

| Path | Contents |
|---|---|
| `<Ramp>/Vp<v>mmps/PL_<needle>/<ink> - *.png` + `<ink>_data.txt` | Per-ink legacy 4-panel PNGs (velocity / shear-rate / shear-stress / system overview) + raw radial-profile text dump, for the Power-Law solver |
| `<Ramp>/Vp<v>mmps/PL_<needle>/all_samples_velocity_<orientation>.png` | **Cross-ink overlay** for that subset — needle + syringe radial velocity profiles, one curve per ink that has parameters for this ramp |
| `<Ramp>/Vp<v>mmps/Cross_<needle>/...` | Same set of outputs for the Cross solver (called with `Rn_in = Rn_out` so it degenerates to a straight cylinder, matching the PL geometry) |
| `slicer_lookup/<Ramp>/slicer_lookup_<ink>_<needle>.csv` | Per-(ink × needle) slicer CSV, PL & Cross side-by-side, swept over the full `Vp` vector |
| `slicer_lookup/<Ramp>/plots_<ink>_<needle>.png` | 4-panel summary (ΔP, wall shear rate, head speed, k_flow) |
| `master_summary_v4.csv` | **Long format**: one row per `{ramp, ink, needle, Vp, model}` |

Master summary columns: `ramp, ink, needle, Vp_mm_s, model, Q_mm3_s,
u_avg_needle_mm_s, u_max_needle_mm_s, tau_wall_needle_Pa,
gamma_w_true_needle_invs, dP_total_kPa, Re_needle`.

**PL vs Cross row asymmetry (read before filtering by column):**

| Column | PL row source | Cross row source |
|---|---|---|
| `u_avg_needle_mm_s` | `results.flow.u_avg_needle` (constant cross-section) | `results.flow.u_avg_needle_exit` (Cross is tapered-aware; exit value) |
| `u_max_needle_mm_s` | `results.flow.u_max_needle` (centreline value from analytical PL profile) | **`NaN`** — Cross algorithm doesn't expose a single `u_max` field |
| `gamma_w_true_needle_invs` | Rabinowitsch-corrected true wall shear rate (`results.shear.gamma_dot_true_needle`) | `results.shear.gamma_dot_wall_needle` — max of the discrete radial-profile shear-rate solution. Equivalent to the wall value but no Rabinowitsch correction is applied (Cross integrates the profile directly) |
| `Re_needle` | Generalised power-law Reynolds | Newtonian-form Reynolds with effective `η_wall = τ_wall / γ̇_wall` |

When comparing PL vs Cross rows side-by-side, treat `gamma_w_true_needle_invs` and `Re_needle` as "wall metrics with different definitions" rather than identical quantities. The `Q_mm3_s`, `tau_wall_needle_Pa`, and `dP_total_kPa` columns are directly comparable.

## Inputs

**Edit the `inks(...)` struct at the top of the script**. Each ink carries:

| Field | From | Notes |
|---|---|---|
| `name` | — | Canonical sample name (e.g. `C15`, `C15 Gira 5.5`, `Bozzano's Hair Gel`) |
| `rho` | — | Density (kg/m³) — used by the legacy algorithms for generalised Reynolds |
| `Rrec_pct` | `extract_recovery_v2.py` (or `Recovery_v1.py`) at deposition shear (typ. 150 s⁻¹) | Used by the v3 slicer layer for `k_flow` |
| `Ramp1` | `Fit_Muitos_Modelos_v5.py` (or v4) — **no-pre-shear** fit | Struct of `K_PL, n_PL, eta0, etaInf, lambda, m_Cross` |
| `Ramp2` | Same script — **post-pre-shear** fit (200 s⁻¹ × 300 s) | Same fields; thixotropy/rheopexy comparison vs Ramp 1. **Optional** — see below |

### Inks with only one ramp

Not every ink has both a no-pre-shear and a post-pre-shear flow curve. To
declare an ink that has **only Ramp 1** (or only Ramp 2), simply omit the
missing field from the struct, or set it to an empty array:

```matlab
inks(4).name     = 'NewInk_X';
inks(4).rho      = 980;
inks(4).Rrec_pct = 70;
inks(4).Ramp1    = struct('K_PL', ..., 'n_PL', ..., 'eta0', ..., ...);
% (no Ramp2 line — this ink has no post-pre-shear data)
```

The driver checks `has_ramp(ink, ramp)` before every per-ink call (legacy
PL, legacy Cross, and slicer layer). Missing slices are skipped silently:

- the per-sample 4-panel PNGs and `_data.txt` are not produced for that
  `{ink, ramp}` combination,
- the ink is omitted from the cross-ink overlay for that ramp,
- no row is written to `master_summary_v4.csv` for that combination,
- the slicer CSV layer skips the ink.

If **no** inks declare a particular ramp, that whole `Ramp1/` or `Ramp2/`
subtree is simply never created. The script will not produce empty
folders.

### Per-subset cross-ink overlay

After the inner ink loop for each `{ramp, Vp, needle, model}` slice closes,
the driver emits one overlay PNG:

```
<Ramp>/Vp_<v>_mmps/PL_<needle>/all_samples_velocity_<orientation>.png
<Ramp>/Vp_<v>_mmps/Cross_<needle>/all_samples_velocity_<orientation>.png
```

Each overlay is a 1×2 figure (needle + syringe velocity profiles), one
coloured curve per ink that had parameters for this slice. This mirrors
the `create_multi_sample_plots` nested function from
`run_solver_improved.m` / `run_solver_Cross_v2.m`, with three deliberate
differences:

- **deterministic filename** (no timestamp suffix) — overwrites on rerun
  so the folder tree stays clean,
- **invalid-profile filter** — if the Cross solver returns NaN for an ink
  at very high Vp (Q above the asymptotic Q_max), that curve is dropped
  from the overlay instead of aborting the whole figure,
- **render at `cfg.fig_dpi`** (default 300) — the legacy default of 1200
  produced ~5–8 MB per overlay; 300 DPI is print-quality at a fraction of
  the size.

**Geometries** (lines 76–82): `Rs`, `R_n`, `L_n`, `Ls`. Defaults: BD 10 mL
syringe + 21G/22G blunt needles. Edit if your hardware differs.

**Configuration block** (lines 26–35):

| Knob | Default | Purpose |
|---|---|---|
| `cfg.output_root` | `'output_v4'` | Where everything is written |
| `cfg.orientation` | `'downward'` | `horizontal` / `upward` / `downward` (affects hydrostatic term) |
| `cfg.include_hydrostat` | `true` | Toggle ρgh contribution to ΔP |
| `cfg.save_legacy_plots` | `true` | Toggle the per-sample 4 PNGs + `_data.txt` (set `false` for fast production runs) |
| `cfg.save_slicer` | `true` | Toggle the v3-style slicer CSVs + summary PNGs |
| `cfg.ramps_to_run` | `{'Ramp1','Ramp2'}` | Run one or both ramps |
| `cfg.legacy_Vp_mode` | `false` | If `true`, uses the old `[3..25] mm/s` grid — that grid is on a mechanically unreachable operating point for the CMC/NE rig and is provided **only** to audit-reproduce pre-2026-05 outputs. New inks with a different printer envelope may legitimately want this range; in that case prefer to set `cfg.Vp_mm_s` directly |
| `cfg.show_figures` | `false` | When `false` (v4.2 default), every figure created during the run is built headless and saved straight to PNG — no GUI windows pop up. Set `true` only for interactive debugging |
| `cfg.fig_dpi` | `300` | **Currently parsed but not threaded** into the legacy algorithms (which hard-code 1200 DPI internally). See "Known TODO" below |

## Folder-name convention

Each `Vp` value becomes a folder tag of the form `Vp<value>mmps`, with any
decimal point in the value replaced by the letter `p`:

| `cfg.Vp_mm_s(k)` | folder tag |
|---|---|
| `10` | `Vp10mmps` |
| `3` | `Vp3mmps` |
| `0.01` | `Vp0p01mmps` |
| `0.003` | `Vp0p003mmps` |

**The number in the folder name is the value of `Vp` in mm/s, not a
truncated form of it.** `Vp0p01mmps` is Vp = 0.01 mm/s, not Vp = 10 mm/s.
The bench-realistic default sweep (`[0.003 … 0.04] mm/s`) and the legacy
sweep (`[3 … 25] mm/s`) differ by roughly three orders of magnitude — see
the project memory `v3_vs_old_simulation_audit` for why.

## How to run

```matlab
cd '<PROJECT ROOT>'
addpath('Export/02_MATLAB')
run_solver_v4
```

## Runtime estimate

Defaults: 2 ramps × 8 Vps × 3 inks × 2 needles × 2 models = **192 legacy calls**
+ 12 slicer calls. Cross is the bottleneck (`fzero` over the
Weissenberg–Rabinowitsch integral). Expect ~3–6 minutes on a typical laptop.

Quick smoke-test: set `cfg.ramps_to_run = {'Ramp1'}` and shorten `cfg.Vp_mm_s`
to e.g. `[0.01]` — this drops the total to ~12 calls (<30 s).

## Differences vs `run_solver_v3.m`

| Aspect | `run_solver_v3.m` | `run_solver_v4.m` |
|---|---|---|
| Per-sample radial-profile PNGs + `_data.txt` | ✗ | ✓ |
| Cross-ink overlay (`all_samples_velocity_*.png`) per subset | ✗ | ✓ (v4.1) |
| Slicer CSV + summary PNG | ✓ | ✓ |
| Master summary format | Wide (PL & Cross in one row) | Long (one row per model) |
| Ramp coverage in one run | Manual edit between runs | Both in one execution |
| Inks with only one ramp | All-or-nothing per script edit | Per-ink `has_ramp` guard — drop the field, the ink is skipped where it has no data |
| Cross algorithm called | `bioprinting_algorithm_v3.m` (vectorised) | `bioprinting_algorithm_cross_v2.m` (per-Vp, with full profile outputs) plus `bioprinting_algorithm_v3.m` for the slicer CSV |
| Figure RAM behaviour | N/A (single PNG) | `close all` per inner iteration to prevent accumulation |
| Ink naming convention | Underscored aliases (`C15_SF_5_5`) | Canonical names (`C15 Gira 5.5`) — both refer to the same material |

## Dependencies

MATLAB R2020a+. Calls the existing solver functions unchanged:

- `bioprinting_algorithm_3.m` (PL, straight cylinder)
- `bioprinting_algorithm_cross_v2.m` (Cross, called with `Rn_in = Rn_out`)
- `bioprinting_algorithm_v3.m` (PL + Cross slicer layer)

All four must be on the path.

## Caveats

- Same flow assumptions as the underlying algorithms: fully developed laminar
  flow, no entrance/exit losses, no wall slip, no yield stress in the solver
  (yield handled separately by `extract_hmax_v3.py`).
- The `cfg.fig_dpi` knob is currently **not** threaded into the legacy
  algorithms — they call `exportgraphics(..., 'Resolution', 1200)` internally,
  which produces large PNGs (~5–8 MB each). Promoting `cfg.fig_dpi` to actually
  control output resolution requires adding a `Resolution` name-value to
  `bioprinting_algorithm_3.m` and `bioprinting_algorithm_cross_v2.m`. Left as
  an explicit follow-up to avoid touching algorithms that have audited output
  history.
- `cfg.legacy_Vp_mode = true` reproduces the pre-2026-05 Vp grid for audit
  only. Do **not** report results at those Vps as if they were achievable on
  the current rig — they are not. See `kflow_inversion_convention.md` and
  related project memories.

## When to promote to `active`

Promote `run_solver_v4.m` (and demote `run_solver_v3.m` to `deprecated`) once
all of the following are true:

1. A smoke-test run with `cfg.ramps_to_run = {'Ramp1'}` and the v3 default Vp
   grid produces `slicer_lookup_*` CSVs and `plots_*` PNGs that are
   **numerically identical** to the `output_v3/` ones, modulo file-naming.
2. The legacy 4-panel PNGs for at least one (ink, needle, Vp) match the
   per-sample outputs from a stand-alone `run_solver_improved.m` run at the
   same Vp.
3. The `master_summary_v4.csv` opens cleanly in pandas/Excel and pivots over
   `(ink, needle, Vp, model)` without column-name collisions.

Update `SCRIPT_REGISTRY.md` and the manuals index in
`docs/manuals/README.md` when promoting.
