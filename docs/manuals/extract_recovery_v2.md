# `extract_recovery_v2.py` — 3iTT Structural Recovery

> **STATUS:** `active`. Replaces `Recovery_v1.py` for Anton Paar
> Rheocompass data. The protocol itself also changed — old code consumed
> a steady-shear ramp battery (compare η across multiple shear rates);
> v2 consumes the modern **3-Interval Thixotropy Test** (3iTT-Osc-Rot-Osc)
> which is one shot per sample and gives both modulus and complex-
> viscosity recovery ratios.

## Purpose

For each 3iTT set (three CSV files: `_1.csv`, `_2.csv`, `_3.csv` sharing
the same basename), compute:

| Quantity | Symbol | Definition |
|---|---|---|
| Modulus recovery | `Rec_G' (%)` | mean G'(interval 3 tail) / mean G'(interval 1 tail) × 100 |
| Complex viscosity recovery | `Rec_|η*| (%)` | mean \|η*\|(int 3 tail) / mean \|η*\|(int 1 tail) × 100 |

Both ratios use the **last N points** of each oscillation interval as
"settled" representative values. Errors propagate from per-tail std.

The rotation interval (interval 2) is parsed and its end-viscosity is
reported as a diagnostic only — it does not enter the ratios.

## Protocol the script expects

```
Interval 1:  Oscillation  at low strain (in LVR)   →  30 points, ~2 min
                                                       baseline G', |η*|
                                                                ↓
Interval 2:  Rotation at high shear rate           →  30 points, ~2 min
                                                       γ̇ ≈ 100–500 s⁻¹
                                                                ↓
Interval 3:  Oscillation at low strain (in LVR)    →  150 points, ~10 min
                                                       recovery trajectory
```

Other 3iTT variants (Rot-Rot-Rot, Osc-Rot-Rot, etc.) are **not yet
supported** — the script errors out cleanly if interval 1 or 3 isn't
oscillation. The underlying parser in `antpar_io` reads them just fine;
only this extractor's recovery math is hard-wired to Osc-…-Osc.

## Inputs

| Item | Source | Required |
|---|---|---|
| Recovery folder containing `_1.csv` / `_2.csv` / `_3.csv` trios | one folder per project | yes |
| Sample include list | edit `SAMPLE_INCLUDE` to filter | optional |

The parser groups files by the part of the filename BEFORE `_1.csv` /
`_2.csv` / `_3.csv` and requires all three to be present per set.
Incomplete sets are warned about and skipped.

## Outputs

The script writes **two pairs** of outputs on every run, all to
`./Analises/Python/Results/`:

| File | Filter | Contents |
|---|---|---|
| `Recovery_v2.txt` / `.csv` | WMa excluded (default) | The trusted scientific values |
| `Recovery_v2_allpoints.txt` / `.csv` | No Status filter | Same algorithm with every measured point, regardless of WMa |

Set `ALSO_WRITE_ALLPOINTS = False` at the top of the script to skip the
second pass. The CSVs differ by the `drop_status_flags` column so two
rows from different files can never be confused.

**Why two outputs?** In practice the WMa filter is usually cosmetic for
recovery, because WMa flags the *first* few points of an interval (motor
settling), while the recovery calculation averages the *last* `N_TAIL_POINTS`
points of each interval — which is rarely WMa-flagged. The all-points
variant lets you confirm this empirically: for Carla's 6 sets the
filtered and all-points modulus recovery values are byte-identical for
every healthy sample. The diff only ever shows up at the front-of-interval
tail-of-data overlap, which is rare.

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/extract_recovery_v2.py
```

Edit `INPUT_FOLDER` at the top of the script.

## Dependencies

`numpy`, plus the local `antpar_io.py` module.

## Tail-averaging strategy

Both baseline (interval 1) and recovered (interval 3) values are
computed as the **mean of the last `N_TAIL_POINTS` valid points** of
their interval, defaulting to 5. Rationale:

- Interval 1 is short (~2 min) and the material is still slowly
  structuring throughout it. The END value is the closest thing to an
  "at-rest" baseline.
- Interval 3 is long (~10 min). The recovery trajectory is monotonically
  rising for healthy samples; the tail is the most-recovered state
  measured.

Reported error per ratio is from first-order Gaussian propagation
through the division.

## Interval-3 cutoff (crosslinking-during-test handler)

Some materials *crosslink* while the recovery interval is running
(e.g., alginate ↔ Ca²⁺ phosphate). G' then shoots far past the interval-1
baseline — producing nonsense "super-recovery" values like 1000 %.

Two cutoff knobs:

```python
# Global cutoff — applied to every sample unless overridden.
INTERVAL3_MAX_POINT: int | None = None    # None = use all points

# Per-sample overrides — looked up by sample short-name (the part of
# the Teste: header BEFORE " - ", e.g. "Amostra PC1").
INTERVAL3_MAX_POINT_OVERRIDES: dict[str, int | None] = {
    "Amostra A":   30,
    "Amostra PC1": 20,    # this sample crosslinked earlier
}
```

If a cutoff is set, only interval-3 points whose `point` value is
**≤ cutoff** enter the tail average. The TXT report shows the cutoff
applied per sample (`int3 cut`), the number of points actually used
(`n_int3_used / n_int3_total`), and the time range of the tail —
so you can see whether the cutoff caught the right window.

**Diagnostic signal**: if the std on the recovered values is huge (>>
the std on baseline) the cutoff is probably too generous and the
crosslinking burst is bleeding in. Tighten the override and re-run.

## Status filtering

Rheocompass flags the first few points of an interval with `WMa`
("Wall Maximum adjusting") while the geometry/torque adjusts. Those
points are usually 2–10 % off the settled value. **By default the
script drops rows containing `WMa`** (the filtered pair) **AND emits
an `_allpoints` pair** with every measured point for diagnostic
comparison.

Two knobs:

```python
DROP_STATUS_FLAGS = ("WMa",)        # what the filtered pair uses
ALSO_WRITE_ALLPOINTS = True         # write the second (no-filter) pair

# Other useful overrides — applied to the FILTERED pair:
DROP_STATUS_FLAGS = ()              # filtered == allpoints (skip the dual mode)
DROP_STATUS_FLAGS = ("WMa", "M-")   # also drop sub-resolution torque
```

## Interpretation cheatsheet

| Pattern | What it usually means |
|---|---|
| `Rec_G' ≈ Rec_|η*|` (both 60–100 %) | Healthy thixotropic recovery |
| `Rec_G' < Rec_|η*|` (e.g. 50 vs 80 %) | Partial recovery — viscous components return faster than elastic network (typical) |
| `Rec_G' < 20 %` with low Rec_\|η*\| too | Network destroyed, no real recovery on this timescale |
| Rec_G' > 100 % with huge std | Crosslinking during the recovery interval — **tighten the cutoff** |
| Rec_G' > 100 % with low std | Sample is still structuring (thixotropic build-up beyond original baseline) |

## Caveats

- The disruption shear rate (interval 2 γ̇) is reported but not used in
  the ratios. Different protocols may use different γ̇ — verify it
  matches your expectation by reading the "DISRUPTION" line in the
  report.
- The script assumes the **TAIL** of interval 1 is the "true" baseline.
  For samples that haven't reached steady state in 2 minutes (rare but
  possible for very slow-structuring materials), the baseline will be
  under-estimated and the recovery % will look artificially high.
  Inspect the baseline std — if it's > 5 % of the mean, interval 1
  was too short.
- Duplicate sample sets (e.g. two runs of PC1 with slightly different
  filenames) keep the first one under the short-name key; the second
  is keyed by its full basename. Both appear in the report.

## Hand-off

Recovery is consumed by:

| Downstream | What it pulls | Future |
|---|---|---|
| `run_solver_v4.m` (slicer layer of `bioprinting_algorithm_v4`) | `inks(i).Rrec_pct` — pick one of the two ratios | manual paste from CSV — pick whichever matches your physical interpretation |
| Manuscript tables | All columns of the CSV | direct LaTeX/Markdown rendering |
| `extract_hmax_v3` (pending) | Not used — h_max draws from SAOS, not 3iTT | — |

Decide which recovery to feed the MATLAB solver based on whether
your printability problem is **elasticity-limited** (use Rec_G') or
**bulk-flow-limited** (use Rec_\|η*\|). For DIW post-extrusion sag,
G' is usually the right choice.

## Versioning note

This is **v2** — supersedes `Recovery_v1.py` which is now legacy. The
protocol math is different (3iTT vs steady-shear ramps), so v1 numbers
are NOT comparable to v2 numbers. Do not mix them across the same
manuscript without an explicit conversion / acknowledgment.
