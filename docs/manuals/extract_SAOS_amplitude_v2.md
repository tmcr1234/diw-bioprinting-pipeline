# `extract_SAOS_amplitude_v2.py` — Amplitude Sweep LVR Extraction

> **STATUS:** `active`. The amplitude-sweep portion of the legacy
> `extract_SAOS_values.py` ported to read Anton Paar Rheocompass CSVs via
> `antpar_io`. The frequency-sweep portion still needs porting (see
> migration roadmap).

## Purpose

For each amplitude (strain) sweep CSV in a folder, compute:

| Quantity | Symbol | Meaning |
|---|---|---|
| Storage modulus plateau | `G'_LVR` | mean G' over the linear region |
| Loss modulus plateau | `G''_LVR` | mean G'' over the linear region |
| Loss factor plateau | `tan δ_LVR` | mean of G''/G' over the linear region |
| Yield strain | `γ_LVR` (%) | strain at which G' first drops below **90 % of plateau**, log-interpolated |
| Yield stress estimate | `σ_y_LVR` (Pa) | `G'_LVR × γ_LVR` — the LVR-endpoint yield stress used downstream by `h_max` |

## Inputs

| Item | Source | Required |
|---|---|---|
| Raw amplitude-sweep CSVs | folder of `*.csv` (skips `*_Analysis.csv`, `*_Analise.csv`, PDFs) | yes |
| Sample include list | edit `SAMPLE_INCLUDE` to filter to a subset | optional |

**Expected file format** — Anton Paar Rheocompass raw export (UTF-16 LE,
tab-separated, comma decimal). Required columns: point, **Deformação de
Cisalhamento** `[1]`, Tensão de Cisalhamento `[Pa]`, **Módulo de
Armazenamento** `[Pa]`, **Módulo de Perda** `[Pa]`, **Fator de Perda**
`[1]`, Torque `[µN·m]`, Status. The parser converts fractional strain to
percent automatically and ships both columns.

## Outputs

The script writes **two pairs** of outputs on every run, all to
`./Analises/Python/Results/`:

| File | Filter | Contents |
|---|---|---|
| `SAOS_amplitude_v2.txt` / `.csv` | M- excluded (default) | The trusted scientific values — what manuscript tables quote and what `extract_hmax_v3` consumes |
| `SAOS_amplitude_v2_allpoints.txt` / `.csv` | No Status filter | Same algorithm with every measured point, regardless of M- flag |

Set `ALSO_WRITE_ALLPOINTS = False` at the top of the script to skip the
second pass. The CSVs differ by the `drop_M_minus` column (`True` vs
`False`) so two rows from different files can never be confused.

**Why two outputs?** The all-points variant lets you see at a glance how
much the M- exclusion is moving your numbers. If the LVR values agree
between the two files, your plateau detector landed in a clean region
anyway and the filter is just cosmetic. If they disagree, the filter
actually matters — and you should report the filtered values (they're
the defensible ones).

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/extract_SAOS_amplitude_v2.py
```

Edit `INPUT_FOLDER` (line ~58) to point at your Strain Sweep folder.

## Dependencies

`numpy`, plus the local `antpar_io.py` module (ships in the same folder).
No `xlrd`, no pandas in the analysis layer.

## LVR detection logic

Identical to the legacy v1 algorithm, kept on purpose so values match the
preprint when re-fitted on equivalent Anton Paar data:

1. **Drop invalid points**: NaN, non-positive G' / G''. **By default,
   also drop rows where Status contains `M-`** (Anton Paar's
   sub-resolution-torque flag — see *Status filtering* below).
2. **Primary detector — sliding-window CV**: walk a 5-point window over
   G'(γ); compute the local coefficient of variation; the longest
   contiguous run with CV < 5 % is the plateau.
3. **Fallback** (used only when the CV detector finds < 3 stable points):
   median ± 15 % over the lower half of the sweep.
4. **γ_LVR** = strain where G' first drops below `0.90 × plateau`,
   **log-interpolated** between bracketing points.

## Status filtering

Anton Paar flags points where the motor torque is below the rheometer's
sensitivity with `M-` in the Status column. Those points are noise floor.
The plateau detector misbehaves on them — it can latch onto a fake
plateau in the noise.

This script **excludes M- by default** (the filtered pair) AND **also
emits an `_allpoints` pair** with every measured point (no Status
filter) for diagnostic comparison.

Two knobs:

```python
DROP_M_MINUS = True            # what the filtered pair uses
ALSO_WRITE_ALLPOINTS = True    # write the second (no-filter) pair
```

For the Carla samples the M- filter typically removes 8–23
sub-resolution points per sweep (out of 41) and shifts G'_LVR by < 2 %
when the plateau region is far from the noise floor. For weaker samples
the effect can be larger — eyeball the two CSVs side by side.

## Caveats

- **High γ_LVR values** (30–70 %) for these samples are real, not a bug.
  Pre-gelled / strongly entangled formulations have shallow elastic
  drops — the 10 % G' drop threshold isn't reached until very large
  strains. If you expect γ_LVR ~ 0.1–5 % (typical crosslinked hydrogel
  range) and you don't see it, your sample isn't behaving like a gel.
- **`tan δ > 1` across the LVR** means G'' > G' — viscous-dominant, NOT
  a true gel. The summary report annotates this interpretation per
  sample. The numerical LVR values are still extracted faithfully — the
  *interpretation* changes (the "yield strain" is the entanglement-
  disruption strain, not a network rupture strain).
- **No temperature column** in Rheocompass amplitude-sweep exports. The
  setpoint temperature is in the filename / `Teste:` header. If you need
  T per point for thermal-history analysis, log it via a separate
  temperature ramp test.

## Hand-off

`extract_hmax_v3.py` reads `SAOS_amplitude_v2.csv` directly (the
**filtered** variant — not the all-points one). If you want to feed
the all-points results into h_max for a sensitivity check, point
`extract_hmax_v3.AMPLITUDE_CSV` at `SAOS_amplitude_v2_allpoints.csv`
manually before running it.

## Versioning note

This is **v2** — supersedes the amplitude-sweep half of
`extract_SAOS_values.py` (the `extract_amplitude_sweep` function). The
frequency-sweep half of the legacy module is still active and not
ported yet.
