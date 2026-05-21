# `extract_SAOS_frequency_v2.py` — Frequency Sweep Value Extraction

> **STATUS:** `active`. The frequency-sweep portion of the legacy
> `extract_SAOS_values.py` ported to read Anton Paar Rheocompass CSVs via
> `antpar_io`. Companion to `extract_SAOS_amplitude_v2.py`. Together they
> replace v1 entirely — `extract_SAOS_values.py` is now legacy.

## Purpose

For each frequency-sweep CSV in a folder, compute:

| Quantity | Symbol | Meaning |
|---|---|---|
| Moduli at SAOS targets | `G'(ω)`, `G''(ω)`, `tan δ(ω)` | picked at ω ≈ 1, 10, 100 rad/s (actual data point closest in log-ω) |
| Low-ω power-law | `G' = G₀ · ωᵝ` | fit over the lowest 6 valid points |
| Extrapolated long-dwell modulus | `G'(ω = 0.01 rad/s)` | from the power-law fit |

The β exponent itself is a useful **liquid-vs-solid descriptor**:
β → 0 means solid-like (G' independent of ω); β → 1 means viscous-like
(G' linear in ω, Maxwell-fluid regime).

## Inputs

| Item | Source | Required |
|---|---|---|
| Raw frequency-sweep CSVs | folder of `*.csv` | yes |
| Sample include list | edit `SAMPLE_INCLUDE` to filter | optional |

**Expected file format** — Anton Paar Rheocompass raw export. Required
columns: point, **Frequencia Angular** `[rad/s]`, **Módulo de
Armazenamento** `[Pa]` (G'), **Módulo de Perda** `[Pa]` (G''), **Fator
de Perda** `[1]` (tan δ), Re(η*) `[mPa·s]`, Im(η*) `[mPa·s]`, |G*| `[Pa]`,
Deformação `[%]`, Tensão `[Pa]`, Torque `[mN·m]`, Status. The parser
converts complex viscosity to Pa·s automatically.

## Outputs

The script writes **two pairs** of outputs on every run, all to
`./Analises/Python/Results/`:

| File | Filter | Contents |
|---|---|---|
| `SAOS_frequency_v2.txt` / `.csv` | M- and ME- excluded (default) | The trusted scientific values — what `extract_hmax_v3` consumes |
| `SAOS_frequency_v2_allpoints.txt` / `.csv` | No Status filter | Same algorithm with every measured point regardless of M-/ME-/taD |

Set `ALSO_WRITE_ALLPOINTS = False` at the top of the script to skip the
second pass. The CSVs differ by the `drop_status_flags` column so two
rows from different files can never be confused.

**Why two outputs?** Most samples produce nearly identical results both
ways — but for samples near the rheometer sensitivity floor the
difference is dramatic. Sample A in the Carla dataset, for instance,
is "unrecoverable" in the filtered pass (the v2 correctly refuses to
fit when only 0 points survive) but produces a β = 7.1, R² = 0.55 fit
when all points are used (visibly garbage but numerically present).
Having both lets you defend the choice in a manuscript.

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/extract_SAOS_frequency_v2.py
```

Edit `INPUT_FOLDER` (line ~58) to point at your Frequency Sweep folder.

## Dependencies

`numpy`, plus the local `antpar_io.py` module.

## Status filtering

At low ω the rheometer hits sub-resolution torque (`M-`) and modulus
estimate errors (`ME-`). Those points have G' on the order of `1e-5 Pa`
— numerical noise. They destroy the low-ω power-law fit if kept.

This script **drops rows whose Status contains `M-` or `ME-` by default**
(the filtered pair) **AND emits an `_allpoints` pair** with no Status
filter for diagnostic comparison. The `taD` flag (tan-delta drift,
usually high-ω inertial corrections) is kept in both passes because the
G' value at those frequencies is still physically meaningful — it's the
*tan δ* that's drifting, not G'.

Two knobs:

```python
DROP_STATUS_FLAGS = ("M-", "ME-")    # what the filtered pair uses
ALSO_WRITE_ALLPOINTS = True          # write the second (no-filter) pair

# Other useful overrides — applied to the FILTERED pair:
DROP_STATUS_FLAGS = ("M-",)          # keep ME- points
DROP_STATUS_FLAGS = ()               # filtered pair becomes the same as
                                     # the all-points pair (don't do this)
```

## ω-target picking strategy

For each target ω ∈ {1, 10, 100} rad/s, the script picks the **valid
data point closest in log-ω**. If a target falls outside the *valid*
ω range by more than a factor of 2, the pick is marked `out-of-range`
in the report. This prevents extrapolating manuscript-table values from
the noisy region of the sweep.

For Carla's data, the target ω = 100 rad/s maps to the actual data
point at ω ≈ 94 rad/s — close enough to be reported. Target ω = 1
sometimes lands at ω ≈ 1.27 or 1.55 when the lowest valid point isn't
exactly 1 rad/s.

## Edge cases handled

| Situation | Behaviour |
|---|---|
| Sample with all 41 points flagged | Reports `"Only 0 usable points after Status filtering — cannot extract."` and still writes a row in the report |
| Fewer than 6 valid points for low-ω fit | Fit uses whatever's available, reports `n_used` so you can spot the warning |
| Target ω outside the valid sweep | Reported as `out-of-range`; not extrapolated into the noise |
| Low-ω fit R² below ~0.95 | Reported numerically — visual inspection recommended (planned: per-sample fit plots in extract_hmax_v3) |

## Caveats

- **The β exponent is best interpreted ALONGSIDE the LVR results.** A
  sample with high G'_LVR and low β is a strong gel; a sample with high
  G'_LVR but high β is more like a structured liquid that happens to be
  stiff at the LVR strain.
- **G'(ω=0.01) is an extrapolation**, not a measurement. The lowest
  measured ω is typically 0.1 rad/s — extrapolating one decade further
  down is reasonable when R² > 0.99, less so when the fit struggles.
  Use it as a "long-dwell modulus estimate" not a hard number.
- **The complex viscosity columns** (`eta_complex_real_Pa_s`,
  `eta_complex_imag_Pa_s`) are parsed and stored in the DataFrame but
  not used by this extractor — they're available for downstream code
  that wants them (e.g., Cox–Merz comparisons against the flow curve).

## Hand-off

`extract_hmax_v3.py` reads `SAOS_frequency_v2.csv` directly (the
**filtered** variant — not the all-points one). If you want to feed
the all-points results into h_max for a sensitivity check, point
`extract_hmax_v3.FREQUENCY_CSV` at `SAOS_frequency_v2_allpoints.csv`
manually before running it.

## Versioning note

This is **v2** of the frequency-sweep extraction — supersedes the
`extract_frequency_sweep()` function in the legacy `extract_SAOS_values.py`.
The amplitude-sweep half of v1 is also superseded (by
`extract_SAOS_amplitude_v2.py`). At this point **the entire legacy
`extract_SAOS_values.py` is superseded** — kept active only because
`extract_hmax_v2.py` still imports it. Once `extract_hmax_v3` lands,
`extract_SAOS_values.py` and `extract_hmax_v2.py` can both be moved to
deprecated.
