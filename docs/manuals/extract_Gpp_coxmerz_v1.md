# `extract_Gpp_coxmerz_v1.py` — Loss-Modulus (G″) Analyses

> **STATUS:** `active` (new, added 2026-06-11). Reads legacy TA-style `.xls`
> via the parser in `extract_SAOS_values.py`. **Not yet ported** to the
> Anton Paar CSV pipeline — see *Caveats* and the migration roadmap in
> `SCRIPT_REGISTRY.md`.

## Purpose

The G′-only SAOS pipeline (`extract_SAOS_values.py` / `extract_hmax_v*.py`)
extracts storage modulus and yield-strain metrics but never uses the loss
modulus G″. This script adds the three G″-based analyses that loss modulus
unlocks:

- **(A) tan δ(ω) + Winter–Chambon gel-point test.** At a critical gel,
  `G' ~ G'' ~ ω^β` and `tan δ = tan(βπ/2)` is *frequency-independent*. The
  script fits `β'` from `G'(low ω)` and `β''` from `G''(low ω)`, reports
  whether G′ and G″ are parallel power laws, the measured `tan δ` and its
  flatness across ω, and the predicted `tan(β'π/2)`. This is the rigorous
  version of the gel-point claim (a flat tan δ is the real test, not just
  `G' > G''`).
- **(B) Cox–Merz comparison.** Steady-shear viscosity `η(γ̇)` (flow curve,
  Ramp 1) vs complex viscosity `|η*(ω)| = |G*|/ω` (frequency sweep). Simple
  liquids obey `η(γ̇) ≈ |η*(ω)|` at `γ̇ = ω`; structured / associating fluids
  deviate (`|η*| > η`), and the deviation scales with at-rest structure.
- **(C) Weak-strain-overshoot / Payne test.** From the amplitude sweep:
  does G″ rise above its LVR plateau before the G′/G″ crossover? A G″
  overshoot is the LAOS type-III fingerprint of an attractive / associating
  network.

These analyses back Figures F8–F10 of the maintainer's manuscript and the
"depletion-dominant network" interpretation.

## Inputs

| Item | Path / format |
|---|---|
| Frequency sweeps | `./Reologia/Frequency/<sample> - Freq CP50 0,1mm.xls` |
| Flow curves (for Cox–Merz) | `./Reologia/Viscosity/<sample> ... .xls` |
| Amplitude sweeps (for overshoot) | `./Reologia/Strain/<sample> - Deformacao CP50 0,1mm.xls` |

File paths are **hard-coded** in the `FREQ_FILES`, `FLOW_FILES`, and amplitude
dicts near the top of the script, and `PROJECT_ROOT` is an absolute path. For a
new ink / new machine, **edit those dicts and `PROJECT_ROOT` before running**
(same convention as `extract_SAOS_values.py` and `Calculo da Pressao de Extrusao.py`).

## Outputs

- `./Analises/Python/Results/SAOS_Gpp_coxmerz_v1.txt` — full text report
- `./Analises/Python/Results/tandelta_coxmerz_per_ink.csv` — machine-readable
- `./Latex/Figs/F8_tandelta.pdf` (+ `.png` preview in `Results/`)
- `./Latex/Figs/F9_coxmerz.pdf` (+ `.png` preview)
- `./Latex/Figs/F10_Gpp_overshoot.pdf` (+ `.png` preview)

## How to run

```bash
cd "<PROJECT ROOT>"        # folder containing Export/ and Reologia/
python Export/01_Python/extract_Gpp_coxmerz_v1.py
```

## Dependencies

`numpy`, `matplotlib`, `xlrd` (via the imported `extract_SAOS_values` parser).
`matplotlib` uses the `Agg` backend (no display needed). No SciPy.

## Caveats

- **Legacy ingestion.** This script reads `.xls` through
  `extract_SAOS_values.read_antpar_xls` / `find_col` / `extract_amplitude_sweep`.
  It does **not** use `antpar_io.py`, so it will not read Anton Paar Rheocompass
  CSVs as-is. To port it, swap those calls for the Anton Paar readers
  (`antpar_io` + the `_v2` SAOS extractors) and source the flow curve from the
  v5 reader for the Cox–Merz panel. Tracked as an open item in `SCRIPT_REGISTRY.md`.
- **Hard-coded absolute `PROJECT_ROOT`.** Edit it (and the file dicts) for any
  project other than the maintainer's, or the script will write figures back
  into the original project tree.
- The Winter–Chambon test is only meaningful when the low-ω fit window actually
  sits in the terminal/gel region; inspect the fitted `β'`, `β''`, and the
  R² before trusting the reported gel-point verdict.
