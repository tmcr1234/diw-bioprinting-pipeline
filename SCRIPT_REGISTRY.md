# Script Version Registry

Single source of truth for **which version of each script is currently active**.
Read this **before** running or modifying any script.

## Status definitions

| Status | Meaning | Can it be run? | Can it be edited? |
|---|---|---|---|
| `active` | Current production version | ✅ yes | ✅ yes (carefully) |
| `draft` | Work-in-progress, not yet promoted | ⚠️ for testing only | ✅ yes |
| `deprecated` | Superseded; kept for audit context | ❌ no | ❌ no |
| `audited` | Output appears in a submitted/published manuscript | ❌ no — frozen | ❌ no — frozen |

If a script's status is unclear, **ask the maintainer before running it**.

---

## Python — `01_Python/`

| Script | Version | Status | Source format | Notes / Output location |
|---|---|---|---|---|
| `antpar_io.py` | — | `active` | Anton Paar Rheocompass CSV (UTF-16 LE) | Library — flow + amplitude + frequency + recovery readers; imported by v5 / v2 |
| `Fit_Muitos_Modelos_v5.py` | v5 | `active` (current) | Anton Paar Rheocompass CSV | `./Analises/Python/Results/FitAll-AntPar-v5.{txt,csv}` |
| `extract_SAOS_amplitude_v2.py` | v2 | `active` (current) | Anton Paar Rheocompass CSV | `./Analises/Python/Results/SAOS_amplitude_v2.{txt,csv}` |
| `extract_SAOS_frequency_v2.py` | v2 | `active` (current) | Anton Paar Rheocompass CSV | `./Analises/Python/Results/SAOS_frequency_v2.{txt,csv}` |
| `extract_recovery_v2.py` | v2 | `active` (current) | Anton Paar Rheocompass CSV (3iTT trio) | `./Analises/Python/Results/Recovery_v2.{txt,csv}` |
| `extract_hmax_v3.py` | v3 | `active` (current) | the two v2 SAOS CSVs (no xlrd) | `./Analises/Python/Results/{SAOS_hmax_v3.txt, printing_parameters_per_ink_v3.csv, Gprime_extrap_*_v3.png}` |
| `Fit_Muitos_Modelos_v4.py` | v4 | `active` (legacy) | TA-style `.xls` (xlrd) | Use only on the pre-migration archive |
| `Recovery_v1.py` | v1 | `active` (legacy) | TA-style `.xls` (xlrd) | **SUPERSEDED** by extract_recovery_v2 — old steady-shear protocol; v2 uses the 3iTT-Osc-Rot-Osc protocol with different math |
| `extract_SAOS_values.py` | — | `active` (legacy) | TA-style `.xls` (xlrd) | **FULLY SUPERSEDED** by extract_SAOS_amplitude_v2 + extract_SAOS_frequency_v2. No longer imported by anything — keep only for archival reference |
| `extract_hmax_v2.py` | v2 | `active` (legacy) | inherits from extract_SAOS_values | **SUPERSEDED** by extract_hmax_v3. Keep until v3 has been used end-to-end on a published dataset, then move to deprecated |
| `Calculo da Pressao de Extrusao.py` | — | `active` (UNDER_REVIEW) | TA-style `.xls` (xlrd) | Custom path set inside the script; **edit before first use** |

### Migration roadmap (which extractors still need an Anton Paar CSV port)

- [x] Flow / viscosity curves → done 2026-05-20 in v5 (`antpar_io.read_flow_curve_csv`)
- [x] Amplitude sweep (SAOS) → done 2026-05-20 in `extract_SAOS_amplitude_v2` (`antpar_io.read_amplitude_sweep_csv`)
- [x] Frequency sweep (SAOS) → done 2026-05-20 in `extract_SAOS_frequency_v2` (`antpar_io.read_frequency_sweep_csv`)
- [x] Recovery test → done 2026-05-20 in `extract_recovery_v2` (`antpar_io.read_recovery_set` — 3iTT-Osc-Rot-Osc protocol; computes both modulus G' and complex viscosity |η*| recovery ratios)
- [x] `extract_hmax_v2.py` → done 2026-05-20 in `extract_hmax_v3.py`. Reads the v2 SAOS CSVs directly (single source of truth — no re-parsing, no re-fitting). Matplotlib lazy-loaded so plots are optional. **The migration is complete.**

**The entire active pipeline is now `xlrd`-free.** The only scripts that
still touch `.xls` are the explicitly legacy ones (Fit_v4, Recovery_v1,
extract_SAOS_values, extract_hmax_v2, Calculo da Pressao de Extrusao),
kept around for pre-migration archive reproducibility.

### Deprecated / audited Python (NOT shipped in this export)

| Script | Version | Status | Why |
|---|---|---|---|
| `Recovery_Prediction_v1.py` | v1 | `deprecated` | Superseded by `Recovery_v1.py` |
| `extract_hmax_v1.py` | v1 | `deprecated` | Superseded by `extract_hmax_v2.py` |
| `Fit_Muitos_Modelos_v1–v3.py` | v1–v3 | `audited` | Earlier fitting versions; outputs frozen |

> **Note on `Fit_Muitos_Modelos_v4.py`**: kept as `active (legacy)` rather
> than deprecated because the pre-migration `.xls` archive is still
> referenced by the published preprint. Anyone re-running those fits must
> use v4 to reproduce the audited numbers; v5 would re-fit on different
> ingestion paths and could give slightly different parameters at the
> last significant figure.

---

## MATLAB — `02_MATLAB/`

| Script | Version | Status | Notes / Output location |
|---|---|---|---|
| `run_solver_v4.m` | v4 | `draft` (2026-05-21) | **One-shot master driver** — joins legacy per-sample radial profiles (PL + Cross) with v3 slicer lookup and a long-format master summary. Output: `output_v4/`. Smoke-test pending; will supersede `run_solver_v3.m` once validated. See [`run_solver_v4.md`](docs/manuals/run_solver_v4.md) |
| `run_solver_v3.m` | v3 | `active` | **Current primary driver.** `output_v3/slicer_lookup_*.csv`, `output_v3/master_summary_v3.csv` |
| `bioprinting_algorithm_v3.m` | v3 | `active` | Called by `run_solver_v3.m` and (for the slicer layer) by `run_solver_v4.m` |
| `run_solver_improved.m` | — | `active` | Legacy PL-only driver, radial profiles |
| `bioprinting_algorithm_3.m` | — | `active` | Called by `run_solver_improved.m` and (for the per-sample PL layer) by `run_solver_v4.m` |
| `run_solver_Cross_v2.m` | v2 | `active` | Cross-only driver, tapered nozzle |
| `bioprinting_algorithm_cross_v2.m` | v2 | `active` | Called by `run_solver_Cross_v2.m` and (for the per-sample Cross layer, degenerate to straight cylinder) by `run_solver_v4.m` |
| `bioprinting_algorithm_conical.m` | — | `active` | Conical-tip PL solver — only needed for conical nozzles |

### Deprecated / audited MATLAB (NOT shipped)

| Script | Status | Why |
|---|---|---|
| `bioprinting_algorithm.m` | `deprecated` | Pre-v3 fixed-Vp Power-Law solver |
| `bioprinting_algorithm_2.m` | `deprecated` | Intermediate revision |
| `bioprinting_algorithm_cross.m` | `deprecated` | Replaced by `_v2` (better convergence brackets) |

---

## Adding a new version

When you produce a new script version, before promoting it to `active`:

1. **Increment the version suffix** in the filename (`_v5`, not in-place edit).
2. Add a row to this registry with status `draft`.
3. Run on at least one known-good dataset and compare against the previous
   `active` version. Document any differences in the script's docstring.
4. Promote the previous version to `deprecated` (or `audited`, if its output
   appears in a submitted manuscript).
5. Flip the new version's status to `active`.

**Never delete a deprecated or audited script.** Keep them for paper-trail
purposes; just remove them from `Export/` and document the removal here.
