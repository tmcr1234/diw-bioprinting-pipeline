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
| `antpar_io.py` | — | `active` | Anton Paar Rheocompass CSV (UTF-16 LE) **+ force-augmented TA LIMS `.txt` (latin-1)** | Library — flow + amplitude + frequency + recovery readers; imported by v5 / v2. **2026-08-26:** adds `read_flow_curve_with_force`, `normal_stress_from_force`, `baseline_correct`, `pair_audited_files`, and the `instrument` field on FlowCurve/AmplitudeSweep/FrequencySweep. The force reader is **specified, not validated** — written from the format contract, not against a specimen file; validate it once on a real export before quoting numbers |
| `Fit_Muitos_Modelos_v5.py` | v5 | `active` (current) | Anton Paar Rheocompass CSV | `./Analises/Python/Results/FitAll-AntPar-v5.{txt,csv}` |
| `extract_SAOS_amplitude_v2.py` | v2 | `active` (current) | Anton Paar Rheocompass CSV | `./Analises/Python/Results/SAOS_amplitude_v2.{txt,csv}` |
| `extract_SAOS_frequency_v2.py` | v2 | `active` (current) | Anton Paar Rheocompass CSV | `./Analises/Python/Results/SAOS_frequency_v2.{txt,csv}` |
| `extract_recovery_v2.py` | v2 | `active` (current) | Anton Paar Rheocompass CSV (3iTT trio) | `./Analises/Python/Results/Recovery_v2.{txt,csv}` |
| `extract_hmax_v3.py` | v3 | `active` (current) | the two v2 SAOS CSVs (no xlrd) | `./Analises/Python/Results/{SAOS_hmax_v3.txt, printing_parameters_per_ink_v3.csv, Gprime_extrap_*_v3.png}` |
| `Fit_Muitos_Modelos_v4.py` | v4 | `active` (legacy) | TA-style `.xls` (xlrd) | Use only on the pre-migration archive |
| `Recovery_v1.py` | v1 | `active` (legacy) | TA-style `.xls` (xlrd) | **SUPERSEDED** by extract_recovery_v2 — old steady-shear protocol; v2 uses the 3iTT-Osc-Rot-Osc protocol with different math |
| `extract_SAOS_values.py` | — | `active` (legacy) | TA-style `.xls` (xlrd) | **SUPERSEDED for SAOS extraction** by extract_SAOS_amplitude_v2 + extract_SAOS_frequency_v2. Still imported by `extract_Gpp_coxmerz_v1.py` for its `.xls` parser (`read_antpar_xls`, `find_col`, `extract_amplitude_sweep`) — keep |
| `extract_Gpp_coxmerz_v1.py` | v1 | `active` (new) | TA-style `.xls` (xlrd), via extract_SAOS_values | **G″ / loss-modulus analyses** the G′-only pipeline omitted: tan δ(ω) + Winter–Chambon gel-point test, Cox–Merz (η vs \|η*\|), G″ weak-strain-overshoot. Hard-coded paths in header dicts — **edit before reuse**. Outputs `Results/SAOS_Gpp_coxmerz_v1.txt`, `tandelta_coxmerz_per_ink.csv`, figs `F8_tandelta`/`F9_coxmerz`/`F10_Gpp_overshoot` (.pdf+.png) |
| `extract_hmax_v2.py` | v2 | `active` (legacy) | inherits from extract_SAOS_values | **SUPERSEDED** by extract_hmax_v3. Keep until v3 has been used end-to-end on a published dataset, then move to deprecated |
| `Calculo da Pressao de Extrusao.py` | — | `active` (UNDER_REVIEW) | TA-style `.xls` (xlrd) | Custom path set inside the script; **edit before first use** |
| `data_config.py` | — | `active` (new 2026-08-26) | — | Library — project-local folder/sample/geometry configuration, loaded from `data_config_local.py` in the project root. The Python counterpart of `inks_local.m`. `filter_files()` raises on a file naming no measuring system |
| `data_config_local.template.py` | — | `template` | — | Copy to the PROJECT ROOT as `data_config_local.py`. Never edited in place |
| `flow_diagnostics.py` | — | `active` (new 2026-08-26) | — | Library — `detect_startup_transient` (A7) and `detect_edge_fracture` (A8). Reports only; never trims data. Consumed by `Fit_Muitos_Modelos_v5.py` |
| `validate_fits.py` | — | `active` (maintainer tool, new 2026-08-26) | both | Cross-checks `Fit_Muitos_Modelos_v4` (audited `.xls`) against `v5` (active CSV) on shared samples and splits disagreement into **fitter-attributable** vs **ingestion-attributable**. Requires the pre-migration `.xls` archive; exits 2 with an explanation when it is absent. Sandboxes the v4 import — v4 is a flat script with no `main()` guard and executes on import |
| `check_flow_rate_consistency.py` | — | `active` (new 2026-08-26) | G-code + numbers | Compares solver `Q`, G-code `Q` and deposit-implied `Q` (`w·h·v_print`) side by side. `--self-test` reproduces the documented C20 discrepancy (0.71 / 0.20 / 0.46–0.51 mm³/s) with no files present. **Standard pre-submission step** |
| `extract_N1_tanner.py` | — | `active` (new 2026-08-26) | force-augmented `.txt` | **Produces `N1_wall_Pa` for `run_solver_v5`.** Reads a force export, excludes the edge-fracture window, log-log interpolates N1 and tau to the needle wall shear rate (read from `master_summary_v4.csv` or given), computes the Tanner beta and prints a paste-ready `inks_local.m` line. **Refuses to extrapolate past the measured shear-rate range** unless `--allow-extrapolation`, which stamps every downstream number. `--self-test` checks the closure maths with no data |
| `tests/test_antpar_force.py` | — | `active` (tests, new 2026-08-26) | synthetic fixtures | 12 acceptance tests for the force layer. Run `python3 tests/test_antpar_force.py` (no pytest needed). Synthetic by necessity — the audited C10–C25 force exports are not in every checkout, so the "reproduces the five audited floors" criterion remains **open** |

### Migration roadmap (which extractors still need an Anton Paar CSV port)

- [x] Flow / viscosity curves → done 2026-05-20 in v5 (`antpar_io.read_flow_curve_csv`)
- [x] Amplitude sweep (SAOS) → done 2026-05-20 in `extract_SAOS_amplitude_v2` (`antpar_io.read_amplitude_sweep_csv`)
- [x] Frequency sweep (SAOS) → done 2026-05-20 in `extract_SAOS_frequency_v2` (`antpar_io.read_frequency_sweep_csv`)
- [x] Recovery test → done 2026-05-20 in `extract_recovery_v2` (`antpar_io.read_recovery_set` — 3iTT-Osc-Rot-Osc protocol; computes both modulus G' and complex viscosity |η*| recovery ratios)
- [x] `extract_hmax_v2.py` → done 2026-05-20 in `extract_hmax_v3.py`. Reads the v2 SAOS CSVs directly (single source of truth — no re-parsing, no re-fitting). Matplotlib lazy-loaded so plots are optional. **The five core test types are migrated.**
- [ ] `extract_Gpp_coxmerz_v1.py` (G″ / tan δ / Cox–Merz / gel-point, added 2026-06-11) → **not yet ported.** Still imports the legacy `.xls` parser from `extract_SAOS_values.py`. To port: swap its `read_antpar_xls`/`extract_amplitude_sweep` calls for the equivalents in `antpar_io.py` + `extract_SAOS_amplitude_v2.py`/`extract_SAOS_frequency_v2.py`, and feed it the flow curve from the v5 Anton Paar reader for the Cox–Merz comparison.

**The core migrated pipeline (Stages 1–4) is `xlrd`-free.** The scripts that
still touch `.xls` are the legacy ingestion ones kept for pre-migration archive
reproducibility (Fit_v4, Recovery_v1, extract_SAOS_values, extract_hmax_v2,
Calculo da Pressao de Extrusao) **plus** `extract_Gpp_coxmerz_v1.py` — a newer
G″ analysis that still reads `.xls` via the legacy parser and has not been
ported to the Anton Paar CSV reader yet (see roadmap below).

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

As of the 2026-06-01 rewrite, the MATLAB layer is a **single active solver
plus its driver**. `bioprinting_algorithm_v4.m` is a verified superset of the
three earlier solvers; `run_solver_v4.m` is the sole driver.

| Script | Version | Status | Notes / Output location |
|---|---|---|---|
| `run_solver_v4.m` | v4 | `active` | **Sole driver.** One `bioprinting_algorithm_v4` call per `(ramp × ink × needle)`: combined `_data.txt`, per-Vp + summary PNGs, slicer CSVs, long-format `master_summary_v4.csv`. Output: `output_v4/`. See [`run_solver_v4.md`](docs/manuals/run_solver_v4.md) |
| `bioprinting_algorithm_v4.m` | v4 | `active` | **Unified superset solver (straight needle).** *2026-08-26: gained inline disclosure of the two heuristic closures, a single `crit_reynolds()` definition replacing two inline copies, and four OPTIONAL closure-override parameters used only by v5. All default-path numbers are unchanged — asserted by `validate_v4` (rel = 0) and by 20 identity checks in `validate_v5`.* Per Vp, both PL + Cross: full syringe+needle physics (pressure decomposition, hydrostatic, Reynolds/flow-regime, radial profiles) + slicer layer (`v_print`, `w_line`, `k_flow`). Validated against the legacy solvers via `validate_v4.m`. See [`bioprinting_algorithm_v4.md`](docs/manuals/bioprinting_algorithm_v4.md) |
| `bioprinting_algorithm_conical.m` | — | `active` | Conical-tip PL solver — only needed for conical nozzles; outside v4's straight-needle scope |
| `validate_v4.m` | — | `active` | Regression check: cross-validates `bioprinting_algorithm_v4` against the legacy PL/Cross solvers (in `archive/scripts/`) + v3 slicer formulas. Maintainer tool, not part of the student workflow. See [`validate_v4.md`](docs/manuals/validate_v4.md) |
| `bioprinting_algorithm_v5.m` | v5 | `draft` | **Physically anchored closures.** Replaces v4's two unsourced heuristics: die swell via the Tanner (1970) closure keyed on measured `N1` (falling back to the v4 heuristic, loudly, when `N1` is absent), and critical Reynolds via Ryan & Johnson (1959). **Does not fork the physics** — resolves the closures and delegates to v4 through its `BetaPL`/`BetaCross`/`RecritFcn` overrides. `draft` because it has never been run against a real measured `N1`; promote to `active` only after it has |
| `run_solver_v5.m` | v5 | `draft` | Driver for `bioprinting_algorithm_v5`. Same sweep as `run_solver_v4`; reads the optional `inks(i).N1_wall_Pa` and runs each ink on Tanner or, without N1, on the tagged heuristic fallback — mixing both in one sweep is allowed and reported. Writes `master_summary_v5.csv` (26 cols: adds `closure`, `Re_crit`, `N1_wall_Pa`) and `closure_report_v5.txt`, a generated disclosure block. **Duplicates only the orchestration loop — no physics, no slicer formula** |
| `validate_v5.m` | — | `active` (new 2026-08-26) | Regression check for v5. Asserts (1) v5 with no `N1` reproduces v4 on 20 numeric identities, (2) `tanner_beta(N1=0) = 0.13` and `ryan_johnson(n=1) = 2099.4`, (3) the closure swap is reported in `beta_report`. Fully synthetic — runs with no rheology data present |

### Retained-for-audit MATLAB — `02_MATLAB/archive/scripts/`

Two deprecated solvers are kept **only** because `validate_v4.m` cross-checks
against them. Their physics is fully folded into `bioprinting_algorithm_v4.m`.
**Do not call them for new work.**

| Script | Version | Status | Why retained |
|---|---|---|---|
| `bioprinting_algorithm_3.m` | — | `deprecated` | Legacy Power-Law solver — ground truth for the PL half of `validate_v4` |
| `bioprinting_algorithm_cross_v2.m` | v2 | `deprecated` | Legacy Cross solver — ground truth for the Cross half of `validate_v4` |

### Removed from this export (2026-06-01)

These were shipped in the 2026-05-21 release but are superseded by the v4
superset and no longer carried. They remain in the maintainer's project tree
(`Analises/MatLab/archive/scripts/`) for the paper trail.

| Script | Status | Superseded by |
|---|---|---|
| `run_solver_v3.m` | `deprecated` | `run_solver_v4.m` |
| `run_solver_improved.m` | `deprecated` | `run_solver_v4.m` |
| `run_solver_Cross_v2.m` | `deprecated` | `run_solver_v4.m` |
| `bioprinting_algorithm_v3.m` | `deprecated` | `bioprinting_algorithm_v4.m` (slicer layer) |
| `bioprinting_algorithm.m`, `bioprinting_algorithm_2.m`, `bioprinting_algorithm_cross.m` | `deprecated` | Pre-v3 solvers; never shipped |

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
