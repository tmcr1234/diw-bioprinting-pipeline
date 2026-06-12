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
| `extract_SAOS_values.py` | — | `active` (legacy) | TA-style `.xls` (xlrd) | **SUPERSEDED for SAOS extraction** by extract_SAOS_amplitude_v2 + extract_SAOS_frequency_v2. Still imported by `extract_Gpp_coxmerz_v1.py` for its `.xls` parser (`read_antpar_xls`, `find_col`, `extract_amplitude_sweep`) — keep |
| `extract_Gpp_coxmerz_v1.py` | v1 | `active` (new) | TA-style `.xls` (xlrd), via extract_SAOS_values | **G″ / loss-modulus analyses** the G′-only pipeline omitted: tan δ(ω) + Winter–Chambon gel-point test, Cox–Merz (η vs \|η*\|), G″ weak-strain-overshoot. Hard-coded paths in header dicts — **edit before reuse**. Outputs `Results/SAOS_Gpp_coxmerz_v1.txt`, `tandelta_coxmerz_per_ink.csv`, figs `F8_tandelta`/`F9_coxmerz`/`F10_Gpp_overshoot` (.pdf+.png) |
| `extract_hmax_v2.py` | v2 | `active` (legacy) | inherits from extract_SAOS_values | **SUPERSEDED** by extract_hmax_v3. Keep until v3 has been used end-to-end on a published dataset, then move to deprecated |
| `Calculo da Pressao de Extrusao.py` | — | `active` (UNDER_REVIEW) | TA-style `.xls` (xlrd) | Custom path set inside the script; **edit before first use** |

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
| `bioprinting_algorithm_v4.m` | v4 | `active` | **Unified superset solver (straight needle).** Per Vp, both PL + Cross: full syringe+needle physics (pressure decomposition, hydrostatic, Reynolds/flow-regime, radial profiles) + slicer layer (`v_print`, `w_line`, `k_flow`). Validated against the legacy solvers via `validate_v4.m`. See [`bioprinting_algorithm_v4.md`](docs/manuals/bioprinting_algorithm_v4.md) |
| `bioprinting_algorithm_conical.m` | — | `active` | Conical-tip PL solver — only needed for conical nozzles; outside v4's straight-needle scope |
| `validate_v4.m` | — | `active` | Regression check: cross-validates `bioprinting_algorithm_v4` against the legacy PL/Cross solvers (in `archive/scripts/`) + v3 slicer formulas. Maintainer tool, not part of the student workflow. See [`validate_v4.md`](docs/manuals/validate_v4.md) |

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
