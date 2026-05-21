# Rheology Data Ingestion — Two Supported Formats

> **READ THIS BEFORE RUNNING ANY PYTHON SCRIPT.**

This pipeline now supports **two** input formats for rheology data, reflecting
the lab's migration from one rheometer workflow to another:

| Workflow | Format | Reader | Status |
|---|---|---|---|
| **Anton Paar Rheocompass** (current) | UTF-16 CSV, tab-separated, decimal comma | `antpar_io.py` | ✅ active (flow curves only — see migration roadmap) |
| Legacy TA-style `.xls` | binary `.xls`, `xlrd`-readable | `pandas + xlrd` | ⚠ legacy — kept for pre-migration archive reproducibility |

**The flow-curve / viscosity-curve path is fully migrated.** Use
`Fit_Muitos_Modelos_v5.py` for new Rheocompass data. The other tests
(recovery, frequency sweep, amplitude sweep) still consume the legacy
format — see the migration roadmap in `SCRIPT_REGISTRY.md`.

---

## A. Anton Paar Rheocompass CSV (current) — flow curves

### Format contract

The parser in `antpar_io.py` is locked to this layout (Rheocompass V1.36.x
raw export, Portuguese locale):

- **Encoding**: UTF-16 LE with BOM
- **Field separator**: TAB
- **Decimal separator**: comma (`,`)
- **Line endings**: CRLF
- **Header**: 8 lines of metadata + 1 units row containing `[1/s]`
- **Columns** (after the units row):
  `point | shear_rate [1/s] | shear_stress [Pa] | viscosity [mPa·s] | T [°C] | torque [mN·m] | status`

**The `_Analysis.csv` siblings are Rheocompass's built-in regression
outputs.** The parser skips them automatically.

### Critical unit conversion

`antpar_io.read_flow_curve_csv` **divides the viscosity column by 1000**
on the way in (mPa·s → Pa·s). Downstream code (the fitter, the MATLAB
solver) all assumes Pa·s. If you ever bypass `antpar_io` and read the
CSV directly, **do not forget this conversion**.

### Quick smoke test on new data

```bash
python Export/01_Python/antpar_io.py /path/to/your/Reologia/Viscosity/
```

Expected output: one summary line per sample, with the shear-rate range
roughly `[0.1, 1000] s⁻¹` and viscosities in Pa·s (typically 0.1–500 Pa·s
for DIW gels). If you see viscosities in the 10²–10⁵ range, you're
probably reading mPa·s — check the unit conversion.

### If Rheocompass changes its export format

The parser is robust against minor drift because it locates the data
block by keyword (`[1/s]`) rather than by fixed row index. Things that
*would* break it:

| What changed | Where to fix |
|---|---|
| Different units row keyword (e.g., `s⁻¹` instead of `[1/s]`) | `antpar_io.py`, search-loop line near `"[1/s]"` |
| Column order changed | `antpar_io.py`, the `rows.append(...)` block — indices 1–5 |
| `.xlsx` export instead of `.csv` | rewrite using `openpyxl` — keep the same `FlowCurve` return type so the fitter doesn't notice |
| Locale change to decimal `.` | `_parse_decimal_comma` already handles it gracefully — but verify on a smoke test |

---

## B. TA-style `.xls` (legacy) — recovery, frequency, amplitude sweeps

### Format contract

These tests are still consumed by the v4-era scripts:

| Script | Sheet name pattern | Skip rows | Critical columns (0-based, post-skip) |
|---|---|---|---|
| `Fit_Muitos_Modelos_v4.py` | `Ramp - 2` (configurable via `SHEET_PREFIX`) | 3 | stress=2, shear_rate=3, viscosity=4 |
| `Recovery_v1.py` | `Ramp*` (≥ 3 sheets) | 3 | viscosity=4 |
| `Calculo da Pressao de Extrusao.py` | `Ramp*` | 3 | stress=2, viscosity=4 |
| `extract_SAOS_values.py` | any (keyword-scanned) | auto | header-based (most tolerant) |
| `extract_hmax_v2.py` | inherits | inherits | inherits |

### Validation checklist for the legacy format

Before running any v4-era script on new `.xls` data:

1. **Open one workbook by hand**. Note: sheet names, header offset, column order.
2. **Compare against the table above.** If anything differs, edit the
   corresponding `skiprows`, `COL_*`, or `SHEET_PREFIX` constants.
3. **Sanity-check the output**: Power-Law `K` 10–1000 Pa·sⁿ for typical DIW
   gels, `n` 0.1–0.6 (shear-thinning), recovery 30–100%.

### Why the legacy path is still here

The pre-migration thesis data lives in `.xls` workbooks and is referenced
by `Latex/Preprint_revised_v17.tex` (and the published preprint at the time
of the migration). If you need to reproduce any of those numbers, you must
re-run **v4** against the original `.xls` files. Re-fitting them through v5
on a new export would give slightly different last-digit answers and would
break the audit trail.

### How to tell which path a script is on

Open the script and look at the docstring banner. Scripts that still need
to be ported to Anton Paar carry a `UNDER REVIEW` block in their docstring.
Scripts already ported (currently only `Fit_Muitos_Modelos_v5`) carry no
banner.

---

## Migration roadmap (rest of the pipeline)

| Test | Current reader | Target | Status |
|---|---|---|---|
| ✅ Flow curve | `Fit_Muitos_Modelos_v5` + `antpar_io.read_flow_curve_csv` | — | done 2026-05-20 |
| ✅ Amplitude (strain) sweep | `extract_SAOS_amplitude_v2` + `antpar_io.read_amplitude_sweep_csv` | — | done 2026-05-20 |
| ✅ Frequency sweep | `extract_SAOS_frequency_v2` + `antpar_io.read_frequency_sweep_csv` | — | done 2026-05-20 |
| ✅ Recovery (3iTT) | `extract_recovery_v2` + `antpar_io.read_recovery_set` | — | done 2026-05-20 (protocol changed: legacy v1 was steady-shear ramps; v2 is 3iTT-Osc-Rot-Osc) |
| ✅ h_max | `extract_hmax_v3` reading the two v2 SAOS CSVs directly | — | done 2026-05-20 — **migration complete**, no more `xlrd` in the active pipeline |

To advance the roadmap, drop one Rheocompass export per test type into a
known folder and ask for the next migration step.

### Format quirks discovered during the amplitude-sweep port (2026-05-20)

- The "analysis" sibling files are named `_Analise.csv` (**Portuguese**)
  for the strain sweep, not `_Analysis.csv` (English) like for the flow
  curves. `list_raw_csvs` now filters both spellings.
- Strain sweeps carry **no temperature column** — Rheocompass omits it.
  Use the filename / `Teste:` header for setpoint T.
- Torque unit differs between test types: **mN·m** in flow curves,
  **µN·m** in amplitude sweeps. The parser preserves the native unit
  in differently-named columns (`torque_mNm`, `torque_uNm`) to make this
  visible — don't sum them across test types without converting first.
- Strain is exported **fractional** `[1]` (0.01 = 1%). The parser keeps
  it fractional and adds a derived `strain_pct` convenience column,
  because γ_LVR is conventionally reported in percent.
- Anton Paar's `M-` Status flag = "torque below sensitivity floor" =
  effectively a noise-floor point. `extract_SAOS_amplitude_v2` drops
  these by default before LVR detection.

### Format quirks discovered during the recovery (3iTT) port (2026-05-20)

- **One 3iTT measurement = THREE separate CSV files** (`_1.csv`, `_2.csv`,
  `_3.csv`), one per interval, sharing a basename. The parser provides
  `read_recovery_interval_csv` for one file and `read_recovery_set` to
  stitch the trio. `list_recovery_basenames` groups files; only
  *complete* trios are returned (incomplete sets are warned).
- The "analysis" sibling is named `_analise.csv` (Portuguese spelling
  again, as with the strain sweep), already filtered by `list_raw_csvs`.
- **Interval kind varies per file** in the same measurement.
  Osc-Rot-Osc means interval 1 and 3 share one schema (G', G'', tan δ,
  |η*|, torque [µN·m]) while interval 2 has a different one (shear rate,
  stress, viscosity, torque). The parser auto-detects via column-header
  keywords (`Armazenamento` → oscillation; `Taxa de Cisalhamento` →
  rotation) and selects the right column mapping.
- **In oscillation intervals, Rheocompass leaves the "Viscosidade" column
  BLANK** but the column still exists in the schema. Parser stores it
  as NaN.
- **Time column** is present in every interval, in seconds. Useful for
  plotting the recovery trajectory.
- **Crosslinking during the recovery interval** (interval 3) is the
  major data-quality risk for chemically-active inks like alginate +
  Ca²⁺. G' shoots above the interval-1 baseline, producing false
  "super-recovery" values. `extract_recovery_v2` handles this with a
  per-sample `INTERVAL3_MAX_POINT_OVERRIDES` cutoff; the huge std on
  the recovered tail flags samples whose cutoff is still too generous.

### Format quirks discovered during the frequency-sweep port (2026-05-20)

- The frequency-sweep CSV has **12 data columns** — more than flow (7)
  or amplitude (8). Includes complex viscosity Re/Im, |G*|, and the
  *controlled* strain amplitude in addition to the moduli.
- Strain in this test is the **controlled SAOS amplitude in percent
  `[%]`**, not fractional `[1]` like the amplitude sweep. The parser
  stores it as-is in `strain_pct`.
- Complex viscosity Re/Im is in `[mPa·s]` (same as flow-curve
  viscosity). The parser converts both to Pa·s.
- Torque is back to `[mN·m]` (like the flow curve), not `[µN·m]` like
  the amplitude sweep. Three test types, two torque units, no other
  pattern — the parser names the column `torque_mNm` here, same as the
  flow curve, so the discrepancy is visible.
- A new Status flag, `ME-`, appears in frequency sweeps at low ω. It
  means "modulus estimate error" — G' has collapsed to numerical noise
  (1e-5 Pa region). `extract_SAOS_frequency_v2` drops both `M-` and
  `ME-` by default before the low-ω power-law fit.
- The `taD` flag (tan-delta drift) is kept by default because it
  usually accompanies high-ω inertial corrections where G' is still
  trustworthy; only tan δ may be unreliable.
