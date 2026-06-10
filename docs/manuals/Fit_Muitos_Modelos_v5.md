# `Fit_Muitos_Modelos_v5.py` — Rheological Model Fitting (Anton Paar CSV)

> **STATUS:** `active`. **Use this for any data collected after the lab moved
> to the Rheocompass workflow.** For legacy TA-style `.xls` archives, use
> `Fit_Muitos_Modelos_v4.py` instead.

## Purpose

Fits the same **seven constitutive models** as v4 (Herschel-Bulkley, Power
Law, Bingham, Casson, Cross, Carreau, Ellis) to a steady-state flow curve and
ranks them by AIC/BIC in log space. **The fitting math is unchanged from v4.**

The only difference is the ingestion layer — v5 reads Anton Paar Rheocompass
raw CSVs via the `antpar_io` module instead of legacy `.xls` workbooks.

## Inputs

| Item | Path / format | Required |
|---|---|---|
| Raw flow-curve CSVs | Any folder of `*.csv` files exported by Rheocompass. `*_Analysis.csv` are skipped automatically. | yes |
| Sample include list | edit `SAMPLE_INCLUDE` (line ~165) — set to `None` to fit every CSV | optional |

**Expected file format** (locked to Rheocompass V1.36.x raw export):

- Encoding: **UTF-16 LE** with BOM
- Field separator: **TAB**
- Decimal separator: **comma** (Portuguese locale)
- Columns (after 8 header lines + units row): No. Ponto, Taxa de
  Cisalhamento [1/s], Tensão de Cisalhamento [Pa], **Viscosidade [mPa·s]**,
  Temperatura [°C], Torque [mN·m], Status
- **Viscosity is in mPa·s — `antpar_io` converts to Pa·s automatically.**

If your Rheocompass build changes the export schema, see *Caveats* below.

## Outputs

Written to `./Analises/Python/Results/`:

| File | Contents |
|---|---|
| `FitAll-AntPar-v5.txt` | Per-sample fit table, R² in log space, AIC, BIC, best-by-domain pick |
| `FitAll-AntPar-v5.csv` | Machine-readable long format — one row per (sample × model), parameters as columns |

## How to run

```bash
cd "<PROJECT ROOT>"
python Export/01_Python/Fit_Muitos_Modelos_v5.py
```

Must be run from the project root (the script uses relative paths). Edit
`INPUT_FOLDER` (line ~161) to point at your Rheocompass output folder.

## Dependencies

`numpy`, `pandas`, `scipy`, plus the local `antpar_io.py` module (ships in
the same folder). **No `xlrd` required** — that's v4's dependency.

## How v5 differs from v4

| Aspect | v4 | v5 |
|---|---|---|
| Source format | TA-style `.xls` workbook | Anton Paar Rheocompass raw CSV |
| Reader | `pandas.ExcelFile(engine='xlrd')` | `antpar_io.read_flow_curve_csv()` |
| Sheet selection | `SHEET_PREFIX = "Ramp - 2"` | n/a (one curve per file) |
| Column mapping | hard-coded indices (`COL_STRESS=2`, …) | from the file's own header row |
| Encoding handling | n/a | transparent UTF-16 LE → str |
| Decimal separator | `.` (TA exported it that way) | `,` (parser auto-converts) |
| Unit conversion | n/a | mPa·s → Pa·s automatic |
| Sample discovery | manual list + suffix string | folder scan, skip `*_Analysis.csv` |
| Fitting math | identical | identical (byte-for-byte) |

## Caveats

- **Ellis tends to win the viscosity-domain AIC** for our Rheocompass data —
  often by < 5 AIC units over Cross. That's the model's 3-parameter
  parsimony bonus, not a physical superiority. If the downstream MATLAB
  solver needs Cross parameters (`run_solver_v4.m` computes the Cross model
  alongside Power-Law), use the Cross row from the CSV regardless of who wins
  the AIC race.
- **HB usually collapses to PL** (`tau0 ~ 1e-16`) when the flow curve doesn't
  reach the LVR plateau at low γ̇. This is correct behaviour, not a bug —
  the data simply doesn't constrain `tau0`. Use `extract_hmax_v3.py` to
  get the yield stress from the amplitude sweep instead.
- **`Status` column is preserved** in the parsed DataFrame but **not used**
  by the fitter. Anton Paar emits flags like `M-` (motor limit hit),
  `Dy_auto` (dynamic equilibration). If you need to drop M- points before
  fitting, filter the DataFrame manually before calling `fit_one_sample`.

## Format drift — what to do if Rheocompass changes its export

The parser locates the data block by searching for the units row containing
`[1/s]` — that string has been stable across Rheocompass versions for many
years. If a future export changes this:

1. Open the CSV in a text editor (any editor that handles UTF-16).
2. Note the actual line containing `1/s` or the shear-rate unit string.
3. Edit `antpar_io.py:read_flow_curve_csv` — the search loop is one line
   (`if "[1/s]" in ln:`). Adjust the substring.

If the column order changes, the parser will still find the data block but
will mis-label columns. The fix is one block in `antpar_io.py` (the
`rows.append(...)` line where columns are pulled by index).

## Hand-off to MATLAB

After running, copy the **Cross row** (eta0, eta_inf, K, m) and the
**PowerLaw row** (K, n) for each ink from `FitAll-AntPar-v5.csv` into
the `inks(...)` struct at the top of `02_MATLAB/run_solver_v4.m` (fields live
on the `Ramp1` / `Ramp2` sub-struct for the no-pre-shear / pre-shear fit).
Field mapping (for a no-pre-shear fit):

| CSV column | MATLAB field |
|---|---|
| `K` (PowerLaw row) | `inks(i).Ramp1.K_PL` |
| `n` (PowerLaw row) | `inks(i).Ramp1.n_PL` |
| `eta0` (Cross row) | `inks(i).Ramp1.eta0` |
| `eta_inf` (Cross row) | `inks(i).Ramp1.etaInf` |
| `K` (Cross row) — Cross uses K for `1/λ` | `inks(i).Ramp1.lambda = 1/K` |
| `m` (Cross row) | `inks(i).Ramp1.m_Cross` |

> ⚠ Note: v5's Cross uses `K` as the inverse of the Cross time constant λ
> (this matches the curve_fit-friendly form `1 + (K·γ̇)ᵐ`). The MATLAB
> solver expects `lambda` directly. **Take `lambda = K_csv` if the CSV
> already named it `K` for the same `(K·γ̇)ᵐ` form** — keeping these
> consistent is on the user. Cross-check by computing
> `η(γ̇=10) = η∞ + (η₀−η∞) / (1 + (λ·10)^m)` in both layers; they should
> agree to within the fit's residual.
