# `antpar_io.py` — Anton Paar Rheocompass CSV Parser

> **STATUS:** `active`. Library module — imported by every Anton-Paar-aware
> Python script (`Fit_Muitos_Modelos_v5`, `extract_SAOS_amplitude_v2`,
> and future ports). Not meant to be run as a pipeline stage, but has a
> CLI smoke test for diagnostics.

## Purpose

Reads raw CSV exports from Anton Paar Rheocompass V1.36.x. Handles the
peculiarities of the export format that trip naive parsers up:

1. **UTF-16 LE encoding** with BOM (most CSV parsers assume UTF-8)
2. **Decimal comma** (Portuguese / European locale)
3. **Inconsistent torque units** across test types (mN·m for flow,
   µN·m for SAOS — careful when mixing)
4. **Viscosity in mPa·s** for flow curves (auto-converted to Pa·s)
5. **Fractional strain `[1]`** for amplitude sweeps (kept fractional,
   plus a convenience `strain_pct` column)

## Supported test types

| Function | Returns | Test type | Status |
|---|---|---|---|
| `read_flow_curve_csv` | `FlowCurve` | Steady-state flow / viscosity curve | ✅ done |
| `read_amplitude_sweep_csv` | `AmplitudeSweep` | SAOS amplitude (strain) sweep | ✅ done |
| `read_frequency_sweep_csv` | `FrequencySweep` | SAOS frequency sweep | ✅ done |
| `read_recovery_interval_csv` | `RecoveryInterval` | One interval of a 3iTT recovery test | ✅ done |
| `read_recovery_set(folder, basename)` | `RecoverySet` | Stitches the 3 intervals into one object | ✅ done |

## Public API

```python
from antpar_io import (
    FlowCurve, AmplitudeSweep, FrequencySweep,
    RecoveryInterval, RecoverySet,
    read_flow_curve_csv,
    read_amplitude_sweep_csv,
    read_frequency_sweep_csv,
    read_recovery_interval_csv,
    read_recovery_set,
    list_raw_csvs,                         # generic — every test type
    list_flow_curves, list_amplitude_sweeps, list_frequency_sweeps,   # back-compat aliases
    list_recovery_basenames,               # groups _1/_2/_3 trios
    read_flow_curves_in_folder,
    read_amplitude_sweeps_in_folder,
    read_frequency_sweeps_in_folder,
    read_recovery_sets_in_folder,
)

# Flow curve
fc = read_flow_curve_csv("path/to/Visco.csv")
# fc.data columns: point, shear_rate_1_s, shear_stress_Pa, viscosity_Pa_s,
#                  temperature_C, torque_mNm, status

# Amplitude sweep
sweep = read_amplitude_sweep_csv("path/to/Strain.csv")
# sweep.data columns: point, strain_frac, strain_pct,
#                     shear_stress_Pa, G_storage_Pa, G_loss_Pa, tan_delta,
#                     torque_uNm, status

# Frequency sweep
fsw = read_frequency_sweep_csv("path/to/Freq.csv")
# fsw.data columns: point, omega_rad_s,
#                   G_storage_Pa, G_loss_Pa, tan_delta,
#                   eta_complex_real_Pa_s, eta_complex_imag_Pa_s, G_complex_Pa,
#                   strain_pct, shear_stress_Pa, torque_mNm, status

# Recovery — 3iTT stored as THREE separate CSVs per sample (_1, _2, _3)
rs = read_recovery_set("Reologia/Recovery - 3D Bioprinting/",
                        "Carla_Amostra A - CP50 25C_12_05_2026")
# rs.interval1 / .interval2 / .interval3 each contain a DataFrame plus
# .kind ("oscillation" or "rotation").
rsets = read_recovery_sets_in_folder("Reologia/Recovery - 3D Bioprinting/")
# returns dict {short_name: RecoverySet}, e.g. {"Amostra A": <RecoverySet>}

# Folder discovery — skips *_Analysis.csv (EN) and *_Analise.csv (PT)
csvs = list_raw_csvs("Reologia/Strain Sweep - 3D Bioprinting/")
sweeps = read_amplitude_sweeps_in_folder("Reologia/Strain Sweep - 3D Bioprinting/")
fsweeps = read_frequency_sweeps_in_folder("Reologia/Frequency Sweep - 3D Bioprinting/")
# returns dict {short_name: <type>}, e.g. {"Amostra A": <FrequencySweep>}
```

Short-name extraction: the dict key is the part of the `Teste:` header
line before the first ` - ` (e.g. `"Amostra A - PP25 1mm 25C"` →
`"Amostra A"`), not the filename. That makes keys consistent across
test types even when filenames vary (`AMOSTRA A - Visco …` vs
`Carla_Amostra A - PP25 …`).

## CLI smoke test

```bash
# Default mode: flow curves
python Export/01_Python/antpar_io.py /path/to/Reologia/Viscosity

# Amplitude-sweep mode (add --amp flag)
python Export/01_Python/antpar_io.py /path/to/Reologia/Strain\ Sweep --amp

# Frequency-sweep mode (add --freq flag)
python Export/01_Python/antpar_io.py /path/to/Reologia/Frequency\ Sweep --freq

# Recovery mode (folder must contain _1/_2/_3 trios per sample)
python Export/01_Python/antpar_io.py /path/to/Reologia/Recovery\ -\ 3D\ Bioprinting --rec

# Inspect one file in detail
python Export/01_Python/antpar_io.py /path/to/Amostra\ A.csv --freq
```

Folder mode prints a one-line per-sample summary (test-type appropriate:
γ̇ / η range for flow; γ / G' / tan δ for amplitude; ω / G' / γ_set
for frequency; intervals + kinds for recovery). File mode dumps the
metadata plus first 5 / last 3 rows.

## What the parser explicitly handles

| Quirk | Where in the file | Handling |
|---|---|---|
| BOM `0xFF 0xFE` | byte 0–1 | stripped after decode |
| UTF-16 LE | whole file | `bytes.decode("utf-16-le")` |
| Decimal `,` | every numeric cell | `str.replace(",", ".")` then `float()` |
| Header offset | lines 1–10 | parsed by label, not index |
| Units row (flow) | line containing `[1/s]` | located by substring search |
| Units row (amplitude) | line containing `[µN·m]` (fallback: `Deforma` header) | located by substring search |
| Units row (frequency) | line containing `[rad/s]` | located by substring search — unique to freq sweep |
| Recovery interval splits across 3 files | filenames end in `_1/_2/_3` | grouped by `list_recovery_basenames`; stitched by `read_recovery_set` |
| Osc vs rotation distinction (recovery) | column header keywords | `Armazenamento`/`G'` → oscillation; `Taxa de Cisalhamento` → rotation |
| Leading tab in data rows | every data line | dropped (cell-zero == `""` check) |
| Scientific notation `9,9E-05` (with comma decimal) | first/last points | comma → dot, then `float()` |
| `*_Analysis.csv` (EN) and `*_Analise.csv` (PT) | filename | both filtered by `list_raw_csvs` |
| Report PDFs co-located in the folder | filename | filtered by extension |
| Viscosity in mPa·s (flow only) | column index 3 | **divided by 1000** before returning |
| Strain in fractional `[1]` (amplitude only) | column index 1 | kept fractional; `strain_pct` derived column added |
| Complex viscosity Re/Im in mPa·s (frequency only) | column indices 5–6 | **divided by 1000** before returning |
| Strain in percent `[%]` (frequency only) | column index 8 | kept as-is — already in percent, unlike amplitude `[1]` |
| Torque unit varies by test | various | flow: `[mN·m]`, amplitude: `[µN·m]`, frequency: `[mN·m]` — preserved in differently-named columns (`torque_mNm`, `torque_uNm`) |

## What the parser does NOT handle (yet)

- **Frequency sweeps** and **amplitude sweeps** — those have different
  column sets. A sister parser will be added when we migrate those
  scripts. For now, frequency/strain ingestion still uses `extract_SAOS_values.py`
  on `.xls` files.
- **Recovery tests** with multiple intervals in one CSV. The current parser
  reads a single interval. Recovery scripts will need a multi-interval
  variant.
- **Excel `.xlsx` exports** from a future Rheocompass build. The decoder
  will raise `ValueError` rather than silently mangle.

## Format contract (the bit that may need updating someday)

The parser is locked to this contract:

```
Line 1:  Projeto:<TAB>Viscosity Curve - 3D Printing
Line 3:  Teste:<TAB>AMOSTRA <X> - Visco <geom> <temp> <gap>
Line 5:  Resultado:<TAB>Viscosity curve <N>
Line 7:  Intervalo e ponto de dados:<TAB>1<TAB><n_points>
Line 8:  Dados do intervalo:<TAB>No. do Ponto<TAB>Taxa de Cisalhamento<TAB>
         Tensão de Cisalhamento<TAB>Viscosidade<TAB>Temperatura<TAB>
         Torque<TAB>Status
Line 10: <TAB><TAB>[1/s]<TAB>[Pa]<TAB>[mPa·s]<TAB>[°C]<TAB>[mN·m]<TAB>
Line 11+: <TAB><point>,<TAB><γ̇>,<TAB><τ>,<TAB><η_mPa·s>,<TAB><T>,<TAB><M>,<TAB><status>
```

If Rheocompass ever changes any of the labelled keys (Projeto, Teste,
Resultado, the units row, column order), update the corresponding section
in `_HEADER_LABELS` and the `rows.append(...)` block.

## Author

T.M.C. Rodrigues (PEMM/COPPE/UFRJ) — 2026-05-20.
