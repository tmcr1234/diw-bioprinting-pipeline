# Manuals Index

One manual per active script. Click through for inputs, outputs, dependencies,
and known caveats. **Always read the manual before running a script for the
first time on a new ink.**

## Python

| Script | Source format | Purpose | Manual |
|---|---|---|---|
| `antpar_io.py` | Anton Paar CSV | Library: parse Rheocompass flow + amplitude + frequency + recovery (3iTT) | [antpar_io.md](antpar_io.md) |
| **`Fit_Muitos_Modelos_v5.py`** | Anton Paar CSV | **Current** — fit 7 rheological models to a flow curve | [Fit_Muitos_Modelos_v5.md](Fit_Muitos_Modelos_v5.md) |
| **`extract_SAOS_amplitude_v2.py`** | Anton Paar CSV | **Current** — LVR extraction from amplitude (strain) sweeps | [extract_SAOS_amplitude_v2.md](extract_SAOS_amplitude_v2.md) |
| **`extract_SAOS_frequency_v2.py`** | Anton Paar CSV | **Current** — G'(ω) picks + low-ω power-law extrapolation | [extract_SAOS_frequency_v2.md](extract_SAOS_frequency_v2.md) |
| **`extract_recovery_v2.py`** | Anton Paar CSV | **Current** — 3iTT modulus recovery + complex viscosity recovery | [extract_recovery_v2.md](extract_recovery_v2.md) |
| **`extract_hmax_v3.py`** | the two v2 SAOS CSVs | **Current** — h_max + N_max under 4 yield criteria (no xlrd) | [extract_hmax_v3.md](extract_hmax_v3.md) |
| `Fit_Muitos_Modelos_v4.py` | TA-style `.xls` | Legacy — for the pre-migration archive | [Fit_Muitos_Modelos_v4.md](Fit_Muitos_Modelos_v4.md) |
| `Recovery_v1.py` | TA-style `.xls` | Legacy steady-shear-ramp recovery — superseded by 3iTT-based extract_recovery_v2 | [Recovery_v1.md](Recovery_v1.md) |
| `extract_SAOS_values.py` | TA-style `.xls` | Legacy: fully superseded by the two v2 scripts; archive reference only | [extract_SAOS_values.md](extract_SAOS_values.md) |
| `extract_hmax_v2.py` | inherits | Legacy — superseded by extract_hmax_v3 | [extract_hmax_v2.md](extract_hmax_v2.md) |
| `Calculo da Pressao de Extrusao.py` | TA-style `.xls` | Quick Hagen–Poiseuille pressure estimate | [Calculo_da_Pressao_de_Extrusao.md](Calculo_da_Pressao_de_Extrusao.md) |

## MATLAB — drivers (run these)

| Script | Purpose | Manual |
|---|---|---|
| `run_solver_v4.m` | **Draft** (2026-05-21) — one-shot master: legacy radial profiles + v3 slicer CSVs + long-format master summary in one run | [run_solver_v4.md](run_solver_v4.md) |
| `run_solver_v3.m` | **Current primary** — PL + Cross, straight cylinder, slicer CSVs | [run_solver_v3.md](run_solver_v3.md) |
| `run_solver_improved.m` | Legacy PL-only driver, radial profiles | [run_solver_improved.md](run_solver_improved.md) |
| `run_solver_Cross_v2.m` | Cross-only driver, tapered nozzle | [run_solver_Cross_v2.md](run_solver_Cross_v2.md) |

## MATLAB — solver functions (called by drivers, not run directly)

| Function | Geometry | Model | Manual |
|---|---|---|---|
| `bioprinting_algorithm_v3.m` | straight cylinder | PL + Cross | [bioprinting_algorithm_v3.md](bioprinting_algorithm_v3.md) |
| `bioprinting_algorithm_3.m` | straight cylinder | PL only | [bioprinting_algorithm_3.md](bioprinting_algorithm_3.md) |
| `bioprinting_algorithm_cross_v2.m` | tapered nozzle | Cross | [bioprinting_algorithm_cross_v2.md](bioprinting_algorithm_cross_v2.md) |
| `bioprinting_algorithm_conical.m` | conical tip | PL | [bioprinting_algorithm_conical.md](bioprinting_algorithm_conical.md) |

## Standalone module

- **Printability pipeline** — `Export/04_Printability/Printability/` has its own
  `pyproject.toml`, test suite, and CLI. See its own README and the SOP at
  `Export/03_SOPs/Printing_Parameters_SOP_v4_EN.pdf`.
