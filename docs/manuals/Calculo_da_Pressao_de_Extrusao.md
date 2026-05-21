# `Calculo da Pressao de Extrusao.py` — Quick Extrusion Pressure Estimate

> **STATUS:** `active` (legacy helper). Superseded for full simulations by the
> MATLAB solver, but useful as a sanity-check / back-of-envelope before running
> the simulator.

## Purpose

Reads the **maximum viscosity** from a flow curve and applies the Hagen–Poiseuille
analogue for a Power-Law fluid to estimate the minimum extrusion pressure
through a 21G needle, 22G needle, and the syringe barrel:

`ΔP = 4 · (L/R) · τ_w`

The script does **not** fit a Power-Law model — it pulls the stress at the
highest viscosity datapoint directly from the workbook. Use it for quick
feasibility numbers, not for publication-grade pressure values.

## Inputs

- Flow-curve `.xls` workbooks in the path set by `file_path` (line 28).
- Sample list in `samples` (line 25).

**You must edit these two paths before running** — the defaults still point to
an old project directory.

## Outputs

- `Pressão.txt` written to the path set by `save_path` (line 29).

Each line:
```
Sample: <name>, Max Viscosidade Row: <i>, Stress: <τ> Pa,
Pressão 21G: <kPa>, Pressão 22G: <kPa>, Pressão Seringa: <kPa>,
Force 21G: <N>, Force 22G: <N>, Force Syringe: <N>
```

## How to run

```bash
cd "<PROJECT ROOT>"
python "Export/01_Python/Calculo da Pressao de Extrusao.py"
```

## Dependencies

`pandas`, `xlrd`, `numpy`.

## Hardcoded geometry (BD 10 mL syringe + standard needles)

| Component | Symbol | Value |
|---|---|---|
| Syringe ID | `d_seringa` | 14.3 mm |
| Syringe length | `l_seringa` | 90 mm |
| 21G needle ID | `d_int_21_g` | 0.515 mm |
| 21G needle length | `l_1_25_pol` | 31.75 mm (1¼″) |
| 22G needle ID | `d_int_22_g` | 0.413 mm |
| 22G needle length | `l_1_pol` | 25.4 mm (1″) |

Update lines 8–13 if your hardware differs.

## Caveats

- Treats the entire ink as a single Newtonian fluid at its **peak viscosity**.
  This is conservative (over-estimates pressure) and is intentionally so for a
  feasibility check.
- Same `.xls` ingestion assumptions as the rest of the Python suite — file
  format is **UNDER REVIEW** for new rheometer exports.
- The legacy default sample names and paths inside this script point to an
  unrelated project ("Artigo Carol") — edit them before first use.
