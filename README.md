# DIW Bioprinting — Portable Analysis Pipeline

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Cite this](https://img.shields.io/badge/Cite-CITATION.cff-brightgreen.svg)](./CITATION.cff)

Self-contained export of the Python + MATLAB tooling used to characterise
DIW (Direct Ink Writing) bioinks, predict their printability from rheology,
and validate prints on the bench. Designed so a new intern or master's student
can pick this up on day one with their own ink and run the full pipeline.

> **License note:** This repository is licensed under [CC BY-NC 4.0](./LICENSE).
> GitHub's automatic license detector does not recognise CC BY-NC 4.0 (it
> only auto-tags 13 OSI-style licenses), so the sidebar may display
> "Other" — the licence is nonetheless the canonical CC BY-NC 4.0 legal
> code from creativecommons.org and is fully binding.

---

## Folder layout

```
Export/
├── README.md                        ← you are here
├── CLAUDE.md                        ← instructions for AI assistants
├── SCRIPT_REGISTRY.md               ← which version of each script is active
├── requirements.txt                 ← Python dependencies
│
├── 01_Python/                       ← rheology fitting + printability math
├── 02_MATLAB/                       ← extrusion simulation
├── 03_SOPs/                         ← printed-protocol PDFs + LaTeX sources
├── 04_Printability/                 ← image-based print-quality pipeline
│
├── docs/
│   ├── 00_workflow.md               ← the four-stage pipeline overview
│   ├── 01_rheology_data_formats.md  ← Anton Paar CSV vs legacy .xls contracts
│   └── manuals/                     ← one .md per script (inputs/outputs/caveats)
│
└── examples/                        ← (empty — place a representative .xls here for smoke tests)
```

---

## Quick start

### 1. Set up Python (one-time)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r Export/requirements.txt
```

### 2. Validate the data format

> The pipeline supports two input formats:
>
> - **Anton Paar Rheocompass CSV** (current, UTF-16, decimal comma) — for
>   flow curves, use the v5 scripts.
> - **Legacy TA-style `.xls`** — for the pre-migration archive only.
>
> Read **`docs/01_rheology_data_formats.md`** for the contract each parser
> expects and how to tell which path applies to your data.

### 3. Run the pipeline on one ink

```bash
cd "<YOUR PROJECT ROOT>"      # not Export/, the folder above it

# Stage 2 — fit the flow curve (Anton Paar CSV)
python Export/01_Python/Fit_Muitos_Modelos_v5.py

# Stage 2 — recovery % (still on legacy .xls — migration pending)
python Export/01_Python/Recovery_v1.py

# Stage 2 — h_max prediction (still on legacy .xls — migration pending)
python Export/01_Python/extract_hmax_v2.py
```

Open `Analises/Python/Results/FitAll-AntPar-v5.txt`. Find the row with the
best AIC for the Power-Law and Cross models. Copy the numbers into the
`samples(...)` struct at the top of `Export/02_MATLAB/run_solver_v3.m`.

```matlab
% In MATLAB:
cd '<YOUR PROJECT ROOT>'
addpath('Export/02_MATLAB')
run_solver_v3
```

This produces `output_v3/slicer_lookup_<sample>_<needle>.csv` — the table you
hand to whoever runs the printer.

### 4. Print

Follow `03_SOPs/Printing_Parameters_SOP_v4_EN.pdf`. The key rule:

> **Slicer Extrusion Multiplier = 1 / k_flow** (from the slicer-lookup CSV)

---

## The pipeline at a glance

```
RHEOLOGY  →  FITTING (Python)  →  SIMULATION (MATLAB)  →  PRINT (SOP)
   .xls         K, n, η₀, λ, m       v_print, ΔP, k_flow      validated
```

Full details and a parameter hand-off table: **`docs/00_workflow.md`**.

---

## What's where

| What you need to do | File to open |
|---|---|
| Understand the four stages | `docs/00_workflow.md` |
| Verify your `.xls` files are compatible | `docs/01_rheology_data_format_UNDER_REVIEW.md` |
| Learn a specific script (inputs, outputs, caveats) | `docs/manuals/<script>.md` |
| See which version of each script is current | `SCRIPT_REGISTRY.md` |
| Read the bench-side printing SOP | `03_SOPs/Printing_Parameters_SOP_v4_EN.pdf` |
| Read the long-form printing-parameters textbook | `03_SOPs/Printing_Parameters_Textbook_v2.tex` (+ PDF) |
| Run image-based scaffold quality analysis | `04_Printability/Printability/README.md` |
| Brief an AI assistant about this codebase | `CLAUDE.md` |

---

## Known limitations

- **Geometry defaults** are BD 10 mL syringe + 21G/22G blunt needles. If you
  print with different hardware, edit `geom_21G` / `geom_22G` in
  `02_MATLAB/run_solver_v3.m`.
- **No yield stress** in the extrusion solvers. Yield is handled separately via
  the `h_max` calculation in `01_Python/extract_hmax_v2.py`.
- **No wall slip, no viscoelasticity, no thixotropy** in the steady-state
  solvers. For pre-shear effects, fit the `Visco_Artur` (pre-sheared) data
  branch separately and compare K, n values.
- **Image-based printability pipeline (`04_Printability/`)** is a complete but
  more recent module. Its CLI and segmentation methods are documented inside
  the package — start with its own README and its `tests/` folder.

---

## Citations and credit

If outputs of these scripts appear in a thesis or publication, cite the
relevant references listed in each script's docstring (e.g., Smay 2002 for the
buckling correction, daSilva 2018 for CMC MW). Cross-check every numerical
output against literature ranges before trusting it.

---

## Contact

Maintainer: **Thiago Matos Custodio Rodrigues** — thiago.mcrodrigues@gmail.com
PEMM/COPPE, Universidade Federal do Rio de Janeiro (UFRJ)
Supervisor: Prof. Rossana Thiré
