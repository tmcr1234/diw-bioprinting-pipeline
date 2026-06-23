# CLAUDE.md — DIW Bioprinting Analysis Pipeline (Portable Export)

Instructions for AI assistants (Claude Code, Copilot, etc.) working with a new
student or researcher who has inherited this analysis pipeline.

> **This file is the portable, student-onboarding version.** The maintainer's
> full project-specific CLAUDE.md lives in the parent project root and
> contains lab-specific details (ink identity, sample registry, ongoing
> manuscripts). Do not assume those details apply here — this pipeline is
> the equipment-and-method-agnostic core.

---

## 1. The pipeline (read first)

All work in this lab follows the same strictly sequential four-stage pipeline.
**Stages are never skipped or reordered.**

```
[1] RHEOLOGY        Anton Paar (or equivalent) — .xls workbooks
       ↓
[2] MODELING        Python — fit K, n, η₀, η∞, λ, m, Rrec
       ↓
[3] SIMULATION      MATLAB — extrusion solver consuming Python parameters
       ↓
[4] PRINT VALIDATION   DIW printer — slicer-lookup CSV from Stage 3
```

Full description: `docs/00_workflow.md`.

---

## 2. Standard test battery (Stage 1)

| Test ID | Type | Protocol | Key output |
|---|---|---|---|
| Flow | Steady-state ramp | Shear rate ramp, equilibrium per point | η(γ̇), τ(γ̇) |
| Flow Sheared | Pre-sheared flow | 200 s⁻¹ × 300 s → ramp | Thixotropy / rheopexy diagnosis |
| Freq | SAOS frequency sweep | Fixed strain in LVR, sweep ω | G'(ω), G''(ω) |
| Strain | SAOS amplitude sweep | Fixed ω, sweep strain | LVR boundary, G', G'' |
| Recovery (50–300 s⁻¹) | 3-step recovery | Low-shear → disruption → low-shear | Structural recovery % |

### Thixotropy vs rheopexy interpretation

Compare Power-Law `K` from `Flow` vs `Flow Sheared`:

- `K` decreases after pre-shear → **thixotropic** (structure broken, recovers)
- `K` increases after pre-shear → **rheopexic** (shear-induced structuring)

This must be discussed explicitly in any manuscript that uses pre-shear data.

---

## 3. Python layer — critical rules

**Dependencies:** `numpy`, `pandas`, `scipy`, `matplotlib`. The current
Anton Paar Rheocompass pipeline is `xlrd`-free. `xlrd` is only needed if
you re-run the legacy `.xls` scripts (Fit_v4, Recovery_v1,
extract_SAOS_values, extract_hmax_v2) against the pre-migration archive;
do not substitute `openpyxl` for `xlrd` (that handles `.xlsx` only).

**All scripts must be run from the project root** (the folder containing
`Export/` and `Reologia/`), not from `Export/01_Python/`. They use relative
paths like `./Reologia/Viscosity`.

**Log-space fitting (Fit_Muitos_Modelos_v4):**
- Objective minimises `Σ(log ŷ − log y)²`
- Equivalent to WLS with `w = 1/y²` — gives each decade equal weight
- AIC/BIC are also computed in log space
- **Stress-domain models (HB, PL, Bingham, Casson) and viscosity-domain models
  (Cross, Carreau, Ellis) have incommensurable AIC values.** Compare within
  each domain only.

**Versioning rule:** every script is registered in `SCRIPT_REGISTRY.md` with
status `active`, `deprecated`, or `audited`. **Never run a deprecated script.
Never modify an audited script** (audited = output appeared in a submitted /
published manuscript and is frozen).

---

## 4. MATLAB layer — critical rules

**Flow assumptions (all solvers):**
- Fully developed laminar flow
- No entrance/exit losses
- No wall slip
- No yield stress (PL & Cross solvers — yield handled separately by
  `extract_hmax_v3.py`)

If any assumption is materially violated for your ink, flag it in the
manuscript.

**Standard hardware geometry** (defaults — override in `geom_*` structs if
your hardware differs):

| Component | ID (mm) | Length (mm) |
|---|---|---|
| Syringe (BD 10 mL) | 14.3 | 90 |
| Needle 21G (blunt) | 0.515 | 31.75 |
| Needle 22G (blunt) | 0.413 | 25.4 |

**Parameter sync (the rule that breaks pipelines):** every time Python re-fits
a model, **update the corresponding fields in the `inks(...)` struct in
`inks_local.m` in your project root on the same day** (each ink carries `Ramp1`
/ `Ramp2` sub-structs of `K_PL, n_PL, eta0, etaInf, lambda, m_Cross`). Stale
MATLAB parameters give silently wrong slicer CSVs.

**Shared Export / `inks_local.m` (project-local parameters):** `Export/` is a
single shared clone, symlinked into each project, so the fitted parameters live
**outside** it. `run_solver_v4` calls `inks_local()` from the current folder
(your project root); each project keeps its own `inks_local.m` there. Copy
`02_MATLAB/inks_local.template.m` to start. Do **not** hard-code parameters back
into the shared `run_solver_v4.m`. (Geometry `geom_21G`/`geom_22G` remains in
`run_solver_v4.m` — it is a shared hardware default; override per project only
if your hardware differs.)

**One solver, one driver:** the MATLAB layer is `bioprinting_algorithm_v4.m`
(unified PL + Cross superset, straight needle) driven by `run_solver_v4.m`.
Use `bioprinting_algorithm_conical.m` only for conical tips. The deprecated
solvers under `archive/scripts/` exist only for the `validate_v4.m` regression
check — never call them for new work.

**Slicer convention:** `k_flow` is **deposition efficiency**, not a slicer
input. The slicer's **Extrusion Multiplier = 1 / k_flow**.

---

## 5. Rheology data ingestion — two supported formats

The pipeline supports both, and **the Anton Paar Rheocompass migration is
complete (2026-05-20)**:

- **Anton Paar Rheocompass CSV** (current, default) — all four test types
  are fully migrated. Reader: `antpar_io.py`. Active scripts:
  - Flow curve fitting → `Fit_Muitos_Modelos_v5.py`
  - Amplitude sweep → `extract_SAOS_amplitude_v2.py`
  - Frequency sweep → `extract_SAOS_frequency_v2.py`
  - Recovery (3iTT-Osc-Rot-Osc) → `extract_recovery_v2.py`
  - h_max → `extract_hmax_v3.py` (reads the two v2 SAOS CSVs directly — no
    re-parsing, no xlrd)
- **Legacy TA-style `.xls`** — kept active **only** for the pre-migration
  archive (audited preprint data). Reader: `pandas + xlrd`. Use
  `Fit_Muitos_Modelos_v4.py`, `Recovery_v1.py`, `extract_SAOS_values.py`,
  `extract_hmax_v2.py`.

**The entire current pipeline is `xlrd`-free.** `xlrd` is only required if
you need to re-run the legacy scripts against the pre-migration archive.

**Full details and validation contract**: `docs/01_rheology_data_formats.md`.

The math layers (log-space fitting, AIC/BIC, MATLAB solver, h_max) are
equipment-agnostic. **Only the ingestion layer changes between v4 and v5.**

---

## 6. Documentation rules

- All outputs (code comments, manuscripts, analysis notes): **English**
- Academic tone: precise, impersonal, evidence-based
- All claims require citations; evidence hierarchy:
  meta-analyses > RCTs > cohort > case reports > expert opinion
- **Read existing LaTeX documents before writing new analysis** — do not
  contradict or duplicate established results without flagging the discrepancy

---

## 7. AI agent directives

1. **Read before writing.** Before producing new content, read the relevant
   manuals in `docs/manuals/` and any existing analysis in the project's
   `Latex/` folder. Never produce results that contradict established work
   without flagging it.

2. **No hallucinated references.** If literature search is needed, use real
   tools (PubMed, Semantic Scholar, ScienceDirect, Web of Science). Never
   fabricate citations.

3. **Parameter sync.** When the Python fit changes, immediately update the
   MATLAB `samples` struct. Never let them drift.

4. **Script-version enforcement.** Before running or modifying any script,
   check `SCRIPT_REGISTRY.md`. Only operate on `active` scripts. Never re-run
   or alter `audited` versions — read their outputs instead.

5. **File integrity.** Never overwrite completed analysis files. For
   extensions, create new files with clear version labels (`_v2`, `_draft2`).

6. **Match the model to the task complexity.** Use the lightest model that
   can do the work; reserve heavy reasoning for genuine multi-step synthesis
   (manuscript discussion, methodology decisions).

7. **Be concise.** Prefer tables and bullet points to dense prose. Don't
   re-print whole files for one-line edits; show diffs.

---

## 8. What you will not find here

This portable export deliberately excludes:

- Raw rheology data (`.xls` workbooks)
- Project-specific manuscripts
- Ink-specific sample identities, statistical-analysis outputs, FTIR data
- Print images and printability test results
- The maintainer's personal `.venv`

The new student supplies all of those for their own ink. This Export ships
**only the tooling and the methodology**.
