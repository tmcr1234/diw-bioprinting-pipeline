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

**Project-local configuration — `data_config_local.py`:** folder names, the
sample list, exclude patterns and the geometry filter live in a
`data_config_local.py` in the **project root**, next to `inks_local.m`. The
shared scripts read it through `data_config.load()`. Copy
`01_Python/data_config_local.template.py` to start. **Do not edit folder
names back into a shared script** — that is a fork, and a fork that differs
only in a folder name is the hardest kind to notice. With no local config
every script resolves to exactly the literals it used before this mechanism
existed, so doing nothing is safe.

**Geometry filtering is not optional bookkeeping.** `data_config.filter_files`
drops files naming a different measuring system and **raises** on a file
naming none. Silently mixing a PP50 parallel-plate run into a CP50 analysis
produces an answer that looks self-consistent and is physically wrong. Set
`GEOMETRY_FILTER` in every project whose filenames carry the geometry.

**Reuse the legacy `.xls` reader — do not re-derive the sheet layout.**
`extract_SAOS_values.read_antpar_xls()` already knows that sheet 0 is
"Details", the data lives on "Ramp - 1", headers are on row 1, units on row 2
and data starts at row 3. Import it for any new one-off audit or verification
script. Re-implementing equivalent logic from scratch has happened at least
once and is pure waste.

**Diagnose the curve before fitting it.** `flow_diagnostics.py` flags the two
artefacts that no fit statistic reveals: a **startup transient** (viscosity
rising at low shear — moved fitted `eta0` by 19–31 % on the CMC series) and
**edge fracture** (stress falling while shear rate rises — the reason one
ink's Power-Law fit sits at R² = 0.505). `Fit_Muitos_Modelos_v5.py` runs both
automatically and, when a transient is found, reports the full-range and
transient-excluded fits **side by side**. Neither function trims data;
silently trimming points is how a fit becomes irreproducible.

**Normal stress / N₁.** `antpar_io.read_flow_curve_with_force()` is the single
entry point for force-augmented exports. It handles the trap that the file is
**latin-1, not the UTF-16 the rest of this module assumes** (assuming UTF-16
does not raise — it silently mojibakes), applies `N1 = 2F/(πR²)`, and tares
with `baseline_correct(..., method="min")`. The rejected
`method="median_low25"` convention **raises with an explanation** so it is not
retried: it gave ~460 Pa where the min convention gives ~802 Pa on the same
run. Filename pairing between audited `.xls` and force `.txt` is trap-laden
and non-systematic — always go through `pair_audited_files()`, never glob.

**Instrument identity travels with the number.** `antpar_io` now extracts the
instrument name from the export header into `FlowCurve.instrument` (and the
two sweep containers). Downstream scripts must echo it into their report
headers. An ARES-G2 dataset sat mislabelled as an MCR in a project CLAUDE.md
for months precisely because nothing ever read that field.

**Cross-version checks.** `validate_fits.py` is the Python counterpart of
`validate_v4.m`: it runs the audited `.xls` fitter and the active CSV fitter
on the same samples and splits any disagreement into a **fitter-attributable**
and an **ingestion-attributable** part. Reporting only "v4 says X, v5 says Y"
leaves you unable to say which changed, and that is the question an editor
asks.

**Before submitting anything that quotes a flow rate**, run
`check_flow_rate_consistency.py`. It compares the solver `Q`, the `Q` read
back out of the G-code that was actually sent to the printer, and the `Q`
implied by the measured deposit (`w·h·v_print`). Three disagreeing values for
one print went unreconciled through a whole manuscript cycle.

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

**Slicer convention (state it exactly this way — a manuscript got it wrong
for three drafts):** three quantities come out of the slicer layer and they
are NOT interchangeable.

```
w_line  = 2*Rn*(1 + beta)                          deposited road width (m)
v_print = Q / (w_line * h_layer)                   head speed (m/s)
k_flow  = (1 + beta)^2 * f_slip * sqrt(Rrec/100)   deposition efficiency (-)
```

- `v_print` is **pure mass conservation of the road cross-section. It carries
  no efficiency term.** Writing `v_print = Q*k_flow/(w*h)` double-counts the
  swell. That exact error survived three manuscript drafts undetected; the
  code was never wrong, the missing documentation was.
- `k_flow` is a **separate slicer input**, applied as
  **Extrusion Multiplier = 1 / k_flow**.

Copying this block into a Methods section cannot get `k_flow`'s placement
wrong. Re-deriving it from memory can.

**Two closures in the solver are in-house heuristics, not literature.** Both
are flagged inline in `bioprinting_algorithm_v4.m` and both are replaced in
`bioprinting_algorithm_v5.m`. If you report a number that depends on either,
disclose it:

| quantity | v4 (heuristic) | v5 (anchored) | what it affects |
|---|---|---|---|
| die swell `beta` | `0.30*(1-n)` — no source; makes swell GROW as the ink gets more shear-thinning, which is backwards, and over-predicts the one measured case by ~2.5x | Tanner (1970) closure from measured `N1`; `N1 -> 0` gives the 0.13 inelastic floor (Nickell 1974) | `w_line`, `v_print`, `k_flow` → **the Extrusion Multiplier used at the bench**. Moving `beta` 0.228 → 0.107 on C15-SF5.5 moves EM 1.102 → 1.356 (+23 %). Touches no pressure or shear output. |
| critical Reynolds | `2100*n^0.75` — no source; collapses monotonically to 222 by n = 0.05 | Ryan & Johnson (1959), non-monotonic, stays in 992–2397 | the laminar verdict string only. `Re_gen ~ 1e-4` here, so the verdict is "laminar" either way. Heuristic under-predicts 3.13x at n = 0.24, 4.31x at n = 0.088. |

**Running the anchored (Tanner) path — the full chain:**

```
1. run_solver_v4                      -> master_summary_v4.csv (gives gamma_w at the wall)
2. extract_N1_tanner.py               -> reads the force export, tares, excludes
     --master-summary ... --needle 21G   edge fracture, interpolates N1 to gamma_w,
     --vp 0.01                           prints a paste-ready inks_local.m line
3. paste  inks(i).N1_wall_Pa = ...    -> into inks_local.m in the project root
4. run_solver_v5                      -> master_summary_v5.csv + closure_report_v5.txt
```

Step 2 **refuses to extrapolate** past the measured shear-rate range. The
needle reaches ~165–218 s⁻¹ at the nominal operating point; a flow curve that
stops at 100 s⁻¹ cannot speak about it, and `--allow-extrapolation` stamps
every number it produces. An ink with no `N1_wall_Pa` still runs — on the
tagged heuristic fallback — so a partial measurement campaign is usable, but
the two closures must not be quoted as if they were the same.

Note that the Tanner `beta` is **needle-dependent** (it goes through `τ_w`),
where the heuristic depends only on `n`. On the C15 smoke run: 0.2134 at 21G
vs 0.1831 at 22G. That is a real effect the heuristic cannot express.

`bioprinting_algorithm_v5.m` does **not** fork the physics — it resolves the
two closures and delegates to v4 through optional override parameters.
`validate_v5.m` asserts that v5 with no `N1` reproduces v4 bit-for-bit.

**Never fork a shared parameter into a project-local script.** `inks_local.m`
in the project root is the single source of truth for `rho`, `K_PL`, `n_PL`,
`eta0`, `etaInf`, `lambda`, `m_Cross` and `Rrec_pct`. A project-local
extrusion or plotting script must call `inks_local()` and read the fields it
needs — it must **never** define its own `rho` or inline `sample` struct.
A local copy does not stay in sync; it silently forks, and the fork is
discovered by a number that no longer matches.

The same rule applies to the solver itself. If a project keeps a private copy
of `bioprinting_algorithm_v4.m` beside its own analysis, that copy will drift.
Use the `Export/` symlink.

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
