# The Pipeline — End to End

Every ink in this lab goes through the same four stages. **Stages are
sequential, never reordered, never skipped.**

```
┌──────────────────────────────────────────────────────────────────┐
│  [1] RHEOLOGY  ─  Anton Paar rheometer + Rheocompass             │
│      Tests: Flow / viscosity curve, Frequency sweep,             │
│             Amplitude sweep, 3iTT-Osc-Rot-Osc recovery           │
│      Output: UTF-16 .csv files  →  Reologia/{Viscosity, Strain   │
│              Sweep, Frequency Sweep, Recovery}                   │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│  [2] MODELING (Python — all read via antpar_io.py, no xlrd)      │
│                                                                   │
│      Fit_Muitos_Modelos_v5.py       →  K, n  /  η₀, η∞, λ, m       │
│      extract_SAOS_amplitude_v2.py   →  G'_LVR, γ_LVR, σ_y_LVR     │
│      extract_SAOS_frequency_v2.py   →  G'(ω), low-ω fit, β        │
│      extract_recovery_v2.py         →  Rec_G' %, Rec_|η*| %       │
│      extract_hmax_v3.py             →  h_max, N_max (from above)  │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│  [3] EXTRUSION SIMULATION (MATLAB)                                │
│      Paste Python params into the inks(…) struct, then run        │
│      run_solver_v4.m  →  ΔP, γ̇_w, v_print, w_line, k_flow,         │
│                          per (ramp × ink × needle × Vp)           │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│  [4] PRINTING VALIDATION                                          │
│      Use slicer_lookup_*.csv to set machine speeds.               │
│      Slicer Extrusion Multiplier  =  1 / k_flow                   │
│      Print the per-shape benchmark suite                          │
│      (SOP_v4 in 03_SOPs/)                                          │
└──────────────────────────────────────────────────────────────────┘
```

## The data hand-off that matters

Stage 2 → Stage 3 is the one parameter sync that breaks pipelines if it falls
out of date. **Whenever you re-fit in Python, update the MATLAB
`inks(...)` struct on the same day.**

In `run_solver_v4.m` each ink carries `Ramp1` (no pre-shear) and optional
`Ramp2` (post pre-shear) sub-structs. Fields below live on the relevant ramp
sub-struct, except `rho` / `Rrec_pct` which are per-ink.

| Python output | MATLAB field (in `run_solver_v4.m`) |
|---|---|
| `Fit_Muitos_Modelos_v5` → Power-Law `K` | `inks(i).Ramp1.K_PL` |
| `Fit_Muitos_Modelos_v5` → Power-Law `n` | `inks(i).Ramp1.n_PL` |
| `Fit_Muitos_Modelos_v5` → Cross `eta0` | `inks(i).Ramp1.eta0` |
| `Fit_Muitos_Modelos_v5` → Cross `etaInf` | `inks(i).Ramp1.etaInf` |
| `Fit_Muitos_Modelos_v5` → Cross `lambda` | `inks(i).Ramp1.lambda` |
| `Fit_Muitos_Modelos_v5` → Cross `m` | `inks(i).Ramp1.m_Cross` |
| `extract_recovery_v2` → Rec_G' or Rec_\|η*\| % (pick one) | `inks(i).Rrec_pct` |
| Pre-sheared (`Visco_Artur`) re-fit of the same fields | `inks(i).Ramp2.*` |

## Common mistakes (the ones a new student will make)

1. **Running the Python script from `Export/01_Python/` directly.** All scripts
   use relative paths like `./Reologia/Viscosity` and assume the project root
   is the current working directory. Run them with `cd <PROJECT ROOT>` first.
2. **Forgetting to update the MATLAB struct after re-fitting.** The MATLAB
   defaults are hardcoded values for *this project's three inks*. Replace them
   for your new ink or the simulation lies.
3. **Trusting the simulation without checking the Python fit's AIC/BIC.** A
   bad fit cascades into a bad simulation cascades into a useless slicer CSV.
4. **Picking a piston velocity from outside the calibrated sweep
   (`[0.003, 0.04] mm/s`).** Anything above 0.05 mm/s implies a head speed the
   lab DIW rig cannot actually reach.
