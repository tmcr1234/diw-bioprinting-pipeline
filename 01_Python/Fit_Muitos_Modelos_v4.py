"""
Fit_Muitos_Modelos_v4.py
========================

╔══════════════════════════════════════════════════════════════════════════╗
║  ⚠ LEGACY SCRIPT — for the pre-migration TA-style .xls archive only.    ║
║  For data exported by Anton Paar Rheocompass (UTF-16 CSV), use          ║
║  Fit_Muitos_Modelos_v5.py instead. v4 is kept active because the         ║
║  published preprint's fits were produced by this version; re-running     ║
║  through v5 would re-fit on different ingestion paths.                   ║
║  See Export/docs/01_rheology_data_formats.md                             ║
╚══════════════════════════════════════════════════════════════════════════╝

Rheology model fitting with log-space least squares (equivalent to minimizing
relative error), plus AIC/BIC for model selection.

WHY LOG-FIT?
    - Rheological data spans orders of magnitude in both stress and viscosity.
    - OLS in linear space lets high-value points dominate; fitting in log space
      weights each decade equally (equivalent to WLS with w_i = 1/y_i^2).
    - For any physically positive output (stress, viscosity) this is the
      statistically correct baseline approach.

AIC / BIC:
    - Computed from the residuals in LOG space (consistent with the fitting objective).
    - Lower is better. BIC penalises extra parameters more strongly than AIC.
    - Use these to rank models per sample objectively.

DOMAIN CONVENTION:
    - HB / Power Law / Bingham / Casson → fit shear stress vs shear rate.
    - Cross / Carreau / Ellis             → fit viscosity vs shear rate
      (parameters are defined in the viscosity domain; converting to stress
       introduces correlated errors from γ̇ measurement).

USAGE:
    1. Edit the "Configuration" block below (samples, file_path, save_path).
    2. Run:  python Fit_Muitos_Modelos_v4.py
    3. Results → ./Results/FitAll-<sheet>-v4.txt  (human-readable)
                  ./Results/FitAll-<sheet>-v4.csv  (machine-readable, easy for plotting)

CHANGES FROM v3:
    - Cross, Carreau: eta_inf lower bound raised 0.0 → 0.001 Pa·s
      (physical floor = water viscosity at 25 °C; zero is thermodynamically impossible).
    - Casson: bounds [-inf,-inf]/[+inf,+inf] → [0,0]/[+inf,+inf]
      (enforce principal square roots a=sqrt(tau0)>=0, b=sqrt(Kc)>=0;
       eliminates symmetric local minimum at (-a,-b)).
    - AIC best-model summary groups by domain before ranking;
      stress-domain and viscosity-domain AIC values are incommensurable.
"""

import os
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# ──────────────────────────────────────────────────────────────────────────────
# 1.  MODEL DEFINITIONS
# ──────────────────────────────────────────────────────────────────────────────

# ── Stress-based models (τ vs γ̇) ─────────────────────────────────────────────

def herschel_bulkley(gd, tau0, K, n):
    """τ = τ₀ + K·γ̇ⁿ"""
    return tau0 + K * gd**n

def power_law(gd, K, n):
    """τ = K·γ̇ⁿ"""
    return K * gd**n

def bingham(gd, tau0, mu_p):
    """τ = τ₀ + μ_p·γ̇"""
    return tau0 + mu_p * gd

def casson(gd, a, b):
    """τ = (a + b·√γ̇)²  with a=√τ₀, b=√Kc"""
    return (a + b * np.sqrt(np.maximum(gd, 0.0)))**2

# ── Viscosity-based models (η vs γ̇) ──────────────────────────────────────────

def cross(gd, eta0, eta_inf, K, m):
    """η = η∞ + (η₀–η∞) / (1 + (K·γ̇)ᵐ)"""
    return eta_inf + (eta0 - eta_inf) / (1.0 + (K * gd)**m)

def carreau(gd, eta0, eta_inf, lam, n):
    """η = η∞ + (η₀–η∞)·[1+(λγ̇)²]^((n-1)/2)"""
    return eta_inf + (eta0 - eta_inf) * (1.0 + (lam * gd)**2)**((n - 1.0) / 2.0)

def ellis(gd, eta0, K_E, alpha):
    """η = η₀ / [1+(K_E·γ̇)^(α-1)]"""
    return eta0 / (1.0 + (K_E * gd)**(alpha - 1.0))


# ──────────────────────────────────────────────────────────────────────────────
# 2.  INITIAL GUESSES & BOUNDS
# ──────────────────────────────────────────────────────────────────────────────

MODEL_STRESS = {
    # name: (function, p0, bounds, param_names)
    "HB": (
        herschel_bulkley,
        [0.1, 1.0, 0.8],
        ([0.0,  0.0,  0.0], [np.inf, np.inf, 2.0]),
        ["tau0", "K", "n"],
    ),
    "PowerLaw": (
        power_law,
        [1.0, 0.8],
        ([0.0,  0.0], [np.inf, 2.0]),
        ["K", "n"],
    ),
    "Bingham": (
        bingham,
        [0.1, 1.0],
        ([0.0,  0.0], [np.inf, np.inf]),
        ["tau0", "mu_p"],
    ),
    "Casson": (
        casson,
        [0.1, 0.5],
        ([0.0, 0.0], [np.inf, np.inf]),   # principal roots: a=sqrt(tau0)>=0, b=sqrt(Kc)>=0
        ["a_sqrt_tau0", "b_sqrt_Kc"],
    ),
}

MODEL_VISCOSITY = {
    "Cross": (
        cross,
        [10.0, 0.1, 1.0, 1.0],
        ([0.0, 0.001, 0.0, 0.0], [np.inf, np.inf, np.inf, 10.0]),
        ["eta0", "eta_inf", "K", "m"],
    ),
    "Carreau": (
        carreau,
        [10.0, 0.1, 1.0, 0.6],
        ([0.0, 0.001, 0.0, 0.0], [np.inf, np.inf, np.inf,  2.0]),
        ["eta0", "eta_inf", "lambda", "n"],
    ),
    "Ellis": (
        ellis,
        [10.0, 0.1, 1.5],
        ([0.0, 0.0,  1.0], [np.inf, np.inf, 10.0]),
        ["eta0", "K_E", "alpha"],
    ),
}


# ──────────────────────────────────────────────────────────────────────────────
# 3.  STATISTICS
# ──────────────────────────────────────────────────────────────────────────────

def r_squared_log(y_true, y_pred):
    """R² computed in log space (consistent with log-fit objective)."""
    ly = np.log(y_true)
    lyh = np.log(np.maximum(y_pred, 1e-300))
    ss_res = np.sum((ly - lyh) ** 2)
    ss_tot = np.sum((ly - np.mean(ly)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan


def aic_bic_log(y_true, y_pred, n_params):
    """
    AIC and BIC from log-space residuals.

    Assumes normal errors in log space (log-normal errors in linear space).
    n_params = number of fitted parameters (NOT counting σ²).
    """
    n = len(y_true)
    residuals = np.log(y_true) - np.log(np.maximum(y_pred, 1e-300))
    sse_log = np.sum(residuals**2)
    # MLE variance estimate
    sigma2 = sse_log / n
    # Log-likelihood under Gaussian assumption in log space
    log_lik = -n / 2.0 * np.log(2 * np.pi * sigma2) - sse_log / (2.0 * sigma2)
    k = n_params + 1   # +1 for σ²
    aic = 2 * k - 2 * log_lik
    bic = k * np.log(n) - 2 * log_lik
    return aic, bic


# ──────────────────────────────────────────────────────────────────────────────
# 4.  LOG-SPACE FITTING WRAPPER
# ──────────────────────────────────────────────────────────────────────────────

def fit_log(model_fn, X, Y, p0, bounds):
    """
    Fit model_fn to (X, Y) by minimising sum of squared LOG residuals.

    Equivalent to WLS with w_i = 1/y_i², which gives equal weight to each
    decade of variation — the physically appropriate weighting for rheology.

    Returns
    -------
    params : array   Best-fit parameters.
    pcov   : 2-D array  Parameter covariance (in log space).
    """
    log_Y = np.log(Y)   # requires Y > 0 (enforced by caller)

    def log_model(x, *p):
        pred = model_fn(x, *p)
        # Clamp to avoid log(0) or log(negative) during iteration
        pred = np.where(np.isfinite(pred) & (pred > 0), pred, 1e-300)
        return np.log(pred)

    params, pcov = curve_fit(
        log_model, X, log_Y,
        p0=p0, bounds=bounds,
        maxfev=50000,
        method="trf",   # Trust Region Reflective — handles bounds correctly
    )
    return params, pcov


# ──────────────────────────────────────────────────────────────────────────────
# 5.  CONFIGURATION  ← edit here
# ──────────────────────────────────────────────────────────────────────────────

SAMPLES = ['C15', 'C15 Gira 5.5', 'Gel de Cabelo Bozano']

FILE_PATH = "./Reologia/Viscosity"
FILE_SUFFIX = " - Visco_Artur CP50 0,1mm.xls"

SAVE_PATH = "./Analises/Python/Results"
SHEET_PREFIX = "Ramp - 2"

_ramp_tag = SHEET_PREFIX.replace(" ", "").replace("-", "")  # e.g. "Ramp2"
OUTPUT_TXT = os.path.join(SAVE_PATH, f"FitAll-{_ramp_tag}-v4.txt")
OUTPUT_CSV = os.path.join(SAVE_PATH, f"FitAll-{_ramp_tag}-v4.csv")

# Column indices in your .xls sheets (0-based, after skiprows=3)
COL_STRESS    = 2
COL_SHEAR_RATE = 3
COL_VISCOSITY  = 4


# ──────────────────────────────────────────────────────────────────────────────
# 6.  MAIN LOOP
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(SAVE_PATH, exist_ok=True)

csv_rows = []   # collect for DataFrame at end

def _safe_val(v):
    return f"{v:.6g}" if np.isfinite(v) else str(v)

with open(OUTPUT_TXT, "w", encoding="utf-8") as fout:

    fout.write("=" * 70 + "\n")
    fout.write("  RHEOLOGY LOG-FIT RESULTS  (v4 — log-space OLS + physical bounds)\n")
    fout.write("=" * 70 + "\n\n")

    for sample in SAMPLES:
        file_name = f"{sample}{FILE_SUFFIX}"
        full_path = os.path.join(FILE_PATH, file_name)

        fout.write(f"{'─'*60}\n")
        fout.write(f"SAMPLE: {sample}\n")
        fout.write(f"{'─'*60}\n")

        if not os.path.exists(full_path):
            msg = f"  [!] File not found: {file_name}\n\n"
            fout.write(msg)
            print(msg.strip())
            continue

        try:
            xls = pd.ExcelFile(full_path, engine="xlrd")
        except Exception as e:
            fout.write(f"  [!] Could not open file: {e}\n\n")
            continue

        for sheet_name in xls.sheet_names:
            if not sheet_name.startswith(SHEET_PREFIX):
                continue

            fout.write(f"\n  Sheet: {sheet_name}\n")

            try:
                df = pd.read_excel(
                    full_path, sheet_name=sheet_name,
                    header=None, skiprows=3, engine="xlrd"
                )
            except Exception as e:
                fout.write(f"    [!] Read error: {e}\n")
                continue

            if df.empty:
                fout.write("    [!] Empty sheet.\n")
                continue

            # ── Extract columns ───────────────────────────────────────────
            try:
                shear_rate = df.iloc[:, COL_SHEAR_RATE].to_numpy(dtype=float)
                shear_stress = df.iloc[:, COL_STRESS].to_numpy(dtype=float)
                viscosity = df.iloc[:, COL_VISCOSITY].to_numpy(dtype=float)
            except Exception as e:
                fout.write(f"    [!] Column extraction error: {e}\n")
                continue

            # ── Clean data ────────────────────────────────────────────────
            # Log-fit requires strictly positive Y; also require γ̇ > 0.
            mask_stress = (
                np.isfinite(shear_rate) & np.isfinite(shear_stress)
                & (shear_rate > 0) & (shear_stress > 0)
            )
            mask_visc = (
                np.isfinite(shear_rate) & np.isfinite(viscosity)
                & (shear_rate > 0) & (viscosity > 0)
            )

            Xs = shear_rate[mask_stress];  Ys = shear_stress[mask_stress]
            Xv = shear_rate[mask_visc];    Yv = viscosity[mask_visc]

            # Casson also needs non-negative sqrt argument (already γ̇>0 here)
            fout.write(f"    Data points (stress-fit): {len(Xs)}  |  (visc-fit): {len(Xv)}\n\n")

            if len(Xs) < 4:
                fout.write("    [!] Too few stress points for fitting.\n")
            if len(Xv) < 4:
                fout.write("    [!] Too few viscosity points for fitting.\n")

            # ── Stress models ─────────────────────────────────────────────
            fout.write("    [Stress-domain fits]\n")

            for model_name, (fn, p0, bounds, pnames) in MODEL_STRESS.items():
                if len(Xs) < len(p0) + 1:
                    fout.write(f"      {model_name}: skipped (not enough points)\n")
                    continue
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        params, _ = fit_log(fn, Xs, Ys, p0, bounds)

                    yhat = fn(Xs, *params)
                    r2 = r_squared_log(Ys, yhat)
                    aic, bic = aic_bic_log(Ys, yhat, len(params))

                    param_str = ", ".join(
                        f"{name}={_safe_val(v)}"
                        for name, v in zip(pnames, params)
                    )
                    fout.write(
                        f"      {model_name:12s}  R²={r2:.5f}  AIC={aic:.2f}  BIC={bic:.2f}\n"
                        f"        params: {param_str}\n"
                    )

                    # Casson extra: back-calculate τ₀ and Kc
                    if model_name == "Casson":
                        a, b = params
                        fout.write(f"        derived: tau0={a**2:.6g}, Kc={b**2:.6g}\n")

                    csv_rows.append({
                        "sample": sample, "sheet": sheet_name,
                        "model": model_name, "domain": "stress",
                        "R2_log": r2, "AIC": aic, "BIC": bic,
                        **{f"p_{n}": v for n, v in zip(pnames, params)},
                    })

                except Exception as e:
                    fout.write(f"      {model_name}: FAILED — {e}\n")

            # ── Viscosity models ──────────────────────────────────────────
            fout.write("\n    [Viscosity-domain fits]\n")

            for model_name, (fn, p0, bounds, pnames) in MODEL_VISCOSITY.items():
                if len(Xv) < len(p0) + 1:
                    fout.write(f"      {model_name}: skipped (not enough points)\n")
                    continue
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        params, _ = fit_log(fn, Xv, Yv, p0, bounds)

                    yhat = fn(Xv, *params)
                    r2 = r_squared_log(Yv, yhat)
                    aic, bic = aic_bic_log(Yv, yhat, len(params))

                    param_str = ", ".join(
                        f"{name}={_safe_val(v)}"
                        for name, v in zip(pnames, params)
                    )
                    fout.write(
                        f"      {model_name:12s}  R²={r2:.5f}  AIC={aic:.2f}  BIC={bic:.2f}\n"
                        f"        params: {param_str}\n"
                    )

                    csv_rows.append({
                        "sample": sample, "sheet": sheet_name,
                        "model": model_name, "domain": "viscosity",
                        "R2_log": r2, "AIC": aic, "BIC": bic,
                        **{f"p_{n}": v for n, v in zip(pnames, params)},
                    })

                except Exception as e:
                    fout.write(f"      {model_name}: FAILED — {e}\n")

        fout.write("\n")

    # ── Summary table: best model per sample ──────────────────────────────
    fout.write("\n" + "=" * 70 + "\n")
    fout.write("  BEST MODEL PER SAMPLE (lowest AIC)\n")
    fout.write("=" * 70 + "\n")

    if csv_rows:
        df_results = pd.DataFrame(csv_rows)
        df_results.to_csv(OUTPUT_CSV, index=False)

        for sample, grp in df_results.groupby("sample"):
            for domain, dgrp in grp.groupby("domain"):
                best = dgrp.loc[dgrp["AIC"].idxmin()]
                fout.write(
                    f"  {sample:<28s} [{domain:9s}] → {best['model']:12s}"
                    f"  (R²={best['R2_log']:.4f}, AIC={best['AIC']:.1f}, BIC={best['BIC']:.1f})\n"
                )

fout_closed = True
print(f"Results saved to:\n  {OUTPUT_TXT}\n  {OUTPUT_CSV}")
