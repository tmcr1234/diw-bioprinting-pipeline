"""
Fit_Muitos_Modelos_v5.py
========================

Rheological model fitting from **Anton Paar Rheocompass CSV exports**.

This is the v5 successor to ``Fit_Muitos_Modelos_v4.py``. The fitting math
(log-space OLS, AIC/BIC for model selection, stress-vs-viscosity domain
separation, physical bounds) is **byte-for-byte identical** to v4. The
only thing that changed is the data-ingestion layer:

    v4 → reads legacy TA-Instruments-style ``.xls`` via ``xlrd``
    v5 → reads Anton Paar Rheocompass raw ``.csv`` via ``antpar_io``

Use v5 for any data collected after the lab moved to the Rheocompass workflow.
Use v4 only to re-run audited fits on the pre-migration ``.xls`` archive.

Models fitted (same as v4)
--------------------------
Stress domain      (fit τ vs γ̇):  Herschel-Bulkley, Power Law, Bingham, Casson
Viscosity domain   (fit η vs γ̇):  Cross, Carreau, Ellis

Outputs
-------
``./Analises/Python/Results/FitAll-AntPar-v5.txt``  human-readable
``./Analises/Python/Results/FitAll-AntPar-v5.csv``  machine-readable

USAGE
-----
1. Edit the ``CONFIGURATION`` block below (input folder, save path,
   optionally a sample include/exclude list).
2. Run from the project root:
       python Export/01_Python/Fit_Muitos_Modelos_v5.py
3. Inspect the .txt report. Best-AIC rows feed directly into the MATLAB
   ``inks(...)`` struct in ``02_MATLAB/run_solver_v4.m``.

Author: T.M.C. Rodrigues - 2026-05-20
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Load the Anton Paar CSV parser from the same folder
HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))
from antpar_io import read_flow_curves_in_folder, FlowCurve  # noqa: E402


# ──────────────────────────────────────────────────────────────────────────────
# 1.  MODEL DEFINITIONS  (unchanged from v4)
# ──────────────────────────────────────────────────────────────────────────────

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
# 2.  INITIAL GUESSES & BOUNDS  (unchanged from v4)
# ──────────────────────────────────────────────────────────────────────────────

MODEL_STRESS = {
    "HB":       (herschel_bulkley, [0.1, 1.0, 0.8],
                 ([0.0, 0.0, 0.0], [np.inf, np.inf, 2.0]),
                 ["tau0", "K", "n"]),
    "PowerLaw": (power_law, [1.0, 0.8],
                 ([0.0, 0.0], [np.inf, 2.0]),
                 ["K", "n"]),
    "Bingham":  (bingham, [0.1, 1.0],
                 ([0.0, 0.0], [np.inf, np.inf]),
                 ["tau0", "mu_p"]),
    "Casson":   (casson, [0.1, 0.5],
                 ([0.0, 0.0], [np.inf, np.inf]),
                 ["a_sqrt_tau0", "b_sqrt_Kc"]),
}

MODEL_VISCOSITY = {
    "Cross":   (cross, [10.0, 0.1, 1.0, 1.0],
                ([0.0, 0.001, 0.0, 0.0], [np.inf, np.inf, np.inf, 10.0]),
                ["eta0", "eta_inf", "K", "m"]),
    "Carreau": (carreau, [10.0, 0.1, 1.0, 0.6],
                ([0.0, 0.001, 0.0, 0.0], [np.inf, np.inf, np.inf, 2.0]),
                ["eta0", "eta_inf", "lambda", "n"]),
    "Ellis":   (ellis, [10.0, 0.1, 1.5],
                ([0.0, 0.0, 1.0], [np.inf, np.inf, 10.0]),
                ["eta0", "K_E", "alpha"]),
}


# ──────────────────────────────────────────────────────────────────────────────
# 3.  STATISTICS  (unchanged from v4)
# ──────────────────────────────────────────────────────────────────────────────

def r_squared_log(y_true, y_pred):
    ly = np.log(y_true)
    lyh = np.log(np.maximum(y_pred, 1e-300))
    ss_res = np.sum((ly - lyh) ** 2)
    ss_tot = np.sum((ly - np.mean(ly)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan


def aic_bic_log(y_true, y_pred, n_params):
    n = len(y_true)
    residuals = np.log(y_true) - np.log(np.maximum(y_pred, 1e-300))
    sse_log = np.sum(residuals**2)
    sigma2 = sse_log / n
    log_lik = -n / 2.0 * np.log(2 * np.pi * sigma2) - sse_log / (2.0 * sigma2)
    k = n_params + 1
    aic = 2 * k - 2 * log_lik
    bic = k * np.log(n) - 2 * log_lik
    return aic, bic


# ──────────────────────────────────────────────────────────────────────────────
# 4.  LOG-SPACE FITTING WRAPPER  (unchanged from v4)
# ──────────────────────────────────────────────────────────────────────────────

def fit_log(model_fn, X, Y, p0, bounds):
    """Minimise Σ(log ŷ − log y)² ; equivalent to WLS with w = 1/y²."""
    log_Y = np.log(Y)

    def log_model(x, *p):
        pred = model_fn(x, *p)
        pred = np.where(np.isfinite(pred) & (pred > 0), pred, 1e-300)
        return np.log(pred)

    params, pcov = curve_fit(
        log_model, X, log_Y,
        p0=p0, bounds=bounds, maxfev=50000, method="trf",
    )
    return params, pcov


# ──────────────────────────────────────────────────────────────────────────────
# 5.  CONFIGURATION  ← edit here
# ──────────────────────────────────────────────────────────────────────────────

# Folder containing raw Rheocompass flow-curve CSVs (the _Analysis.csv
# files are skipped automatically by antpar_io.list_flow_curves).
INPUT_FOLDER = "./Reologia/Viscosity"

# Optional: restrict to a subset of samples. Leave as None to fit every
# raw CSV in INPUT_FOLDER. Keys must match the short filename prefix
# (text before " - "), e.g. ["AMOSTRA A", "AMOSTRA C"].
SAMPLE_INCLUDE: list[str] | None = None

# Where to write the report and CSV.
SAVE_PATH = "./Analises/Python/Results"
TAG = "AntPar-v5"   # appears in output filenames


# ──────────────────────────────────────────────────────────────────────────────
# 6.  FIT ONE SAMPLE
# ──────────────────────────────────────────────────────────────────────────────

def _safe_val(v):
    return f"{v:.6g}" if np.isfinite(v) else str(v)


def fit_one_sample(short_name: str, fc: FlowCurve, fout, csv_rows: list):
    """Fit all 7 models for one FlowCurve, write to fout and csv_rows."""
    df = fc.data

    # Drop non-positive points (log-fit can't tolerate them) and any NaNs.
    gd = df["shear_rate_1_s"].to_numpy()
    tau = df["shear_stress_Pa"].to_numpy()
    eta = df["viscosity_Pa_s"].to_numpy()
    mask = (gd > 0) & (tau > 0) & (eta > 0) & np.isfinite(gd) & np.isfinite(tau) & np.isfinite(eta)
    gd, tau, eta = gd[mask], tau[mask], eta[mask]

    fout.write(f"{'─'*70}\n")
    fout.write(f"SAMPLE: {short_name}   ({fc.sample_name})\n")
    fout.write(f"  source : {fc.source_path.name}\n")
    fout.write(f"  n_pts  : {len(gd)} (after positivity filter)\n")
    fout.write(f"  γ̇ range: [{gd.min():.3g}, {gd.max():.3g}] 1/s\n")
    fout.write(f"  η range: [{eta.min():.3g}, {eta.max():.3g}] Pa·s\n")
    fout.write(f"  T̄      : {df['temperature_C'].mean():.2f} °C\n")
    fout.write(f"{'─'*70}\n\n")

    if len(gd) < 4:
        fout.write("  [!] Too few positive data points to fit. Skipping.\n\n")
        return

    # ---- Stress-domain fits ----
    fout.write("  STRESS-DOMAIN MODELS  (fit τ vs γ̇)\n")
    fout.write(f"  {'model':<10} | {'params':<55} | {'R²(log)':>8} | {'AIC':>9} | {'BIC':>9}\n")
    fout.write("  " + "-" * 100 + "\n")
    stress_results = {}
    for name, (fn, p0, bnds, pnames) in MODEL_STRESS.items():
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                params, _ = fit_log(fn, gd, tau, p0, bnds)
            tau_hat = fn(gd, *params)
            r2 = r_squared_log(tau, tau_hat)
            aic, bic = aic_bic_log(tau, tau_hat, len(params))
            pstr = ", ".join(f"{k}={_safe_val(v)}" for k, v in zip(pnames, params))
            fout.write(f"  {name:<10} | {pstr:<55} | {r2:>8.4f} | {aic:>9.2f} | {bic:>9.2f}\n")
            stress_results[name] = dict(params=params, pnames=pnames,
                                        r2=r2, aic=aic, bic=bic, domain="stress")
        except Exception as e:
            fout.write(f"  {name:<10} | FAILED: {e}\n")

    # ---- Viscosity-domain fits ----
    fout.write("\n  VISCOSITY-DOMAIN MODELS  (fit η vs γ̇)\n")
    fout.write(f"  {'model':<10} | {'params':<55} | {'R²(log)':>8} | {'AIC':>9} | {'BIC':>9}\n")
    fout.write("  " + "-" * 100 + "\n")
    visc_results = {}
    for name, (fn, p0, bnds, pnames) in MODEL_VISCOSITY.items():
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                params, _ = fit_log(fn, gd, eta, p0, bnds)
            eta_hat = fn(gd, *params)
            r2 = r_squared_log(eta, eta_hat)
            aic, bic = aic_bic_log(eta, eta_hat, len(params))
            pstr = ", ".join(f"{k}={_safe_val(v)}" for k, v in zip(pnames, params))
            fout.write(f"  {name:<10} | {pstr:<55} | {r2:>8.4f} | {aic:>9.2f} | {bic:>9.2f}\n")
            visc_results[name] = dict(params=params, pnames=pnames,
                                      r2=r2, aic=aic, bic=bic, domain="viscosity")
        except Exception as e:
            fout.write(f"  {name:<10} | FAILED: {e}\n")

    # ---- Best per domain ----
    fout.write("\n  BEST FIT BY AIC (within each domain — DO NOT cross-compare across domains)\n")
    if stress_results:
        best_s = min(stress_results.items(), key=lambda kv: kv[1]["aic"])
        fout.write(f"    stress    →  {best_s[0]}   (AIC={best_s[1]['aic']:.2f})\n")
    if visc_results:
        best_v = min(visc_results.items(), key=lambda kv: kv[1]["aic"])
        fout.write(f"    viscosity →  {best_v[0]}   (AIC={best_v[1]['aic']:.2f})\n")
    fout.write("\n")

    # ---- Accumulate CSV rows ----
    for name, r in {**stress_results, **visc_results}.items():
        row = {"sample": short_name, "model": name, "domain": r["domain"],
               "R2_log": r["r2"], "AIC": r["aic"], "BIC": r["bic"]}
        for pname, pval in zip(r["pnames"], r["params"]):
            row[pname] = pval
        csv_rows.append(row)


# ──────────────────────────────────────────────────────────────────────────────
# 7.  MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main() -> int:
    os.makedirs(SAVE_PATH, exist_ok=True)
    out_txt = os.path.join(SAVE_PATH, f"FitAll-{TAG}.txt")
    out_csv = os.path.join(SAVE_PATH, f"FitAll-{TAG}.csv")

    print(f"Scanning {INPUT_FOLDER} for raw Rheocompass CSVs ...")
    curves = read_flow_curves_in_folder(INPUT_FOLDER)
    if SAMPLE_INCLUDE is not None:
        curves = {k: v for k, v in curves.items() if k in SAMPLE_INCLUDE}
    if not curves:
        print(f"  [!] No raw flow-curve CSVs found in {INPUT_FOLDER}.")
        return 1
    print(f"  parsed {len(curves)} sample(s): {', '.join(curves.keys())}")

    csv_rows: list[dict] = []
    with open(out_txt, "w", encoding="utf-8") as fout:
        fout.write("=" * 78 + "\n")
        fout.write(f"  RHEOLOGY LOG-FIT RESULTS  ({TAG})  -  Anton Paar Rheocompass CSV\n")
        fout.write("=" * 78 + "\n\n")
        for short_name, fc in curves.items():
            fit_one_sample(short_name, fc, fout, csv_rows)

    pd.DataFrame(csv_rows).to_csv(out_csv, index=False)
    print(f"[Saved] {out_txt}")
    print(f"[Saved] {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
