"""
extract_SAOS_frequency_v2.py
============================

Frequency-sweep value extraction from **Anton Paar Rheocompass** CSVs.

This is the v2 successor to the frequency-sweep portion of the legacy
``extract_SAOS_values.py``. The extraction logic (ω-target picks +
log-log low-ω power-law fit for extrapolation) is preserved from v1.
Only the data-ingestion layer changed:

    v1 → reads TA-style ``.xls`` via ``xlrd``
    v2 → reads Anton Paar Rheocompass raw ``.csv`` via ``antpar_io``

Per-sample outputs
------------------
* G', G'', tan δ at the three canonical SAOS targets ω ∈ {1, 10, 100} rad/s,
  with the actual data point used (closest to each target in log-ω).
* Low-ω power-law fit  G' = G0 · ω^β  over the n_low lowest valid points.
* Extrapolated G'(ω = 0.01 rad/s) — relevant for long-dwell h_max criteria.

Status filtering
----------------
At low ω the rheometer hits sub-resolution torque (``M-``) and modulus-
estimate errors (``ME-``). These points have G' essentially equal to
numerical noise (1e-5 Pa region), and including them in the low-ω
power-law fit destroys the extrapolation. By default this script
**drops rows whose Status contains M- or ME-** before fitting. The
``taD`` flag (tan-delta drift) is left alone — it usually accompanies
high-ω inertial corrections and the G' value is still meaningful.

Outputs
-------
``./Analises/Python/Results/SAOS_frequency_v2.txt``   human-readable
``./Analises/Python/Results/SAOS_frequency_v2.csv``   machine-readable

USAGE
-----
1. Edit the ``CONFIGURATION`` block below.
2. From the project root:
       python Export/01_Python/extract_SAOS_frequency_v2.py
3. The CSV columns ``Gp_omega1_Pa``, ``beta_lowomega``, etc. are exactly
   the fields ``extract_hmax_v3`` (pending) will consume.

Author: T.M.C. Rodrigues - 2026-05-20
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))
from antpar_io import (   # noqa: E402
    read_frequency_sweeps_in_folder, FrequencySweep,

)
from data_config import load as _load_data_config  # noqa: E402

_CFG = _load_data_config()   # data_config_local.py in the project root,
                             # else the built-in defaults below.


# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION  ← edit here
# ──────────────────────────────────────────────────────────────────────────────

# A12: the folder name is no longer a literal in this shared script.
# It comes from data_config_local.py in the project root (the Python
# counterpart of inks_local.m), so a project with different folder
# names configures instead of forking. With no local config, this
# resolves to exactly the literal that used to be hard-coded here:
#     "./Reologia/Frequency Sweep - 3D Bioprinting"
INPUT_FOLDER = _CFG.folder("frequency")

# Optional include list — None to process every CSV in INPUT_FOLDER.
SAMPLE_INCLUDE = _CFG.sample_include   # A12: set it in data_config_local.py

SAVE_PATH = _CFG.results_dir     # A12: was "./Analises/Python/Results"
TAG = "SAOS_frequency_v2"

# Targets (rad/s) for picking representative G', G'', tan δ.
OMEGA_TARGETS = (1.0, 10.0, 100.0)

# Low-ω power-law fit: take the LOWEST n_low *valid* points and fit
#     log G' = log G0 + β · log ω
N_LOW = 6

# Frequency to extrapolate G' to (rad/s). 0.01 ~ quasi-static / long dwell.
OMEGA_EXTRAP = 0.01

# Drop rows whose Status contains any of these substrings before fitting.
# 'M-'  = torque below sensitivity
# 'ME-' = modulus estimate error (G' collapsed to numerical noise)
# 'taD' = tan-delta drift (high-ω inertial corrections) — kept by default
DROP_STATUS_FLAGS = ("M-", "ME-")

# Also write an "_allpoints" variant of the report and CSV, using EVERY
# measured point (no Status-based filtering). Useful as a diagnostic
# comparison to see how much the M-/ME- exclusion shifts the picks and
# the low-ω power-law fit. Set to False to skip the second pass.
ALSO_WRITE_ALLPOINTS = True


# ──────────────────────────────────────────────────────────────────────────────
# CORE EXTRACTION
# ──────────────────────────────────────────────────────────────────────────────

def _filter_valid(df, drop_flags: tuple[str, ...]):
    """Drop NaN / non-positive points and (optionally) flagged Status rows."""
    mask = (
        np.isfinite(df["omega_rad_s"])
        & np.isfinite(df["G_storage_Pa"])
        & np.isfinite(df["G_loss_Pa"])
        & (df["omega_rad_s"] > 0)
        & (df["G_storage_Pa"] > 0)
        & (df["G_loss_Pa"] > 0)
    )
    for flag in drop_flags:
        mask &= ~df["status"].str.contains(flag, regex=False, na=False)
    return df.loc[mask].sort_values("omega_rad_s").reset_index(drop=True)


def _pick_at_target(omega: np.ndarray, target: float) -> int:
    """Index of the data point closest to ``target`` in log-ω."""
    return int(np.argmin(np.abs(np.log10(omega) - np.log10(target))))


def _fit_low_omega_powerlaw(omega: np.ndarray, gp: np.ndarray, n_low: int):
    """
    Fit  log G' = log G0 + β · log ω  over the lowest ``n_low`` points.

    Returns (G0, beta, R²_log, n_used). If fewer than 3 valid points,
    returns (nan, nan, nan, n_used).
    """
    n = min(n_low, len(omega))
    if n < 3:
        return float("nan"), float("nan"), float("nan"), n
    lw = np.log10(omega[:n])
    lg = np.log10(gp[:n])
    A = np.vstack([lw, np.ones_like(lw)]).T
    coef, *_ = np.linalg.lstsq(A, lg, rcond=None)
    beta, log_G0 = coef
    G0 = 10.0 ** log_G0
    pred = log_G0 + beta * lw
    ss_res = float(np.sum((lg - pred) ** 2))
    ss_tot = float(np.sum((lg - np.mean(lg)) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(G0), float(beta), float(R2), n


def analyse_one(sweep: FrequencySweep,
                drop_flags: tuple[str, ...] = DROP_STATUS_FLAGS) -> dict:
    """Pick ω-target values + low-ω fit for one frequency sweep.

    ``drop_flags``: substrings; rows whose Status contains ANY of them
    are excluded before fitting / picking. Pass an empty tuple ``()`` to
    use every measured point regardless of status (diagnostic mode).
    """
    df = _filter_valid(sweep.data, drop_flags)
    info: dict = {
        "sample": sweep.sample_name.split(" - ", 1)[0].strip(),
        "full_name": sweep.sample_name,
        "source": sweep.source_path.name,
        "n_total": len(sweep.data),
        "n_used": len(df),
        "drop_status_flags": ",".join(drop_flags),
    }
    if len(df) < 3:
        info["error"] = (
            f"Only {len(df)} usable points after Status filtering — "
            "cannot extract."
        )
        return info

    omega = df["omega_rad_s"].to_numpy()
    gp = df["G_storage_Pa"].to_numpy()
    gpp = df["G_loss_Pa"].to_numpy()
    tan_d = df["tan_delta"].to_numpy()

    # ω-target picks
    info["omega_range_used"] = (float(omega[0]), float(omega[-1]))
    picks = {}
    for tgt in OMEGA_TARGETS:
        # Skip targets outside the swept range — would be a wild extrapolation.
        if tgt < omega[0] * 0.5 or tgt > omega[-1] * 2.0:
            picks[tgt] = None
            continue
        i = _pick_at_target(omega, tgt)
        picks[tgt] = {
            "omega_actual": float(omega[i]),
            "Gp": float(gp[i]),
            "Gpp": float(gpp[i]),
            "tan": float(tan_d[i]),
        }
    info["picks"] = picks

    # Low-ω power-law fit + extrapolation
    G0, beta, R2, n_used_fit = _fit_low_omega_powerlaw(omega, gp, N_LOW)
    info["lowomega_G0_Pa"] = G0
    info["lowomega_beta"] = beta
    info["lowomega_R2"] = R2
    info["lowomega_n_used"] = n_used_fit
    info["lowomega_fit_strain_range"] = (
        float(omega[0]), float(omega[min(n_used_fit, len(omega)) - 1])
    ) if n_used_fit >= 1 else None
    info["Gp_extrap_Pa"] = (
        float(G0 * (OMEGA_EXTRAP ** beta)) if np.isfinite(G0) and np.isfinite(beta)
        else float("nan")
    )
    return info


# ──────────────────────────────────────────────────────────────────────────────
# REPORTING
# ──────────────────────────────────────────────────────────────────────────────

def _write_report(infos: list[dict], path: str, drop_flags: tuple[str, ...]):
    flags = ", ".join(drop_flags) or "(none)"
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write(f"  SAOS FREQUENCY-SWEEP EXTRACTION  ({TAG})\n")
        f.write(f"  Anton Paar Rheocompass CSV  via antpar_io.read_frequency_sweep_csv\n")
        if drop_flags:
            f.write(f"  Status flags EXCLUDED before fit/picks: {flags}\n")
        else:
            f.write("  ALL POINTS used — no Status-flag filtering (diagnostic mode).\n")
            f.write("  Compare against the filtered report to see the impact of M-/ME- exclusion.\n")
        f.write(f"  Low-ω power-law: G' = G0 · ω^β over lowest {N_LOW} valid points\n")
        f.write(f"  Extrapolation target: ω = {OMEGA_EXTRAP} rad/s\n")
        f.write(f"  ω targets for picks: {OMEGA_TARGETS}\n")
        f.write("=" * 78 + "\n\n")
        for info in infos:
            f.write("─" * 70 + "\n")
            f.write(f"SAMPLE: {info['sample']}   ({info['full_name']})\n")
            f.write(f"  source : {info['source']}\n")
            f.write(f"  points : {info['n_used']} used / {info['n_total']} total\n")
            f.write("─" * 70 + "\n")
            if "error" in info:
                f.write(f"  [!] {info['error']}\n\n")
                continue

            wlo, whi = info["omega_range_used"]
            f.write(f"  Valid ω range : [{wlo:.3g}, {whi:.3g}] rad/s\n\n")

            gp_h = "G' (Pa)"
            gpp_h = "G'' (Pa)"
            f.write("  ω-target picks:\n")
            f.write(f"    {'target':>8} | {'ω_actual':>10} | {gp_h:>9} | "
                    f"{gpp_h:>9} | {'tan δ':>7}\n")
            f.write("    " + "-" * 60 + "\n")
            for tgt in OMEGA_TARGETS:
                p = info["picks"].get(tgt)
                if p is None:
                    f.write(f"    {tgt:>8.0f} | {'out-of-range':>10} |  -        |  -        | -\n")
                else:
                    f.write(f"    {tgt:>8.0f} | {p['omega_actual']:>10.4g} | "
                            f"{p['Gp']:>9.3f} | {p['Gpp']:>9.3f} | {p['tan']:>7.3f}\n")
            f.write("\n")

            G0, beta, R2 = info["lowomega_G0_Pa"], info["lowomega_beta"], info["lowomega_R2"]
            f.write("  Low-ω power-law fit:\n")
            f.write(f"    G' = {G0:.4g} · ω^{beta:.4f}     "
                    f"(R²_log = {R2:.4f}, n_used = {info['lowomega_n_used']})\n")
            f.write(f"    Extrapolated  G'(ω = {OMEGA_EXTRAP}) = "
                    f"{info['Gp_extrap_Pa']:.4g} Pa\n\n")

        # Summary table
        rows = [i for i in infos if "error" not in i]
        if rows:
            gp1_h = "G'(ω=1)"
            gp10_h = "G'(ω=10)"
            gp100_h = "G'(ω=100)"
            gpex_h = "G'(ω=0.01)"
            f.write("=" * 78 + "\n")
            f.write("  SUMMARY TABLE\n")
            f.write("=" * 78 + "\n")
            f.write(f"  {'Sample':<14} | {gp1_h:>10} | {gp10_h:>10} | "
                    f"{gp100_h:>10} | {'β':>7} | {gpex_h:>12}\n")
            f.write("  " + "-" * 80 + "\n")
            for i in rows:
                def _gp_at(t):
                    p = i["picks"].get(t)
                    return f"{p['Gp']:>10.3f}" if p else f"{'out-of-range':>10}"
                gp_extrap = i["Gp_extrap_Pa"]
                gp_extrap_s = f"{gp_extrap:>12.4g}" if np.isfinite(gp_extrap) else f"{'nan':>12}"
                beta_s = (f"{i['lowomega_beta']:>7.3f}"
                          if np.isfinite(i["lowomega_beta"]) else f"{'nan':>7}")
                f.write(f"  {i['sample']:<14} | {_gp_at(1.0)} | {_gp_at(10.0)} | "
                        f"{_gp_at(100.0)} | {beta_s} | {gp_extrap_s}\n")


def _write_csv(infos: list[dict], path: str):
    rows = [i for i in infos if "error" not in i]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "sample", "full_name", "source",
            "n_total", "n_used",
            "omega_used_lo", "omega_used_hi",
            "Gp_omega1_Pa",   "Gpp_omega1_Pa",   "tan_omega1",   "omega_actual_1",
            "Gp_omega10_Pa",  "Gpp_omega10_Pa",  "tan_omega10",  "omega_actual_10",
            "Gp_omega100_Pa", "Gpp_omega100_Pa", "tan_omega100", "omega_actual_100",
            "lowomega_G0_Pa", "lowomega_beta", "lowomega_R2", "lowomega_n_used",
            "omega_extrap_rad_s", "Gp_extrap_Pa",
            "drop_status_flags",
        ])
        for i in rows:
            def _cells(t):
                p = i["picks"].get(t)
                if p is None:
                    return ["", "", "", ""]
                return [f"{p['Gp']:.6g}", f"{p['Gpp']:.6g}",
                        f"{p['tan']:.6g}", f"{p['omega_actual']:.6g}"]
            wlo, whi = i["omega_range_used"]
            w.writerow([
                i["sample"], i["full_name"], i["source"],
                i["n_total"], i["n_used"],
                f"{wlo:.6g}", f"{whi:.6g}",
                *_cells(1.0),
                *_cells(10.0),
                *_cells(100.0),
                f"{i['lowomega_G0_Pa']:.6g}", f"{i['lowomega_beta']:.6g}",
                f"{i['lowomega_R2']:.6g}", i["lowomega_n_used"],
                f"{OMEGA_EXTRAP}", f"{i['Gp_extrap_Pa']:.6g}",
                i["drop_status_flags"],
            ])


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def _run_one_pass(sweeps: dict, drop_flags: tuple[str, ...],
                  suffix: str) -> tuple[str, str]:
    """Run one analysis pass with the given drop_flags, write outputs."""
    out_txt = os.path.join(SAVE_PATH, f"{TAG}{suffix}.txt")
    out_csv = os.path.join(SAVE_PATH, f"{TAG}{suffix}.csv")
    infos = [analyse_one(s, drop_flags=drop_flags) for s in sweeps.values()]
    _write_report(infos, out_txt, drop_flags=drop_flags)
    _write_csv(infos, out_csv)
    return out_txt, out_csv


def main() -> int:
    os.makedirs(SAVE_PATH, exist_ok=True)

    print(f"Scanning {INPUT_FOLDER} for raw Rheocompass frequency-sweep CSVs ...")
    sweeps = read_frequency_sweeps_in_folder(INPUT_FOLDER)
    if SAMPLE_INCLUDE is not None:
        sweeps = {k: v for k, v in sweeps.items() if k in SAMPLE_INCLUDE}
    if not sweeps:
        print(f"  [!] No raw frequency-sweep CSVs found in {INPUT_FOLDER}.")
        return 1
    print(f"  parsed {len(sweeps)} sample(s): {', '.join(sweeps.keys())}")

    # Pass 1: default (filtered) — what the scientific tables use.
    out_txt, out_csv = _run_one_pass(sweeps, drop_flags=DROP_STATUS_FLAGS,
                                     suffix="")
    print(f"[Saved] {out_txt}")
    print(f"[Saved] {out_csv}")

    # Pass 2: diagnostic (no Status-based filtering).
    if ALSO_WRITE_ALLPOINTS:
        out_txt2, out_csv2 = _run_one_pass(sweeps, drop_flags=(),
                                           suffix="_allpoints")
        print(f"[Saved] {out_txt2}   (no Status filtering — diagnostic)")
        print(f"[Saved] {out_csv2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
