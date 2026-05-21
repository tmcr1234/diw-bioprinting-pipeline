"""
extract_SAOS_amplitude_v2.py
============================

LVR extraction from **Anton Paar Rheocompass** amplitude (strain) sweeps.

This is the v2 successor to the amplitude-sweep portion of
``extract_SAOS_values.py``. The LVR detection logic (sliding-window CV
plateau, log-interpolated γ_LVR at the 10% G' drop) is preserved from v1.
Only the data-ingestion layer changed:

    v1 → reads TA-style ``.xls`` via ``xlrd`` (multiple test types in one module)
    v2 → reads Anton Paar Rheocompass raw ``.csv`` via ``antpar_io``
         (one module per test type)

Per-sample outputs
------------------
* G'_LVR, G''_LVR, tan δ_LVR (means over the plateau)
* γ_LVR (%) — strain at which G' first falls below 90 % of plateau
* σ_y_LVR (Pa) = G'_LVR × γ_LVR — yield-stress estimate (LVR endpoint)
* Plateau window indices (so you can audit the picked plateau)

Status filtering
----------------
Anton Paar flags points where motor torque is below sensitivity with
``M-`` in the Status column. Those points are sub-resolution noise and
distort plateau detection. By default this script **excludes M- points**
before fitting; set ``DROP_M_MINUS = False`` to keep them.

Outputs
-------
``./Analises/Python/Results/SAOS_amplitude_v2.txt``  human-readable
``./Analises/Python/Results/SAOS_amplitude_v2.csv``  machine-readable

USAGE
-----
1. Edit the ``CONFIGURATION`` block below.
2. From the project root:
       python Export/01_Python/extract_SAOS_amplitude_v2.py
3. The CSV columns ``Gp_LVR_Pa`` and ``gamma_LVR_pct`` feed the h_max
   calculation. Once the frequency-sweep port lands, ``extract_hmax_v3``
   will pull from this CSV directly.

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
    read_amplitude_sweeps_in_folder, AmplitudeSweep,
)


# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION  ← edit here
# ──────────────────────────────────────────────────────────────────────────────

INPUT_FOLDER = "./Reologia/Strain Sweep - 3D Bioprinting"

# Optional include list — None to process every CSV in INPUT_FOLDER.
SAMPLE_INCLUDE: list[str] | None = None

SAVE_PATH = "./Analises/Python/Results"
TAG = "SAOS_amplitude_v2"

# Plateau-detection parameters (identical defaults to the legacy v1 logic).
CV_WINDOW = 5            # sliding window for local CV(G')
CV_THRESHOLD = 0.05      # 5 % local coefficient of variation = plateau
MIN_PLATEAU_LEN = 3      # require at least this many consecutive stable points
G_DROP_FRACTION = 0.90   # γ_LVR = first strain where G' < 0.90 × plateau

# Filter Anton Paar sub-resolution flag ('M-') before LVR detection.
DROP_M_MINUS = True

# Also write an "_allpoints" variant of the report and CSV, using EVERY
# measured point (no Status-based filtering). Useful as a diagnostic
# comparison to see how much the M- exclusion shifts the LVR result.
# Set to False to skip the second pass.
ALSO_WRITE_ALLPOINTS = True


# ──────────────────────────────────────────────────────────────────────────────
# LVR DETECTION
# ──────────────────────────────────────────────────────────────────────────────

def _filter_valid(df, drop_m_minus: bool):
    """Drop NaN / non-positive points and (optionally) M- status rows."""
    mask = (
        np.isfinite(df["strain_frac"])
        & np.isfinite(df["G_storage_Pa"])
        & np.isfinite(df["G_loss_Pa"])
        & (df["strain_frac"] > 0)
        & (df["G_storage_Pa"] > 0)
        & (df["G_loss_Pa"] > 0)
    )
    if drop_m_minus:
        # Status field is free-form text; substring match is the safe check.
        mask &= ~df["status"].str.contains("M-", regex=False, na=False)
    return df.loc[mask].sort_values("strain_frac").reset_index(drop=True)


def _longest_stable_run(gp: np.ndarray, win: int, cv_threshold: float):
    """Return (start, stop_inclusive) index range of the longest run of
    points where the local CV(G') over a window of ``win`` is below
    ``cv_threshold``. Returns ``None`` if no such run exists.
    """
    n = len(gp)
    if n < win:
        return None
    cv = np.full(n, np.inf)
    for k in range(n - win + 1):
        seg = gp[k:k + win]
        m = float(np.mean(seg))
        if m > 0:
            cv[k + win // 2] = np.std(seg) / m
    stable = cv < cv_threshold
    best = None
    cur_start = None
    for k in range(n):
        if stable[k]:
            if cur_start is None:
                cur_start = k
            if best is None or k - cur_start > best[1] - best[0]:
                best = (cur_start, k)
        else:
            cur_start = None
    return best


def _gamma_lvr_log_interp(strain_frac, gp, plateau_end_idx, plateau_val, drop_frac):
    """Strain at which G' first drops below ``drop_frac × plateau``,
    interpolated in log-strain space (rheologists do it this way to
    average out the geometric sweep step)."""
    threshold = drop_frac * plateau_val
    last = plateau_end_idx
    for k in range(last + 1, len(gp)):
        if gp[k] < threshold:
            x0, x1 = np.log10(strain_frac[k - 1]), np.log10(strain_frac[k])
            y0, y1 = gp[k - 1], gp[k]
            if y1 == y0:
                return float(strain_frac[k])
            frac = (threshold - y0) / (y1 - y0)
            return float(10 ** (x0 + frac * (x1 - x0)))
    # No drop detected within the swept range — return the last point as a lower bound.
    return float(strain_frac[-1])


def analyse_one(sweep: AmplitudeSweep, drop_m_minus: bool = DROP_M_MINUS) -> dict:
    """Run LVR extraction on a single AmplitudeSweep.

    ``drop_m_minus``: when True, the Anton Paar ``M-`` (sub-resolution
    torque) status rows are dropped before plateau detection; when False,
    every measured point is used regardless of status flag.
    """
    df = _filter_valid(sweep.data, drop_m_minus)
    info: dict = {
        "sample": sweep.sample_name.split(" - ", 1)[0].strip(),
        "full_name": sweep.sample_name,
        "source": sweep.source_path.name,
        "n_total": len(sweep.data),
        "n_used": len(df),
        "drop_M_minus": drop_m_minus,
    }
    if len(df) < max(CV_WINDOW, MIN_PLATEAU_LEN + 2):
        info["error"] = (
            f"Only {len(df)} usable points after filtering — too few for "
            "robust plateau detection."
        )
        return info

    strain = df["strain_frac"].to_numpy()
    gp = df["G_storage_Pa"].to_numpy()
    gpp = df["G_loss_Pa"].to_numpy()
    tan_d = df["tan_delta"].to_numpy()

    # Primary detection: sliding-window CV.
    best = _longest_stable_run(gp, CV_WINDOW, CV_THRESHOLD)
    if best and (best[1] - best[0] + 1) >= MIN_PLATEAU_LEN:
        sel = np.arange(best[0], best[1] + 1)
        plateau_method = f"CV<{CV_THRESHOLD:.0%} window={CV_WINDOW}"
    else:
        # Fallback: median ± 15% over the lower half.
        n_low = max(5, len(gp) // 2)
        median_low = np.median(gp[:n_low])
        in_plateau = np.abs(gp[:n_low] - median_low) <= 0.15 * median_low
        sel = np.where(in_plateau)[0]
        if len(sel) < MIN_PLATEAU_LEN:
            sel = np.arange(min(5, len(gp)))
        plateau_method = "fallback: median±15% (lower half)"

    Gp_plateau = float(np.mean(gp[sel]))
    Gpp_plateau = float(np.mean(gpp[sel]))
    tan_plateau = float(np.mean(tan_d[sel]))
    gamma_lvr_frac = _gamma_lvr_log_interp(
        strain, gp, plateau_end_idx=sel[-1],
        plateau_val=Gp_plateau, drop_frac=G_DROP_FRACTION,
    )

    info.update({
        "plateau_method": plateau_method,
        "plateau_indices": [int(i) for i in sel],
        "plateau_strain_range_pct": (float(strain[sel[0]] * 100),
                                     float(strain[sel[-1]] * 100)),
        "Gp_LVR_Pa": Gp_plateau,
        "Gpp_LVR_Pa": Gpp_plateau,
        "tan_LVR": tan_plateau,
        "gamma_LVR_pct": gamma_lvr_frac * 100.0,
        "sigma_y_LVR_Pa": Gp_plateau * gamma_lvr_frac,
        "n_plateau_points": int(len(sel)),
    })
    return info


# ──────────────────────────────────────────────────────────────────────────────
# REPORTING
# ──────────────────────────────────────────────────────────────────────────────

def _write_report(infos: list[dict], path: str, drop_m_minus: bool):
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write(f"  SAOS AMPLITUDE-SWEEP LVR EXTRACTION  ({TAG})\n")
        f.write(f"  Anton Paar Rheocompass CSV  via antpar_io.read_amplitude_sweep_csv\n")
        if drop_m_minus:
            f.write("  M- (sub-resolution) points EXCLUDED before plateau detection.\n")
        else:
            f.write("  ALL POINTS used — no Status-flag filtering (diagnostic mode).\n")
            f.write("  Compare against the filtered report to see the impact of M- exclusion.\n")
        f.write("=" * 78 + "\n\n")
        for info in infos:
            f.write("─" * 70 + "\n")
            f.write(f"SAMPLE: {info['sample']}   ({info['full_name']})\n")
            f.write(f"  source       : {info['source']}\n")
            f.write(f"  points       : {info['n_used']} used / {info['n_total']} total\n")
            f.write("─" * 70 + "\n")
            if "error" in info:
                f.write(f"  [!] {info['error']}\n\n")
                continue
            f.write(f"  plateau method  : {info['plateau_method']}\n")
            ps = info["plateau_strain_range_pct"]
            f.write(f"  plateau strain  : {ps[0]:.4g} – {ps[1]:.4g} %  "
                    f"({info['n_plateau_points']} points)\n")
            f.write(f"\n")
            f.write(f"  G'_LVR    = {info['Gp_LVR_Pa']:9.3f}  Pa\n")
            f.write(f"  G''_LVR   = {info['Gpp_LVR_Pa']:9.3f}  Pa\n")
            f.write(f"  tan δ_LVR = {info['tan_LVR']:9.4f}     ({_interpret_tan(info['tan_LVR'])})\n")
            f.write(f"  γ_LVR     = {info['gamma_LVR_pct']:9.4f}  %\n")
            f.write(f"  σ_y_LVR   = {info['sigma_y_LVR_Pa']:9.3f}  Pa   "
                    "(= G'_LVR × γ_LVR; LVR-endpoint yield-stress estimate)\n\n")

        # Summary table at the end — directly usable in a manuscript.
        rows = [i for i in infos if "error" not in i]
        if rows:
            f.write("=" * 78 + "\n")
            f.write("  SUMMARY TABLE\n")
            f.write("=" * 78 + "\n")
            gp_col = "G'_LVR (Pa)"
            gpp_col = "G''_LVR (Pa)"
            f.write(f"  {'Sample':<14} | {gp_col:>12} | {gpp_col:>13} | "
                    f"{'tan δ':>7} | {'γ_LVR (%)':>10} | {'σ_y (Pa)':>10}\n")
            f.write("  " + "-" * 80 + "\n")
            for i in rows:
                f.write(f"  {i['sample']:<14} | {i['Gp_LVR_Pa']:>12.2f} | "
                        f"{i['Gpp_LVR_Pa']:>13.2f} | {i['tan_LVR']:>7.3f} | "
                        f"{i['gamma_LVR_pct']:>10.4f} | {i['sigma_y_LVR_Pa']:>10.3f}\n")


def _interpret_tan(t: float) -> str:
    if t < 0.5:
        return "elastic-dominant gel"
    if t < 1.0:
        return "elastic-dominant"
    if t < 2.0:
        return "viscous-dominant (weak gel / entangled solution)"
    return "viscous-dominant (sol / liquid)"


def _write_csv(infos: list[dict], path: str):
    rows = [i for i in infos if "error" not in i]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "sample", "full_name", "source",
            "n_total", "n_used", "n_plateau_points",
            "plateau_strain_pct_lo", "plateau_strain_pct_hi",
            "Gp_LVR_Pa", "Gpp_LVR_Pa", "tan_LVR",
            "gamma_LVR_pct", "sigma_y_LVR_Pa",
            "plateau_method", "drop_M_minus",
        ])
        for i in rows:
            ps = i["plateau_strain_range_pct"]
            w.writerow([
                i["sample"], i["full_name"], i["source"],
                i["n_total"], i["n_used"], i["n_plateau_points"],
                f"{ps[0]:.6g}", f"{ps[1]:.6g}",
                f"{i['Gp_LVR_Pa']:.6g}", f"{i['Gpp_LVR_Pa']:.6g}",
                f"{i['tan_LVR']:.6g}", f"{i['gamma_LVR_pct']:.6g}",
                f"{i['sigma_y_LVR_Pa']:.6g}",
                i["plateau_method"], i["drop_M_minus"],
            ])


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def _run_one_pass(sweeps: dict, drop_m_minus: bool, suffix: str) -> tuple[str, str]:
    """Run one analysis pass with the given filter setting, write outputs."""
    out_txt = os.path.join(SAVE_PATH, f"{TAG}{suffix}.txt")
    out_csv = os.path.join(SAVE_PATH, f"{TAG}{suffix}.csv")
    infos = [analyse_one(s, drop_m_minus=drop_m_minus) for s in sweeps.values()]
    _write_report(infos, out_txt, drop_m_minus=drop_m_minus)
    _write_csv(infos, out_csv)
    return out_txt, out_csv


def main() -> int:
    os.makedirs(SAVE_PATH, exist_ok=True)

    print(f"Scanning {INPUT_FOLDER} for raw Rheocompass amplitude-sweep CSVs ...")
    sweeps = read_amplitude_sweeps_in_folder(INPUT_FOLDER)
    if SAMPLE_INCLUDE is not None:
        sweeps = {k: v for k, v in sweeps.items() if k in SAMPLE_INCLUDE}
    if not sweeps:
        print(f"  [!] No raw amplitude-sweep CSVs found in {INPUT_FOLDER}.")
        return 1
    print(f"  parsed {len(sweeps)} sample(s): {', '.join(sweeps.keys())}")

    # Pass 1: default (filtered) — what the scientific tables use.
    out_txt, out_csv = _run_one_pass(sweeps, drop_m_minus=DROP_M_MINUS, suffix="")
    print(f"[Saved] {out_txt}")
    print(f"[Saved] {out_csv}")

    # Pass 2: diagnostic (no Status-based filtering).
    if ALSO_WRITE_ALLPOINTS:
        out_txt2, out_csv2 = _run_one_pass(sweeps, drop_m_minus=False,
                                           suffix="_allpoints")
        print(f"[Saved] {out_txt2}   (no Status filtering — diagnostic)")
        print(f"[Saved] {out_csv2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
