"""
extract_recovery_v2.py
======================

Structural recovery from **Anton Paar Rheocompass 3iTT** (3-Interval
Thixotropy Test) data.

This is the v2 successor to ``Recovery_v1.py``. The old script computed
recovery from a steady-shear ramp protocol (compare η at the END low-shear
ramp vs the START low-shear ramp). The new 3iTT protocol replaces that
with an oscillation-rotation-oscillation cycle, which gives:

  * **Modulus recovery (%)**  — based on G' (storage modulus), the
    elastic-structure indicator
  * **Complex viscosity recovery (%)** — based on |η*|, an
    integrated elasticity+dissipation indicator

Both ratios compare the *recovered* oscillation interval (#3) to the
*baseline* oscillation interval (#1). The rotation interval (#2) is the
disruption phase and is not used in the ratio — its data is still parsed
and the breakdown is reported as a sanity check.

Reference / recovered point selection
-------------------------------------
By default both the reference (interval 1) and recovered (interval 3)
values are the **mean of the last ``N_TAIL_POINTS`` valid points** of
their respective intervals — the "settled" state at the end of each
oscillation. ``N_TAIL_POINTS`` defaults to 5; tweak via config.

Interval-3 cutoff (for crosslinking-during-test cases)
------------------------------------------------------
Some samples crosslink while the recovery interval is running, which
makes G' shoot above the interval-1 baseline (false "super-recovery").
For those samples set ``INTERVAL3_MAX_POINT`` to the point number BEYOND
WHICH data should be discarded. You can set it globally (applies to
every sample) or per-sample via ``INTERVAL3_MAX_POINT_OVERRIDES``.

Outputs
-------
``./Analises/Python/Results/Recovery_v2.txt``  human-readable
``./Analises/Python/Results/Recovery_v2.csv``  machine-readable

USAGE
-----
1. Edit the ``CONFIGURATION`` block below.
2. From the project root:
       python Export/01_Python/extract_recovery_v2.py

Author: T.M.C. Rodrigues - 2026-05-20
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))
from antpar_io import (   # noqa: E402
    read_recovery_sets_in_folder, RecoverySet,
)


# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION  ← edit here
# ──────────────────────────────────────────────────────────────────────────────

INPUT_FOLDER = "./Reologia/Recovery - 3D Bioprinting"

# Optional include list — None to process every set in INPUT_FOLDER.
SAMPLE_INCLUDE: list[str] | None = None

SAVE_PATH = "./Analises/Python/Results"
TAG = "Recovery_v2"

# How many points at the END of each oscillation interval to average for
# the "settled" reference / recovered value.
N_TAIL_POINTS = 5

# Global interval-3 cutoff (last point number to KEEP). None = use all.
# For Carla's calcium-phosphate-alginate dataset (crosslinking during the
# recovery interval), set this to 30 — beyond point 30 the data reflects
# crosslinking, not structural recovery.
INTERVAL3_MAX_POINT: int | None = None

# Per-sample overrides. Key is the short sample name (text before " - "
# in the Teste: header). These override INTERVAL3_MAX_POINT for the
# matching sample, leaving the global value to apply to all other samples.
INTERVAL3_MAX_POINT_OVERRIDES: dict[str, int | None] = {
    # Example for Carla's dataset:
    # "Amostra A":   30,
    # "Amostra B":   30,
    # "Amostra C":   30,
    # "Amostra D":   30,
    # "Amostra PC1": 30,
}

# Drop rows whose Status contains any of these flags BEFORE computing means.
# 'WMa' = wall slip / motor adjustment warning at the very start of an
#         interval; the next few points are usually settled. Default: drop.
DROP_STATUS_FLAGS: tuple[str, ...] = ("WMa",)

# Also write an "_allpoints" variant of the report and CSV, using EVERY
# measured point (no Status-based filtering). Useful as a diagnostic
# comparison to see how much the WMa exclusion shifts the recovery
# ratios. Set to False to skip the second pass.
ALSO_WRITE_ALLPOINTS = True


# ──────────────────────────────────────────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def _cutoff_for(sample_short: str) -> int | None:
    """Return the interval-3 cutoff to apply to this sample."""
    if sample_short in INTERVAL3_MAX_POINT_OVERRIDES:
        return INTERVAL3_MAX_POINT_OVERRIDES[sample_short]
    return INTERVAL3_MAX_POINT


def _filter_osc(df, drop_flags: tuple[str, ...] = DROP_STATUS_FLAGS,
                point_max: int | None = None):
    """Drop NaN / non-positive G' rows; optionally cap at ``point_max``
    and drop Status-flagged rows."""
    mask = (
        np.isfinite(df["G_storage_Pa"])
        & np.isfinite(df["eta_complex_Pa_s"])
        & (df["G_storage_Pa"] > 0)
        & (df["eta_complex_Pa_s"] > 0)
    )
    for flag in drop_flags:
        mask &= ~df["status"].str.contains(flag, regex=False, na=False)
    if point_max is not None:
        mask &= df["point"] <= point_max
    return df.loc[mask].sort_values("point").reset_index(drop=True)


def _tail_stats(arr: np.ndarray, n_tail: int) -> tuple[float, float]:
    """Return (mean, std) of the last ``n_tail`` finite, positive values."""
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr) & (arr > 0)]
    if len(arr) == 0:
        return float("nan"), float("nan")
    tail = arr[-min(n_tail, len(arr)):]
    return float(np.mean(tail)), float(np.std(tail, ddof=1) if len(tail) > 1 else 0.0)


def _ratio_recovery(recovered: float, baseline: float,
                    rec_std: float, base_std: float) -> tuple[float, float]:
    """Recovery (%) = recovered / baseline × 100 with propagated error."""
    if not (np.isfinite(recovered) and np.isfinite(baseline) and baseline > 0):
        return float("nan"), float("nan")
    pct = recovered / baseline * 100.0
    rel_err = math.sqrt(
        (rec_std / recovered) ** 2 + (base_std / baseline) ** 2
    ) if recovered > 0 else float("nan")
    err_pct = pct * rel_err if math.isfinite(rel_err) else float("nan")
    return pct, err_pct


# ──────────────────────────────────────────────────────────────────────────────
# CORE EXTRACTION
# ──────────────────────────────────────────────────────────────────────────────

def analyse_one(rs: RecoverySet,
                drop_flags: tuple[str, ...] = DROP_STATUS_FLAGS) -> dict:
    """Compute modulus + complex-viscosity recovery for one 3iTT set.

    ``drop_flags``: substrings; rows whose Status contains ANY of them
    are excluded from the per-interval baseline / recovered tails. Pass
    an empty tuple ``()`` to use every measured point regardless of
    Status (diagnostic mode).
    """
    short = rs.sample_name.split(" - ", 1)[0].strip()
    cutoff = _cutoff_for(short)

    info: dict = {
        "sample": short,
        "full_name": rs.sample_name,
        "basename": rs.source_basename,
        "interval_label": rs.interval_label,
        "interval1_kind": rs.interval1.kind,
        "interval2_kind": rs.interval2.kind,
        "interval3_kind": rs.interval3.kind,
        "interval3_cutoff": cutoff,
        "n_tail_points": N_TAIL_POINTS,
        "drop_status_flags": ",".join(drop_flags),
    }

    # Sanity: this extractor is built for Osc-Rot-Osc. If a future 3iTT
    # variant ships with different interval kinds, we surface a clear error.
    if rs.interval1.kind != "oscillation" or rs.interval3.kind != "oscillation":
        info["error"] = (
            f"This extractor expects oscillation intervals 1 and 3, "
            f"but got int1={rs.interval1.kind}, int3={rs.interval3.kind}. "
            f"Update extract_recovery_v2 if you're using a non-Osc-Rot-Osc 3iTT."
        )
        return info

    # ---- Baseline (interval 1) ----
    df1 = _filter_osc(rs.interval1.data, drop_flags=drop_flags)
    info["n_int1_used"] = len(df1)
    info["n_int1_total"] = len(rs.interval1.data)
    if len(df1) < N_TAIL_POINTS:
        info["error"] = f"Only {len(df1)} valid baseline points (need ≥ {N_TAIL_POINTS})."
        return info
    gp1_mean, gp1_std = _tail_stats(df1["G_storage_Pa"].to_numpy(), N_TAIL_POINTS)
    eta1_mean, eta1_std = _tail_stats(df1["eta_complex_Pa_s"].to_numpy(), N_TAIL_POINTS)
    info["Gp_baseline_Pa"] = gp1_mean
    info["Gp_baseline_std_Pa"] = gp1_std
    info["eta_baseline_Pa_s"] = eta1_mean
    info["eta_baseline_std_Pa_s"] = eta1_std
    info["baseline_tail_time_s"] = (
        float(df1["time_s"].iloc[-N_TAIL_POINTS]),
        float(df1["time_s"].iloc[-1]),
    )

    # ---- Disruption (interval 2) — diagnostic only ----
    df2 = rs.interval2.data
    visc_disruption_end = float(df2["viscosity_Pa_s"].iloc[-N_TAIL_POINTS:].mean())
    info["disruption_shear_rate_1_s"] = float(df2["shear_rate_1_s"].median())
    info["viscosity_disruption_end_Pa_s"] = visc_disruption_end

    # ---- Recovered (interval 3, with cutoff if any) ----
    df3 = _filter_osc(rs.interval3.data, drop_flags=drop_flags, point_max=cutoff)
    info["n_int3_used"] = len(df3)
    info["n_int3_total"] = len(rs.interval3.data)
    if len(df3) < N_TAIL_POINTS:
        info["error"] = (
            f"Only {len(df3)} valid recovery points within cutoff "
            f"(cutoff={cutoff}). Need ≥ {N_TAIL_POINTS}."
        )
        return info
    gp3_mean, gp3_std = _tail_stats(df3["G_storage_Pa"].to_numpy(), N_TAIL_POINTS)
    eta3_mean, eta3_std = _tail_stats(df3["eta_complex_Pa_s"].to_numpy(), N_TAIL_POINTS)
    info["Gp_recovered_Pa"] = gp3_mean
    info["Gp_recovered_std_Pa"] = gp3_std
    info["eta_recovered_Pa_s"] = eta3_mean
    info["eta_recovered_std_Pa_s"] = eta3_std
    info["recovered_tail_time_s"] = (
        float(df3["time_s"].iloc[-N_TAIL_POINTS]),
        float(df3["time_s"].iloc[-1]),
    )
    info["recovered_tail_point_range"] = (
        int(df3["point"].iloc[-N_TAIL_POINTS]),
        int(df3["point"].iloc[-1]),
    )

    # ---- Recovery ratios ----
    rec_G_pct, err_G_pct = _ratio_recovery(gp3_mean, gp1_mean, gp3_std, gp1_std)
    rec_eta_pct, err_eta_pct = _ratio_recovery(eta3_mean, eta1_mean, eta3_std, eta1_std)
    info["modulus_recovery_pct"] = rec_G_pct
    info["modulus_recovery_err_pct"] = err_G_pct
    info["eta_recovery_pct"] = rec_eta_pct
    info["eta_recovery_err_pct"] = err_eta_pct
    return info


# ──────────────────────────────────────────────────────────────────────────────
# REPORTING
# ──────────────────────────────────────────────────────────────────────────────

def _write_report(infos: list[dict], path: str, drop_flags: tuple[str, ...]):
    flags = ", ".join(drop_flags) or "(none)"
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write(f"  3iTT RECOVERY EXTRACTION  ({TAG})\n")
        f.write(f"  Anton Paar Rheocompass CSV  via antpar_io.read_recovery_set\n")
        f.write(f"  Reference / recovered = mean of last {N_TAIL_POINTS} valid points per interval\n")
        if drop_flags:
            f.write(f"  Status flags EXCLUDED before averaging: {flags}\n")
        else:
            f.write("  ALL POINTS used — no Status-flag filtering (diagnostic mode).\n")
            f.write("  Compare against the filtered report to see the impact of WMa exclusion.\n")
        f.write(f"  Global interval-3 cutoff: "
                f"{INTERVAL3_MAX_POINT if INTERVAL3_MAX_POINT is not None else 'none (use all)'}\n")
        if INTERVAL3_MAX_POINT_OVERRIDES:
            f.write(f"  Per-sample overrides: {INTERVAL3_MAX_POINT_OVERRIDES}\n")
        f.write("=" * 78 + "\n\n")

        for info in infos:
            f.write("─" * 70 + "\n")
            f.write(f"SAMPLE: {info['sample']}   ({info['full_name']})\n")
            f.write(f"  basename : {info['basename']}\n")
            f.write(f"  label    : {info['interval_label']}\n")
            f.write(f"  intervals: 1={info['interval1_kind']}, "
                    f"2={info['interval2_kind']}, 3={info['interval3_kind']}\n")
            cutoff = info["interval3_cutoff"]
            f.write(f"  int3 cut : {cutoff if cutoff is not None else 'none'}\n")
            f.write("─" * 70 + "\n")
            if "error" in info:
                f.write(f"  [!] {info['error']}\n\n")
                continue

            t_lo, t_hi = info["baseline_tail_time_s"]
            f.write(f"\n  BASELINE  (interval 1, points used = {info['n_int1_used']}/"
                    f"{info['n_int1_total']}, tail t ∈ [{t_lo:.1f}, {t_hi:.1f}] s)\n")
            f.write(f"    G'        = {info['Gp_baseline_Pa']:9.3f} ± "
                    f"{info['Gp_baseline_std_Pa']:.3f} Pa\n")
            f.write(f"    |η*|      = {info['eta_baseline_Pa_s']:9.3f} ± "
                    f"{info['eta_baseline_std_Pa_s']:.3f} Pa·s\n")

            f.write(f"\n  DISRUPTION  (interval 2, γ̇ = "
                    f"{info['disruption_shear_rate_1_s']:.1f} s⁻¹)\n")
            f.write(f"    η_end     = {info['viscosity_disruption_end_Pa_s']:9.3f} Pa·s"
                    f"   (diagnostic only)\n")

            t_lo, t_hi = info["recovered_tail_time_s"]
            p_lo, p_hi = info["recovered_tail_point_range"]
            f.write(f"\n  RECOVERED  (interval 3, points used = {info['n_int3_used']}/"
                    f"{info['n_int3_total']}, "
                    f"tail = points {p_lo}–{p_hi}, t ∈ [{t_lo:.1f}, {t_hi:.1f}] s)\n")
            f.write(f"    G'        = {info['Gp_recovered_Pa']:9.3f} ± "
                    f"{info['Gp_recovered_std_Pa']:.3f} Pa\n")
            f.write(f"    |η*|      = {info['eta_recovered_Pa_s']:9.3f} ± "
                    f"{info['eta_recovered_std_Pa_s']:.3f} Pa·s\n")

            f.write(f"\n  RECOVERY RATIOS\n")
            f.write(f"    Modulus recovery (G')          = "
                    f"{info['modulus_recovery_pct']:6.2f} ± "
                    f"{info['modulus_recovery_err_pct']:.2f} %\n")
            f.write(f"    Complex viscosity recovery |η*| = "
                    f"{info['eta_recovery_pct']:6.2f} ± "
                    f"{info['eta_recovery_err_pct']:.2f} %\n\n")

        # Summary table
        rows = [i for i in infos if "error" not in i]
        if rows:
            f.write("=" * 78 + "\n")
            f.write("  SUMMARY TABLE\n")
            f.write("=" * 78 + "\n")
            f.write(f"  {'Sample':<18} | {'G''_base':>10} | {'G''_rec':>10} | "
                    f"{'Rec_G (%)':>11} | {'η*_base':>10} | {'η*_rec':>10} | "
                    f"{'Rec_η* (%)':>11} | {'cutoff':>7}\n")
            f.write("  " + "-" * 110 + "\n")
            for i in rows:
                cutoff = i["interval3_cutoff"]
                cutoff_s = f"{cutoff}" if cutoff is not None else "—"
                f.write(f"  {i['sample']:<18} | {i['Gp_baseline_Pa']:>10.2f} | "
                        f"{i['Gp_recovered_Pa']:>10.2f} | "
                        f"{i['modulus_recovery_pct']:>11.2f} | "
                        f"{i['eta_baseline_Pa_s']:>10.2f} | "
                        f"{i['eta_recovered_Pa_s']:>10.2f} | "
                        f"{i['eta_recovery_pct']:>11.2f} | "
                        f"{cutoff_s:>7}\n")


def _write_csv(infos: list[dict], path: str):
    rows = [i for i in infos if "error" not in i]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "sample", "full_name", "basename", "interval_label",
            "interval1_kind", "interval2_kind", "interval3_kind",
            "interval3_cutoff", "n_tail_points", "drop_status_flags",
            "n_int1_used", "n_int1_total",
            "n_int3_used", "n_int3_total",
            "Gp_baseline_Pa", "Gp_baseline_std_Pa",
            "eta_baseline_Pa_s", "eta_baseline_std_Pa_s",
            "disruption_shear_rate_1_s", "viscosity_disruption_end_Pa_s",
            "Gp_recovered_Pa", "Gp_recovered_std_Pa",
            "eta_recovered_Pa_s", "eta_recovered_std_Pa_s",
            "modulus_recovery_pct", "modulus_recovery_err_pct",
            "eta_recovery_pct", "eta_recovery_err_pct",
        ])
        for i in rows:
            cutoff = i["interval3_cutoff"]
            w.writerow([
                i["sample"], i["full_name"], i["basename"], i["interval_label"],
                i["interval1_kind"], i["interval2_kind"], i["interval3_kind"],
                cutoff if cutoff is not None else "",
                i["n_tail_points"], i.get("drop_status_flags", ""),
                i["n_int1_used"], i["n_int1_total"],
                i["n_int3_used"], i["n_int3_total"],
                f"{i['Gp_baseline_Pa']:.6g}", f"{i['Gp_baseline_std_Pa']:.6g}",
                f"{i['eta_baseline_Pa_s']:.6g}", f"{i['eta_baseline_std_Pa_s']:.6g}",
                f"{i['disruption_shear_rate_1_s']:.6g}",
                f"{i['viscosity_disruption_end_Pa_s']:.6g}",
                f"{i['Gp_recovered_Pa']:.6g}", f"{i['Gp_recovered_std_Pa']:.6g}",
                f"{i['eta_recovered_Pa_s']:.6g}", f"{i['eta_recovered_std_Pa_s']:.6g}",
                f"{i['modulus_recovery_pct']:.6g}", f"{i['modulus_recovery_err_pct']:.6g}",
                f"{i['eta_recovery_pct']:.6g}", f"{i['eta_recovery_err_pct']:.6g}",
            ])


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def _run_one_pass(sets: dict, drop_flags: tuple[str, ...],
                  suffix: str) -> tuple[str, str]:
    """Run one analysis pass with the given drop_flags, write outputs."""
    out_txt = os.path.join(SAVE_PATH, f"{TAG}{suffix}.txt")
    out_csv = os.path.join(SAVE_PATH, f"{TAG}{suffix}.csv")
    infos = [analyse_one(rs, drop_flags=drop_flags) for rs in sets.values()]
    _write_report(infos, out_txt, drop_flags=drop_flags)
    _write_csv(infos, out_csv)
    return out_txt, out_csv


def main() -> int:
    os.makedirs(SAVE_PATH, exist_ok=True)

    print(f"Scanning {INPUT_FOLDER} for raw Rheocompass 3iTT recovery sets ...")
    sets = read_recovery_sets_in_folder(INPUT_FOLDER)
    if SAMPLE_INCLUDE is not None:
        sets = {k: v for k, v in sets.items() if k in SAMPLE_INCLUDE}
    if not sets:
        print(f"  [!] No complete recovery sets found in {INPUT_FOLDER}.")
        return 1
    print(f"  loaded {len(sets)} set(s): {', '.join(sets.keys())}")

    # Pass 1: default (filtered) — what the scientific tables use.
    out_txt, out_csv = _run_one_pass(sets, drop_flags=DROP_STATUS_FLAGS, suffix="")
    print(f"[Saved] {out_txt}")
    print(f"[Saved] {out_csv}")

    # Pass 2: diagnostic (no Status-based filtering).
    if ALSO_WRITE_ALLPOINTS:
        out_txt2, out_csv2 = _run_one_pass(sets, drop_flags=(),
                                           suffix="_allpoints")
        print(f"[Saved] {out_txt2}   (no Status filtering — diagnostic)")
        print(f"[Saved] {out_csv2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
