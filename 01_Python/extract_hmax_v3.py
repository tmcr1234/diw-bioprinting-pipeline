"""
extract_hmax_v3.py
==================

Maximum self-supported scaffold height (h_max) and layer-count predictions
per ink, computed directly from the **v2 SAOS CSVs** (no .xls anywhere).

This is the v3 successor to ``extract_hmax_v2.py``. The physics is
identical; only the data inputs changed:

    v2 → reads .xls amplitude + frequency sweeps via xlrd
         (re-parses raw data + re-fits the low-ω power law every run)
    v3 → reads SAOS_amplitude_v2.csv + SAOS_frequency_v2.csv directly
         (single source of truth — the v2 extractors did the fits once)

Why ``v3`` matters
------------------
Removes the last `xlrd` dependency from the pipeline. Re-runs are fast
(no re-parsing) and deterministic (uses the same fits the manuscript
tables will quote). For audit trail: regenerate the upstream v2 CSVs,
then re-run v3 — no risk of v3 producing a number that drifts from the
SAOS reports.

Four h_max criteria  (with ρ = 1000 kg/m³, g = 9.81 m/s²)
---------------------------------------------------------
(A) LVR endpoint:           σ_max = G'_LVR · γ_LVR (fractional)
                             → h_max = σ_max / (ρ g)
(B) Practical:              σ_max = G'(ω = 1 rad/s)  (measured)
                             → h_max = σ_max / (ρ g)
(C) Quasi-static:           σ_max = G'(ω = 0.01 rad/s)  (extrapolated)
                             → h_max = σ_max / (ρ g)
(D) Buckling-corrected:     Smay 2002 — h_max(L) = √(1.94 G' / ρ g) · √L
                             where L is the unsupported-span layer length

Layer-count prediction
----------------------
For each criterion's h_max and each typical DIW layer height
``h_layer ∈ {0.3, 0.4, 0.5} mm``:
    N_max = floor(h_max / h_layer)

Inputs
------
``SAOS_amplitude_v2.csv``  ← from ``extract_SAOS_amplitude_v2.py``
    must contain columns: sample, Gp_LVR_Pa, gamma_LVR_pct
``SAOS_frequency_v2.csv``  ← from ``extract_SAOS_frequency_v2.py``
    must contain columns: sample, Gp_omega1_Pa, lowomega_G0_Pa,
                          lowomega_beta, Gp_extrap_Pa, omega_extrap_rad_s

The two CSVs are joined on `sample`. Inks present in only one of them
are skipped with a warning.

Outputs
-------
``./Analises/Python/Results/SAOS_hmax_v3.txt``                    (report)
``./Analises/Python/Results/printing_parameters_per_ink_v3.csv``  (table)
``./Analises/Python/Results/Gprime_extrap_<ink>_v3.png``          (optional)

Optional plotting
-----------------
If ``FREQUENCY_RAW_FOLDER`` is set to the folder containing the raw
frequency-sweep CSVs, per-ink low-ω extrapolation plots are generated
using ``antpar_io.read_frequency_sweep_csv``. The fit overlaid on the
plot uses the EXACT (G0, β) from the CSV — no re-fitting, no risk of
drift from the CSV values. Set to ``None`` to skip plotting.

Author: T.M.C. Rodrigues - 2026-05-20
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# matplotlib is lazy-imported inside _plot_extrapolation() — see notes in the
# plotting helper. Core h_max math has zero plotting dependency.

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))


# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION  ← edit here
# ──────────────────────────────────────────────────────────────────────────────

# Inputs — produced upstream by the v2 SAOS scripts.
AMPLITUDE_CSV = "./Analises/Python/Results/SAOS_amplitude_v2.csv"
FREQUENCY_CSV = "./Analises/Python/Results/SAOS_frequency_v2.csv"

# Optional — set to a folder of raw Rheocompass frequency-sweep CSVs to
# enable per-ink low-ω extrapolation plots. None = skip plotting.
FREQUENCY_RAW_FOLDER: str | None = None

# Outputs
SAVE_PATH = "./Analises/Python/Results"
TAG = "SAOS_hmax_v3"

# Physical constants
RHO = 1000.0   # kg/m³ — adjust if your ink is significantly denser
G_ACC = 9.81   # m/s²

# Reporting parameters
H_LAYER_TARGETS_MM = [0.30, 0.40, 0.50]         # for N_max layer counts
LAYER_HEIGHTS_FOR_BUCKLE_MM = [0.5, 1.0, 2.0]   # unsupported-span lengths

# Sanity-check the freq v2 extrapolation target. The v2 CSV records its
# `omega_extrap_rad_s` per row; this is the expected default. If the v2
# CSV used a different ω_extrap, the report will note it per-ink.
OMEGA_EXTRAP_EXPECTED = 0.01


# ──────────────────────────────────────────────────────────────────────────────
# CORE FORMULAE
# ──────────────────────────────────────────────────────────────────────────────

def h_max_stacked_mm(sigma_max_Pa: float) -> float:
    """Static stacked height limit:  h = σ_max / (ρ g)   →  mm."""
    if not np.isfinite(sigma_max_Pa) or sigma_max_Pa <= 0:
        return float("nan")
    return 1000.0 * sigma_max_Pa / (RHO * G_ACC)


def h_max_buckle_mm(Gp_Pa: float, L_mm: float) -> float:
    """Smay-2002 buckling-corrected unsupported-span height (mm)."""
    if not np.isfinite(Gp_Pa) or Gp_Pa <= 0 or L_mm <= 0:
        return float("nan")
    L_m = L_mm * 1e-3
    return 1000.0 * np.sqrt(1.94 * Gp_Pa / (RHO * G_ACC)) * np.sqrt(L_m)


def n_layers_max(h_max_mm: float, h_layer_mm: float) -> int:
    if not np.isfinite(h_max_mm) or h_max_mm <= 0 or h_layer_mm <= 0:
        return 0
    return int(np.floor(h_max_mm / h_layer_mm))


# ──────────────────────────────────────────────────────────────────────────────
# JOIN THE TWO V2 CSVs ON `sample`
# ──────────────────────────────────────────────────────────────────────────────

def _load_amplitude(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"sample", "Gp_LVR_Pa", "gamma_LVR_pct"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: missing required columns {missing}. "
            f"Regenerate via extract_SAOS_amplitude_v2.py."
        )
    return df[["sample", "Gp_LVR_Pa", "gamma_LVR_pct"]].copy()


def _load_frequency(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"sample", "Gp_omega1_Pa", "lowomega_G0_Pa",
              "lowomega_beta", "Gp_extrap_Pa", "omega_extrap_rad_s"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: missing required columns {missing}. "
            f"Regenerate via extract_SAOS_frequency_v2.py."
        )
    return df[["sample", "Gp_omega1_Pa", "lowomega_G0_Pa", "lowomega_beta",
               "Gp_extrap_Pa", "omega_extrap_rad_s"]].copy()


def _merge(amp: pd.DataFrame, freq: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Inner join + log which samples were dropped from each side."""
    merged = pd.merge(amp, freq, on="sample", how="inner")
    only_amp = set(amp["sample"]) - set(freq["sample"])
    only_freq = set(freq["sample"]) - set(amp["sample"])
    warnings = []
    for s in sorted(only_amp):
        warnings.append(f"  [!] {s}: in amplitude CSV but missing from frequency CSV — skipped")
    for s in sorted(only_freq):
        warnings.append(f"  [!] {s}: in frequency CSV but missing from amplitude CSV — skipped")
    return merged, warnings


# ──────────────────────────────────────────────────────────────────────────────
# PER-INK COMPUTATION
# ──────────────────────────────────────────────────────────────────────────────

def compute_one(row: pd.Series) -> dict:
    sample = str(row["sample"])
    Gp_LVR = float(row["Gp_LVR_Pa"])
    gamma_LVR_frac = float(row["gamma_LVR_pct"]) / 100.0
    Gp_omega1 = float(row["Gp_omega1_Pa"])
    G0 = float(row["lowomega_G0_Pa"])
    beta = float(row["lowomega_beta"])
    Gp_extrap = float(row["Gp_extrap_Pa"])
    omega_extrap = float(row["omega_extrap_rad_s"])

    sig_lvr = Gp_LVR * gamma_LVR_frac    # criterion (A)
    sig_g1 = Gp_omega1                   # criterion (B)
    sig_gex = Gp_extrap                  # criterion (C)

    info = {
        "sample": sample,
        # Inputs (pass-through for traceability)
        "Gp_LVR_Pa": Gp_LVR,
        "gamma_LVR_pct": gamma_LVR_frac * 100,
        "Gp_omega1_Pa": Gp_omega1,
        "G0_lowomega": G0,
        "beta_lowomega": beta,
        "omega_extrap_rad_s": omega_extrap,
        "Gp_extrap_Pa": Gp_extrap,
        # Yield-stress estimates per criterion
        "sigma_LVR_Pa": sig_lvr,
        "sigma_G1_Pa":  sig_g1,
        "sigma_G01_Pa": sig_gex,
        # h_max per criterion (mm)
        "h_LVR_mm": h_max_stacked_mm(sig_lvr),
        "h_G1_mm":  h_max_stacked_mm(sig_g1),
        "h_G01_mm": h_max_stacked_mm(sig_gex),
    }
    # Buckling-corrected at each unsupported-span length
    for L in LAYER_HEIGHTS_FOR_BUCKLE_MM:
        info[f"h_buckle_L{L}mm_mm"] = h_max_buckle_mm(Gp_omega1, L)
    # N_max layer counts under each criterion / layer height
    for crit_name, h in (("LVR", info["h_LVR_mm"]),
                         ("G1",  info["h_G1_mm"]),
                         ("G01", info["h_G01_mm"])):
        for hL in H_LAYER_TARGETS_MM:
            info[f"Nmax_{crit_name}_hL{hL:.2f}"] = n_layers_max(h, hL)
    return info


# ──────────────────────────────────────────────────────────────────────────────
# PLOTTING (optional)
# ──────────────────────────────────────────────────────────────────────────────

def _plot_extrapolation(sample: str, raw_folder: Path, info: dict,
                        output_dir: Path) -> Path | None:
    """Per-ink low-ω extrapolation plot using the freq v2 fit parameters.

    matplotlib is imported lazily here — if the user doesn't ask for plots
    (FREQUENCY_RAW_FOLDER = None) the core h_max calculation never touches it,
    so a minimal venv without matplotlib still runs the pipeline successfully.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"    [skip] matplotlib not installed — install it to enable plots")
        return None

    from antpar_io import read_frequency_sweep_csv, list_raw_csvs   # noqa: E402

    # Find the raw freq CSV whose Teste: header matches this short sample name.
    # Cheapest robust match: parse each CSV in the folder until short-name match.
    target = None
    for p in list_raw_csvs(raw_folder):
        try:
            fsw = read_frequency_sweep_csv(p)
        except Exception:
            continue
        short = fsw.sample_name.split(" - ", 1)[0].strip()
        if short == sample:
            target = (p, fsw)
            break
    if target is None:
        return None
    _path, fsw = target
    df = fsw.data
    valid = (df["G_storage_Pa"] > 0) & (df["omega_rad_s"] > 0)
    w = df.loc[valid, "omega_rad_s"].to_numpy()
    Gp = df.loc[valid, "G_storage_Pa"].to_numpy()
    if len(w) < 2:
        return None
    order = np.argsort(w)
    w, Gp = w[order], Gp[order]

    G0 = info["G0_lowomega"]
    beta = info["beta_lowomega"]
    omega_extrap = info["omega_extrap_rad_s"]

    plt.figure(figsize=(6.0, 4.5))
    plt.loglog(w, Gp, "o", label="measured")
    w_ext = np.logspace(np.log10(omega_extrap * 0.5),
                        np.log10(w.max() * 1.1), 60)
    Gp_ext = G0 * w_ext ** beta
    plt.loglog(w_ext, Gp_ext, "--",
               label=f"v2 fit:  G'={G0:.2g}·ω^{beta:.3f}")
    plt.axvline(omega_extrap, color="gray", ls=":", lw=0.8)
    plt.text(omega_extrap * 1.15, 1.5 * np.nanmin(Gp[Gp > 0]),
             f"ω = {omega_extrap}", color="gray", fontsize=9)
    plt.xlabel("Angular frequency  ω  (rad/s)")
    plt.ylabel("Storage modulus  G'  (Pa)")
    plt.title(f"{sample}: low-ω extrapolation (v3, fit from SAOS_frequency_v2)")
    plt.legend()
    plt.grid(True, which="both", ls=":", alpha=0.5)
    plt.tight_layout()
    out_path = output_dir / f"Gprime_extrap_{sample}_v3.png"
    plt.savefig(out_path, dpi=160)
    plt.close()
    return out_path


# ──────────────────────────────────────────────────────────────────────────────
# REPORTING
# ──────────────────────────────────────────────────────────────────────────────

def _write_report(infos: list[dict], warnings: list[str], path: str):
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write(f"  SAOS-DERIVED MAX SCAFFOLD HEIGHT  ({TAG})\n")
        f.write(f"  Inputs: {AMPLITUDE_CSV} + {FREQUENCY_CSV}\n")
        f.write(f"  ρ = {RHO} kg/m³ | g = {G_ACC} m/s²\n")
        f.write(f"  h_layer targets (mm): {H_LAYER_TARGETS_MM}\n")
        f.write(f"  Unsupported-span L (mm) for buckling: {LAYER_HEIGHTS_FOR_BUCKLE_MM}\n")
        f.write("=" * 78 + "\n\n")
        if warnings:
            f.write("WARNINGS\n")
            for w in warnings:
                f.write(w + "\n")
            f.write("\n")

        # Per-ink detail
        for info in infos:
            sample = info["sample"]
            f.write("─" * 70 + "\n")
            f.write(f"SAMPLE: {sample}\n")
            f.write("─" * 70 + "\n")
            f.write(f"  Inputs:\n")
            f.write(f"    G'_LVR        = {info['Gp_LVR_Pa']:9.3f}  Pa\n")
            f.write(f"    γ_LVR         = {info['gamma_LVR_pct']:9.4f}  %\n")
            f.write(f"    G'(ω=1)       = {info['Gp_omega1_Pa']:9.3f}  Pa\n")
            f.write(f"    low-ω fit     = G' = {info['G0_lowomega']:.4g} · ω^{info['beta_lowomega']:.4f}\n")
            f.write(f"    G'(ω={info['omega_extrap_rad_s']}) = {info['Gp_extrap_Pa']:9.4g}  Pa  (extrapolated)\n")
            f.write(f"\n  Yield-stress estimates:\n")
            f.write(f"    σ_LVR (A) = G'_LVR · γ_LVR = {info['sigma_LVR_Pa']:9.3f} Pa\n")
            f.write(f"    σ_G1  (B) = G'(ω=1)        = {info['sigma_G1_Pa']:9.3f} Pa\n")
            f.write(f"    σ_G01 (C) = G'(ω={info['omega_extrap_rad_s']}) "
                    f"= {info['sigma_G01_Pa']:9.4g} Pa\n")
            f.write(f"\n  h_max  (stacked, full sub-layer support):\n")
            f.write(f"    (A) LVR     :  h = {info['h_LVR_mm']:7.3f}  mm\n")
            f.write(f"    (B) G'(ω=1) :  h = {info['h_G1_mm']:7.3f}  mm\n")
            f.write(f"    (C) G'(ω=0.01) (extrap): h = {info['h_G01_mm']:7.4g}  mm\n")
            f.write(f"\n  h_max  (Smay 2002 buckling-corrected, G' = G'(ω=1)):\n")
            for L in LAYER_HEIGHTS_FOR_BUCKLE_MM:
                f.write(f"    L = {L:.1f} mm  →  h = "
                        f"{info[f'h_buckle_L{L}mm_mm']:7.3f}  mm\n")
            f.write(f"\n  N_max layer count  (criterion × h_layer):\n")
            f.write(f"    {'criterion':<6} | " +
                    " ".join(f"{f'h_layer={hL:.2f}mm':>14}" for hL in H_LAYER_TARGETS_MM) + "\n")
            for crit_name, crit_label in (("LVR", "(A) LVR"),
                                          ("G1",  "(B) G'(ω=1)"),
                                          ("G01", "(C) G'(ω=.01)")):
                cells = " ".join(
                    f"{info[f'Nmax_{crit_name}_hL{hL:.2f}']:>14d}"
                    for hL in H_LAYER_TARGETS_MM
                )
                f.write(f"    {crit_label:<13} | {cells}\n")
            f.write("\n")

        # Summary table
        if infos:
            f.write("=" * 78 + "\n")
            f.write("  SUMMARY  —  h_max (mm) per ink per criterion\n")
            f.write("=" * 78 + "\n")
            f.write(f"  {'Sample':<14} | {'h_LVR':>8} | {'h_G(ω=1)':>9} | "
                    f"{'h_G(ω=.01)':>11} | {'h_buckle(L=1mm)':>16}\n")
            f.write("  " + "-" * 70 + "\n")
            for info in infos:
                f.write(f"  {info['sample']:<14} | "
                        f"{info['h_LVR_mm']:>8.3f} | "
                        f"{info['h_G1_mm']:>9.3f} | "
                        f"{info['h_G01_mm']:>11.4g} | "
                        f"{info['h_buckle_L1.0mm_mm']:>16.3f}\n")


def _write_csv(infos: list[dict], path: str):
    if not infos:
        return
    columns = list(infos[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for info in infos:
            w.writerow([
                f"{v:.6g}" if isinstance(v, float) else v
                for v in (info[c] for c in columns)
            ])


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main() -> int:
    os.makedirs(SAVE_PATH, exist_ok=True)
    out_txt = os.path.join(SAVE_PATH, f"{TAG}.txt")
    out_csv = os.path.join(SAVE_PATH, "printing_parameters_per_ink_v3.csv")

    print(f"Loading inputs:")
    print(f"  amplitude → {AMPLITUDE_CSV}")
    print(f"  frequency → {FREQUENCY_CSV}")
    if not Path(AMPLITUDE_CSV).exists() or not Path(FREQUENCY_CSV).exists():
        print("  [!] One of the input CSVs is missing. Run the v2 SAOS extractors first.")
        return 1

    amp = _load_amplitude(AMPLITUDE_CSV)
    freq = _load_frequency(FREQUENCY_CSV)
    merged, warnings = _merge(amp, freq)
    print(f"  joined → {len(merged)} ink(s): {', '.join(merged['sample'])}")
    for w in warnings:
        print(w)

    infos = [compute_one(row) for _, row in merged.iterrows()]
    _write_report(infos, warnings, out_txt)
    _write_csv(infos, out_csv)
    print(f"[Saved] {out_txt}")
    print(f"[Saved] {out_csv}")

    if FREQUENCY_RAW_FOLDER:
        raw_folder = Path(FREQUENCY_RAW_FOLDER)
        if not raw_folder.is_dir():
            print(f"  [!] FREQUENCY_RAW_FOLDER not a directory; skipping plots.")
        else:
            print(f"  generating per-ink plots from {raw_folder} ...")
            for info in infos:
                out = _plot_extrapolation(info["sample"], raw_folder, info,
                                          Path(SAVE_PATH))
                if out is not None:
                    print(f"    [plot] {out}")
                else:
                    print(f"    [skip] no raw freq CSV matched '{info['sample']}'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
