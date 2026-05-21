"""
extract_hmax_v2.py
==================

╔══════════════════════════════════════════════════════════════════════════╗
║  ⚠ LEGACY — superseded by extract_hmax_v3.py.                            ║
║  v3 reads the v2 SAOS CSVs (SAOS_amplitude_v2.csv + SAOS_frequency_v2.csv)║
║  directly — no xlrd, no .xls re-parsing, single source of truth for the  ║
║  upstream fits. v2 is kept here only as the audit reference for any      ║
║  previously published h_max numbers.                                     ║
║  See Export/docs/01_rheology_data_formats.md and the v3 manual.          ║
╚══════════════════════════════════════════════════════════════════════════╝

Maximum self-supported scaffold height (h_max) per ink, with:

  * G'(omega) plateau from frequency sweep (measured at omega in [1, 240] rad/s)
  * Low-omega power-law extrapolation  G' ~ G0 * omega^beta
    -> evaluated at omega = 0.01 rad/s (quasi-static, long-dwell relevant)
  * gamma_LVR endpoint from amplitude sweep
    -> yield stress estimate tau_y = G'_LVR * gamma_LVR

Four criteria for h_max (rho = 1000 kg/m^3, g = 9.81 m/s^2):

  (A) LVR-endpoint:        sigma_max = G'_LVR * gamma_LVR
  (B) Practical G'(omega=1 rad/s):   sigma_max = G'(omega=1)
  (C) Quasi-static G'(omega=0.01 rad/s) (extrapolated)
  (D) Buckling-corrected (Smay 2002): h_max(L) = sqrt(1.94 G' / (rho g) ) * sqrt(L)

Also reports N_max (number of stacked layers) for h_layer in {0.3, 0.4, 0.5} mm.

Outputs:
  - stdout summary
  - Analises/Python/Results/SAOS_hmax_v2.txt          (full report)
  - Analises/Python/Results/printing_parameters_per_ink.csv  (machine-readable)
  - Analises/Python/Results/Gprime_extrap_<ink>.png   (per-ink plot, low-omega fit)

Status registry (Section 6.9):
  - Script: extract_hmax_v2.py
  - Version: v2
  - Status:  active
  - Supersedes: extract_hmax_v1.py (kept for audit only)

Author: T.M.C. Rodrigues -- 2026-05-14
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract_SAOS_values import (
    read_antpar_xls, find_col, extract_amplitude_sweep
)

PROJECT_ROOT = "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE"
RHO = 1000.0
G_ACC = 9.81
OMEGA_EXTRAP = 0.01
H_LAYER_TARGETS_MM = [0.30, 0.40, 0.50]   # for N_max layer-count predictions
LAYER_HEIGHTS_FOR_BUCKLE_MM = [0.5, 1.0, 2.0]

FREQ_FILES = {
    "C15":    "Reologia/Frequency/C15 - Freq CP50 0,1mm.xls",
    "NE":     "Reologia/Frequency/C15 Gira 5.5  - Freq CP50 0,1mm.xls",
    "Bozano": "Reologia/Frequency/Gel de Cabelo Bozano - Freq CP50 0,1mm.xls",
}
STRAIN_FILES = {
    "C15":    "Reologia/Strain/C15 - Deformacao CP50 0,1mm.xls",
    "NE":     "Reologia/Strain/C15 Gira 5.5  - Deformacao CP50 0,1mm.xls",
    "Bozano": "Reologia/Strain/Gel de Cabelo Bozano - Deformacao CP50 0,1mm.xls",
}


def load_freq_curve(filepath):
    headers, units, arr, _ = read_antpar_xls(filepath)
    i_w  = find_col(headers, "Angular frequency", "Angular Frequency",
                    "Omega", "omega")
    i_Gp = find_col(headers, "Storage modulus", "G'")
    if i_w is None or i_Gp is None:
        raise RuntimeError(f"Could not locate columns in {filepath}")
    w  = arr[:, i_w]
    Gp = arr[:, i_Gp]
    valid = ~np.isnan(w) & ~np.isnan(Gp) & (w > 0) & (Gp > 0)
    w, Gp = w[valid], Gp[valid]
    order = np.argsort(w)
    return w[order], Gp[order]


def fit_low_omega_power_law(w, Gp, n_low=6):
    n = min(n_low, len(w))
    lw, lg = np.log10(w[:n]), np.log10(Gp[:n])
    A = np.vstack([lw, np.ones_like(lw)]).T
    coef, *_ = np.linalg.lstsq(A, lg, rcond=None)
    beta, log_G0 = coef
    G0 = 10.0 ** log_G0
    pred = log_G0 + beta * lw
    ss_res = float(np.sum((lg - pred) ** 2))
    ss_tot = float(np.sum((lg - np.mean(lg)) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(G0), float(beta), R2, n


def h_max_stacked_mm(sigma_max_Pa):
    return 1000.0 * sigma_max_Pa / (RHO * G_ACC)


def h_max_buckle_mm(Gp_Pa, L_mm):
    L_m = L_mm * 1e-3
    return 1000.0 * np.sqrt(1.94 * Gp_Pa / (RHO * G_ACC)) * np.sqrt(L_m)


def n_layers_max(h_max_mm, h_layer_mm):
    if h_layer_mm <= 0:
        return 0
    return int(np.floor(h_max_mm / h_layer_mm))


def plot_Gprime_extrap(ink, w, Gp, G0, beta, output_dir):
    plt.figure(figsize=(6.0, 4.5))
    plt.loglog(w, Gp, "o", label="measured")
    w_ext = np.logspace(np.log10(1e-3), np.log10(w.max() * 1.1), 60)
    Gp_ext = G0 * w_ext ** beta
    plt.loglog(w_ext, Gp_ext, "--",
               label=f"fit:  G'={G0:.1f} ω^{beta:.3f}")
    plt.axvline(OMEGA_EXTRAP, color="gray", ls=":", lw=0.8)
    plt.text(OMEGA_EXTRAP*1.15, 1.5*min(Gp),
             f"ω = {OMEGA_EXTRAP}", color="gray", fontsize=9)
    plt.xlabel("Angular frequency ω (rad/s)")
    plt.ylabel("Storage modulus G' (Pa)")
    plt.title(f"{ink}: low-ω extrapolation")
    plt.legend()
    plt.grid(True, which="both", ls=":", alpha=0.5)
    plt.tight_layout()
    fig_path = os.path.join(output_dir, f"Gprime_extrap_{ink}.png")
    plt.savefig(fig_path, dpi=160)
    plt.close()
    return fig_path


def main():
    os.chdir(PROJECT_ROOT)

    out_dir = "Analises/Python/Results"
    os.makedirs(out_dir, exist_ok=True)

    out_lines = []
    def emit(s=""):
        print(s)
        out_lines.append(s)

    emit("=" * 78)
    emit("SAOS-DERIVED MAX SCAFFOLD HEIGHT — extract_hmax_v2.py")
    emit(f"rho={RHO} kg/m^3 | g={G_ACC} m/s^2 | extrapolation target = {OMEGA_EXTRAP} rad/s")
    emit(f"h_layer targets for N_max (mm): {H_LAYER_TARGETS_MM}")
    emit("=" * 78)

    # ---- (1) Amplitude sweep -> gamma_LVR + G'_LVR ----
    gamma_lvr = {}
    Gp_lvr = {}
    emit("\n[1] Amplitude sweep — γ_LVR (10% drop) and G'_LVR")
    for name, path in STRAIN_FILES.items():
        if not os.path.exists(path):
            emit(f"   {name}: FILE NOT FOUND ({path})")
            continue
        info, _, _ = extract_amplitude_sweep(path)
        if info is None or "error" in info:
            emit(f"   {name}: parse error -> {info.get('error') if info else 'no data'}")
            continue
        gamma_lvr[name] = info["gamma_LVR_pct"] / 100.0
        Gp_lvr[name]    = info["Gp_LVR"]
        emit(f"   {name:<8}  G'_LVR={info['Gp_LVR']:8.1f} Pa  "
             f"γ_LVR={info['gamma_LVR_pct']:6.3f}%")

    # ---- (2) Frequency sweep -> measured + extrapolation ----
    Gp_meas, Gp_ext = {}, {}
    fit_params = {}
    emit("\n[2] Frequency sweep — G'(ω=1) measured, G'(ω=0.01) extrapolated")
    for name, path in FREQ_FILES.items():
        if not os.path.exists(path):
            emit(f"   {name}: FILE NOT FOUND ({path})"); continue
        w, gp = load_freq_curve(path)
        idx1 = int(np.argmin(np.abs(np.log10(w) - np.log10(1.0))))
        Gp_meas[name] = float(gp[idx1])
        G0, beta, R2, n_used = fit_low_omega_power_law(w, gp, n_low=6)
        Gp_ext[name] = float(G0 * (OMEGA_EXTRAP ** beta))
        fit_params[name] = (G0, beta, R2, n_used)
        emit(f"   {name:<8} ω∈[{w[0]:.2f},{w[-1]:.1f}] rad/s n={len(w)}")
        emit(f"            G'(ω=1)   = {Gp_meas[name]:8.1f} Pa  (measured)")
        emit(f"            G' = {G0:.2f} ω^{beta:.3f}  (R²={R2:.3f}, n_low={n_used})")
        emit(f"            G'(ω=0.01)= {Gp_ext[name]:8.1f} Pa  (extrapolated)")
        plot_path = plot_Gprime_extrap(name, w, gp, G0, beta, out_dir)
        emit(f"            plot -> {plot_path}")

    # ---- (3) h_max under 3 stacked + buckling ----
    emit("\n[3] h_max [mm] — stacked, full sub-layer support (h = σ_max / (ρ g))")
    emit("")
    header = (f"  {'Ink':<8} | {'σ_LVR':>9} {'h_LVR':>8} | "
              f"{'G(ω=1)':>8} {'h_G1':>8} | "
              f"{'G(ω=.01)':>10} {'h_G01':>8}")
    emit(header); emit("  " + "-"*len(header))
    rows = {}
    for ink in FREQ_FILES.keys():
        if ink not in Gp_meas: continue
        g1, gex = Gp_meas[ink], Gp_ext[ink]
        if ink in gamma_lvr and ink in Gp_lvr:
            sig_lvr = Gp_lvr[ink] * gamma_lvr[ink]
        else:
            sig_lvr = float("nan")
        h_lvr  = h_max_stacked_mm(sig_lvr)
        h_g1   = h_max_stacked_mm(g1)
        h_gex  = h_max_stacked_mm(gex)
        rows[ink] = dict(sig_lvr=sig_lvr, h_lvr=h_lvr,
                         sig_g1=g1, h_g1=h_g1,
                         sig_gex=gex, h_gex=h_gex)
        emit(f"  {ink:<8} | {sig_lvr:>9.1f} {h_lvr:>8.2f} | "
             f"{g1:>8.1f} {h_g1:>8.2f} | "
             f"{gex:>10.1f} {h_gex:>8.2f}")

    # ---- (4) Buckling correction (multiple spans) ----
    emit("\n[4] Buckling-corrected unsupported-span (Smay 2002): "
         "h_max(L) = √(1.94 G' / ρg) · √L")
    emit("")
    bh = "  " + f"{'Ink':<8} | {'G(ω=1) Pa':>10} | " + \
         " ".join(f"h(L={L:.1f}mm)".rjust(10) for L in LAYER_HEIGHTS_FOR_BUCKLE_MM)
    emit(bh); emit("  " + "-"*len(bh))
    for ink, g1 in Gp_meas.items():
        hs = "  " + f"{ink:<8} | {g1:>10.1f} | "
        for L in LAYER_HEIGHTS_FOR_BUCKLE_MM:
            hs += f"{h_max_buckle_mm(g1, L):>10.2f} "
        emit(hs)

    # ---- (5) N_max layer count predictions ----
    emit("\n[5] Predicted maximum layer count N_max  (= floor(h_max / h_layer))")
    emit(f"    For typical DIW layer heights {H_LAYER_TARGETS_MM} mm:")
    emit("")
    nh = "  " + f"{'Ink':<8} | {'criterion':<10} | " + \
         " ".join(f"h_layer={hL:.2f}mm".rjust(15) for hL in H_LAYER_TARGETS_MM)
    emit(nh); emit("  " + "-"*len(nh))
    for ink, r in rows.items():
        for crit, h in (("LVR", r["h_lvr"]), ("G(ω=1)", r["h_g1"]),
                        ("G(ω=.01)", r["h_gex"])):
            cells = " ".join(f"{n_layers_max(h, hL):>15d}"
                             for hL in H_LAYER_TARGETS_MM)
            emit(f"  {ink:<8} | {crit:<10} | {cells}")
        emit("")

    # ---- (6) Save machine-readable CSV ----
    csv_path = os.path.join(out_dir, "printing_parameters_per_ink.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "ink",
            "Gp_LVR_Pa", "gamma_LVR_strain",
            "Gp_omega1_Pa", "beta_lowomega", "R2_fit",
            "Gp_omega001_Pa",
            "sigma_max_LVR_Pa",     "h_max_LVR_mm",
            "sigma_max_G1_Pa",      "h_max_G1_mm",
            "sigma_max_G001_Pa",    "h_max_G001_mm",
            "h_buckle_L0p5mm",      "h_buckle_L1mm",      "h_buckle_L2mm",
            "Nmax_LVR_hL0p3", "Nmax_LVR_hL0p4", "Nmax_LVR_hL0p5",
            "Nmax_G1_hL0p3",  "Nmax_G1_hL0p4",  "Nmax_G1_hL0p5",
            "Nmax_G001_hL0p3","Nmax_G001_hL0p4","Nmax_G001_hL0p5",
        ])
        for ink in FREQ_FILES.keys():
            if ink not in rows: continue
            r  = rows[ink]
            g1 = Gp_meas[ink]
            G0, beta, R2, _ = fit_params[ink]
            row = [
                ink,
                Gp_lvr.get(ink, float("nan")),
                gamma_lvr.get(ink, float("nan")),
                g1, beta, R2,
                Gp_ext[ink],
                r["sig_lvr"], r["h_lvr"],
                r["sig_g1"],  r["h_g1"],
                r["sig_gex"], r["h_gex"],
                h_max_buckle_mm(g1, 0.5),
                h_max_buckle_mm(g1, 1.0),
                h_max_buckle_mm(g1, 2.0),
                n_layers_max(r["h_lvr"], 0.3),
                n_layers_max(r["h_lvr"], 0.4),
                n_layers_max(r["h_lvr"], 0.5),
                n_layers_max(r["h_g1"],  0.3),
                n_layers_max(r["h_g1"],  0.4),
                n_layers_max(r["h_g1"],  0.5),
                n_layers_max(r["h_gex"], 0.3),
                n_layers_max(r["h_gex"], 0.4),
                n_layers_max(r["h_gex"], 0.5),
            ]
            w.writerow([round(x, 4) if isinstance(x, float) else x for x in row])
    emit(f"\n[Saved CSV] {csv_path}")

    # ---- (7) Save report ----
    out_path = os.path.join(out_dir, "SAOS_hmax_v2.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"[Saved report] {out_path}")


if __name__ == "__main__":
    main()
