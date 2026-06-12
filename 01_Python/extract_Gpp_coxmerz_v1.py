"""
extract_Gpp_coxmerz_v1.py
=========================
Loss-modulus (G'') analyses that were missing from the G'-only SAOS pipeline
(extract_hmax_v2.py / extract_SAOS_values.py). Three corrections:

  (A) tan(delta)(omega) over the frequency sweep + rigorous Winter-Chambon
      gel-point test. At a critical gel  G' ~ G'' ~ omega^beta  and
      tan(delta) = tan(beta*pi/2)  is FREQUENCY-INDEPENDENT.
      -> fit beta' from G'(low-omega) and beta'' from G''(low-omega);
         report whether G' and G'' are parallel power laws, the measured
         tan(delta) and its flatness, and the predicted tan(beta'*pi/2).

  (B) Cox-Merz comparison:  eta(gammadot) [flow curve, Ramp 1] vs
      |eta*(omega)| = |G*|/omega [frequency sweep].  Liquids obey
      eta(gammadot) ~= |eta*(omega)| at gammadot = omega; structured /
      associating fluids deviate (|eta*| > eta). Deviation scales with
      the strength of the at-rest structure.

  (C) Weak-strain-overshoot / Payne test (amplitude sweep): does G'' rise
      above its LVR plateau before the G'/G'' crossover? A G'' overshoot is
      the LAOS type-III fingerprint of an attractive / associating network.

Reuses the validated AresG2/Anton-Paar .xls parser from extract_SAOS_values.py.

Outputs:
  - Analises/Python/Results/SAOS_Gpp_coxmerz_v1.txt      (full report)
  - Analises/Python/Results/tandelta_coxmerz_per_ink.csv (machine-readable)
  - Latex/Figs/F8_tandelta.pdf   (+ .png preview in Results/)
  - Latex/Figs/F9_coxmerz.pdf    (+ .png preview)
  - Latex/Figs/F10_Gpp_overshoot.pdf (+ .png preview)

Status registry (Section 6.9): extract_Gpp_coxmerz_v1.py | v1 | active | new

Author: T.M.C. Rodrigues -- 2026-06-11
"""

import os
import sys
import csv
import numpy as np
import xlrd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract_SAOS_values import read_antpar_xls, find_col, extract_amplitude_sweep

PROJECT_ROOT = "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE"
FIGS = os.path.join(PROJECT_ROOT, "Latex", "Figs")
RESULTS = os.path.join(PROJECT_ROOT, "Analises", "Python", "Results")

FREQ_FILES = {
    "C15":    "Reologia/Frequency/C15 - Freq CP50 0,1mm.xls",
    "NE":     "Reologia/Frequency/C15 Gira 5.5  - Freq CP50 0,1mm.xls",
    "Bozano": "Reologia/Frequency/Gel de Cabelo Bozano - Freq CP50 0,1mm.xls",
}
FLOW_FILES = {
    "C15":    "Reologia/Viscosity/C15 - Visco CP50 0,1mm.xls",
    "NE":     "Reologia/Viscosity/C15 Gira 5.5 - Visco CP50 0,1mm.xls",
    "Bozano": "Reologia/Viscosity/Gel de Cabelo Bozano - Visco CP50 0,1mm.xls",
}
STRAIN_FILES = {
    "C15":    "Reologia/Strain/C15 - Deformacao CP50 0,1mm.xls",
    "NE":     "Reologia/Strain/C15 Gira 5.5  - Deformacao CP50 0,1mm.xls",
    "Bozano": "Reologia/Strain/Gel de Cabelo Bozano - Deformacao CP50 0,1mm.xls",
}

COLORS = {"C15": "#1f77b4", "NE": "#ff7f0e", "Bozano": "#2ca02c"}
LABELS = {"C15": "C15", "NE": "C15 + Gira 5.5", "Bozano": "Bozano"}
N_LOW = 6  # low-omega / low-strain points used for the power-law fits


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------
def load_freq_full(path):
    """omega, G', G'', tan_delta, |eta*| from a frequency sweep .xls."""
    headers, units, arr, _ = read_antpar_xls(path)
    i_w = find_col(headers, "Angular frequency", "Angular Frequency", "Omega", "omega")
    i_gp = find_col(headers, "Storage modulus", "G'")
    i_gpp = find_col(headers, "Loss modulus", "G''", "G\"")
    if None in (i_w, i_gp, i_gpp):
        raise RuntimeError(f"freq columns not found in {path}: {headers}")
    w, gp, gpp = arr[:, i_w], arr[:, i_gp], arr[:, i_gpp]
    v = (~np.isnan(w) & ~np.isnan(gp) & ~np.isnan(gpp) & (w > 0) & (gp > 0) & (gpp > 0))
    w, gp, gpp = w[v], gp[v], gpp[v]
    o = np.argsort(w)
    w, gp, gpp = w[o], gp[o], gpp[o]
    tand = gpp / gp
    etastar = np.sqrt(gp ** 2 + gpp ** 2) / w
    return w, gp, gpp, tand, etastar


def load_flow(path):
    """shear_rate, viscosity from a flow-ramp .xls (Ramp 1, standard)."""
    wb = xlrd.open_workbook(path)
    for si in range(wb.nsheets):
        sh = wb.sheet_by_index(si)
        for ri in range(sh.nrows):
            row = [str(sh.cell_value(ri, ci)).strip() for ci in range(sh.ncols)]
            has_sr = any("Shear rate" in c or "Shear Rate" in c for c in row)
            has_vi = any("Viscosity" in c for c in row)
            if has_sr and has_vi:
                i_sr = next(k for k, c in enumerate(row) if "Shear rate" in c or "Shear Rate" in c)
                i_vi = next(k for k, c in enumerate(row) if "Viscosity" in c)
                data = []
                for rr in range(ri + 1, sh.nrows):
                    try:
                        sr = float(sh.cell_value(rr, i_sr))
                        vi = float(sh.cell_value(rr, i_vi))
                        data.append((sr, vi))
                    except (ValueError, TypeError):
                        continue
                a = np.array(data)
                m = (a[:, 0] > 0) & (a[:, 1] > 0) & np.isfinite(a[:, 0]) & np.isfinite(a[:, 1])
                a = a[m]
                o = np.argsort(a[:, 0])
                return a[o, 0], a[o, 1]
    raise RuntimeError(f"flow columns not found in {path}")


def load_strain_full(path):
    """strain(%), G', G'' from an amplitude sweep .xls."""
    headers, units, arr, _ = read_antpar_xls(path)
    i_s = find_col(headers, "Oscillation strain", "Deformation", "Deforma", "Strain")
    i_gp = find_col(headers, "Storage modulus", "G'")
    i_gpp = find_col(headers, "Loss modulus", "G''", "G\"")
    if None in (i_s, i_gp, i_gpp):
        raise RuntimeError(f"strain columns not found in {path}: {headers}")
    s, gp, gpp = arr[:, i_s], arr[:, i_gp], arr[:, i_gpp]
    v = (~np.isnan(s) & ~np.isnan(gp) & ~np.isnan(gpp) & (s > 0) & (gp > 0) & (gpp > 0))
    s, gp, gpp = s[v], gp[v], gpp[v]
    o = np.argsort(s)
    return s[o], gp[o], gpp[o]


def fit_power_low(x, y, n=N_LOW):
    """log-log slope (exponent) + prefactor over the n lowest-x points."""
    n = min(n, len(x))
    lx, ly = np.log10(x[:n]), np.log10(y[:n])
    A = np.vstack([lx, np.ones_like(lx)]).T
    (b, c), *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = c + b * lx
    ssr = float(np.sum((ly - pred) ** 2))
    sst = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 - ssr / sst if sst > 0 else float("nan")
    return float(b), float(10 ** c), float(r2)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    os.chdir(PROJECT_ROOT)
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)

    out = []
    def emit(s=""):
        print(s)
        out.append(s)

    emit("=" * 78)
    emit("  G'' (LOSS MODULUS) ANALYSES  --  extract_Gpp_coxmerz_v1.py")
    emit("  tan(delta) / Winter-Chambon  |  Cox-Merz  |  G'' weak-strain-overshoot")
    emit("=" * 78)

    rows = {}  # per-ink summary for CSV

    # ============================ (A) tan(delta) / Winter-Chambon ============
    emit("\n[A] FREQUENCY SWEEP: tan(delta)(omega) and Winter-Chambon gel-point test")
    emit("    Critical gel  =>  G' ~ G'' ~ w^beta  with tan(delta)=tan(beta*pi/2) "
         "frequency-independent.\n")
    emit(f"  {'Ink':<14} {'beta_Gp':>8} {'beta_Gpp':>9} {'tand(w=1)':>10} "
         f"{'tand_low_avg':>13} {'tan(b*pi/2)':>12} {'tand slope':>11}")
    emit("  " + "-" * 80)
    freq = {}
    for ink, path in FREQ_FILES.items():
        w, gp, gpp, tand, etastar = load_freq_full(path)
        freq[ink] = (w, gp, gpp, tand, etastar)
        b_gp, _, _ = fit_power_low(w, gp)
        b_gpp, _, _ = fit_power_low(w, gpp)
        # tan delta at w~1 and averaged over the low-omega window
        i1 = int(np.argmin(np.abs(np.log10(w) - 0.0)))
        tand_w1 = float(tand[i1])
        tand_low = float(np.mean(tand[:N_LOW]))
        tand_pred = float(np.tan(b_gp * np.pi / 2.0))
        # flatness: slope of tan delta vs log10(omega) over low-omega window
        slope, _, _ = fit_power_low(w, tand) if np.all(tand[:N_LOW] > 0) else (float("nan"),) * 3
        # (fit_power_low returns log-log slope; for tan-delta-vs-logw use linear-in-logw)
        ll = np.log10(w[:N_LOW])
        slope_lin = float(np.polyfit(ll, tand[:N_LOW], 1)[0])
        emit(f"  {LABELS[ink]:<14} {b_gp:>8.3f} {b_gpp:>9.3f} {tand_w1:>10.3f} "
             f"{tand_low:>13.3f} {tand_pred:>12.3f} {slope_lin:>11.3f}")
        rows.setdefault(ink, {}).update(
            beta_Gp=b_gp, beta_Gpp=b_gpp, tand_w1=tand_w1,
            tand_low_avg=tand_low, tand_pred_WC=tand_pred, tand_slope=slope_lin)
    emit("\n  Reading: tan(delta) flat in omega AND beta_Gp ~ beta_Gpp  => critical-gel-like.")
    emit("  C15/composite expected tan(delta)~0.85 (soft transient gel); Bozano <<1 (permanent).")

    # ============================ (B) Cox-Merz ==============================
    emit("\n[B] COX-MERZ:  eta(gammadot) [flow]  vs  |eta*(omega)| [SAOS]")
    emit("    Ratio |eta*|/eta over the overlapping decade(s); >1 = structured fluid.\n")
    emit(f"  {'Ink':<14} {'overlap (1/s)':>16} {'mean |eta*|/eta':>16} {'max |eta*|/eta':>15}")
    emit("  " + "-" * 64)
    flow = {}
    for ink, path in FLOW_FILES.items():
        try:
            sr, eta = load_flow(path)
        except Exception as e:
            emit(f"  {LABELS[ink]:<14}  FLOW LOAD FAILED: {e}")
            continue
        flow[ink] = (sr, eta)
        w, gp, gpp, tand, etastar = freq[ink]
        lo = max(sr.min(), w.min())
        hi = min(sr.max(), w.max())
        msk = (w >= lo) & (w <= hi)
        wo = w[msk]
        eta_at_w = 10 ** np.interp(np.log10(wo), np.log10(sr), np.log10(eta))
        ratio = etastar[msk] / eta_at_w
        rmean, rmax = float(np.mean(ratio)), float(np.max(ratio))
        emit(f"  {LABELS[ink]:<14} {f'{lo:.2g}-{hi:.2g}':>16} {rmean:>16.2f} {rmax:>15.2f}")
        rows.setdefault(ink, {}).update(
            coxmerz_overlap_lo=lo, coxmerz_overlap_hi=hi,
            coxmerz_ratio_mean=rmean, coxmerz_ratio_max=rmax)
    emit("\n  Reading: a larger |eta*|/eta => stronger at-rest structure not probed by")
    emit("  steady shear (depletion network). Compare composite vs C15.")

    # ============================ (C) G'' overshoot =========================
    emit("\n[C] AMPLITUDE SWEEP: G'' weak-strain-overshoot (Payne / LAOS type III)")
    emit("    Peak G'' above its LVR plateau before the G'/G'' crossover = attractive network.\n")
    emit(f"  {'Ink':<14} {'Gpp_LVR(Pa)':>12} {'Gpp_peak(Pa)':>13} {'peak/LVR':>9} "
         f"{'gamma_peak(%)':>14} {'cross gamma(%)':>15}")
    emit("  " + "-" * 82)
    strain = {}
    for ink, path in STRAIN_FILES.items():
        try:
            s, gp, gpp = load_strain_full(path)
            info, _, _ = extract_amplitude_sweep(path)  # robust LVR plateau
        except Exception as e:
            emit(f"  {LABELS[ink]:<14}  STRAIN LOAD FAILED: {e}")
            continue
        strain[ink] = (s, gp, gpp)
        # robust LVR baseline (sliding-window plateau, low-torque noise excluded)
        gpp_lvr = float(info["Gpp_LVR"])
        g_lvr = float(info["gamma_LVR_pct"])
        # the physically meaningful G'' overshoot is in the YIELDING region
        # (gamma >= gamma_LVR); the low-strain region is low-torque/inertial noise
        yld = s >= g_lvr
        if np.any(yld):
            sy, gppy = s[yld], gpp[yld]
            j = int(np.argmax(gppy))
            gpp_peak = float(gppy[j]); gamma_peak = float(sy[j])
        else:
            j = int(np.argmax(gpp)); gpp_peak = float(gpp[j]); gamma_peak = float(s[j])
        ratio = gpp_peak / gpp_lvr
        # G'/G'' crossover (flow point) searched only at/above the LVR
        cross = float("nan")
        for k in range(1, len(s)):
            if s[k] >= g_lvr and gpp[k] >= gp[k] and gpp[k - 1] < gp[k - 1]:
                x0, x1 = np.log10(s[k - 1]), np.log10(s[k])
                d0, d1 = gp[k - 1] - gpp[k - 1], gp[k] - gpp[k]
                frac = d0 / (d0 - d1) if (d0 - d1) != 0 else 0.0
                cross = float(10 ** (x0 + frac * (x1 - x0)))
                break
        emit(f"  {LABELS[ink]:<14} {gpp_lvr:>12.1f} {gpp_peak:>13.1f} {ratio:>9.3f} "
             f"{gamma_peak:>14.2f} {cross:>15.2f}")
        rows.setdefault(ink, {}).update(
            Gpp_LVR=gpp_lvr, Gpp_peak=gpp_peak, Gpp_overshoot_ratio=ratio,
            gamma_Gpp_peak_pct=gamma_peak, crossover_gamma_pct=cross)
    emit("\n  Reading: peak/LVR > ~1.05 with the peak located near the end of the LVR")
    emit("  is a weak-strain-overshoot, the signature of a yielding attractive network.")

    # ============================ CSV =======================================
    csv_path = os.path.join(RESULTS, "tandelta_coxmerz_per_ink.csv")
    fields = ["ink", "beta_Gp", "beta_Gpp", "tand_w1", "tand_low_avg",
              "tand_pred_WC", "tand_slope",
              "coxmerz_overlap_lo", "coxmerz_overlap_hi",
              "coxmerz_ratio_mean", "coxmerz_ratio_max",
              "Gpp_LVR", "Gpp_peak", "Gpp_overshoot_ratio",
              "gamma_Gpp_peak_pct", "crossover_gamma_pct"]
    with open(csv_path, "w", newline="") as f:
        wtr = csv.writer(f)
        wtr.writerow(fields)
        for ink in FREQ_FILES:
            r = rows.get(ink, {})
            wtr.writerow([ink] + [round(r.get(k, float("nan")), 5)
                                  if isinstance(r.get(k), float) else r.get(k, "")
                                  for k in fields[1:]])
    emit(f"\n[Saved CSV] {csv_path}")

    # ============================ FIGURES ===================================
    # F8: tan(delta)(omega) + Winter-Chambon reference lines
    plt.figure(figsize=(6.4, 4.6))
    for ink in FREQ_FILES:
        w, gp, gpp, tand, etastar = freq[ink]
        plt.semilogx(w, tand, "o-", ms=3.5, color=COLORS[ink], label=LABELS[ink])
        pred = rows[ink]["tand_pred_WC"]
        plt.axhline(pred, color=COLORS[ink], ls="--", lw=0.9, alpha=0.7)
    plt.xlabel(r"Angular frequency $\omega$ (rad s$^{-1}$)")
    plt.ylabel(r"$\tan\delta = G''/G'$")
    plt.title(r"Loss tangent vs frequency (dashed = $\tan(\beta\pi/2)$ Winter--Chambon)")
    plt.grid(True, which="both", ls=":", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    for ext, d in ((".pdf", FIGS), (".png", RESULTS)):
        plt.savefig(os.path.join(d, "F8_tandelta" + ext), dpi=170)
    plt.close()

    # F9: Cox-Merz
    plt.figure(figsize=(6.4, 4.6))
    for ink in FREQ_FILES:
        if ink in flow:
            sr, eta = flow[ink]
            plt.loglog(sr, eta, "-", color=COLORS[ink], lw=1.6,
                       label=f"{LABELS[ink]}  $\\eta(\\dot\\gamma)$")
        w, gp, gpp, tand, etastar = freq[ink]
        plt.loglog(w, etastar, "o", ms=4, mfc="none", color=COLORS[ink],
                   label=f"{LABELS[ink]}  $|\\eta^*(\\omega)|$")
    plt.xlabel(r"$\dot\gamma$ (s$^{-1}$)  /  $\omega$ (rad s$^{-1}$)")
    plt.ylabel(r"$\eta,\ |\eta^*|$ (Pa s)")
    plt.title(r"Cox--Merz: steady $\eta(\dot\gamma)$ (lines) vs $|\eta^*(\omega)|$ (open)")
    plt.grid(True, which="both", ls=":", alpha=0.5)
    plt.legend(fontsize=7, ncol=1)
    plt.tight_layout()
    for ext, d in ((".pdf", FIGS), (".png", RESULTS)):
        plt.savefig(os.path.join(d, "F9_coxmerz" + ext), dpi=170)
    plt.close()

    # F10: amplitude sweep G', G'' (overshoot)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=False)
    for ax, ink in zip(axes, FREQ_FILES):
        if ink not in strain:
            continue
        s, gp, gpp = strain[ink]
        ax.loglog(s, gp, "o-", ms=3, color=COLORS[ink], label="G'")
        ax.loglog(s, gpp, "s--", ms=3, color=COLORS[ink], alpha=0.6, label="G''")
        gpp_lvr = rows[ink]["Gpp_LVR"]
        ax.axhline(gpp_lvr, color="gray", ls=":", lw=0.8)
        gp_pk = rows[ink]["gamma_Gpp_peak_pct"]
        ax.axvline(gp_pk, color=COLORS[ink], ls=":", lw=0.8, alpha=0.6)
        ax.set_title(f"{LABELS[ink]}  (G'' peak/LVR={rows[ink]['Gpp_overshoot_ratio']:.2f})",
                     fontsize=9)
        ax.set_xlabel(r"strain $\gamma$ (%)")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("G', G'' (Pa)")
    fig.suptitle("Amplitude sweep: G'' weak-strain-overshoot test", y=1.02)
    fig.tight_layout()
    for ext, d in ((".pdf", FIGS), (".png", RESULTS)):
        fig.savefig(os.path.join(d, "F10_Gpp_overshoot" + ext), dpi=170, bbox_inches="tight")
    plt.close(fig)

    emit(f"[Saved figures] {FIGS}/F8_tandelta.pdf, F9_coxmerz.pdf, F10_Gpp_overshoot.pdf")

    # ============================ report ====================================
    rep = os.path.join(RESULTS, "SAOS_Gpp_coxmerz_v1.txt")
    with open(rep, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"[Saved report] {rep}")


if __name__ == "__main__":
    main()
