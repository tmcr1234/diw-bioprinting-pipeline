"""
flow_diagnostics.py
===================

Shape diagnostics for steady-shear flow curves, run BEFORE fitting.

Two artefacts have each cost this group a round of external questions, and
neither was detectable from a fit statistic alone:

A7  STARTUP TRANSIENT. The first stretch of a shear ramp is not equilibrium
    data — the instrument is still ramping and the structure is still
    breaking down, so viscosity RISES with shear rate there instead of
    falling. Including that stretch shifted fitted eta0 by 19-31% on this
    lab's CMC series. It is not noise and it is not reproducible run to
    run; it must be identified and reported, not averaged over.

A8  EDGE FRACTURE. At high shear the sample can tear away from the fixture
    edge. Shear stress then DROPS discontinuously while shear rate keeps
    rising — the measurement has stopped being a measurement. On the C20
    ink this happens at 9-15 1/s (tau 510 -> 415 Pa, N1 2527 -> 59 Pa) and
    it is why that ink's Power-Law fit sits at R^2 = 0.505. Nothing about
    a low R^2 tells you WHICH points are bad.

Neither function edits data. They return a report; the caller decides.
That is deliberate: silently trimming points is how a fit becomes
irreproducible.

Author: T.M.C. Rodrigues - 2026-08-26
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class TransientReport:
    detected: bool
    n_excluded: int          # points to drop from the head of the curve
    gamma_dot_cut: float     # shear rate at the cut, NaN if none
    rise_fraction: float     # how much eta rises across the suspect head
    message: str


@dataclass
class FractureReport:
    detected: bool
    index_range: tuple       # (i0, i1) inclusive, or (-1,-1)
    gamma_dot_range: tuple   # (lo, hi), or (nan, nan)
    max_drop_fraction: float # largest single-step relative stress drop
    message: str


def detect_startup_transient(
    shear_rate,
    viscosity,
    min_rise: float = 0.02,
    max_head_fraction: float = 0.5,
) -> TransientReport:
    """
    Flag a non-monotonic low-shear head on a flow curve.

    A shear-thinning ink must have eta falling monotonically with gamma_dot.
    Any leading stretch where eta RISES is instrument startup, not rheology.

    The cut is placed at the LAST index of the rising head — i.e. at the
    viscosity maximum — so that everything at or below it is excluded.

    Parameters
    ----------
    min_rise : relative rise required before the head is called a transient.
        Defaults to 2%, comfortably above point-to-point scatter.
    max_head_fraction : refuse to call more than this fraction of the curve
        a transient. If the "head" is most of the curve, the ink is not
        simply shear-thinning and the diagnosis belongs to a human.

    Returns
    -------
    TransientReport. ``n_excluded`` is a SUGGESTION; the caller reports both
    fits rather than substituting one for the other.
    """
    gd = np.asarray(shear_rate, dtype=float)
    eta = np.asarray(viscosity, dtype=float)
    ok = np.isfinite(gd) & np.isfinite(eta) & (gd > 0) & (eta > 0)
    gd, eta = gd[ok], eta[ok]
    order = np.argsort(gd)
    gd, eta = gd[order], eta[order]

    if gd.size < 8:
        return TransientReport(False, 0, float("nan"), 0.0,
                               "too few points to diagnose a transient")

    imax = int(np.argmax(eta))
    if imax == 0:
        return TransientReport(False, 0, float("nan"), 0.0,
                               "viscosity is monotonically decreasing — no startup head")

    rise = (eta[imax] - eta[0]) / eta[0]
    head_frac = (imax + 1) / gd.size

    if rise < min_rise:
        return TransientReport(
            False, 0, float("nan"), rise,
            f"low-shear head rises only {rise*100:.1f}% (< {min_rise*100:.0f}%) — "
            f"treating as scatter, not a transient")

    if head_frac > max_head_fraction:
        return TransientReport(
            True, imax + 1, float(gd[imax]), rise,
            f"[!] The rising head is {head_frac*100:.0f}% of the curve "
            f"({imax+1}/{gd.size} points, +{rise*100:.1f}% in eta). That is too "
            f"much of the measurement to call startup. Do NOT auto-exclude — "
            f"this ink may be rheopexic or shear-structuring, which is a "
            f"result, not an artefact. Compare Flow vs Flow-Sheared first.")

    return TransientReport(
        True, imax + 1, float(gd[imax]), rise,
        f"[!] Startup transient: eta rises {rise*100:.1f}% over the first "
        f"{imax+1} of {gd.size} points, up to gamma_dot = {gd[imax]:.4g} 1/s. "
        f"Those points are not equilibrium data. Both fits are reported below; "
        f"state which one the manuscript uses.")


def detect_edge_fracture(
    shear_rate,
    shear_stress,
    normal_stress=None,
    min_drop: float = 0.05,
) -> FractureReport:
    """
    Flag a discontinuous stress DROP at rising shear rate — edge fracture.

    In steady shear, tau must increase monotonically with gamma_dot for any
    inelastic or shear-thinning fluid. A fall means the sample is no longer
    fully in the gap.

    If ``normal_stress`` (N1) is supplied it is used as corroboration only:
    a simultaneous collapse of N1 makes edge fracture near-certain, and is
    said so in the message. Its absence never suppresses the flag.

    Parameters
    ----------
    min_drop : relative single-step fall in tau that counts. 5% is well
        above instrument scatter on these gels and well below the ~19% drop
        seen on the C20 fracture.
    """
    gd = np.asarray(shear_rate, dtype=float)
    tau = np.asarray(shear_stress, dtype=float)
    ok = np.isfinite(gd) & np.isfinite(tau) & (gd > 0) & (tau > 0)
    idx = np.flatnonzero(ok)
    gd_o, tau_o = gd[ok], tau[ok]
    order = np.argsort(gd_o)
    gd_o, tau_o, idx = gd_o[order], tau_o[order], idx[order]

    if gd_o.size < 5:
        return FractureReport(False, (-1, -1), (float("nan"), float("nan")), 0.0,
                              "too few points to diagnose edge fracture")

    rel_step = np.diff(tau_o) / tau_o[:-1]
    bad = np.flatnonzero(rel_step < -min_drop)
    if bad.size == 0:
        largest_fall = float(max(0.0, -rel_step.min()))
        return FractureReport(
            False, (-1, -1), (float("nan"), float("nan")), largest_fall,
            f"shear stress rises monotonically within tolerance "
            f"(largest fall {largest_fall*100:.1f}% < {min_drop*100:.0f}%)")

    i0, i1 = int(bad[0]), int(bad[-1] + 1)
    worst = float(-rel_step[bad].min())

    corroboration = ""
    if normal_stress is not None:
        n1 = np.asarray(normal_stress, dtype=float)[ok][order]
        if np.isfinite(n1[i0]) and np.isfinite(n1[i1]) and n1[i0] > 0:
            n1_drop = (n1[i0] - n1[i1]) / n1[i0]
            if n1_drop > 0.5:
                corroboration = (
                    f" N1 collapses across the same window "
                    f"({n1[i0]:.4g} -> {n1[i1]:.4g} Pa, -{n1_drop*100:.0f}%), "
                    f"which makes edge fracture near-certain rather than suspected.")
            else:
                corroboration = (
                    f" N1 does NOT collapse here ({n1[i0]:.4g} -> {n1[i1]:.4g} Pa), "
                    f"so consider slip or a genuine flow instability too.")

    return FractureReport(
        True, (i0, i1), (float(gd_o[i0]), float(gd_o[i1])), worst,
        f"[!] Possible edge fracture near gamma_dot = {gd_o[i0]:.4g}-{gd_o[i1]:.4g} 1/s: "
        f"shear stress falls {worst*100:.1f}% while shear rate rises "
        f"({tau_o[i0]:.4g} -> {tau_o[i1]:.4g} Pa). Points in that window are not "
        f"a measurement of the material.{corroboration} Exclude them by hand and "
        f"say so — a low R^2 on the full range does not tell a reader WHICH "
        f"points were bad.")
