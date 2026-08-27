#!/usr/bin/env python3
"""
extract_N1_tanner.py
====================

Turn a force-augmented flow curve into the one number the v5 solver needs:
``N1_wall_Pa``, the first normal-stress difference at the needle wall shear
rate — and report the die-swell closure it implies.

WHERE THIS SITS
---------------
    force-augmented export (.txt)
        |  antpar_io.read_flow_curve_with_force   (latin-1, N1 = 2F/piR^2,
        |                                          min-baseline tare)
        v
    N1(gamma_dot), tared
        |  THIS SCRIPT: exclude edge fracture, interpolate to gamma_dot_wall
        v
    N1_wall_Pa  ->  inks_local.m
        |
        v
    run_solver_v5  ->  bioprinting_algorithm_v5  ->  Tanner beta  ->  k_flow

Without this step ``bioprinting_algorithm_v5`` has no N1 and falls back to
the v4 heuristic every time, which defeats the point of having it.

THE TWO GUARDS THAT MATTER
--------------------------
1. NO SILENT EXTRAPOLATION. The needle wall shear rate (~165-218 1/s at the
   nominal operating point) is often ABOVE the top of the rheometer flow
   curve. Extrapolating N1 past the measured range is how a die-swell number
   becomes fiction. This script REFUSES by default and tells you by how much
   you are short. ``--allow-extrapolation`` overrides it, loudly, and stamps
   every downstream number as extrapolated.

2. THE TARE IS BIASED DOWNWARD, AND BY A KNOWN AMOUNT. The min-baseline
   convention treats the quietest point of the run as zero normal stress.
   On a ramp that point is the lowest shear rate, where N1 is small but not
   zero, so every N1 comes out short by exactly N1(gamma_dot_min). On a
   synthetic power-law check this was 26.17 Pa out of 1590.75 — a 1.6%
   under-estimate that propagates straight into beta, always downward.
   The reported ``baseline_floor_Pa`` lets you see how big the subtraction
   was. Start the flow curve as low as the torque resolution allows.

3. EDGE FRACTURE IS EXCLUDED BEFORE INTERPOLATING. Once the sample tears
   away from the fixture, N1 collapses (2527 -> 59 Pa on the known C20 case)
   and both tau and N1 stop describing the material. Interpolating through
   that region silently drags the Tanner ratio toward zero.
   ``flow_diagnostics.detect_edge_fracture`` finds it, using the N1 collapse
   itself as corroboration.

WHICH tau THE RATIO USES
------------------------
Tanner needs N1_w/(2*tau_w), both at the same shear rate. There are two
defensible choices and they differ:

  * RHEOMETER tau, interpolated from the same curve as N1. This is the
    material's recoverable strain, self-consistently measured.
  * SOLVER tau_w at the needle wall, which carries the Rabinowitsch
    correction and the actual flow geometry.

``bioprinting_algorithm_v5`` uses the SOLVER tau_w, because the swell it is
predicting happens in the needle, not in the rheometer. This script reports
BOTH so the gap between them is visible rather than assumed away. If they
disagree badly, that is a finding about your flow-curve fit, not a rounding
issue.

USAGE
-----
    # gamma_dot_wall taken straight from a solver run
    python3 Export/01_Python/extract_N1_tanner.py \\
        --force-folder "./Reologia/Flow with force" \\
        --master-summary output_v4/master_summary_v4.csv \\
        --needle 21G --vp 0.01

    # or state the wall shear rate yourself
    python3 Export/01_Python/extract_N1_tanner.py \\
        --force-folder "./Reologia/Flow with force" --gamma-w 200

    # sanity-check the closure maths with no data present
    python3 Export/01_Python/extract_N1_tanner.py --self-test

Author: T.M.C. Rodrigues (PEMM/COPPE/UFRJ) - 2026-08-26
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))

from antpar_io import read_flow_curve_with_force, GEOMETRY_RADIUS_M  # noqa: E402
from flow_diagnostics import detect_edge_fracture  # noqa: E402

try:
    from data_config import load as _load_data_config
except Exception:                                    # pragma: no cover
    _load_data_config = None


# ----------------------------------------------------------------------
# The closure itself
# ----------------------------------------------------------------------

NICKELL_FLOOR = 0.13
"""Inelastic die-swell floor. Nickell, Tanner & Yamada 1974,
J Fluid Mech 65(1):189-206, doi:10.1017/S0022112074001339. This is the value
Tanner's closure returns as N1 -> 0, which is why the two are one closure
rather than two competing ones."""


def tanner_beta(N1_w: float, tau_w: float) -> float:
    """
    Die-swell ratio minus one, from the measured first normal-stress
    difference. Tanner 1970, J Polym Sci A-2 8(12):2067-2078,
    doi:10.1002/pol.1970.160081203.

        d/D  = 0.13 + [ 1 + 0.5*( N1_w / (2*tau_w) )^2 ]^(1/6)
        beta = d/D - 1

    Kept byte-identical to the MATLAB local function of the same name in
    bioprinting_algorithm_v5.m. If you change one, change both, and re-run
    validate_v5.

    FREE-JET DEFINITION. Tanner describes a jet leaving a long capillary
    into air. A DIW road is deposited against a substrate at a finite
    stand-off and is still swelling when it lands, so the measured road
    width is expected to exceed this estimate. Report both; do not tune the
    0.13 or the exponent to close the gap.
    """
    if not math.isfinite(tau_w) or tau_w <= 0:
        raise ValueError(f"tau_w must be finite and positive, got {tau_w}")
    if N1_w < 0:
        raise ValueError(
            f"N1_w = {N1_w:g} Pa is negative. A negative first normal-stress "
            f"difference in steady shear on these gels is essentially always "
            f"an un-tared transducer, not a measurement. The tare is applied "
            f"by antpar_io.baseline_correct(method='min') — check it ran.")
    ratio = N1_w / (2.0 * tau_w)
    return NICKELL_FLOOR + (1.0 + 0.5 * ratio ** 2) ** (1.0 / 6.0) - 1.0


def heuristic_beta(n: float) -> float:
    """The v4 closure, for comparison only: beta = 0.30*(1-n). No source."""
    return max(0.0, 0.30 * (1.0 - n))


def k_flow(beta: float, Rrec_pct: float, f_slip: float = 1.0) -> float:
    """Deposition efficiency. Extrusion Multiplier = 1/k_flow."""
    return (1.0 + beta) ** 2 * f_slip * math.sqrt(Rrec_pct / 100.0)


# ----------------------------------------------------------------------
# Interpolation with the two guards
# ----------------------------------------------------------------------

class ExtrapolationRefused(RuntimeError):
    """The requested wall shear rate lies outside the measured range."""


def interpolate_at(gamma_dot, values, gamma_target: float,
                   allow_extrapolation: bool = False, label: str = "value"):
    """
    Log-log interpolate one channel to a target shear rate.

    Log-log because both N1 and tau are power-law-like over the working
    range; linear interpolation across a decade of shear rate is a much
    worse approximation than it looks on a linear plot.

    Returns (value, was_extrapolated).
    """
    gd = np.asarray(gamma_dot, dtype=float)
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(gd) & np.isfinite(v) & (gd > 0) & (v > 0)
    gd, v = gd[ok], v[ok]
    if gd.size < 2:
        raise ValueError(f"{label}: fewer than two usable points to interpolate")
    order = np.argsort(gd)
    gd, v = gd[order], v[order]

    lo, hi = gd[0], gd[-1]
    outside = gamma_target < lo or gamma_target > hi
    if outside and not allow_extrapolation:
        short_by = (gamma_target / hi) if gamma_target > hi else (lo / gamma_target)
        raise ExtrapolationRefused(
            f"{label}: the target wall shear rate {gamma_target:.4g} 1/s lies "
            f"outside the measured range [{lo:.4g}, {hi:.4g}] 1/s "
            f"({short_by:.2f}x short). Extrapolating N1 past the measurement is "
            f"how a die-swell number becomes fiction.\n"
            f"  The right fix is to re-measure the flow curve out to the needle "
            f"wall shear rate. The needle reaches ~165-218 1/s at the nominal "
            f"operating point, and a flow curve that stops at 100 1/s cannot "
            f"speak about it.\n"
            f"  If you must proceed anyway, pass --allow-extrapolation. Every "
            f"number downstream will be stamped EXTRAPOLATED and must be "
            f"reported as such.")

    out = float(10 ** np.interp(np.log10(gamma_target), np.log10(gd), np.log10(v)))
    return out, bool(outside)


def n1_at_wall(fc, gamma_w: float, allow_extrapolation: bool = False) -> dict:
    """
    Extract N1 and tau at the wall shear rate from one force-augmented curve,
    excluding any edge-fracture window first.
    """
    df = fc.data
    gd = df["shear_rate_1_s"].to_numpy(dtype=float)
    tau = df["shear_stress_Pa"].to_numpy(dtype=float)
    n1 = df["N1_Pa"].to_numpy(dtype=float)

    frac = detect_edge_fracture(gd, tau, n1)
    keep = np.ones(gd.size, dtype=bool)
    if frac.detected:
        lo, hi = frac.gamma_dot_range
        # Everything at or beyond the onset is suspect: once the sample has
        # torn, later points do not recover.
        keep &= gd < lo

    if keep.sum() < 2:
        raise ValueError(
            f"{fc.sample_name}: excluding the edge-fracture window leaves "
            f"{keep.sum()} usable point(s). This curve cannot support a "
            f"Tanner estimate; re-measure with a smaller gap or a serrated "
            f"fixture.")

    N1_w, n1_extrap = interpolate_at(gd[keep], n1[keep], gamma_w,
                                     allow_extrapolation, "N1")
    tau_rheo, tau_extrap = interpolate_at(gd[keep], tau[keep], gamma_w,
                                          allow_extrapolation, "tau (rheometer)")
    return {
        "sample": fc.sample_name,
        "instrument": fc.instrument or "(not recorded)",
        "geometry": fc.geometry,
        "baseline_floor_Pa": fc.baseline_floor_Pa,
        "n_points_total": int(gd.size),
        "n_points_used": int(keep.sum()),
        "gamma_measured_range": (float(np.nanmin(gd)), float(np.nanmax(gd))),
        "fracture": frac,
        "gamma_w": gamma_w,
        "N1_wall_Pa": N1_w,
        "tau_wall_rheometer_Pa": tau_rheo,
        "extrapolated": n1_extrap or tau_extrap,
    }


# ----------------------------------------------------------------------
# Pulling the wall shear rate out of a solver run
# ----------------------------------------------------------------------

def gamma_w_from_master_summary(path: Path, ink: str | None = None,
                                needle: str = "21G", vp_mm_s: float = 0.01,
                                model: str = "PowerLaw") -> dict:
    """
    Read the needle wall shear rate per ink out of ``master_summary_v4.csv``.

    Closes the loop: the shear rate at which N1 must be read is the one the
    solver actually predicts for the operating point you print at, not a
    round number someone remembered.
    """
    import pandas as pd
    df = pd.read_csv(path)
    need = {"ink", "needle", "Vp_mm_s", "model", "gamma_w_needle_invs"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(
            f"{path.name} lacks {sorted(missing)}. If this CSV predates "
            f"2026-08-26 it was written before the geometry/beta columns were "
            f"added — re-run run_solver_v4 to regenerate it.")

    sel = df[(df["needle"] == needle) & (df["model"] == model)]
    if sel.empty:
        raise ValueError(f"no rows for needle={needle!r} model={model!r} in {path.name}")
    # nearest Vp, rather than exact equality on a float
    sel = sel.assign(_d=(sel["Vp_mm_s"] - vp_mm_s).abs())
    out = {}
    for name, grp in sel.groupby("ink"):
        row = grp.loc[grp["_d"].idxmin()]
        out[str(name)] = {"gamma_w": float(row["gamma_w_needle_invs"]),
                          "tau_w_solver_Pa": float(row["tau_wall_needle_Pa"]),
                          "Vp_mm_s": float(row["Vp_mm_s"])}
    if ink is not None:
        out = {k: v for k, v in out.items() if k == ink}
    return out


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------

def report_one(res: dict, n_PL: float | None, Rrec_pct: float | None,
               tau_w_solver: float | None) -> dict:
    """Print the closure for one ink and return the numbers for the snippet."""
    print("\n" + "=" * 74)
    print(f"  {res['sample']}")
    print("=" * 74)
    print(f"  instrument           : {res['instrument']}")
    print(f"  geometry             : {res['geometry']} "
          f"(R = {GEOMETRY_RADIUS_M.get(res['geometry'], float('nan'))*1e3:.1f} mm)")
    print(f"  baseline floor tared : {res['baseline_floor_Pa']:.4g} Pa")
    lo, hi = res["gamma_measured_range"]
    print(f"  measured gamma range : {lo:.4g} - {hi:.4g} 1/s "
          f"({res['n_points_used']}/{res['n_points_total']} points used)")
    if res["fracture"].detected:
        print(f"  [!] {res['fracture'].message}")
        print("      -> the fractured window was EXCLUDED before interpolating.")
    print(f"  target wall gamma    : {res['gamma_w']:.4g} 1/s"
          + ("   *** EXTRAPOLATED ***" if res["extrapolated"] else ""))
    print(f"  N1 at wall           : {res['N1_wall_Pa']:.6g} Pa")
    print(f"  tau at wall (rheom.) : {res['tau_wall_rheometer_Pa']:.6g} Pa")

    out = {"N1_wall_Pa": res["N1_wall_Pa"], "extrapolated": res["extrapolated"]}

    # --- the closure, both tau choices side by side ---
    print(f"\n  {'tau source':<22} {'tau_w (Pa)':>12} {'N1/(2tau)':>11} {'beta':>9}")
    print("  " + "-" * 58)
    b_rheo = tanner_beta(res["N1_wall_Pa"], res["tau_wall_rheometer_Pa"])
    print(f"  {'rheometer (measured)':<22} {res['tau_wall_rheometer_Pa']:>12.4g} "
          f"{res['N1_wall_Pa']/(2*res['tau_wall_rheometer_Pa']):>11.4f} {b_rheo:>9.4f}")
    b_solver = None
    if tau_w_solver:
        b_solver = tanner_beta(res["N1_wall_Pa"], tau_w_solver)
        print(f"  {'solver (needle wall)':<22} {tau_w_solver:>12.4g} "
              f"{res['N1_wall_Pa']/(2*tau_w_solver):>11.4f} {b_solver:>9.4f}   <- v5 uses this")
        spread = abs(b_solver - b_rheo) / max(b_rheo, 1e-12) * 100
        print(f"\n  the two tau choices differ by {spread:.1f}% in beta")
        if spread > 20:
            print("  [!] That is a large spread. The solver's tau_w comes from the "
                  "fitted\n      flow curve; if it disagrees this much with the "
                  "measured tau at the same\n      shear rate, check the fit before "
                  "trusting either beta.")
    out["beta_tanner"] = b_solver if b_solver is not None else b_rheo

    # --- what it costs at the bench ---
    if n_PL is not None and Rrec_pct is not None:
        b_h = heuristic_beta(n_PL)
        b_t = out["beta_tanner"]
        print(f"\n  {'closure':<22} {'beta':>9} {'w_line/2Rn':>12} {'k_flow':>9} {'Extr.Mult':>10}")
        print("  " + "-" * 66)
        for label, b in (("v4 heuristic", b_h),
                         ("Nickell floor (N1=0)", NICKELL_FLOOR),
                         ("Tanner (measured N1)", b_t)):
            kf = k_flow(b, Rrec_pct)
            print(f"  {label:<22} {b:>9.4f} {1+b:>12.4f} {kf:>9.4f} {1/kf:>10.4f}")
        em_h, em_t = 1 / k_flow(b_h, Rrec_pct), 1 / k_flow(b_t, Rrec_pct)
        print(f"\n  switching v4 -> Tanner moves the Extrusion Multiplier "
              f"{em_h:.4f} -> {em_t:.4f} ({(em_t/em_h-1)*100:+.1f}%)")
        out["beta_heuristic"] = b_h
        out["EM_heuristic"] = em_h
        out["EM_tanner"] = em_t
    else:
        print("\n  (pass --n-pl and --rrec for this ink to see the slicer consequence)")
    return out


def emit_snippet(results: dict) -> None:
    """Print an inks_local.m-ready block. Pasting beats retyping."""
    print("\n" + "=" * 74)
    print("  PASTE INTO inks_local.m  (one line per ink, next to .Rrec_pct)")
    print("=" * 74)
    for i, (ink, r) in enumerate(sorted(results.items()), start=1):
        stamp = "   % *** EXTRAPOLATED - report as such ***" if r.get("extrapolated") else ""
        print(f"inks({i}).N1_wall_Pa = {r['N1_wall_Pa']:.6g};"
              f"   % Pa, at the needle wall shear rate  [{ink}]{stamp}")
    print("\n  Then run run_solver_v5 instead of run_solver_v4. Any ink without")
    print("  N1_wall_Pa falls back to the v4 heuristic and says so in its output.")


# ----------------------------------------------------------------------

def _self_test() -> int:
    print("SELF-TEST — closure maths only, no data required.\n")
    ok = True

    b0 = tanner_beta(0.0, 500.0)
    print(f"  tanner_beta(N1=0, tau=500)   = {b0:.6f}   (must be exactly the "
          f"Nickell floor {NICKELL_FLOOR})")
    ok &= abs(b0 - NICKELL_FLOOR) < 1e-12

    # monotone in N1: more elasticity, more swell
    bs = [tanner_beta(x, 500.0) for x in (0, 250, 500, 1000, 2000, 4000)]
    print(f"  beta rises with N1           : "
          + " -> ".join(f"{b:.4f}" for b in bs))
    ok &= all(b2 > b1 for b1, b2 in zip(bs, bs[1:]))

    # the ratio is what matters, not the absolute scale
    a = tanner_beta(2000.0, 500.0)
    b = tanner_beta(4000.0, 1000.0)
    print(f"  scale-invariant in N1/tau    : {a:.6f} vs {b:.6f}")
    ok &= abs(a - b) < 1e-12

    try:
        tanner_beta(-1.0, 500.0)
        print("  [FAIL] negative N1 accepted")
        ok = False
    except ValueError:
        print("  negative N1                  : refused (un-tared transducer)")

    gd = np.logspace(-1, 2, 40)
    v = 100 * gd ** 0.5
    got, ex = interpolate_at(gd, v, 10.0)
    print(f"  log-log interp on a power law: {got:.4f} (exact {100*10**0.5:.4f}), "
          f"extrapolated={ex}")
    ok &= abs(got - 100 * 10 ** 0.5) / (100 * 10 ** 0.5) < 1e-9

    try:
        interpolate_at(gd, v, 200.0)
        print("  [FAIL] extrapolation allowed by default")
        ok = False
    except ExtrapolationRefused as e:
        print(f"  target above measured range  : refused "
              f"({str(e).split(chr(10))[0][-40:].strip()})")

    got2, ex2 = interpolate_at(gd, v, 200.0, allow_extrapolation=True)
    print(f"  with --allow-extrapolation   : {got2:.4f}, flagged={ex2}")
    ok &= ex2

    print("\n  SELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--force-folder", type=Path, default=None,
                    help="folder of force-augmented exports (default: the "
                         "'force' entry of data_config_local.py)")
    ap.add_argument("--file", type=Path, action="append", default=None,
                    help="one specific force export; repeatable")
    ap.add_argument("--geometry", default="CP50", help="measuring system (default CP50)")
    ap.add_argument("--force-unit", default="N", choices=["N", "mN", "gf"])
    ap.add_argument("--gamma-w", type=float, default=None,
                    help="needle wall shear rate, 1/s")
    ap.add_argument("--master-summary", type=Path, default=None,
                    help="read gamma_w (and solver tau_w) per ink from a "
                         "master_summary_v4.csv")
    ap.add_argument("--needle", default="21G")
    ap.add_argument("--vp", type=float, default=0.01, help="piston velocity, mm/s")
    ap.add_argument("--n-pl", type=float, default=None,
                    help="Power-Law n, to show the v4 heuristic alongside")
    ap.add_argument("--rrec", type=float, default=None,
                    help="Rrec_pct, to show k_flow and the Extrusion Multiplier")
    ap.add_argument("--allow-extrapolation", action="store_true",
                    help="proceed past the measured shear-rate range (flagged)")
    args = ap.parse_args(argv)

    if args.self_test:
        return _self_test()

    if args.gamma_w is None and args.master_summary is None:
        ap.error("supply --gamma-w or --master-summary; the wall shear rate is "
                 "not guessable and getting it wrong changes beta")

    folder = args.force_folder
    if folder is None and not args.file and _load_data_config is not None:
        folder = Path(_load_data_config().folder("force"))

    files = list(args.file or [])
    if folder is not None:
        if not folder.is_dir():
            print(f"[!] {folder} not found. Force-augmented exports are produced "
                  f"by re-exporting the flow curve with the axial-force channel "
                  f"switched on; they are not part of the standard export.")
            return 2
        files += sorted(p for p in folder.iterdir()
                        if p.suffix.lower() in (".txt", ".csv"))
    if not files:
        print("[!] no force exports found.")
        return 2

    solver = {}
    if args.master_summary:
        try:
            solver = gamma_w_from_master_summary(
                args.master_summary, None, args.needle, args.vp)
            print(f"[gamma_w] from {args.master_summary} "
                  f"(needle {args.needle}, Vp ~ {args.vp} mm/s):")
            for k, v in solver.items():
                print(f"    {k:<16} gamma_w = {v['gamma_w']:.4g} 1/s   "
                      f"tau_w = {v['tau_w_solver_Pa']:.4g} Pa   "
                      f"(Vp {v['Vp_mm_s']} mm/s)")
        except Exception as e:
            print(f"[!] could not read {args.master_summary}: {e}")
            return 2

    results, nfail = {}, 0
    for f in files:
        try:
            fc = read_flow_curve_with_force(f, geometry=args.geometry,
                                            force_unit=args.force_unit)
        except Exception as e:
            print(f"\n[!] {f.name}: {e}")
            nfail += 1
            continue

        match = next((v for k, v in solver.items() if k in fc.sample_name), None)
        gamma_w = args.gamma_w if args.gamma_w is not None else (
            match["gamma_w"] if match else None)
        if gamma_w is None:
            print(f"\n[!] {fc.sample_name}: no gamma_w — not in the master "
                  f"summary and no --gamma-w given. Skipped.")
            nfail += 1
            continue

        try:
            res = n1_at_wall(fc, gamma_w, args.allow_extrapolation)
        except (ExtrapolationRefused, ValueError) as e:
            print(f"\n[!] {fc.sample_name}:\n    {e}")
            nfail += 1
            continue

        results[fc.sample_name] = report_one(
            res, args.n_pl, args.rrec,
            match["tau_w_solver_Pa"] if match else None)

    if results:
        emit_snippet(results)
    if nfail:
        print(f"\n[!] {nfail} curve(s) could not be reduced. Nothing was guessed.")
    return 0 if results and not nfail else (0 if results else 1)


if __name__ == "__main__":
    sys.exit(main())
