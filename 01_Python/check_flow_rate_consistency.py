#!/usr/bin/env python3
"""
check_flow_rate_consistency.py
==============================

Cross-check the three independent estimates of volumetric flow rate that a
DIW print run produces, and report their disagreement.

WHY THIS EXISTS
---------------
Three values of Q for the same print surfaced during a manuscript audit and
were never reconciled:

    solver Q               0.71  mm^3/s   (from the extrusion model)
    G-code Q               ~0.2  mm^3/s   (two independent readings agreed)
    mass-conservation Q    0.46-0.51 mm^3/s   (from the measured deposit)

They bear directly on whether the statement "the computed multiplier was
loaded into the slicer for every print" is defensible. Nothing in the
pipeline computed all three, so nobody noticed until a reviewer's question
forced the comparison by hand. This turns that into a standard
pre-submission step.

THE THREE ROUTES, AND WHAT A MISMATCH MEANS
-------------------------------------------
1. SOLVER Q = pi*Rs^2 * Vp. What the piston displaces. Depends only on
   syringe bore and piston speed, so it is the least likely to be wrong and
   the most likely to be irrelevant — it is what the machine WOULD deliver
   with no slip and no compliance.

2. G-CODE Q = (extruded volume)/(time), read back out of the file that was
   actually sent to the printer. If this disagrees with the solver, the
   computed multiplier was not what got loaded, or the slicer's filament
   geometry does not match the syringe.

3. MASS-CONSERVATION Q = w * h * v_print, from the road you can measure on
   the print. This is ground truth for what was actually deposited. If it
   sits between the other two, you are looking at deposition efficiency
   (that is what k_flow is for). If it exceeds the solver Q, something is
   wrong with the geometry, not the efficiency.

USAGE
-----
    # all three
    python3 Export/01_Python/check_flow_rate_consistency.py \\
        --solver-q 0.71 \\
        --gcode print_C20_21G.gcode --filament-diameter 1.75 \\
        --width 0.63 --height 0.36 --speed 7.33

    # or let the solver Q be computed from the hardware
    ... --syringe-id 14.3 --piston-speed 0.01

    # reproduce the documented C20 discrepancy without any files present
    python3 Export/01_Python/check_flow_rate_consistency.py --self-test

Exit code 0 if the spread is within tolerance, 1 otherwise.

Author: T.M.C. Rodrigues - 2026-08-26
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

_G1 = re.compile(r"^G1\b")
_AXIS = re.compile(r"([XYZEF])\s*(-?\d*\.?\d+)", re.IGNORECASE)


def solver_flow_rate(syringe_id_mm: float, piston_speed_mm_s: float) -> float:
    """Q = pi * Rs^2 * Vp, in mm^3/s. Rs is the syringe INNER RADIUS."""
    Rs = syringe_id_mm / 2.0
    return math.pi * Rs ** 2 * piston_speed_mm_s


def deposit_flow_rate(width_mm: float, height_mm: float, speed_mm_s: float) -> float:
    """Q = w * h * v_print, in mm^3/s — mass conservation of the road."""
    return width_mm * height_mm * speed_mm_s


def gcode_flow_rate(path: Path, filament_diameter_mm: float = 1.75,
                    volumetric_e: bool = False) -> dict:
    """
    Read back the mean volumetric flow rate from a G-code file.

    Only G1 moves with a POSITIVE E delta count as extrusion; travels and
    retractions are excluded from both the volume and the time, because
    including travel time is the single most common way to get a G-code Q
    that is quietly too low.

    Parameters
    ----------
    filament_diameter_mm : used to convert an E axis expressed in mm of
        filament into volume. IGNORED when `volumetric_e` is set. If the
        slicer's filament diameter does not match the syringe this number
        is where the mismatch enters — that is a finding, not a nuisance.
    volumetric_e : True when the firmware/slicer emits E directly in mm^3.

    Returns dict with Q_mm3_s, total_volume_mm3, total_time_s, n_segments.
    """
    area = 1.0 if volumetric_e else math.pi * (filament_diameter_mm / 2.0) ** 2
    x = y = z = e = 0.0
    feed_mm_min = 0.0
    absolute_e = True
    total_v = total_t = 0.0
    nseg = 0

    for raw in Path(path).read_text(errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("M82"):
            absolute_e = True
            continue
        if line.startswith("M83"):
            absolute_e = False
            continue
        if line.startswith("G92"):
            for ax, val in _AXIS.findall(line):
                if ax.upper() == "E":
                    e = float(val)
            continue
        if not _G1.match(line):
            continue

        nx, ny, nz, ne = x, y, z, e
        for ax, val in _AXIS.findall(line):
            ax = ax.upper()
            v = float(val)
            if ax == "X":
                nx = v
            elif ax == "Y":
                ny = v
            elif ax == "Z":
                nz = v
            elif ax == "E":
                ne = v if absolute_e else e + v
            elif ax == "F":
                feed_mm_min = v

        de = ne - e
        dist = math.dist((x, y, z), (nx, ny, nz))
        if de > 0 and feed_mm_min > 0:
            # Time from the commanded feedrate. For a pure-extrusion move
            # with no XYZ travel, fall back to the E axis length.
            travel = dist if dist > 0 else abs(de)
            total_t += travel / (feed_mm_min / 60.0)
            total_v += de * area
            nseg += 1

        x, y, z, e = nx, ny, nz, ne

    if total_t <= 0:
        raise ValueError(
            f"{Path(path).name}: no extruding G1 moves with a feedrate were "
            f"found. Check that the file is G-code and that F is set — a "
            f"G-code Q of zero is a parsing failure, not a slow print.")

    return {"Q_mm3_s": total_v / total_t, "total_volume_mm3": total_v,
            "total_time_s": total_t, "n_segments": nseg}


def report(values: dict, tol_pct: float = 15.0) -> int:
    """Print the three routes side by side and return an exit code."""
    present = {k: v for k, v in values.items() if v is not None}
    print("=" * 72)
    print("  FLOW-RATE CONSISTENCY CHECK")
    print("=" * 72)
    if not present:
        print("  Nothing to compare — supply at least two routes.")
        return 2

    width = max(len(k) for k in present)
    print(f"\n  {'route':<{width}}   {'Q (mm^3/s)':>12}")
    print("  " + "-" * (width + 16))
    for k, v in present.items():
        print(f"  {k:<{width}}   {v:>12.4f}")

    if len(present) < 2:
        print("\n  Only one route available — no cross-check performed.")
        return 2

    print(f"\n  {'pair':<{2*width+5}}   {'disagreement':>13}")
    print("  " + "-" * (2 * width + 21))
    keys = list(present)
    worst = 0.0
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = present[keys[i]], present[keys[j]]
            ref = max(abs(a), abs(b))
            d = abs(a - b) / ref * 100 if ref else 0.0
            worst = max(worst, d)
            mark = "  <-- exceeds tolerance" if d > tol_pct else ""
            print(f"  {keys[i]} vs {keys[j]:<{width}}   {d:>12.1f}%{mark}")

    print(f"\n  worst pairwise disagreement: {worst:.1f}%  (tolerance {tol_pct:.0f}%)")
    if worst > tol_pct:
        print(
            "\n  READ THIS BEFORE SUBMITTING. A spread this size means the three\n"
            "  routes are not describing the same print. Likely causes, in the\n"
            "  order worth checking:\n"
            "    - the computed extrusion multiplier (1/k_flow) was not the one\n"
            "      loaded into the slicer for this run;\n"
            "    - the slicer's filament diameter does not match the syringe, so\n"
            "      the E axis means something different from what it is read as;\n"
            "    - travel moves are being counted in the G-code time;\n"
            "    - the deposit width was measured on a road that had spread after\n"
            "      landing, which inflates the mass-conservation route.\n"
            "  Do not average them.")
        return 1
    print("\n  The three routes agree within tolerance.")
    return 0


def _self_test() -> int:
    """Reproduce the documented C20 discrepancy with no files present."""
    print("SELF-TEST: the unreconciled C20 case, from the manuscript audit.\n"
          "These are the three numbers that were found by hand and never\n"
          "resolved. The tool exists so this comparison happens before a\n"
          "reviewer asks, not after.\n")
    rc = report({
        "solver": 0.71,
        "gcode": 0.20,
        "deposit": 0.485,     # midpoint of the reported 0.46-0.51
    })
    print("\n  Cross-check of the helper functions themselves:")
    q = solver_flow_rate(14.3, 0.01)
    print(f"    solver_flow_rate(ID=14.3 mm, Vp=0.01 mm/s) = {q:.4f} mm^3/s")
    d = deposit_flow_rate(0.63, 0.36, 7.33)
    print(f"    deposit_flow_rate(w=0.63, h=0.36, v=7.33)  = {d:.4f} mm^3/s")
    assert abs(q - math.pi * (14.3 / 2) ** 2 * 0.01) < 1e-12
    assert abs(d - 0.63 * 0.36 * 7.33) < 1e-12
    print("    both analytic — OK")
    return 0 if rc == 1 else 1   # the self-test EXPECTS a flagged discrepancy


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true",
                    help="reproduce the documented C20 discrepancy, no files needed")
    ap.add_argument("--solver-q", type=float, default=None,
                    help="solver volumetric flow rate, mm^3/s")
    ap.add_argument("--syringe-id", type=float, default=None,
                    help="syringe inner diameter, mm (with --piston-speed)")
    ap.add_argument("--piston-speed", type=float, default=None,
                    help="piston velocity, mm/s (with --syringe-id)")
    ap.add_argument("--gcode", type=Path, default=None, help="G-code file to read back")
    ap.add_argument("--filament-diameter", type=float, default=1.75,
                    help="slicer filament diameter, mm (default 1.75)")
    ap.add_argument("--volumetric-e", action="store_true",
                    help="E axis is already in mm^3")
    ap.add_argument("--width", type=float, default=None, help="measured road width, mm")
    ap.add_argument("--height", type=float, default=None, help="layer height, mm")
    ap.add_argument("--speed", type=float, default=None, help="print speed, mm/s")
    ap.add_argument("--tol", type=float, default=15.0,
                    help="pairwise tolerance in percent (default 15)")
    args = ap.parse_args(argv)

    if args.self_test:
        return _self_test()

    values: dict = {}

    if args.solver_q is not None:
        values["solver"] = args.solver_q
    elif args.syringe_id and args.piston_speed:
        values["solver"] = solver_flow_rate(args.syringe_id, args.piston_speed)

    if args.gcode:
        try:
            g = gcode_flow_rate(args.gcode, args.filament_diameter, args.volumetric_e)
            values["gcode"] = g["Q_mm3_s"]
            print(f"  [gcode] {g['n_segments']} extruding segments, "
                  f"{g['total_volume_mm3']:.2f} mm^3 over {g['total_time_s']:.1f} s")
        except Exception as e:
            print(f"  [!] G-code read failed: {e}")

    if args.width and args.height and args.speed:
        values["deposit"] = deposit_flow_rate(args.width, args.height, args.speed)

    return report(values, tol_pct=args.tol)


if __name__ == "__main__":
    sys.exit(main())
