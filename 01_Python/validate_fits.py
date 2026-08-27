#!/usr/bin/env python3
"""
validate_fits.py
================

Python counterpart of ``validate_v4.m``: cross-check the legacy and current
flow-curve fitters against each other on the same data, and say WHERE they
disagree.

WHY
---
``Fit_Muitos_Modelos_v4.py`` (audited, TA-style ``.xls`` via xlrd) and
``Fit_Muitos_Modelos_v5.py`` (active, Anton Paar Rheocompass CSV) are both
registered ``active`` and both get run. They can disagree at the last
significant figure — that is how the still-unreconciled eta0 discrepancy of
2240.9 (one version) vs 2259.2 (the other) on C20 was found, by hand, late.
The MATLAB solvers have had an automated version cross-check since the v4
rewrite; the Python fitters never got one.

WHAT IT SEPARATES
-----------------
A raw "v4 says X, v5 says Y" is not actionable, because the two versions
differ in TWO independent ways: how the file is parsed, and how the fit is
run. This tool runs three combinations so the cause is isolated:

    A. v4 reader + v4 fitter      (the audited number)
    B. v5 reader + v5 fitter      (the current number)
    C. v4 reader + v5 fitter      (the discriminator)

    A vs C differs  ->  the FITTER changed.
    C vs B differs  ->  the INGESTION changed (parsing, unit conversion,
                        point filtering).

Reporting only A vs B, as a naive check would, leaves you unable to say
which, and that is precisely the question an editor asks.

USAGE
-----
    python3 Export/01_Python/validate_fits.py \\
        --xls-folder  "./Reologia/archive-xls" \\
        --csv-folder  "./Reologia/Viscosity" \\
        --tol 0.005

Samples are matched on the short name. Anything present in only one folder
is listed and skipped, never silently dropped.

Exit code 0 if every shared parameter agrees within tolerance, 1 otherwise.

Author: T.M.C. Rodrigues - 2026-08-26
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))


def _load_module(name: str, filename: str):
    """Import a fitter module for its model tables and fit_log, WITHOUT
    letting it run its own analysis.

    ``Fit_Muitos_Modelos_v4.py`` is a flat script: it has no ``main()`` and
    no ``if __name__ == "__main__"`` guard, so importing it executes the
    entire audited analysis and writes report files into whatever directory
    happens to be current. It is registered ``active (legacy)`` and its
    output reproduces numbers in a submitted manuscript, so it must not be
    edited to add a guard. The import is therefore sandboxed instead: run it
    with the working directory pointed at a throwaway folder and its stdout
    swallowed, so any files it writes are discarded with that folder.

    If v4 ever does get a proper main() guard, this wrapper stays correct
    and simply stops mattering.
    """
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod

    cwd = os.getcwd()
    with tempfile.TemporaryDirectory(prefix="validate_fits_sandbox_") as sandbox:
        os.chdir(sandbox)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                spec.loader.exec_module(mod)
        except SystemExit:
            pass
        finally:
            os.chdir(cwd)
    return mod


def _fit_all(mod, gd, y_stress, y_visc) -> dict:
    """Run one module's whole model battery and return {param: value}."""
    out: dict[str, float] = {}
    for label, table, ydata in (("stress", mod.MODEL_STRESS, y_stress),
                                ("visc", mod.MODEL_VISCOSITY, y_visc)):
        for name, (fn, p0, bnds, pnames) in table.items():
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    params, _ = mod.fit_log(fn, gd, ydata, p0, bnds)
            except Exception:
                continue
            for pn, pv in zip(pnames, params):
                out[f"{name}.{pn}"] = float(pv)
    return out


def _positive_mask(gd, tau, eta):
    m = (gd > 0) & (tau > 0) & (eta > 0)
    return m & np.isfinite(gd) & np.isfinite(tau) & np.isfinite(eta)


def compare(a: dict, b: dict, tol: float) -> tuple[list, int]:
    rows, nbad = [], 0
    for key in sorted(set(a) | set(b)):
        va, vb = a.get(key), b.get(key)
        if va is None or vb is None:
            rows.append((key, va, vb, float("nan"), "MISSING"))
            nbad += 1
            continue
        scale = max(abs(va), abs(vb), 1e-30)
        rel = abs(va - vb) / scale
        status = "ok" if rel <= tol else "DIFFERS"
        if status != "ok":
            nbad += 1
        rows.append((key, va, vb, rel, status))
    return rows, nbad


def _print_block(title: str, rows: list, only_bad: bool = False) -> None:
    print(f"\n--- {title} ---")
    print(f"  {'parameter':<26} {'left':>14} {'right':>14} {'rel diff':>11}  status")
    shown = 0
    for key, va, vb, rel, status in rows:
        if only_bad and status == "ok":
            continue
        sa = "  --  " if va is None else f"{va:14.6g}"
        sb = "  --  " if vb is None else f"{vb:14.6g}"
        sr = "   --  " if not np.isfinite(rel) else f"{rel*100:10.4f}%"
        print(f"  {key:<26} {sa} {sb} {sr}  {status}")
        shown += 1
    if shown == 0:
        print("  (all parameters agree within tolerance)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xls-folder", required=True,
                    help="folder of legacy TA-style .xls flow curves (the v4 path)")
    ap.add_argument("--csv-folder", required=True,
                    help="folder of Anton Paar Rheocompass CSVs (the v5 path)")
    ap.add_argument("--tol", type=float, default=5e-3,
                    help="relative tolerance per parameter (default 0.5%%)")
    ap.add_argument("--sample", default=None,
                    help="restrict to one short sample name")
    args = ap.parse_args(argv)

    v4 = _load_module("fit_v4", "Fit_Muitos_Modelos_v4.py")
    v5 = _load_module("fit_v5", "Fit_Muitos_Modelos_v5.py")
    from antpar_io import read_flow_curves_in_folder

    xls_dir, csv_dir = Path(args.xls_folder), Path(args.csv_folder)
    if not xls_dir.is_dir():
        print(f"[!] {xls_dir} not found. The pre-migration .xls archive is not "
              f"present in every checkout of this pipeline; without it the "
              f"legacy half of this comparison cannot run.")
        return 2
    if not csv_dir.is_dir():
        print(f"[!] {csv_dir} not found.")
        return 2

    csv_curves = read_flow_curves_in_folder(str(csv_dir))
    xls_files = sorted(list(xls_dir.glob("*.xls")))
    xls_by_short = {p.stem.split(" - ")[0].strip(): p for p in xls_files}

    shared = sorted(set(csv_curves) & set(xls_by_short))
    only_csv = sorted(set(csv_curves) - set(xls_by_short))
    only_xls = sorted(set(xls_by_short) - set(csv_curves))
    if args.sample:
        shared = [s for s in shared if s == args.sample]

    print("=" * 78)
    print("  validate_fits — Fit_Muitos_Modelos v4 (audited .xls) vs v5 (active CSV)")
    print("=" * 78)
    print(f"  shared samples : {shared or '(none)'}")
    if only_csv:
        print(f"  CSV only       : {only_csv}   (skipped, not silently dropped)")
    if only_xls:
        print(f"  XLS only       : {only_xls}   (skipped, not silently dropped)")
    if not shared:
        print("\n[!] No sample exists in both formats — nothing to cross-check.")
        return 2

    total_bad = 0
    for short in shared:
        print(f"\n{'='*78}\nSAMPLE: {short}\n{'='*78}")

        # --- v5 path: Rheocompass CSV ---
        df5 = csv_curves[short].data
        gd5 = df5["shear_rate_1_s"].to_numpy()
        tau5 = df5["shear_stress_Pa"].to_numpy()
        eta5 = df5["viscosity_Pa_s"].to_numpy()
        m5 = _positive_mask(gd5, tau5, eta5)
        B = _fit_all(v5, gd5[m5], tau5[m5], eta5[m5])

        # --- v4 path: legacy .xls, parsed by the legacy reader ---
        try:
            gd4, tau4, eta4 = _read_legacy_xls(v4, xls_by_short[short])
        except Exception as e:
            print(f"  [!] legacy .xls read failed: {e}")
            total_bad += 1
            continue
        m4 = _positive_mask(gd4, tau4, eta4)
        A = _fit_all(v4, gd4[m4], tau4[m4], eta4[m4])
        C = _fit_all(v5, gd4[m4], tau4[m4], eta4[m4])   # discriminator

        print(f"  points: v4 reader {m4.sum()} | v5 reader {m5.sum()}")

        rows_fit, bad_fit = compare(A, C, args.tol)
        _print_block("FITTER effect   (v4 reader+fitter  vs  v4 reader+v5 fitter)",
                     rows_fit, only_bad=True)
        rows_ing, bad_ing = compare(C, B, args.tol)
        _print_block("INGESTION effect (v4 reader+v5 fitter  vs  v5 reader+fitter)",
                     rows_ing, only_bad=True)
        rows_tot, bad_tot = compare(A, B, args.tol)
        _print_block("TOTAL           (audited v4  vs  active v5)", rows_tot, only_bad=True)

        print(f"\n  verdict: {bad_fit} fitter-attributable, {bad_ing} "
              f"ingestion-attributable, {bad_tot} total parameter(s) outside "
              f"{args.tol*100:.2f}%")
        total_bad += bad_tot

    print("\n" + "=" * 78)
    if total_bad == 0:
        print("  validate_fits: v4 and v5 agree on every shared parameter.")
    else:
        print(f"  validate_fits: {total_bad} parameter(s) DISAGREE. Read the "
              f"FITTER/INGESTION split above before touching either script.")
    print("=" * 78)
    return 0 if total_bad == 0 else 1


def _read_legacy_xls(v4mod, path: Path):
    """Pull gamma_dot / tau / eta out of a legacy TA-style .xls.

    Uses the legacy reader that already exists rather than reimplementing
    the sheet layout (B18): ``extract_SAOS_values.read_antpar_xls`` knows
    that sheet 0 is "Details", the data lives on "Ramp - 1", headers are on
    row 1, units on row 2, and data starts at row 3. Every one-off audit
    script that re-derived that layout from scratch was wasted work.
    """
    from extract_SAOS_values import read_antpar_xls, find_col
    df = read_antpar_xls(str(path))
    gd = df[find_col(df, ["shear rate", "taxa"])].to_numpy(dtype=float)
    tau = df[find_col(df, ["shear stress", "stress", "tensao"])].to_numpy(dtype=float)
    try:
        eta = df[find_col(df, ["viscosity", "viscosidade"])].to_numpy(dtype=float)
    except Exception:
        eta = tau / np.where(gd > 0, gd, np.nan)
    return gd, tau, eta


if __name__ == "__main__":
    sys.exit(main())
