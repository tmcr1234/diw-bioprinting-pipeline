"""
extract_SAOS_values.py
======================

╔══════════════════════════════════════════════════════════════════════════╗
║  ⚠ ANTON PAAR MIGRATION PENDING (frequency + amplitude sweeps)           ║
║  Most format-tolerant of the legacy extractors (header found by          ║
║  keyword scan) but still assumes legacy .xls (xlrd-readable). For        ║
║  Anton Paar Rheocompass CSV data, port to v2 + antpar_io is pending —    ║
║  drop a sample CSV to advance the migration.                             ║
║  See Export/docs/01_rheology_data_formats.md                             ║
╚══════════════════════════════════════════════════════════════════════════╝

Extract exact SAOS (Small-Amplitude Oscillatory Shear) values from
Anton Paar MCR 501 .xls workbooks for Sections 3.1 (amplitude sweep)
and 3.2 (frequency sweep) of the manuscript Results_Discussion_v2.tex.

For each group (C15, C15 Gira 5.5 / NE, Bozano):
  - Amplitude sweep -> G'_LVR, G''_LVR, tan delta_LVR, gamma_LVR (10% drop)
  - Frequency sweep -> G' and tan delta at omega ~ 1, 10, 100 rad/s
"""

import os
import xlrd
import numpy as np


PROJECT_ROOT = "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE"


def dump_sheets(filepath, max_rows=20, max_cols=12):
    """Diagnostic: print first rows of every sheet."""
    wb = xlrd.open_workbook(filepath)
    print(f"--- {os.path.basename(filepath)} ({wb.nsheets} sheets) ---")
    for si in range(wb.nsheets):
        sh = wb.sheet_by_index(si)
        print(f"  Sheet[{si}] '{sh.name}' rows={sh.nrows} cols={sh.ncols}")
        for ri in range(min(sh.nrows, max_rows)):
            row = [str(sh.cell_value(ri, ci)).strip()
                   for ci in range(min(sh.ncols, max_cols))]
            print(f"    R{ri}: {row}")


def read_antpar_xls(filepath):
    """
    Read Anton Paar MCR .xls; returns (headers, units, ndarray).
    Locates a header row containing typical SAOS keywords, then collects
    numeric rows below it.
    """
    wb = xlrd.open_workbook(filepath)
    keywords = ["G'", "G''", "Omega", "omega", "Frequency", "Freq.",
                "Deformation", "Deforma", "Strain", "Tan", "tan(",
                "Storage", "Loss"]
    for si in range(wb.nsheets):
        sh = wb.sheet_by_index(si)
        for ri in range(sh.nrows):
            row = [str(sh.cell_value(ri, ci)).strip()
                   for ci in range(sh.ncols)]
            hits = sum(1 for r in row for kw in keywords if kw in r)
            if hits >= 2:
                headers = row
                # Units row often immediately after header row
                units = []
                start = ri + 1
                if ri + 1 < sh.nrows:
                    u_row = [str(sh.cell_value(ri + 1, ci)).strip()
                             for ci in range(sh.ncols)]
                    # Heuristic: units row has [Pa], [%], [rad/s], [s], [-]
                    if any("[" in v or v in ("Pa", "%", "rad/s", "s", "-",
                                             "1/s", "Hz") for v in u_row):
                        units = u_row
                        start = ri + 2
                data_rows = []
                for rr in range(start, sh.nrows):
                    vals = []
                    for ci in range(sh.ncols):
                        v = sh.cell_value(rr, ci)
                        try:
                            vals.append(float(v))
                        except (ValueError, TypeError):
                            # Empty/strings -> NaN
                            vals.append(np.nan)
                    if not all(np.isnan(v) for v in vals):
                        data_rows.append(vals)
                if data_rows:
                    arr = np.array(data_rows)
                    return headers, units, arr, si
    return [], [], np.array([]), -1


def find_col(headers, *keywords):
    """Return index of first header containing any of the keywords (case-sensitive substring)."""
    for i, h in enumerate(headers):
        for kw in keywords:
            if kw in h:
                return i
    return None


def extract_amplitude_sweep(filepath):
    """Returns dict with G'_LVR, G''_LVR, tan_d_LVR, gamma_LVR."""
    headers, units, arr, sidx = read_antpar_xls(filepath)
    if arr.size == 0:
        return None, headers, units
    # Identify columns (TA Instruments AresG2 export schema)
    i_strain = find_col(headers, "Oscillation strain", "Deformation",
                        "Deforma", "Strain")
    i_Gp = find_col(headers, "Storage modulus", "G'")
    i_Gpp = find_col(headers, "Loss modulus", "G''", "G\"")
    i_tan = find_col(headers, "Tan(delta)", "tan", "Tan")
    info = {"headers": headers, "units": units, "sheet_index": sidx,
            "n_rows": len(arr), "i_strain": i_strain,
            "i_Gp": i_Gp, "i_Gpp": i_Gpp, "i_tan": i_tan}
    if None in (i_strain, i_Gp, i_Gpp):
        info["error"] = "Required columns not found"
        return info, headers, units
    strain = arr[:, i_strain]
    Gp = arr[:, i_Gp]
    Gpp = arr[:, i_Gpp]
    tan_d = arr[:, i_tan] if i_tan is not None else Gpp / Gp
    # Sort ascending by strain
    order = np.argsort(strain)
    strain = strain[order]; Gp = Gp[order]; Gpp = Gpp[order]; tan_d = tan_d[order]
    # Filter: physically valid points only (G'>0, finite)
    valid = (~np.isnan(Gp) & ~np.isnan(strain) & ~np.isnan(Gpp)
             & (Gp > 0) & (Gpp > 0) & np.isfinite(tan_d))
    s = strain[valid]; gp = Gp[valid]; gpp = Gpp[valid]; td = tan_d[valid]

    # Robust LVR detection:
    # 1) Find the longest run of points where the local CV (std/mean) of G'
    #    in a sliding window of 5 is below 5 percent — this is the plateau.
    # 2) Use those indices for the LVR average.
    win = 5
    cv = np.full(len(gp), np.inf)
    for k in range(len(gp) - win + 1):
        seg = gp[k:k + win]
        m = np.mean(seg)
        if m > 0:
            cv[k + win // 2] = np.std(seg) / m
    stable = cv < 0.05  # 5% local coefficient of variation
    # Build the longest contiguous run of stable points
    best_run = (0, 0)
    cur_start = None
    for k in range(len(stable)):
        if stable[k]:
            if cur_start is None:
                cur_start = k
            if k - cur_start + 1 > best_run[1] - best_run[0] + 1:
                best_run = (cur_start, k)
        else:
            cur_start = None
    if best_run[1] - best_run[0] >= 3:
        sel = np.arange(best_run[0], best_run[1] + 1)
    else:
        # Fallback to median ± 15% over the lower half of the data
        n_low = max(5, len(gp) // 2)
        median_low = np.median(gp[:n_low])
        in_plateau = np.abs(gp[:n_low] - median_low) <= 0.15 * median_low
        plateau_idx = np.where(in_plateau)[0]
        sel = plateau_idx if len(plateau_idx) >= 3 else np.arange(min(5, len(gp)))
    Gp_plateau = float(np.mean(gp[sel]))
    Gpp_plateau = float(np.mean(gpp[sel]))
    tan_plateau = float(np.mean(td[sel]))

    # gamma_LVR: smallest strain (above plateau region) at which G' falls
    # below 90% of the plateau value.
    threshold = 0.9 * Gp_plateau
    gamma_lvr = None
    last_in = sel[-1] if len(sel) else 0
    for k in range(last_in + 1, len(gp)):
        if gp[k] < threshold:
            x0, x1 = np.log10(s[k - 1]), np.log10(s[k])
            y0, y1 = gp[k - 1], gp[k]
            if y1 != y0:
                frac = (threshold - y0) / (y1 - y0)
                gamma_lvr = float(10 ** (x0 + frac * (x1 - x0)))
            else:
                gamma_lvr = float(s[k])
            break
    if gamma_lvr is None:
        gamma_lvr = float(s[-1])
    # Replace diagnostic fields
    info["plateau_indices"] = sel.tolist() if hasattr(sel, "tolist") else list(sel)
    info["Gp_first5"] = gp[:5].tolist()
    info["strain_first5"] = s[:5].tolist()
    info["plateau_strain_range"] = (float(s[sel[0]]), float(s[sel[-1]])) if len(sel) else None
    n_avg = len(sel)
    info.update({
        "Gp_LVR": Gp_plateau, "Gpp_LVR": Gpp_plateau,
        "tan_LVR": tan_plateau, "gamma_LVR_pct": gamma_lvr,
        "strain_range": (float(s[0]), float(s[-1])),
        "Gp_first5": gp[:n_avg].tolist(),
        "strain_first5": s[:n_avg].tolist(),
    })
    return info, headers, units


def extract_frequency_sweep(filepath, omega_targets=(1.0, 10.0, 100.0)):
    """Return dict with G', G'', tan delta at omega ~ targets (rad/s)."""
    headers, units, arr, sidx = read_antpar_xls(filepath)
    if arr.size == 0:
        return None, headers, units
    i_w = find_col(headers, "Angular frequency", "Angular Frequency",
                   "Omega", "omega")
    i_freq_hz = None
    if i_w is None:
        i_freq_hz = find_col(headers, "Frequency", "Freq.")
    i_Gp = find_col(headers, "Storage modulus", "G'")
    i_Gpp = find_col(headers, "Loss modulus", "G''", "G\"")
    i_tan = find_col(headers, "Tan(delta)", "tan", "Tan")
    info = {"headers": headers, "units": units, "sheet_index": sidx,
            "n_rows": len(arr), "i_w": i_w, "i_freq_hz": i_freq_hz,
            "i_Gp": i_Gp, "i_Gpp": i_Gpp, "i_tan": i_tan}
    if i_w is None and i_freq_hz is None:
        info["error"] = "No frequency column found"
        return info, headers, units
    if i_Gp is None or i_Gpp is None:
        info["error"] = "G' or G'' column missing"
        return info, headers, units
    if i_w is not None:
        w = arr[:, i_w]
        unit_label = (units[i_w] if i_w < len(units) else "rad/s")
    else:
        f_hz = arr[:, i_freq_hz]
        w = 2 * np.pi * f_hz
        unit_label = "Hz->rad/s"
    Gp = arr[:, i_Gp]; Gpp = arr[:, i_Gpp]
    tan_d = arr[:, i_tan] if i_tan is not None else Gpp / Gp
    valid = ~np.isnan(w) & ~np.isnan(Gp) & (w > 0) & (Gp > 0)
    w = w[valid]; Gp = Gp[valid]; Gpp = Gpp[valid]; tan_d = tan_d[valid]
    order = np.argsort(w)
    w = w[order]; Gp = Gp[order]; Gpp = Gpp[order]; tan_d = tan_d[order]
    info["omega_unit"] = unit_label
    info["omega_range"] = (float(w[0]), float(w[-1]))
    points = {}
    for target in omega_targets:
        # Pick the data point with omega closest (in log) to target
        idx = int(np.argmin(np.abs(np.log10(w) - np.log10(target))))
        points[target] = {
            "omega_actual": float(w[idx]),
            "Gp": float(Gp[idx]),
            "Gpp": float(Gpp[idx]),
            "tan": float(tan_d[idx]),
        }
    info["points"] = points
    return info, headers, units


def main():
    os.chdir(PROJECT_ROOT)

    freq_files = {
        "C15":    "Reologia/Frequency/C15 - Freq CP50 0,1mm.xls",
        "NE":     "Reologia/Frequency/C15 Gira 5.5  - Freq CP50 0,1mm.xls",
        "Bozano": "Reologia/Frequency/Gel de Cabelo Bozano - Freq CP50 0,1mm.xls",
    }
    strain_files = {
        "C15":    "Reologia/Strain/C15 - Deformacao CP50 0,1mm.xls",
        "NE":     "Reologia/Strain/C15 Gira 5.5  - Deformacao CP50 0,1mm.xls",
        "Bozano": "Reologia/Strain/Gel de Cabelo Bozano - Deformacao CP50 0,1mm.xls",
    }

    print("=" * 70)
    print("STRUCTURE DUMP — first .xls of each test")
    print("=" * 70)
    dump_sheets(strain_files["C15"])
    print()
    dump_sheets(freq_files["C15"])

    print("\n" + "=" * 70)
    print("AMPLITUDE (STRAIN) SWEEP — LVR plateau + gamma_LVR (10% drop)")
    print("=" * 70)
    strain_results = {}
    for name, path in strain_files.items():
        print(f"\n>>> {name}: {path}")
        if not os.path.exists(path):
            print("   FILE NOT FOUND")
            continue
        info, headers, units = extract_amplitude_sweep(path)
        if info is None:
            print("   No data parsed.")
            continue
        if "error" in info:
            print(f"   PARSE ERROR: {info['error']}")
            print(f"   headers: {headers}")
            print(f"   units:   {units}")
            continue
        strain_results[name] = info
        print(f"   sheet={info['sheet_index']}  rows={info['n_rows']}")
        print(f"   strain range: {info['strain_range']}")
        print(f"   first 5 strain: {info['strain_first5']}")
        print(f"   first 5 Gp:     {info['Gp_first5']}")
        print(f"   G'_LVR  = {info['Gp_LVR']:.2f} Pa")
        print(f"   G''_LVR = {info['Gpp_LVR']:.2f} Pa")
        print(f"   tan d   = {info['tan_LVR']:.4f}")
        print(f"   gamma_LVR (10% drop) = {info['gamma_LVR_pct']:.3f} %")

    print("\n" + "=" * 70)
    print("FREQUENCY SWEEP — G' and tan d at omega ~ 1, 10, 100 rad/s")
    print("=" * 70)
    freq_results = {}
    for name, path in freq_files.items():
        print(f"\n>>> {name}: {path}")
        if not os.path.exists(path):
            print("   FILE NOT FOUND")
            continue
        info, headers, units = extract_frequency_sweep(path)
        if info is None:
            print("   No data parsed.")
            continue
        if "error" in info:
            print(f"   PARSE ERROR: {info['error']}")
            print(f"   headers: {headers}")
            print(f"   units:   {units}")
            continue
        freq_results[name] = info
        print(f"   sheet={info['sheet_index']}  rows={info['n_rows']}  "
              f"omega unit: {info['omega_unit']}  range={info['omega_range']}")
        for tgt, pt in info["points"].items():
            print(f"   omega ~ {tgt:>5} rad/s  ->  omega_actual={pt['omega_actual']:.3f}, "
                  f"G'={pt['Gp']:.2f} Pa,  G''={pt['Gpp']:.2f} Pa,  tan={pt['tan']:.4f}")

    print("\n" + "=" * 70)
    print("SUMMARY — values for manuscript tables")
    print("=" * 70)
    print("\n[Amplitude sweep — Table tab:strain_summary]")
    print(f"{'Group':<8} | {'Gp(Pa)':>8} | {'Gpp(Pa)':>8} | {'tan d':>7} | {'gLVR(%)':>8}")
    for name, info in strain_results.items():
        print(f"{name:<8} | {info['Gp_LVR']:>8.1f} | {info['Gpp_LVR']:>8.1f} | "
              f"{info['tan_LVR']:>7.3f} | {info['gamma_LVR_pct']:>8.1f}")

    print("\n[Frequency sweep — Table tab:freq_summary]")
    print(f"{'Group':<8} | {'w':>5} | {'w_real':>7} | {'Gp(Pa)':>8} | {'Gpp(Pa)':>8} | {'tan d':>7}")
    for name, info in freq_results.items():
        for tgt, pt in info["points"].items():
            print(f"{name:<8} | {tgt:>5.0f} | {pt['omega_actual']:>7.3f} | "
                  f"{pt['Gp']:>8.1f} | {pt['Gpp']:>8.1f} | {pt['tan']:>7.3f}")


if __name__ == "__main__":
    main()
