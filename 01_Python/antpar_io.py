"""
antpar_io.py
============

I/O module for Anton Paar Rheocompass CSV exports.

This module replaces the legacy ``xlrd``-based ``.xls`` ingestion. The
Rheocompass software exports test results as **UTF-16 LE, tab-separated**
CSVs with a **comma decimal separator** (Portuguese / European locale).

Supported test types
--------------------
* Flow / viscosity curve   →  ``read_flow_curve_csv``       →  ``FlowCurve``
* Amplitude (strain) sweep →  ``read_amplitude_sweep_csv``  →  ``AmplitudeSweep``
* Frequency sweep          →  ``read_frequency_sweep_csv``  →  ``FrequencySweep``
* Recovery interval (3iTT) →  ``read_recovery_interval_csv``→  ``RecoveryInterval``
* Recovery set (3iTT trio) →  ``read_recovery_set(folder, basename)`` →  ``RecoverySet``

Folder helpers
--------------
* ``list_raw_csvs(folder)``                  — every raw CSV, skips
  ``*_Analysis.csv`` (English) and ``*_Analise.csv`` (Portuguese).
* ``list_flow_curves``,  ``list_amplitude_sweeps`` — back-compat aliases.
* ``read_flow_curves_in_folder``             — dict[short_name, FlowCurve]
* ``read_amplitude_sweeps_in_folder``        — dict[short_name, AmplitudeSweep]

Format contract (Rheocompass V1.36.x raw export)
------------------------------------------------
* Encoding:  UTF-16 LE with BOM (``0xFF 0xFE``)
* Separator: TAB (``\\t``)
* Decimal:   COMMA (``,``)
* Line endings: CRLF
* Header lines 1-8 contain Projeto / Teste / Resultado / Intervalo /
  column-header row. Line 10 holds units in brackets. Data starts at
  line 11.

CRITICAL UNIT CONVERSIONS
-------------------------
* Flow curve viscosity column is exported in **mPa·s** → converted to **Pa·s**.
* Amplitude sweep strain column is exported as **fractional** ([1]) → kept
  fractional, but a convenience column ``strain_pct`` (= strain × 100) is
  added because LVR results are conventionally reported in percent.
* Amplitude sweep torque is exported in **µN·m** (≠ flow curve mN·m) → kept
  as-is in the column ``torque_uNm``.

Author: T.M.C. Rodrigues - 2026-05-20
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# Data containers
# ----------------------------------------------------------------------

@dataclass
class FlowCurve:
    """One parsed flow-curve CSV (viscosity vs shear rate)."""
    sample_name: str
    project: str
    interval_label: str
    n_points: int
    source_path: Path
    data: pd.DataFrame = field(repr=False)
    # Columns: point, shear_rate_1_s, shear_stress_Pa,
    #          viscosity_Pa_s (CONVERTED from mPa·s),
    instrument: str = ""     # B20: from the export header, "" if absent
    #          temperature_C, torque_mNm, status


@dataclass
class AmplitudeSweep:
    """One parsed amplitude (strain) sweep CSV — SAOS, fixed frequency."""
    sample_name: str
    project: str
    interval_label: str
    n_points: int
    source_path: Path
    data: pd.DataFrame = field(repr=False)
    # Columns: point, strain_frac, strain_pct (= strain_frac*100),
    #          shear_stress_Pa, G_storage_Pa, G_loss_Pa, tan_delta,
    #          torque_uNm, status
    instrument: str = ""     # B20: from the export header, "" if absent
    # NOTE: no temperature column — Rheocompass amplitude sweeps don't
    # emit T. Use the filename / Teste: header to recover setpoint.


@dataclass
class FrequencySweep:
    """One parsed frequency sweep CSV — SAOS, fixed strain (inside LVR)."""
    sample_name: str
    project: str
    interval_label: str
    n_points: int
    source_path: Path
    data: pd.DataFrame = field(repr=False)
    # Columns: point, omega_rad_s, G_storage_Pa, G_loss_Pa, tan_delta,
    #          eta_complex_real_Pa_s (CONVERTED from mPa·s),
    #          eta_complex_imag_Pa_s (CONVERTED from mPa·s),
    #          G_complex_Pa, strain_pct, shear_stress_Pa, torque_mNm, status
    instrument: str = ""     # B20: from the export header, "" if absent
    # NOTE: no temperature column (same as amplitude sweep). Strain here is
    # the controlled SAOS amplitude in PERCENT (≠ amplitude-sweep [1]
    # fractional convention). Torque is in mN·m (≠ amplitude-sweep µN·m).


@dataclass
class RecoveryInterval:
    """One interval of a 3iTT recovery test (a single ``_<N>.csv`` file)."""
    interval_num: int               # 1, 2, or 3
    kind: str                       # "oscillation" or "rotation"
    sample_name: str
    project: str
    interval_label: str             # e.g. "3ITT-Osc-Rot-Osc 1"
    n_points: int
    source_path: Path
    data: pd.DataFrame = field(repr=False)
    # Oscillation columns: point, time_s, G_storage_Pa, G_loss_Pa, tan_delta,
    #                      eta_complex_Pa_s (from mPa·s), eta_steady_Pa_s (NaN
    #                      — only meaningful in rotation), torque_uNm, status
    # Rotation columns:    point, time_s, shear_rate_1_s, shear_stress_Pa,
    #                      viscosity_Pa_s (from mPa·s), status


@dataclass
class RecoverySet:
    """A complete 3iTT-Osc-Rot-Osc recovery set (three CSV files stitched)."""
    sample_name: str
    project: str
    interval_label: str
    source_basename: str
    interval1: RecoveryInterval     # oscillation — baseline
    interval2: RecoveryInterval     # rotation    — disruption
    interval3: RecoveryInterval     # oscillation — recovery
    # NOTE: For other 3iTT variants (Rot-Rot-Rot, etc.) the kinds will differ;
    # this dataclass holds whatever the underlying files contain. The
    # downstream extractor checks ``kind`` per interval.


# ----------------------------------------------------------------------
# Shared low-level helpers
# ----------------------------------------------------------------------

_HEADER_LABELS = {
    "projeto": "project",
    "teste": "sample_name",
    "resultado": "interval_label",
    # Instrument identity. The header block always carried this field; the
    # parser simply never extracted it, which is how an ARES-G2 dataset sat
    # mislabelled as an MCR in a project CLAUDE.md for months. Downstream
    # scripts must echo `instrument` into their report headers so the
    # machine that produced a number travels with the number.
    "instrumento": "instrument",
    "instrument": "instrument",
    "instrument name": "instrument",
    "dispositivo": "instrument",
    "device": "instrument",
}

# Filename suffixes that mark Rheocompass "analysis" exports (regression
# outputs, NOT raw measurement data). Both spellings encountered in the
# wild: English in newer builds, Portuguese in older / localised ones.
_ANALYSIS_SUFFIXES = ("_analysis.csv", "_analise.csv")


def _decode_utf16(path: Path) -> str:
    """Read an Anton Paar CSV as text, transparently stripping the UTF-16 BOM."""
    raw = Path(path).read_bytes()
    if raw.startswith(b"\xff\xfe"):
        return raw.decode("utf-16-le")[1:]   # drop BOM character
    if raw.startswith(b"\xfe\xff"):
        return raw.decode("utf-16-be")[1:]
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError as e:
        raise ValueError(
            f"{path.name}: cannot decode as UTF-16 LE/BE or UTF-8. "
            f"This parser expects Rheocompass CSV exports."
        ) from e


def _split_tabs(line: str) -> list[str]:
    """Split a tab-separated line and strip whitespace from each field."""
    return [cell.strip() for cell in line.split("\t")]


def _parse_decimal_comma(token: str) -> float:
    """Parse a number that uses ``,`` as the decimal separator. Empty → NaN.

    Also handles scientific notation written with comma (``9,9E-05``) and
    plain integer cells (``133``).
    """
    s = token.strip()
    if not s or s.lower() in {"nan", "n/a", "-"}:
        return float("nan")
    return float(s.replace(",", "."))


def _parse_header_block(lines: list[str]) -> dict:
    """Extract Projeto / Teste / Resultado / n_points from the top of the file."""
    meta: dict = {}
    for ln in lines[:10]:
        if not ln.strip():
            continue
        first, _, rest = ln.partition("\t")
        key = first.strip().rstrip(":").lower()
        if key in _HEADER_LABELS:
            meta[_HEADER_LABELS[key]] = rest.split("\t")[0].strip()
        elif "intervalo e ponto de dados" in key:
            parts = [p for p in rest.split("\t") if p.strip()]
            if len(parts) >= 2:
                try:
                    meta["n_points"] = int(parts[-1])
                except ValueError:
                    pass
    return meta


def _find_units_row(lines: list[str], required_token: str) -> int:
    """Return the index of the units row matching ``required_token`` (e.g. ``[1/s]``)."""
    for i, ln in enumerate(lines):
        if required_token in ln:
            return i
    raise ValueError(
        f"Could not find units row containing {required_token!r}. "
        f"File structure has drifted from the expected Rheocompass layout."
    )


def _iter_data_rows(lines: list[str], data_start: int, n_cols: int):
    """Yield (point, cells_after_leading_tab) for every numeric row."""
    for ln in lines[data_start:]:
        cells = _split_tabs(ln)
        # Rheocompass starts data lines with a leading TAB so the point
        # number lands in column index 1. Drop the leading empty cell.
        if cells and cells[0] == "":
            cells = cells[1:]
        if len(cells) < n_cols:
            continue
        try:
            point = int(cells[0])
        except ValueError:
            continue
        yield point, cells


# ----------------------------------------------------------------------
# Flow curve  (viscosity vs shear rate)
# ----------------------------------------------------------------------

def read_flow_curve_csv(path: str | Path) -> FlowCurve:
    """
    Parse one raw Anton Paar Rheocompass flow-curve CSV.

    Viscosity in the returned DataFrame is in **Pa·s** (converted from the
    mPa·s as exported by Rheocompass).
    """
    path = Path(path)
    if path.name.lower().endswith(_ANALYSIS_SUFFIXES):
        raise ValueError(
            f"{path.name} is a Rheocompass *analysis* (regression) export, "
            "not a raw curve. Skip it."
        )

    text = _decode_utf16(path)
    lines = text.splitlines()
    meta = _parse_header_block(lines)
    data_start = _find_units_row(lines, "[1/s]") + 1

    rows = []
    for point, cells in _iter_data_rows(lines, data_start, n_cols=7):
        try:
            shear_rate = _parse_decimal_comma(cells[1])
            shear_stress = _parse_decimal_comma(cells[2])
            viscosity_mPa_s = _parse_decimal_comma(cells[3])
            temperature = _parse_decimal_comma(cells[4])
            torque = _parse_decimal_comma(cells[5])
        except ValueError as e:
            raise ValueError(f"{path.name} (point {point}): {e}") from e
        status = cells[6] if len(cells) > 6 else ""
        rows.append((point, shear_rate, shear_stress,
                     viscosity_mPa_s / 1000.0,   # mPa·s → Pa·s
                     temperature, torque, status))

    if not rows:
        raise ValueError(f"{path.name}: no numeric data rows parsed.")

    df = pd.DataFrame(rows, columns=[
        "point", "shear_rate_1_s", "shear_stress_Pa",
        "viscosity_Pa_s", "temperature_C", "torque_mNm", "status",
    ])

    return FlowCurve(
        instrument=meta.get("instrument", ""),
        sample_name=meta.get("sample_name", path.stem),
        project=meta.get("project", ""),
        interval_label=meta.get("interval_label", ""),
        n_points=meta["n_points"] if isinstance(meta.get("n_points"), int) else len(df),
        source_path=path,
        data=df,
    )


# ----------------------------------------------------------------------
# Force-augmented flow curves — normal stress / N1
# ----------------------------------------------------------------------
#
# WHY THIS EXISTS
# ---------------
# The die-swell closure in the MATLAB solver needs the first normal-stress
# difference N1. Rheocompass' standard flow-curve export does not carry a
# normal-force channel, so N1 has to come from a re-export with the axial
# force column switched on. That produces a DIFFERENT file format from
# everything else this module reads, and every step of handling it — the
# encoding, the column mapping, the unit conversion, the baseline
# correction — was previously rebuilt ad hoc in throwaway scripts.
#
# FORMAT WARNING (this is the trap)
# ---------------------------------
# The force-augmented export is **latin-1**, NOT the UTF-16 LE this module
# assumes everywhere else. Assuming UTF-16 does not raise — it silently
# mojibakes. `_decode_force_export` therefore never calls `_decode_utf16`.
#
# PROVENANCE OF THIS PARSER
# -------------------------
# Written from the documented format contract, NOT against a specimen file
# (the force exports are not present in every checkout of this pipeline).
# It is deliberately header-driven rather than positional so it survives a
# column-order change, and `read_flow_curve_with_force` raises a loud,
# specific error rather than guessing when it cannot find a channel.
# BEFORE TRUSTING ITS OUTPUT on a new machine, run it once against a file
# whose values you can read by eye and confirm the mapping.

_FORCE_ENCODING = "latin-1"

# Canonical channel names → the header substrings that mark them. Matched
# case-insensitively against the normalised header cell. Order matters:
# the first hit wins, so put the more specific token first.
_FORCE_COLUMN_TOKENS = {
    "shear_rate_1_s":   ("shear rate", "taxa de cisalhamento", "gamma dot"),
    "shear_stress_Pa":  ("shear stress", "tensao de cisalhamento", "stress"),
    "viscosity_Pa_s":   ("viscosity", "viscosidade"),
    "nstress_Pa":       ("normal stress", "tensao normal"),
    "force_N":          ("axial force", "normal force", "forca axial", "force"),
    "time_s":           ("step time", "tempo", "time"),
    "temperature_C":    ("temperature", "temperatura"),
}

# Geometry radii (m). N1 = 2F/(pi R^2) needs the PLATE/CONE radius, not the
# nominal diameter in the geometry's name.
GEOMETRY_RADIUS_M = {
    "CP50": 25.0e-3,   # 50 mm cone   → 25 mm radius
    "PP50": 25.0e-3,   # 50 mm plate  → 25 mm radius
    "CP25": 12.5e-3,
    "PP25": 12.5e-3,
}


def normal_stress_from_force(F_newtons, R_meters: float):
    """
    First normal-stress difference from the measured axial force.

        N1 = 2F / (pi R^2)

    Parameters
    ----------
    F_newtons : float or array-like
        Axial (normal) force, in NEWTONS. If the export is in grams-force
        or millinewtons, convert BEFORE calling — this function does no
        unit sniffing.
    R_meters : float
        Cone/plate radius in metres (25e-3 for a 50 mm geometry — the
        RADIUS, not the diameter that names the fixture). See
        ``GEOMETRY_RADIUS_M``.

    Returns
    -------
    float or numpy.ndarray — N1 in Pa.
    """
    if R_meters <= 0:
        raise ValueError(f"R_meters must be positive, got {R_meters}")
    return 2.0 * np.asarray(F_newtons, dtype=float) / (np.pi * R_meters ** 2)


def baseline_correct(series, method: str = "min"):
    """
    Remove the per-run transducer offset from a force / normal-stress channel.

    WHY IT IS NEEDED
    ----------------
    The normal-force transducer is not necessarily re-zeroed between sample
    loadings, so raw readings carry a large, sign-varying offset that is
    constant within one run. Without a tare the reported N1 is meaningless,
    and a negative N1 (physically impossible in steady shear on these gels)
    is the usual symptom.

    METHOD
    ------
    ``"min"`` (the correct convention): floor = min(series) over the whole
    run; corrected = series - floor. The run's quietest point is taken as
    the zero.

    KNOWN, QUANTIFIED BIAS OF THE MIN CONVENTION. It assumes the quietest
    point has zero TRUE normal stress. On a shear ramp the quietest point is
    the lowest shear rate, where the material usually already has a small but
    non-zero N1, so the tare over-subtracts by exactly N1(gamma_dot_min) and
    every corrected value is low by that amount. Verified on a synthetic
    power-law curve: true N1 at 200 1/s = 1590.75 Pa, recovered = 1564.55 Pa,
    short by 26.17 Pa = N1(0.1 1/s) exactly, a 1.6% under-estimate.

    The bias is always DOWNWARD, so it makes the resulting Tanner die-swell
    estimate conservative rather than optimistic — which is the safe
    direction, but it is not zero. It shrinks as the flow curve starts at a
    lower shear rate. Start the ramp as low as the torque resolution allows,
    and if you need N1 to better than a few percent, measure a genuine
    zero-shear baseline instead of inferring one.

    ``"median_low25"`` — REJECTED, KEPT ONLY SO IT IS NOT RETRIED. Taking
    the median of the 25 lowest points was tried first and gives a
    materially wrong answer: on the C20 run it returned ~460 Pa where the
    min convention returns ~802 Pa, a 74% error, because the low tail of a
    ramp is populated by genuine low-shear data, not by baseline. Selecting
    it raises. If you believe you need it, you have found a different
    problem — probably a drifting rather than offset baseline, which this
    function is the wrong tool for.

    Parameters
    ----------
    series : array-like
    method : {"min"}

    Returns
    -------
    (corrected, floor) : (numpy.ndarray, float)
    """
    arr = np.asarray(series, dtype=float)
    if method == "median_low25":
        raise ValueError(
            "baseline_correct(method='median_low25') is a REJECTED convention. "
            "It was tried and gives ~460 Pa where the correct min convention "
            "gives ~802 Pa on the same C20 run (74% error), because the low "
            "tail of a shear ramp holds real data, not baseline. Use "
            "method='min'."
        )
    if method != "min":
        raise ValueError(f"unknown baseline method {method!r}; only 'min' is supported")
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("baseline_correct: series has no finite values")
    floor = float(np.min(finite))
    return arr - floor, floor


# Filename pairing between the audited .xls workbooks and the force .txt
# re-exports. These are KNOWN ONE-OFF EXCEPTIONS, not a naming scheme, so
# this is a hard-coded lookup rather than fuzzy matching — fuzzy matching
# would silently pick the wrong C10 file, which is the failure this table
# exists to prevent. Extend it by hand when a new ink is re-exported.
_AUDITED_FILE_PAIRS = {
    # ink : (xls basename, force-txt basename, note)
    "C25": (
        "C25",
        "C25",
        "the force export has a DOUBLE space where the .xls has one",
    ),
    "gel": (
        "gel (Ac)",
        "gel ",
        "the force export drops the '(Ac)' and adds a TRAILING space",
    ),
    "C10": (
        "C10 (1)",
        "C10 (1)",
        "TWO candidate files exist; only the one with '(1)' is audited",
    ),
}


def pair_audited_files(ink: str) -> dict:
    """
    Return the audited .xls / force-.txt basenames for one ink.

    The pairing is trap-laden and the traps are not systematic — see
    ``_AUDITED_FILE_PAIRS``. Always go through this function rather than
    globbing, so the exceptions are discovered once and stay discovered.

    Raises KeyError with the list of known inks if `ink` is not registered.
    """
    key = ink.strip()
    if key not in _AUDITED_FILE_PAIRS:
        raise KeyError(
            f"{ink!r} has no audited file pairing registered. Known: "
            f"{sorted(_AUDITED_FILE_PAIRS)}. Add it to _AUDITED_FILE_PAIRS "
            f"by hand after checking the actual filenames — do NOT glob."
        )
    xls, txt, note = _AUDITED_FILE_PAIRS[key]
    return {"ink": key, "xls_basename": xls, "force_txt_basename": txt, "note": note}


@dataclass
class FlowCurveWithForce:
    """A flow curve that additionally carries a tared normal-stress channel."""
    sample_name: str
    project: str
    instrument: str
    geometry: str
    radius_m: float
    baseline_floor_Pa: float
    n_points: int
    source_path: Path
    data: "pd.DataFrame" = field(repr=False)
    # Columns: point, shear_rate_1_s, shear_stress_Pa, viscosity_Pa_s,
    #          nstress_raw_Pa, N1_Pa (tared), plus time_s / temperature_C
    #          when the export carries them.


def _decode_force_export(path: Path) -> list[str]:
    """Read a force-augmented export as latin-1 text, CRLF-tolerant."""
    raw = Path(path).read_bytes()
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        raise ValueError(
            f"{Path(path).name} is UTF-16, not the latin-1 force-export "
            "format. Use read_flow_curve_csv for Rheocompass exports."
        )
    return raw.decode(_FORCE_ENCODING).replace("\r\n", "\n").split("\n")


def _normalise_header_cell(cell: str) -> str:
    return " ".join(cell.strip().strip('"').lower().split())


def _map_force_columns(header_cells: list[str]) -> dict:
    """Map canonical channel names to column indices, header-driven."""
    norm = [_normalise_header_cell(c) for c in header_cells]
    mapping: dict[str, int] = {}
    for canon, tokens in _FORCE_COLUMN_TOKENS.items():
        for tok in tokens:
            hit = next((i for i, h in enumerate(norm)
                        if tok in h and i not in mapping.values()), None)
            if hit is not None:
                mapping[canon] = hit
                break
    return mapping


def _parse_any_decimal(token: str) -> float:
    """Parse a number written with either a comma or a point decimal mark."""
    t = token.strip().strip('"')
    if not t:
        return float("nan")
    if "," in t and "." in t:
        # e.g. 1.234,56 (European thousands) vs 1,234.56 (US thousands)
        t = t.replace(".", "").replace(",", ".") if t.rfind(",") > t.rfind(".") \
            else t.replace(",", "")
    elif "," in t:
        t = t.replace(",", ".")
    try:
        return float(t)
    except ValueError:
        return float("nan")


def read_flow_curve_with_force(
    path: str | Path,
    geometry: str = "CP50",
    radius_m: float | None = None,
    force_unit: str = "N",
) -> FlowCurveWithForce:
    """
    Parse a force-augmented flow-curve export and return a TARED N1 channel.

    One call replaces the whole ad hoc chain: latin-1 decode (B16), header-
    driven column mapping, either decimal convention, N1 = 2F/(pi R^2)
    conversion (B19) when only a force channel is present, and the
    min-baseline tare (A5).

    Parameters
    ----------
    path : file to read.
    geometry : fixture name, used to look up the radius (see
        ``GEOMETRY_RADIUS_M``). Ignored if `radius_m` is given.
    radius_m : explicit radius in metres; overrides `geometry`.
    force_unit : {"N", "mN", "gf"} — unit of the axial-force column, if the
        export carries force rather than normal stress. No sniffing is done;
        state it.

    Returns
    -------
    FlowCurveWithForce, whose ``data`` frame carries ``N1_Pa`` alongside the
    usual flow-curve columns. ``baseline_floor_Pa`` records what was
    subtracted, so the tare is auditable rather than invisible.
    """
    path = Path(path)
    lines = _decode_force_export(path)

    # ---- locate the header row: the first line naming a shear-rate channel
    hdr_idx = None
    for i, line in enumerate(lines[:200]):
        cells = line.split("\t") if "\t" in line else line.split()
        norm = " ".join(_normalise_header_cell(c) for c in cells)
        if any(tok in norm for tok in _FORCE_COLUMN_TOKENS["shear_rate_1_s"]):
            hdr_idx = i
            break
    if hdr_idx is None:
        raise ValueError(
            f"{path.name}: no header row naming a shear-rate channel was found "
            f"in the first 200 lines. This parser is header-driven by design; "
            f"check that the export actually carries column titles."
        )

    header_cells = lines[hdr_idx].split("\t")
    cols = _map_force_columns(header_cells)
    if "shear_rate_1_s" not in cols:
        raise ValueError(f"{path.name}: shear-rate column not identified.")
    if "nstress_Pa" not in cols and "force_N" not in cols:
        raise ValueError(
            f"{path.name}: neither a normal-stress nor an axial-force column "
            f"was found. Re-export from the instrument software with the "
            f"axial force channel switched on — that is the whole point of "
            f"this reader. Columns seen: "
            f"{[_normalise_header_cell(c) for c in header_cells]}"
        )

    # ---- instrument identity (B20), scanned from the pre-header block
    instrument = ""
    for line in lines[:hdr_idx]:
        low = line.lower()
        for label in ("instrument name", "instrument", "instrumento", "device"):
            if label in low:
                parts = line.split("\t") if "\t" in line else line.split(":", 1)
                if len(parts) > 1 and parts[-1].strip():
                    instrument = parts[-1].strip().strip('"')
                    break
        if instrument:
            break

    # ---- data rows
    rows = []
    for line in lines[hdr_idx + 1:]:
        if not line.strip():
            continue
        cells = line.split("\t")
        if len(cells) <= max(cols.values()):
            continue
        vals = {c: _parse_any_decimal(cells[i]) for c, i in cols.items()}
        if not np.isfinite(vals["shear_rate_1_s"]):
            continue        # units row, section break, or trailing text
        rows.append(vals)

    if not rows:
        raise ValueError(f"{path.name}: header found but no numeric data rows parsed.")

    df = pd.DataFrame(rows)
    df.insert(0, "point", np.arange(1, len(df) + 1))

    R = radius_m if radius_m is not None else GEOMETRY_RADIUS_M.get(geometry)
    if R is None:
        raise ValueError(
            f"unknown geometry {geometry!r}; pass radius_m explicitly or add "
            f"it to GEOMETRY_RADIUS_M. Known: {sorted(GEOMETRY_RADIUS_M)}"
        )

    # ---- raw normal stress: measured directly, or from the force channel
    if "nstress_Pa" in df.columns and np.isfinite(df["nstress_Pa"]).any():
        nstress_raw = df["nstress_Pa"].to_numpy(dtype=float)
    else:
        scale = {"N": 1.0, "mN": 1e-3, "gf": 9.80665e-3}[force_unit]
        nstress_raw = normal_stress_from_force(
            df["force_N"].to_numpy(dtype=float) * scale, R)

    N1, floor = baseline_correct(nstress_raw, method="min")
    df["nstress_raw_Pa"] = nstress_raw
    df["N1_Pa"] = N1

    return FlowCurveWithForce(
        sample_name=path.stem,
        project="",
        instrument=instrument,
        geometry=geometry,
        radius_m=R,
        baseline_floor_Pa=floor,
        n_points=len(df),
        source_path=path,
        data=df,
    )


# ----------------------------------------------------------------------
# Amplitude (strain) sweep — SAOS, fixed ω, swept γ
# ----------------------------------------------------------------------

def read_amplitude_sweep_csv(path: str | Path) -> AmplitudeSweep:
    """
    Parse one raw Anton Paar Rheocompass amplitude-sweep CSV.

    Columns (in the returned DataFrame):
        point, strain_frac, strain_pct,
        shear_stress_Pa, G_storage_Pa, G_loss_Pa, tan_delta,
        torque_uNm, status

    Strain is stored both fractional (``strain_frac``, native Rheocompass
    unit ``[1]``) and as percent (``strain_pct`` = strain_frac × 100) for
    convenience — γ_LVR is conventionally reported in %.
    """
    path = Path(path)
    if path.name.lower().endswith(_ANALYSIS_SUFFIXES):
        raise ValueError(
            f"{path.name} is a Rheocompass *analysis* (regression) export, "
            "not a raw sweep. Skip it."
        )

    text = _decode_utf16(path)
    lines = text.splitlines()
    meta = _parse_header_block(lines)

    # The amplitude-sweep units row contains [1] (strain), three [Pa]
    # columns (stress, G', G''), [1] again (tan δ), and [µN·m] (torque).
    # "[µN·m]" is unique to amplitude sweeps — flow curves use [mN·m],
    # frequency sweeps use [mN·m], so the µ prefix is the discriminator.
    try:
        data_start = _find_units_row(lines, "[µN·m]") + 1   # µN·m
    except ValueError:
        # Some Rheocompass builds emit ASCII "uN.m" or "[uN*m]" instead.
        # Fall back to the column-header row.
        try:
            data_start = _find_units_row(lines, "Deforma") + 3   # header + blank + units
        except ValueError as e:
            raise ValueError(
                f"{path.name}: not an amplitude-sweep export. "
                "Expected '[µN·m]' or 'Deforma' anchor in header."
            ) from e

    rows = []
    for point, cells in _iter_data_rows(lines, data_start, n_cols=8):
        try:
            strain_frac = _parse_decimal_comma(cells[1])
            shear_stress = _parse_decimal_comma(cells[2])
            g_storage = _parse_decimal_comma(cells[3])
            g_loss = _parse_decimal_comma(cells[4])
            tan_delta = _parse_decimal_comma(cells[5])
            torque_uNm = _parse_decimal_comma(cells[6])
        except ValueError as e:
            raise ValueError(f"{path.name} (point {point}): {e}") from e
        status = cells[7] if len(cells) > 7 else ""
        rows.append((point, strain_frac, strain_frac * 100.0,
                     shear_stress, g_storage, g_loss, tan_delta,
                     torque_uNm, status))

    if not rows:
        raise ValueError(f"{path.name}: no numeric data rows parsed.")

    df = pd.DataFrame(rows, columns=[
        "point", "strain_frac", "strain_pct",
        "shear_stress_Pa", "G_storage_Pa", "G_loss_Pa", "tan_delta",
        "torque_uNm", "status",
    ])

    return AmplitudeSweep(
        instrument=meta.get("instrument", ""),
        sample_name=meta.get("sample_name", path.stem),
        project=meta.get("project", ""),
        interval_label=meta.get("interval_label", ""),
        n_points=meta["n_points"] if isinstance(meta.get("n_points"), int) else len(df),
        source_path=path,
        data=df,
    )


# ----------------------------------------------------------------------
# Frequency sweep — SAOS, fixed strain (in LVR), swept ω
# ----------------------------------------------------------------------

def read_frequency_sweep_csv(path: str | Path) -> FrequencySweep:
    """
    Parse one raw Anton Paar Rheocompass frequency-sweep CSV.

    Columns (in the returned DataFrame):
        point, omega_rad_s,
        G_storage_Pa, G_loss_Pa, tan_delta,
        eta_complex_real_Pa_s, eta_complex_imag_Pa_s,
        G_complex_Pa,
        strain_pct, shear_stress_Pa,
        torque_mNm, status

    Complex viscosity (Re and Im) is stored in **Pa·s** (converted from
    Rheocompass's native mPa·s). Strain is in **percent** (Rheocompass's
    native unit for this test type — NOT fractional like amplitude sweeps,
    where the unit is [1]).
    """
    path = Path(path)
    if path.name.lower().endswith(_ANALYSIS_SUFFIXES):
        raise ValueError(
            f"{path.name} is a Rheocompass *analysis* (regression) export, "
            "not a raw sweep. Skip it."
        )

    text = _decode_utf16(path)
    lines = text.splitlines()
    meta = _parse_header_block(lines)

    # "[rad/s]" is unique to the frequency sweep (flow curves use [1/s]
    # for shear rate; amplitude sweeps have no rad/s column at all).
    try:
        data_start = _find_units_row(lines, "[rad/s]") + 1
    except ValueError as e:
        raise ValueError(
            f"{path.name}: not a frequency-sweep export. "
            "Expected '[rad/s]' anchor in header."
        ) from e

    rows = []
    for point, cells in _iter_data_rows(lines, data_start, n_cols=12):
        try:
            omega = _parse_decimal_comma(cells[1])
            g_storage = _parse_decimal_comma(cells[2])
            g_loss = _parse_decimal_comma(cells[3])
            tan_delta = _parse_decimal_comma(cells[4])
            eta_re_mPa_s = _parse_decimal_comma(cells[5])
            eta_im_mPa_s = _parse_decimal_comma(cells[6])
            g_complex = _parse_decimal_comma(cells[7])
            strain_pct = _parse_decimal_comma(cells[8])
            shear_stress = _parse_decimal_comma(cells[9])
            torque_mNm = _parse_decimal_comma(cells[10])
        except ValueError as e:
            raise ValueError(f"{path.name} (point {point}): {e}") from e
        status = cells[11] if len(cells) > 11 else ""
        rows.append((point, omega,
                     g_storage, g_loss, tan_delta,
                     eta_re_mPa_s / 1000.0,   # mPa·s → Pa·s
                     eta_im_mPa_s / 1000.0,   # mPa·s → Pa·s
                     g_complex, strain_pct, shear_stress,
                     torque_mNm, status))

    if not rows:
        raise ValueError(f"{path.name}: no numeric data rows parsed.")

    df = pd.DataFrame(rows, columns=[
        "point", "omega_rad_s",
        "G_storage_Pa", "G_loss_Pa", "tan_delta",
        "eta_complex_real_Pa_s", "eta_complex_imag_Pa_s",
        "G_complex_Pa", "strain_pct", "shear_stress_Pa",
        "torque_mNm", "status",
    ])

    return FrequencySweep(
        instrument=meta.get("instrument", ""),
        sample_name=meta.get("sample_name", path.stem),
        project=meta.get("project", ""),
        interval_label=meta.get("interval_label", ""),
        n_points=meta["n_points"] if isinstance(meta.get("n_points"), int) else len(df),
        source_path=path,
        data=df,
    )


# ----------------------------------------------------------------------
# Recovery — 3iTT (3-Interval Thixotropy Test)
# ----------------------------------------------------------------------
#
# Rheocompass exports a 3iTT measurement as THREE separate CSV files,
# one per interval, all sharing the same basename and differing only in
# the trailing ``_<N>.csv`` suffix (N = 1, 2, 3). There is also a
# ``*_analise.csv`` sibling that holds Rheocompass's built-in regression
# / recovery output — filtered out by ``list_raw_csvs``.
#
# For the Osc-Rot-Osc variant (which is what's in the wild so far):
#   * Interval 1 = oscillation at low strain (in LVR) → baseline G', |η*|
#   * Interval 2 = rotation at high shear rate     → disruption phase
#   * Interval 3 = oscillation at low strain again → recovery trajectory
#
# The reader auto-detects oscillation vs rotation per file (oscillation
# files have a G' / "Módulo de Armazenamento" header; rotation files
# have a "Taxa de Cisalhamento" header).


def _parse_interval_number(meta: dict, lines: list[str]) -> int:
    """Find the interval number from the 'Intervalo e ponto de dados:' line."""
    for ln in lines[:10]:
        if "intervalo e ponto de dados" in ln.lower():
            parts = [p.strip() for p in ln.split("\t") if p.strip()]
            # parts[0] = "Intervalo e ponto de dados:", parts[1] = N, parts[2] = n_points
            if len(parts) >= 2:
                try:
                    return int(parts[1])
                except ValueError:
                    pass
    raise ValueError("Could not locate interval number in header.")


def read_recovery_interval_csv(path: str | Path) -> RecoveryInterval:
    """
    Parse one interval of a 3iTT recovery CSV.

    Auto-detects whether the interval is *oscillation* (G', G'', tan δ,
    |η*|) or *rotation* (shear rate, stress, viscosity) by inspecting
    the column-header row.

    Time column is in **seconds** and is preserved as `time_s`. For
    oscillation intervals, complex viscosity is converted from mPa·s
    to Pa·s (``eta_complex_Pa_s``); steady viscosity is also stored
    but is NaN (Rheocompass leaves that column blank during oscillation).
    For rotation intervals, steady viscosity is converted to Pa·s.
    """
    path = Path(path)
    if path.name.lower().endswith(_ANALYSIS_SUFFIXES):
        raise ValueError(
            f"{path.name} is a Rheocompass *analysis* (regression) export, "
            "not a raw recovery interval. Skip it."
        )

    text = _decode_utf16(path)
    lines = text.splitlines()
    meta = _parse_header_block(lines)
    interval_num = _parse_interval_number(meta, lines)

    # Decide oscillation vs rotation from the column-header row keywords.
    # "Mod. Armaz" (Módulo de Armazenamento) ⇒ oscillation.
    # "Taxa de Cis" (Taxa de Cisalhamento)   ⇒ rotation.
    header_text = "\n".join(lines[:12])
    is_oscillation = ("Armazenamento" in header_text
                      or "Storage" in header_text
                      or "G'" in header_text)
    is_rotation = ("Taxa de Cisalhamento" in header_text
                   or "Shear rate" in header_text)

    if is_oscillation and not is_rotation:
        kind = "oscillation"
        data_start = _find_units_row(lines, "[s]") + 1   # time unit anchor
        rows = []
        # Oscillation row layout (post leading-tab strip):
        # [0]=point, [1]=time, [2]=G', [3]=G'', [4]=tan δ,
        # [5]=|η*| (mPa·s), [6]=η_steady (blank in osc), [7]=torque (µN·m),
        # [8]=status
        for point, cells in _iter_data_rows(lines, data_start, n_cols=8):
            try:
                t = _parse_decimal_comma(cells[1])
                gp = _parse_decimal_comma(cells[2])
                gpp = _parse_decimal_comma(cells[3])
                tan_d = _parse_decimal_comma(cells[4])
                eta_c_mPa = _parse_decimal_comma(cells[5])
                eta_steady_mPa = _parse_decimal_comma(cells[6]) if len(cells) > 6 else float("nan")
                torque_uNm = _parse_decimal_comma(cells[7]) if len(cells) > 7 else float("nan")
            except ValueError as e:
                raise ValueError(f"{path.name} (point {point}): {e}") from e
            status = cells[8] if len(cells) > 8 else ""
            rows.append((point, t, gp, gpp, tan_d,
                         eta_c_mPa / 1000.0,         # mPa·s → Pa·s
                         eta_steady_mPa / 1000.0,    # NaN/1000 = NaN; safe
                         torque_uNm, status))
        df = pd.DataFrame(rows, columns=[
            "point", "time_s", "G_storage_Pa", "G_loss_Pa", "tan_delta",
            "eta_complex_Pa_s", "eta_steady_Pa_s", "torque_uNm", "status",
        ])

    elif is_rotation and not is_oscillation:
        kind = "rotation"
        data_start = _find_units_row(lines, "[s]") + 1
        rows = []
        # Rotation row layout:
        # [0]=point, [1]=time, [2]=shear rate, [3]=stress,
        # [4]=viscosity (mPa·s), [5]=status
        for point, cells in _iter_data_rows(lines, data_start, n_cols=5):
            try:
                t = _parse_decimal_comma(cells[1])
                gd = _parse_decimal_comma(cells[2])
                tau = _parse_decimal_comma(cells[3])
                eta_mPa = _parse_decimal_comma(cells[4])
            except ValueError as e:
                raise ValueError(f"{path.name} (point {point}): {e}") from e
            status = cells[5] if len(cells) > 5 else ""
            rows.append((point, t, gd, tau, eta_mPa / 1000.0, status))
        df = pd.DataFrame(rows, columns=[
            "point", "time_s", "shear_rate_1_s", "shear_stress_Pa",
            "viscosity_Pa_s", "status",
        ])

    else:
        raise ValueError(
            f"{path.name}: cannot determine whether interval is oscillation "
            "or rotation. Header did not contain expected keywords."
        )

    if df.empty:
        raise ValueError(f"{path.name}: no numeric data rows parsed.")

    return RecoveryInterval(
        interval_num=interval_num,
        kind=kind,
        sample_name=meta.get("sample_name", path.stem),
        project=meta.get("project", ""),
        interval_label=meta.get("interval_label", ""),
        n_points=meta["n_points"] if isinstance(meta.get("n_points"), int) else len(df),
        source_path=path,
        data=df,
    )


def _recovery_basename(path: Path) -> str | None:
    """Extract the shared basename of a 3iTT set from one of its file paths.

    Returns the basename string (e.g. ``"Carla_Amostra A - CP50 25C_12_05_2026"``)
    or None if the filename doesn't look like a 3iTT interval (i.e. doesn't
    end with ``_1.csv`` / ``_2.csv`` / ``_3.csv``).
    """
    stem = path.stem
    for n in ("_1", "_2", "_3"):
        if stem.endswith(n):
            return stem[: -len(n)]
    return None


def list_recovery_basenames(folder: str | Path) -> list[str]:
    """Return the distinct basenames of every complete 3iTT set in a folder.

    A "complete" set has all three of ``_1.csv``, ``_2.csv``, ``_3.csv``
    present. Incomplete sets are omitted (logged via ``print`` so the user
    notices). Analysis siblings are ignored.
    """
    folder = Path(folder)
    files = list_raw_csvs(folder)
    seen: dict[str, set[int]] = {}
    for p in files:
        base = _recovery_basename(p)
        if base is None:
            continue
        try:
            n = int(p.stem[-1])
        except ValueError:
            continue
        seen.setdefault(base, set()).add(n)
    complete = []
    for base, ns in sorted(seen.items()):
        if ns == {1, 2, 3}:
            complete.append(base)
        else:
            print(f"  [!] incomplete 3iTT set '{base}': intervals present = {sorted(ns)}")
    return complete


def read_recovery_set(folder: str | Path, basename: str) -> RecoverySet:
    """Load all three intervals of a 3iTT set sharing the same basename."""
    folder = Path(folder)
    paths = {n: folder / f"{basename}_{n}.csv" for n in (1, 2, 3)}
    intervals: dict[int, RecoveryInterval] = {}
    for n, p in paths.items():
        if not p.exists():
            raise FileNotFoundError(f"Missing interval {n} for '{basename}': {p}")
        intervals[n] = read_recovery_interval_csv(p)

    # Sanity: the three files should agree on sample_name / project label.
    sample_names = {iv.sample_name for iv in intervals.values()}
    if len(sample_names) > 1:
        # Rheocompass occasionally re-saves with slight metadata drift; warn,
        # don't crash — just take the interval-1 label as canonical.
        print(f"  [!] sample_name disagrees across intervals of '{basename}': {sample_names}")

    return RecoverySet(
        sample_name=intervals[1].sample_name,
        project=intervals[1].project,
        interval_label=intervals[1].interval_label,
        source_basename=basename,
        interval1=intervals[1],
        interval2=intervals[2],
        interval3=intervals[3],
    )


def read_recovery_sets_in_folder(folder: str | Path) -> dict[str, RecoverySet]:
    """Parse every complete 3iTT set in a folder, keyed by short sample name."""
    folder = Path(folder)
    out: dict[str, RecoverySet] = {}
    for base in list_recovery_basenames(folder):
        rs = read_recovery_set(folder, base)
        short = rs.sample_name.split(" - ", 1)[0].strip() or base
        # If two basenames map to the same short name (e.g. PC1 vs "PC1 1"),
        # append the trailing distinguisher so neither is lost.
        if short in out:
            # Fall back to the full basename for the duplicate.
            short = base
        out[short] = rs
    return out


# ----------------------------------------------------------------------
# Folder discovery
# ----------------------------------------------------------------------

def list_raw_csvs(folder: str | Path) -> list[Path]:
    """
    Return every raw Anton Paar CSV in ``folder``.

    Filters out:
        * ``*_Analysis.csv``  (English) — Rheocompass regression export
        * ``*_Analise.csv``   (Portuguese) — same, localised filename
        * Anything not ending in ``.csv`` (e.g. report PDFs)
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise NotADirectoryError(f"Not a directory: {folder}")
    out = []
    for p in sorted(folder.iterdir()):
        if not p.is_file():
            continue
        if p.suffix.lower() != ".csv":
            continue
        if p.name.lower().endswith(_ANALYSIS_SUFFIXES):
            continue
        out.append(p)
    return out


# Back-compatible aliases — all call the same generic filter.
list_flow_curves = list_raw_csvs
list_amplitude_sweeps = list_raw_csvs
list_frequency_sweeps = list_raw_csvs
# (Recovery has its own list_recovery_basenames defined above — that one
# groups the per-interval files into sets, so it's not a simple alias.)


def _short_name_from_metadata(fc_or_as) -> str:
    """Sample short-name = part of ``Teste:`` header before the first ' - '."""
    full = fc_or_as.sample_name or fc_or_as.source_path.stem
    short = full.split(" - ", 1)[0].strip()
    return short or full


def read_flow_curves_in_folder(folder: str | Path) -> dict[str, FlowCurve]:
    """Parse every raw flow-curve CSV in a folder. Keyed by short sample name."""
    out: dict[str, FlowCurve] = {}
    for p in list_raw_csvs(folder):
        fc = read_flow_curve_csv(p)
        out[_short_name_from_metadata(fc)] = fc
    return out


def read_amplitude_sweeps_in_folder(folder: str | Path) -> dict[str, AmplitudeSweep]:
    """Parse every raw amplitude-sweep CSV in a folder. Keyed by short sample name."""
    out: dict[str, AmplitudeSweep] = {}
    for p in list_raw_csvs(folder):
        sweep = read_amplitude_sweep_csv(p)
        out[_short_name_from_metadata(sweep)] = sweep
    return out


def read_frequency_sweeps_in_folder(folder: str | Path) -> dict[str, FrequencySweep]:
    """Parse every raw frequency-sweep CSV in a folder. Keyed by short sample name."""
    out: dict[str, FrequencySweep] = {}
    for p in list_raw_csvs(folder):
        sweep = read_frequency_sweep_csv(p)
        out[_short_name_from_metadata(sweep)] = sweep
    return out


# ----------------------------------------------------------------------
# CLI smoke test:  python antpar_io.py <folder-or-csv> [--amp]
# ----------------------------------------------------------------------

def _cli(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: python antpar_io.py <folder-or-csv> [--amp | --freq | --rec]")
        print("  --amp    parse as amplitude (strain) sweep")
        print("  --freq   parse as frequency sweep")
        print("  --rec    parse 3iTT recovery sets (folder must contain _1/_2/_3 trios)")
        print("  (default: flow / viscosity curve)")
        return 2

    if "--amp" in argv:
        kind, reader, summary = "amplitude", read_amplitude_sweep_csv, _amp_summary
        parse_folder = read_amplitude_sweeps_in_folder
    elif "--freq" in argv:
        kind, reader, summary = "frequency", read_frequency_sweep_csv, _freq_summary
        parse_folder = read_frequency_sweeps_in_folder
    elif "--rec" in argv:
        kind, reader, summary = "recovery", read_recovery_interval_csv, _rec_summary
        parse_folder = read_recovery_sets_in_folder
    else:
        kind, reader, summary = "flow", read_flow_curve_csv, _flow_summary
        parse_folder = read_flow_curves_in_folder

    args = [a for a in argv[1:] if not a.startswith("--")]
    target = Path(args[0])

    if target.is_dir():
        curves = parse_folder(target)
        print(f"Found {len(curves)} raw {kind}-sweep CSVs in {target}:")
        for short, c in curves.items():
            print("  " + summary(short, c))
        return 0

    c = reader(target)
    print(f"Project        : {c.project}")
    print(f"Sample (Teste:): {c.sample_name}")
    print(f"Interval label : {c.interval_label}")
    print(f"Declared rows  : {c.n_points}")
    print(f"Parsed rows    : {len(c.data)}")
    print()
    print(c.data.head(5).to_string(index=False))
    print("...")
    print(c.data.tail(3).to_string(index=False))
    return 0


def _flow_summary(short: str, fc: FlowCurve) -> str:
    df = fc.data
    return (f"{short:<14} "
            f"n={len(df):3d}  "
            f"γ̇=[{df['shear_rate_1_s'].min():.3g}, "
            f"{df['shear_rate_1_s'].max():.3g}] 1/s  "
            f"η=[{df['viscosity_Pa_s'].min():.3g}, "
            f"{df['viscosity_Pa_s'].max():.3g}] Pa·s  "
            f"T̄={df['temperature_C'].mean():.2f}°C")


def _amp_summary(short: str, sweep: AmplitudeSweep) -> str:
    df = sweep.data
    return (f"{short:<14} "
            f"n={len(df):3d}  "
            f"γ=[{df['strain_pct'].min():.3g}, "
            f"{df['strain_pct'].max():.3g}] %  "
            f"G'=[{df['G_storage_Pa'].min():.3g}, "
            f"{df['G_storage_Pa'].max():.3g}] Pa  "
            f"tan δ ≈ {df['tan_delta'].median():.2f}")


def _freq_summary(short: str, sweep: FrequencySweep) -> str:
    df = sweep.data
    return (f"{short:<14} "
            f"n={len(df):3d}  "
            f"ω=[{df['omega_rad_s'].min():.3g}, "
            f"{df['omega_rad_s'].max():.3g}] rad/s  "
            f"G'=[{df['G_storage_Pa'].min():.3g}, "
            f"{df['G_storage_Pa'].max():.3g}] Pa  "
            f"γ_set ≈ {df['strain_pct'].median():.2f}%")


def _rec_summary(short: str, rs: RecoverySet) -> str:
    iv1, iv2, iv3 = rs.interval1, rs.interval2, rs.interval3
    return (f"{short:<14} "
            f"int1({iv1.kind[:3]})={iv1.n_points}  "
            f"int2({iv2.kind[:3]})={iv2.n_points}  "
            f"int3({iv3.kind[:3]})={iv3.n_points}  "
            f"label='{rs.interval_label}'")


if __name__ == "__main__":
    import sys
    sys.exit(_cli(sys.argv))
