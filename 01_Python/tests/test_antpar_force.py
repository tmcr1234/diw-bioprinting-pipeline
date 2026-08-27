"""
test_antpar_force.py
====================

Acceptance tests for the force-augmented ingestion layer in ``antpar_io``
(items A4, A5, B16, B17, B19 of the 2026-08-26 pipeline audit).

WHY THESE ARE SYNTHETIC
-----------------------
The audit's own acceptance criteria are stated against the C10-C25 force
exports ("reproduces the five audited floor values exactly", "~5e-6 relative
difference vs. the audited .xls"). Those files are not present in every
checkout of this pipeline, so those criteria cannot be closed here. What
CAN be closed here is that the logic is right: the conversion is analytic,
the tare convention is the one that was verified correct, the rejected
convention is actively blocked, and the encoding trap raises instead of
mojibaking.

WHEN THE REAL DATA IS AVAILABLE, add a test that reads one real export and
asserts the floor against the lab notebook. Until then, do not describe
this reader as validated — it is *specified*.

Run:  python3 -m pytest tests/test_antpar_force.py -q
      (or: python3 tests/test_antpar_force.py  for a dependency-free run)
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import antpar_io as aio  # noqa: E402


# ----------------------------------------------------------------------
# B19 — N1 = 2F/(pi R^2)
# ----------------------------------------------------------------------

def test_normal_stress_from_force_analytic():
    # A 50 mm fixture is 25 mm RADIUS. Getting this wrong is a 4x error.
    got = aio.normal_stress_from_force(1.0, 25e-3)
    want = 2.0 / (math.pi * 25e-3 ** 2)
    assert abs(got - want) / want < 1e-12, (got, want)
    assert abs(got - 1018.5916) < 1e-3, got


def test_normal_stress_scales_linearly_and_inverse_square():
    a = aio.normal_stress_from_force(1.0, 25e-3)
    assert abs(aio.normal_stress_from_force(3.0, 25e-3) - 3 * a) < 1e-9
    # Halving the radius quadruples the stress.
    assert abs(aio.normal_stress_from_force(1.0, 12.5e-3) - 4 * a) < 1e-9


def test_normal_stress_rejects_bad_radius():
    for bad in (0.0, -1e-3):
        try:
            aio.normal_stress_from_force(1.0, bad)
        except ValueError:
            continue
        raise AssertionError(f"R={bad} should have raised")


# ----------------------------------------------------------------------
# A5 — baseline / tare
# ----------------------------------------------------------------------

def test_baseline_min_recovers_a_known_offset():
    true_signal = np.array([0.0, 12.0, 55.0, 210.0, 802.0, 640.0])
    for offset in (-350.0, 0.0, 1250.0):
        corrected, floor = aio.baseline_correct(true_signal + offset, method="min")
        assert abs(floor - offset) < 1e-12, (floor, offset)
        assert np.allclose(corrected, true_signal), corrected


def test_rejected_median_convention_is_blocked_and_explains_itself():
    try:
        aio.baseline_correct(np.arange(100.0), method="median_low25")
    except ValueError as e:
        msg = str(e)
        assert "REJECTED" in msg
        assert "802" in msg and "460" in msg, "the worked numbers must survive"
        return
    raise AssertionError("median_low25 must raise, not silently compute")


def test_why_median_low25_is_wrong_on_a_ramp_shape():
    # The rejected convention fails because the low tail of a shear ramp is
    # REAL low-shear data, not baseline. Reproduce that shape and show the
    # two conventions disagree materially, which is the whole lesson.
    offset = 500.0
    ramp = np.concatenate([np.linspace(0.0, 400.0, 40), np.linspace(400.0, 2500.0, 60)])
    raw = ramp + offset
    _, floor_min = aio.baseline_correct(raw, method="min")
    floor_median_low25 = float(np.median(np.sort(raw)[:25]))
    assert abs(floor_min - offset) < 1e-9
    # The rejected convention lands well above the true offset.
    assert floor_median_low25 - offset > 100.0, floor_median_low25
    assert (floor_median_low25 - floor_min) / floor_min > 0.2


def test_baseline_ignores_nans_but_needs_some_finite_data():
    corrected, floor = aio.baseline_correct([np.nan, 10.0, 4.0, np.nan], method="min")
    assert floor == 4.0
    assert np.isnan(corrected[0]) and corrected[2] == 0.0
    try:
        aio.baseline_correct([np.nan, np.nan], method="min")
    except ValueError:
        return
    raise AssertionError("all-NaN series must raise")


# ----------------------------------------------------------------------
# B16 — the encoding trap
# ----------------------------------------------------------------------

_FORCE_FIXTURE = (
    "Filename\tSYNTH_C20_flow_force.txt\r\n"
    "Instrument name\tAresG2\r\n"
    "Geometry\t50mm cone\r\n"
    "\r\n"
    "Shear rate\tShear stress\tViscosity\tAxial force\tStep time\tTemperature\r\n"
    "1/s\tPa\tPa.s\tN\ts\t°C\r\n"
    "0,1000\t120,50\t1205,0\t0,4000\t10,0\t25,0\r\n"
    "1,0000\t260,00\t260,00\t0,6000\t20,0\t25,0\r\n"
    "10,000\t420,00\t42,000\t1,2000\t30,0\t25,0\r\n"
    "100,00\t690,00\t6,9000\t2,0000\t40,0\t25,0\r\n"
)


def _write_fixture(tmpdir: Path, text: str = _FORCE_FIXTURE) -> Path:
    f = tmpdir / "SYNTH_C20_flow_force.txt"
    f.write_bytes(text.encode("latin-1"))
    return f


def test_utf16_input_raises_instead_of_mojibaking(tmp_path=None):
    tmp = Path(tmp_path) if tmp_path else Path("./_t_utf16"); tmp.mkdir(exist_ok=True)
    f = tmp / "wrong_encoding.txt"
    f.write_bytes(b"\xff\xfe" + "Shear rate\tPa\n".encode("utf-16-le"))
    try:
        aio.read_flow_curve_with_force(f)
    except ValueError as e:
        assert "UTF-16" in str(e)
        return
    finally:
        f.unlink(missing_ok=True)
        if not tmp_path:
            tmp.rmdir()
    raise AssertionError("a UTF-16 file must raise, not decode as latin-1")


# ----------------------------------------------------------------------
# A4 — one call returns N1_Pa alongside the usual columns
# ----------------------------------------------------------------------

def test_read_flow_curve_with_force_end_to_end(tmp_path=None):
    tmp = Path(tmp_path) if tmp_path else Path("./_t_force"); tmp.mkdir(exist_ok=True)
    f = _write_fixture(tmp)
    try:
        fc = aio.read_flow_curve_with_force(f, geometry="CP50")

        assert fc.n_points == 4, fc.n_points
        assert "N1_Pa" in fc.data.columns
        assert fc.radius_m == 25e-3

        # B20: the instrument identity must survive into the object.
        assert fc.instrument == "AresG2", fc.instrument

        # Decimal commas parsed, shear rate spans four decades.
        np.testing.assert_allclose(
            fc.data["shear_rate_1_s"].to_numpy(), [0.1, 1.0, 10.0, 100.0])
        np.testing.assert_allclose(
            fc.data["shear_stress_Pa"].to_numpy(), [120.5, 260.0, 420.0, 690.0])

        # N1 from force, then tared on the minimum: the first point is the
        # quietest, so it must come out at exactly zero.
        raw = aio.normal_stress_from_force(
            np.array([0.4, 0.6, 1.2, 2.0]), 25e-3)
        want = raw - raw.min()
        np.testing.assert_allclose(fc.data["N1_Pa"].to_numpy(), want, rtol=1e-12)
        assert fc.data["N1_Pa"].iloc[0] == 0.0
        assert abs(fc.baseline_floor_Pa - raw.min()) < 1e-9
    finally:
        f.unlink(missing_ok=True)
        if not tmp_path:
            tmp.rmdir()


def test_missing_force_channel_raises_loudly(tmp_path=None):
    tmp = Path(tmp_path) if tmp_path else Path("./_t_noforce"); tmp.mkdir(exist_ok=True)
    text = (
        "Instrument name\tAresG2\r\n\r\n"
        "Shear rate\tShear stress\tViscosity\r\n"
        "1/s\tPa\tPa.s\r\n"
        "1,0\t260,0\t260,0\r\n"
    )
    f = _write_fixture(tmp, text)
    try:
        aio.read_flow_curve_with_force(f)
    except ValueError as e:
        assert "axial force" in str(e).lower()
        return
    finally:
        f.unlink(missing_ok=True)
        if not tmp_path:
            tmp.rmdir()
    raise AssertionError("a force-less export must raise, not return NaN N1")


# ----------------------------------------------------------------------
# B17 — the filename pairing traps
# ----------------------------------------------------------------------

def test_pair_audited_files_carries_each_known_trap():
    c25 = aio.pair_audited_files("C25")
    assert "double space" in c25["note"].lower()

    gel = aio.pair_audited_files("gel")
    assert gel["xls_basename"] == "gel (Ac)"
    assert gel["force_txt_basename"].endswith(" "), "trailing space must survive"

    c10 = aio.pair_audited_files("C10")
    assert "(1)" in c10["xls_basename"]


def test_unknown_ink_raises_with_the_known_list():
    try:
        aio.pair_audited_files("C99")
    except KeyError as e:
        assert "C25" in str(e) and "do NOT glob" in str(e)
        return
    raise AssertionError("an unregistered ink must raise, never fuzzy-match")


# ----------------------------------------------------------------------

if __name__ == "__main__":
    import tempfile
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        with tempfile.TemporaryDirectory() as td:
            try:
                if "tmp_path" in fn.__code__.co_varnames:
                    fn(td)
                else:
                    fn()
                print(f"  [ok]   {fn.__name__}")
            except Exception as e:                        # noqa: BLE001
                print(f"  [FAIL] {fn.__name__}: {e}")
                failed += 1
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
