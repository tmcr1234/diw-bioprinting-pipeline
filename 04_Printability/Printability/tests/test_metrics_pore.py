import numpy as np
from metrics_pore import compute_pore_report, ideal_pr_from_spec
from gcode_parser import RasterSpec, LineSegment
from registration import Alignment


def test_pr_of_rectangular_pore_is_reasonable():
    H, W = 200, 200
    mask = np.zeros((H, W), dtype=bool)
    mask[40:50, 20:180] = True
    mask[150:160, 20:180] = True
    spec = RasterSpec(
        strands=[LineSegment(0, 0, 16, 0), LineSegment(0, 12, 16, 12)],
        turns=[],
        strand_spacing_mm=12.0,
        strand_length_mm=16.0,
        bounding_box_mm=(0, 0, 16, 12),
        Vp_mm_s=10.0, Flow=1.0,
        source_gcode=None,
    )
    affine = np.array([[10.0, 0, 20], [0, 10.0, 40]])
    alignment = Alignment(affine, iou=1.0, success=True)
    rep = compute_pore_report(mask, spec, alignment)
    assert len(rep.Pr_values) == 1
    # Pore ≈ 100×160 rect (rows 50-150, cols 20-180): Pr = (520)²/(16·16000) ≈ 1.06
    assert 0.9 < rep.Pr_values[0] < 1.3, f"Pr={rep.Pr_values[0]}"


def test_ideal_pr_of_square_is_one():
    spec = RasterSpec(
        strands=[], turns=[],
        strand_spacing_mm=10.0, strand_length_mm=10.0,
        bounding_box_mm=(0, 0, 10, 10),
        Vp_mm_s=10.0, Flow=1.0, source_gcode=None,
    )
    assert abs(ideal_pr_from_spec(spec) - 1.0) < 1e-9
