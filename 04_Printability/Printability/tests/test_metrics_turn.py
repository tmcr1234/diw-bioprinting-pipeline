import numpy as np
from metrics_turn import compute_turn_report
from gcode_parser import RasterSpec, LineSegment, TurnPoint
from registration import Alignment


def test_clean_turn_counts():
    H, W = 200, 200
    mask = np.zeros((H, W), dtype=bool)
    mask[80:90, 50:150] = True
    mask[110:120, 50:150] = True
    mask[80:120, 145:155] = True
    spec = RasterSpec(
        strands=[LineSegment(0, 0, 10, 0), LineSegment(0, 3, 10, 3)],
        turns=[TurnPoint(10, 1.5, 0, 1)],
        strand_spacing_mm=3.0,
        strand_length_mm=10.0,
        bounding_box_mm=(0, 0, 10, 3),
        Vp_mm_s=10.0, Flow=1.0, source_gcode=None,
    )
    affine = np.array([[10.0, 0, 50], [0, 10.0, 80]])
    alignment = Alignment(affine, iou=1.0, success=True)
    rep = compute_turn_report(mask, spec, alignment, median_strand_width_px=10.0)
    assert rep.n_total == 1
    assert rep.fidelity >= 0.0
