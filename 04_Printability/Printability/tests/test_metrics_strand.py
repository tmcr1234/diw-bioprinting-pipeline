import numpy as np
from metrics_strand import compute_strand_report
from gcode_parser import RasterSpec, LineSegment
from registration import Alignment


def test_strand_widths_match_synthetic():
    H, W = 200, 200
    mask = np.zeros((H, W), dtype=bool)
    mask[95:105, 20:180] = True
    mask[125:135, 20:180] = True
    spec = RasterSpec(
        strands=[LineSegment(0, 0, 16, 0), LineSegment(0, 3, 16, 3)],
        turns=[],
        strand_spacing_mm=3.0,
        strand_length_mm=16.0,
        bounding_box_mm=(0, 0, 16, 3),
        Vp_mm_s=10.0, Flow=1.0,
        source_gcode=None,
    )
    affine = np.array([[10.0, 0, 20], [0, 10.0, 100]])
    alignment = Alignment(affine, iou=1.0, success=True)
    rep = compute_strand_report(mask, spec, alignment)
    assert len(rep.widths_px) == 2
    assert all(abs(w - 10) < 2 for w in rep.widths_px), f"widths={rep.widths_px}"
