import numpy as np
from registration import align_gcode_to_image
from gcode_parser import parse_gcode


def _make_synthetic_spec_and_mask(tmp_path):
    g = tmp_path / "g.gcode"
    g.write_text(
        "G1 X0 Y0 E0\n" "G1 X20 Y0 E1\n"
        "G1 X20 Y2 E1.1\n" "G1 X0 Y2 E2.1\n"
        "G1 X0 Y4 E2.2\n" "G1 X20 Y4 E3.2\n"
    )
    spec = parse_gcode(g)
    H, W = 200, 200
    mask = np.zeros((H, W), dtype=bool)
    mask[80:120, 40:160] = True
    return spec, mask


def test_alignment_succeeds_for_axis_aligned_specimen(tmp_path):
    spec, mask = _make_synthetic_spec_and_mask(tmp_path)
    align = align_gcode_to_image(spec, mask)
    assert align.success
    assert align.iou > 0.5
    assert align.affine.shape == (2, 3)
