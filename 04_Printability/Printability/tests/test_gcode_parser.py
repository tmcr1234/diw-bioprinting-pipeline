from pathlib import Path
from gcode_parser import parse_gcode, parse_filename, RasterSpec

REAL_GCODE = Path("/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE/Impressao 3D/S-Test/Vp10-Flow1.gcode")


def test_parse_filename_extracts_Vp_and_Flow():
    Vp, Flow = parse_filename(Path("S-Test_Vp10_Flow1,57_SF5.5_Black_flash.heic"))
    assert Vp == 10.0
    assert Flow == 1.57


def test_parse_filename_handles_Fr_token():
    Vp, Flow = parse_filename(Path("Vp7-Fr1,5.gcode"))
    assert Vp == 7.0
    assert Flow == 1.5


def test_parse_real_gcode_returns_RasterSpec():
    if not REAL_GCODE.exists():
        import pytest
        pytest.skip("Real G-code not available in this env")
    spec = parse_gcode(REAL_GCODE)
    assert isinstance(spec, RasterSpec)
    assert len(spec.strands) >= 2
    assert spec.strand_spacing_mm > 0
    assert spec.Vp_mm_s == 10.0
    assert spec.Flow == 1.0


def test_strands_have_constant_Y(tmp_path):
    g = tmp_path / "tiny.gcode"
    g.write_text(
        "G1 X0 Y0 E0\n"
        "G1 X10 Y0 E1\n"
        "G1 X10 Y1 E1.1\n"
        "G1 X0 Y1 E2.1\n"
    )
    spec = parse_gcode(g)
    assert len(spec.strands) == 2
    assert abs(spec.strand_spacing_mm - 1.0) < 1e-6
