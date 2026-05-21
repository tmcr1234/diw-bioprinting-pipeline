"""Build a synthetic session under tests/fixtures/session/ for end-to-end tests.

Writes:
- session/Impressao 3D/S-Test/<gcode>.gcode
- session/Impressao 3D/S-Test/SF5.5/Black/_reference/*.heic (4)
- session/Impressao 3D/S-Test/SF5.5/Black/*.heic (2 pairs)
"""
from pathlib import Path
import numpy as np
from PIL import Image
import pillow_heif

pillow_heif.register_heif_opener()

ROOT = Path(__file__).parent / "session"


def main():
    photo = ROOT / "Impressao 3D" / "S-Test" / "SF5.5" / "Black"
    photo.mkdir(parents=True, exist_ok=True)
    gcode_dir = ROOT / "Impressao 3D" / "S-Test"
    gcode_dir.mkdir(parents=True, exist_ok=True)

    g = gcode_dir / "S-Test_Vp10_Flow1.gcode"
    g.write_text("\n".join([
        "G1 X10 Y10 E0",
        "G1 X30 Y10 E1",
        "G1 X30 Y13 E1.1",
        "G1 X10 Y13 E2.1",
        "G1 X10 Y16 E2.2",
        "G1 X30 Y16 E3.2",
        "G1 X30 Y19 E3.3",
        "G1 X10 Y19 E4.3",
    ]))

    rng = np.random.default_rng(0)
    ref_dir = photo / "_reference"
    ref_dir.mkdir(exist_ok=True)
    base = rng.integers(0, 5, (200, 200, 3), dtype=np.uint8)
    for name in ["bg_pre_flash", "bg_pre_noflash", "bg_post_flash", "bg_post_noflash"]:
        Image.fromarray(base).save(ref_dir / f"{name}.heic", format="HEIF")

    for flow_label in ["1", "1,57"]:
        img = base.copy()
        for y in [40, 60, 80, 100]:
            img[y:y + 4, 40:160] = [230, 230, 230]
        flash_p = photo / f"S-Test_Vp10_Flow{flow_label}_SF5.5_Black_flash.heic"
        noflash_p = photo / f"S-Test_Vp10_Flow{flow_label}_SF5.5_Black_noflash.heic"
        Image.fromarray(img).save(flash_p, format="HEIF")
        Image.fromarray(img).save(noflash_p, format="HEIF")


if __name__ == "__main__":
    main()
