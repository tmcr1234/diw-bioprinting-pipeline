import numpy as np
from PIL import Image
import pillow_heif
pillow_heif.register_heif_opener()
from background_model import build_background_model
from segmentation.family5_bgdiff import mask_bgdiff
from tests.fixtures.synthetic import make_synthetic


def _write_ref(tmp_path):
    ref = tmp_path / "_reference"
    ref.mkdir()
    bg, _ = make_synthetic("Black", H=200, W=200)
    bg[50:150, 50:150] = [3, 3, 3]
    for name in ["bg_pre_flash", "bg_pre_noflash", "bg_post_flash", "bg_post_noflash"]:
        Image.fromarray(bg).save(ref / f"{name}.heic", format="HEIF")
    return ref


def test_bgdiff_recovers_foreground(tmp_path):
    ref = _write_ref(tmp_path)
    model = build_background_model(ref)
    test_img, gt = make_synthetic("Black", H=200, W=200, seed=42)
    mask = mask_bgdiff(flash=test_img, noflash=test_img, bg="Black", model=model)
    iou = (mask & gt).sum() / (mask | gt).sum()
    assert iou > 0.80
