from pathlib import Path
import numpy as np
import pytest
from io_heic import load_heic, load_pair

FIXTURES = Path(__file__).parent / "fixtures"


def test_load_heic_returns_uint8_rgb():
    arr = load_heic(FIXTURES / "tiny.heic")
    assert arr.dtype == np.uint8
    assert arr.ndim == 3
    assert arr.shape[2] == 3
    assert arr.shape[:2] == (64, 64)


def test_load_pair_returns_two_arrays():
    flash, noflash = load_pair(FIXTURES / "tiny.heic", FIXTURES / "tiny.heic")
    assert flash.shape == noflash.shape


def test_load_pair_rejects_mismatched_shapes(tmp_path):
    from PIL import Image
    import pillow_heif
    pillow_heif.register_heif_opener()
    Image.fromarray(np.zeros((32, 32, 3), dtype=np.uint8)).save(tmp_path / "small.heic", format="HEIF")
    with pytest.raises(ValueError, match="shape"):
        load_pair(FIXTURES / "tiny.heic", tmp_path / "small.heic")
