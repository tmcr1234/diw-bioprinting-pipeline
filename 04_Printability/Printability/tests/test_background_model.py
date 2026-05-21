from pathlib import Path
import numpy as np
from PIL import Image
import pillow_heif
pillow_heif.register_heif_opener()
from background_model import build_background_model, save_model, load_model


def _write_ref(tmp_path: Path) -> Path:
    ref = tmp_path / "_reference"
    ref.mkdir()
    rng = np.random.default_rng(0)
    base = rng.integers(0, 20, (64, 64, 3), dtype=np.uint8)
    for name in ["bg_pre_flash", "bg_pre_noflash", "bg_post_flash", "bg_post_noflash"]:
        noisy = base + rng.integers(-2, 3, base.shape).astype(np.int8)
        noisy = np.clip(noisy.astype(int), 0, 255).astype(np.uint8)
        Image.fromarray(noisy).save(ref / f"{name}.heic", format="HEIF")
    return ref


def test_build_model_returns_correct_shapes(tmp_path):
    ref = _write_ref(tmp_path)
    model = build_background_model(ref)
    assert model.mu.shape == (64, 64, 3)
    assert model.sigma.shape == (64, 64, 3)
    assert model.ghost_mask.shape == (64, 64)
    assert model.ghost_mask.dtype == bool


def test_model_round_trip(tmp_path):
    ref = _write_ref(tmp_path)
    model = build_background_model(ref)
    save_model(model, tmp_path / "bg.npz")
    loaded = load_model(tmp_path / "bg.npz")
    np.testing.assert_array_equal(model.mu, loaded.mu)
    np.testing.assert_array_equal(model.ghost_mask, loaded.ghost_mask)
