import numpy as np
from quality import compute_quality
from background_model import BackgroundModel


def _dummy_bg(H, W):
    return BackgroundModel(
        mu=np.zeros((H, W, 3), dtype=np.float32),
        sigma=np.ones((H, W, 3), dtype=np.float32),
        ghost_mask=np.zeros((H, W), dtype=bool),
        drift_score=0.0,
    )


def test_quality_flags_clipped():
    H, W = 100, 100
    img = np.full((H, W, 3), 255, dtype=np.uint8)
    fg = np.zeros((H, W), dtype=bool)
    fg[40:60, 40:60] = True
    rep = compute_quality(img, img, fg, _dummy_bg(H, W))
    assert rep.flags["clipped"]


def test_quality_flags_blurry():
    H, W = 200, 200
    img = np.full((H, W, 3), 128, dtype=np.uint8)
    fg = np.ones((H, W), dtype=bool)
    rep = compute_quality(img, img, fg, _dummy_bg(H, W))
    assert rep.flags["blurry"]


def test_quality_does_not_flag_clean_image():
    rng = np.random.default_rng(0)
    H, W = 200, 200
    img = rng.integers(80, 180, (H, W, 3), dtype=np.uint8)
    fg = np.zeros((H, W), dtype=bool)
    fg[50:150, 50:150] = True
    img[fg] = 230
    rep = compute_quality(img, img, fg, _dummy_bg(H, W))
    assert not rep.flags["clipped"]
    assert not rep.flags["blurry"]
