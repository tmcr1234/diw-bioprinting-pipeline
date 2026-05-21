"""Build synthetic foreground-on-coloured-background images for unit tests."""
import numpy as np


def make_synthetic(bg: str, H: int = 200, W: int = 200, seed: int = 0):
    """Returns (img, ground_truth_mask). Foreground = a centred 100×100 square,
    coloured white. Background = nearly black or nearly blue with noise."""
    rng = np.random.default_rng(seed)
    img = np.zeros((H, W, 3), dtype=np.uint8)
    if bg == "Black":
        img[:] = rng.integers(0, 5, (H, W, 3))
    elif bg == "Blue":
        img[..., 2] = rng.integers(160, 180, (H, W))
        img[..., 0:2] = rng.integers(0, 15, (H, W, 2))
    else:
        raise ValueError(bg)
    mask = np.zeros((H, W), dtype=bool)
    mask[50:150, 50:150] = True
    img[mask] = [240, 240, 240]
    return img, mask
