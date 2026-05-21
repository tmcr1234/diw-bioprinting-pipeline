"""Method 15: Mahalanobis distance vs the per-session background model."""
from __future__ import annotations
import numpy as np
from background_model import BackgroundModel


def _maha_mask(img: np.ndarray, model: BackgroundModel, k: float) -> np.ndarray:
    d = np.abs(img.astype(np.float32) - model.mu) / model.sigma
    d_max = d.max(axis=2)
    fg = d_max > k
    return fg & (~model.ghost_mask)


def mask_bgdiff(flash, noflash, bg, *, model: BackgroundModel, k: float = 3.0):
    return _maha_mask(flash, model, k) | _maha_mask(noflash, model, k)
