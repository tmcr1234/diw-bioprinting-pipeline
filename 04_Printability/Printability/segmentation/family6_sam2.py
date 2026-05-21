"""Method 14: SAM2 Hiera-Large via the SAM2Service.

When `service` is None (operator passed --no-sam2 or SAM2 failed to load),
return an empty mask so downstream consensus simply has one fewer voter.
"""
from __future__ import annotations
import numpy as np


def mask_sam2(img: np.ndarray, bg: str, *, service=None) -> np.ndarray:
    if service is None:
        return np.zeros(img.shape[:2], dtype=bool)
    return service.predict(img)
