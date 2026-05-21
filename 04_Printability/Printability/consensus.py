"""Per-pixel majority vote across segmentation methods → consensus mask."""
from __future__ import annotations
import math
import numpy as np
from scipy.ndimage import label


def is_degenerate(mask: np.ndarray) -> bool:
    fg = mask.sum()
    if fg == 0 or fg == mask.size:
        return True
    _, n_components = label(mask)
    return n_components > 100


def compute_consensus(masks: dict[str, np.ndarray]) -> tuple[np.ndarray, float]:
    good = [m for m in masks.values() if not is_degenerate(m)]
    if not good:
        H, W = next(iter(masks.values())).shape
        return np.zeros((H, W), dtype=bool), 0.0
    N = len(good)
    threshold = math.ceil(N / 2)
    votes = np.sum(np.stack(good, axis=0), axis=0)
    consensus = votes >= threshold
    if consensus.sum() == 0:
        return consensus, 0.0
    agreement = float(votes[consensus].mean() / N)
    return consensus, agreement
