"""Family 3: local / adaptive threshold methods.

All methods finalize by keeping the largest connected component to suppress
high-frequency local-threshold noise outside the specimen.
"""
from __future__ import annotations
import numpy as np
import cv2
from scipy.ndimage import label
from skimage.filters import threshold_sauvola, threshold_niblack


def _gray(img: np.ndarray) -> np.ndarray:
    return cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)


def _largest_component(mask: np.ndarray) -> np.ndarray:
    labels, n = label(mask)
    if n <= 1:
        return mask
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    return labels == counts.argmax()


def mask_adaptive_gauss(img, bg):
    g = _gray(img)
    raw = cv2.adaptiveThreshold(
        g, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 51, -10
    ).astype(bool)
    return _largest_component(raw)


def mask_sauvola(img, bg):
    g = _gray(img).astype(np.float32) / 255.0
    t = threshold_sauvola(g, window_size=51, k=0.2)
    return _largest_component(g > t)


def mask_niblack(img, bg):
    g = _gray(img).astype(np.float32) / 255.0
    t = threshold_niblack(g, window_size=51, k=-0.2)
    return _largest_component(g > t)
