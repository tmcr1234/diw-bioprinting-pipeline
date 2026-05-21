"""Family 2: global threshold methods on grayscale.

All methods finalize by keeping only the largest connected component to
suppress salt-and-pepper background noise that the global threshold can't
separate from the specimen.
"""
from __future__ import annotations
import numpy as np
import cv2
from scipy.ndimage import label
from skimage.filters import threshold_otsu, threshold_triangle, threshold_li, threshold_yen


def _gray(img: np.ndarray) -> np.ndarray:
    return cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)


def _largest_component(mask: np.ndarray) -> np.ndarray:
    labels, n = label(mask)
    if n <= 1:
        return mask
    counts = np.bincount(labels.ravel())
    counts[0] = 0  # ignore background label
    return labels == counts.argmax()


def _apply(thr_fn, img, bg):
    g = _gray(img)
    t = thr_fn(g)
    return _largest_component(g > t)


def mask_otsu(img, bg):
    return _apply(threshold_otsu, img, bg)


def mask_triangle(img, bg):
    return _apply(threshold_triangle, img, bg)


def mask_li(img, bg):
    return _apply(threshold_li, img, bg)


def mask_yen(img, bg):
    return _apply(threshold_yen, img, bg)
