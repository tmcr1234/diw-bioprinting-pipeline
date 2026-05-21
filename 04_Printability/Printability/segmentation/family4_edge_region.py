"""Family 4: edge and region segmentation methods."""
from __future__ import annotations
import numpy as np
import cv2
from scipy.ndimage import label
from skimage.segmentation import felzenszwalb


def _gray(img):
    return cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)


def _largest_component(mask: np.ndarray) -> np.ndarray:
    labels, n = label(mask)
    if n <= 1:
        return mask
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    return labels == counts.argmax()


def mask_canny(img, bg):
    """Canny edge → morphological close → fill interior via flood from corner."""
    g = _gray(img)
    edges = cv2.Canny(g, 50, 150)
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (7, 7))
    closed = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)
    h, w = closed.shape
    flood = closed.copy()
    cv2.floodFill(flood, np.zeros((h + 2, w + 2), np.uint8), (0, 0), 255)
    interior = (flood == 0)
    return _largest_component(interior | (closed > 0))


def mask_felzenszwalb(img, bg):
    """Region segmentation, foreground = largest segment whose mean luminance
    is sufficiently far from the expected background luminance."""
    segments = felzenszwalb(img, scale=500, sigma=0.8, min_size=500)
    g = _gray(img)
    bg_luminance = {"Black": 5.0, "Blue": 60.0}[bg]
    best_label, best_score = 0, -1.0
    for label_id in np.unique(segments):
        region = segments == label_id
        mean_l = g[region].mean()
        distance_from_bg = abs(mean_l - bg_luminance)
        # Score = region area weighted by colour distance from background
        score = region.sum() * distance_from_bg
        if distance_from_bg > 50 and score > best_score:
            best_score, best_label = score, label_id
    return segments == best_label
