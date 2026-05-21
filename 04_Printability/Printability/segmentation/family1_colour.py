"""Family 1: colour-space segmentation methods.

Method 1 — HSV threshold (bg-aware).
Method 2 — Lab a/b chrominance threshold.
Method 3 — k-means K=2 in Lab.
Method 4 — k-means K=3 in RGB, foreground = brightest cluster.
"""
from __future__ import annotations
import numpy as np
import cv2
from sklearn.cluster import KMeans


def mask_hsv(img: np.ndarray, bg: str) -> np.ndarray:
    hsv = cv2.cvtColor(img, cv2.COLOR_RGB2HSV)
    if bg == "Black":
        return hsv[..., 2] > 60
    elif bg == "Blue":
        H = hsv[..., 0]
        return ~((H >= 95) & (H <= 135))
    raise ValueError(bg)


def mask_lab_ab(img: np.ndarray, bg: str) -> np.ndarray:
    lab = cv2.cvtColor(img, cv2.COLOR_RGB2LAB)
    b = lab[..., 2].astype(int) - 128
    if bg == "Black":
        return lab[..., 0] > 60
    elif bg == "Blue":
        return b > -10
    raise ValueError(bg)


def _kmeans_mask(features: np.ndarray, n_clusters: int, shape: tuple[int, int]) -> np.ndarray:
    """Cluster pixels; return mask of cluster with highest mean brightness."""
    km = KMeans(n_clusters=n_clusters, n_init=4, random_state=0)
    labels = km.fit_predict(features)
    brightness = km.cluster_centers_[:, :3].mean(axis=1)
    fg_label = int(np.argmax(brightness))
    return (labels == fg_label).reshape(shape)


def mask_kmeans_lab(img: np.ndarray, bg: str) -> np.ndarray:
    lab = cv2.cvtColor(img, cv2.COLOR_RGB2LAB).astype(np.float32).reshape(-1, 3)
    return _kmeans_mask(lab, n_clusters=2, shape=img.shape[:2])


def mask_kmeans_rgb(img: np.ndarray, bg: str) -> np.ndarray:
    rgb = img.astype(np.float32).reshape(-1, 3)
    return _kmeans_mask(rgb, n_clusters=3, shape=img.shape[:2])
