"""Bed-colour calibration from a background photo.

Used when reference photos are NOT pose-locked to foreground shots, so the
background subtraction path in `background_model.py` is unreliable. Instead
we learn the HSV distribution of the tape from each BG photo and apply it
as a colour prior to segment the bed independently in every foreground.

Workflow:
  1. Load BG photo (HEIC or JPEG).
  2. Bootstrap a coarse "dark + desaturated" mask (V < 0.45, S < 0.40).
  3. Keep the largest connected component → assumed bed.
  4. Sample HSV stats inside that component.
  5. Persist {H, S, V} percentile bounds + bbox area to JSON.

The same module is later imported by the foreground segmentation step to
build a bed mask without needing pixel-aligned references.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import cv2
import numpy as np
from PIL import Image

import pillow_heif
pillow_heif.register_heif_opener()


# -- I/O ---------------------------------------------------------------------

def load_rgb(path: Path) -> np.ndarray:
    """Load HEIC/JPEG/PNG into uint8 RGB."""
    img = Image.open(path).convert("RGB")
    return np.asarray(img, dtype=np.uint8)


# -- Calibration -------------------------------------------------------------

@dataclass
class TapeCalibration:
    flash: bool
    h_lo: float
    h_hi: float
    s_lo: float
    s_hi: float
    v_lo: float
    v_hi: float
    bed_area_px: int
    bed_area_fraction: float
    image_shape: tuple[int, int]
    source_photo: str


def _select_bed_component(mask: np.ndarray) -> np.ndarray:
    """Pick the connected component most likely to be the bed.

    Scores each candidate by area × centrality. Centrality is the inverse
    distance of the component centroid to image center, normalized to [0,1].
    This biases against off-center clutter (cables, dark boxes near the
    frame) even when they are large in area.
    """
    num, labels, stats, centroids = cv2.connectedComponentsWithStats(
        mask.astype(np.uint8), connectivity=8
    )
    if num <= 1:
        return mask
    h, w = mask.shape
    cy, cx = h / 2.0, w / 2.0
    diag = float(np.hypot(h, w))
    best_score = -1.0
    best_label = 1
    for i in range(1, num):
        area = stats[i, cv2.CC_STAT_AREA]
        if area < 0.02 * h * w:  # ignore tiny blobs
            continue
        ccx, ccy = centroids[i]
        d = float(np.hypot(ccy - cy, ccx - cx))
        centrality = 1.0 - (d / (diag / 2.0))  # 1.0 at center, ~0 at corner
        score = area * max(centrality, 0.05)
        if score > best_score:
            best_score = score
            best_label = i
    return labels == best_label


def _bootstrap_bed_mask(hsv: np.ndarray) -> np.ndarray:
    """Bootstrap a bed mask via GrabCut seeded on the central rectangle.

    GrabCut builds Gaussian mixture models for fg/bg from the user-supplied
    bounding box (here: the central 70% of the image) and iteratively
    refines a pixel-wise segmentation. Robust to texture and uneven
    lighting in a way pure HSV thresholding is not.
    """
    bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)  # GrabCut wants BGR
    h, w = bgr.shape[:2]
    # Central 70% rectangle — the bed dominates the middle of every shot.
    rx, ry = int(0.15 * w), int(0.15 * h)
    rw, rh = int(0.70 * w), int(0.70 * h)
    rect = (rx, ry, rw, rh)

    mask = np.zeros((h, w), np.uint8)
    bgd_model = np.zeros((1, 65), np.float64)
    fgd_model = np.zeros((1, 65), np.float64)
    cv2.grabCut(bgr, mask, rect, bgd_model, fgd_model, 5, cv2.GC_INIT_WITH_RECT)
    bed = (mask == cv2.GC_FGD) | (mask == cv2.GC_PR_FGD)

    # Clean small holes and apply centrality-weighted LCC.
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7))
    bed = cv2.morphologyEx(bed.astype(np.uint8), cv2.MORPH_CLOSE, kernel)
    return _select_bed_component(bed.astype(bool))


def calibrate_from_background(
    bg_path: Path,
    flash: bool,
    pct_lo: float = 2.0,
    pct_hi: float = 98.0,
) -> tuple[TapeCalibration, np.ndarray]:
    """Learn HSV bounds for the tape and return (calibration, bed_mask).

    `pct_lo` / `pct_hi` define which percentile of the sampled-tape HSV
    distribution becomes the lower / upper threshold. 2/98 is robust to
    a few thousand contaminating pixels (specular highlights, glue spots).
    """
    rgb = load_rgb(bg_path)
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    bed = _bootstrap_bed_mask(hsv)

    H, S, V = cv2.split(hsv)
    h_lo, h_hi = np.percentile(H[bed], [pct_lo, pct_hi])
    s_lo, s_hi = np.percentile(S[bed], [pct_lo, pct_hi])
    v_lo, v_hi = np.percentile(V[bed], [pct_lo, pct_hi])

    calib = TapeCalibration(
        flash=flash,
        h_lo=float(h_lo), h_hi=float(h_hi),
        s_lo=float(s_lo), s_hi=float(s_hi),
        v_lo=float(v_lo), v_hi=float(v_hi),
        bed_area_px=int(bed.sum()),
        bed_area_fraction=float(bed.sum() / bed.size),
        image_shape=(int(rgb.shape[0]), int(rgb.shape[1])),
        source_photo=str(bg_path.name),
    )
    return calib, bed


# -- Apply learned calibration to any photo ---------------------------------

def bed_roi_grabcut(rgb: np.ndarray, central_frac: float = 0.70,
                    n_iter: int = 5) -> np.ndarray:
    """Segment the bed ROI via GrabCut with a central-rectangle prior.

    This is the production path for foreground photos. GrabCut builds
    GMM fg/bg colour models from the prior box and refines pixel-wise,
    so it tolerates pose drift, lighting drift, on-bed labels, and the
    printed object itself (all stay inside the bed ROI).

    Returns the convex hull of the segmented region — a clean polygonal
    bed ROI suitable as a downstream mask. Anything outside this hull
    is rejected for further analysis.
    """
    bgr = cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)
    h, w = bgr.shape[:2]
    margin = (1.0 - central_frac) / 2.0
    rect = (int(margin * w), int(margin * h),
            int(central_frac * w), int(central_frac * h))

    mask = np.zeros((h, w), np.uint8)
    bgd = np.zeros((1, 65), np.float64)
    fgd = np.zeros((1, 65), np.float64)
    cv2.grabCut(bgr, mask, rect, bgd, fgd, n_iter, cv2.GC_INIT_WITH_RECT)
    fg = (mask == cv2.GC_FGD) | (mask == cv2.GC_PR_FGD)

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9))
    fg = cv2.morphologyEx(fg.astype(np.uint8), cv2.MORPH_CLOSE, kernel)
    bed = _select_bed_component(fg.astype(bool))

    # Convex hull → bed ROI polygon.
    pts = np.column_stack(np.nonzero(bed))[:, ::-1]  # (x, y)
    if len(pts) < 3:
        return bed
    hull = cv2.convexHull(pts.astype(np.int32))
    hull_mask = np.zeros_like(bed, dtype=np.uint8)
    cv2.fillConvexPoly(hull_mask, hull, 1)
    return hull_mask.astype(bool)


def bed_mask_from_calibration(rgb: np.ndarray, calib: TapeCalibration) -> np.ndarray:
    """Back-compat wrapper. Calibration is now a documentation/QA artifact;
    actual segmentation runs GrabCut per-photo. The `calib` argument is
    accepted but unused — kept so existing callers don't break.
    """
    del calib  # documentation-only
    return bed_roi_grabcut(rgb)


def hsv_consistency_score(rgb: np.ndarray, bed: np.ndarray,
                          calib: TapeCalibration) -> dict[str, float]:
    """QA: how well do bed-pixel HSV stats match the BG calibration?

    Returns the fraction of bed pixels falling within each learned bound.
    Low scores flag substrate / lighting drift between BG and foreground.
    """
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    H, S, V = cv2.split(hsv)
    Hb, Sb, Vb = H[bed], S[bed], V[bed]
    return {
        "h_within": float(((Hb >= calib.h_lo) & (Hb <= calib.h_hi)).mean()),
        "s_within": float(((Sb >= calib.s_lo) & (Sb <= calib.s_hi)).mean()),
        "v_within": float(((Vb >= calib.v_lo) & (Vb <= calib.v_hi)).mean()),
    }


# -- Persistence -------------------------------------------------------------

def save_calibration(calib: TapeCalibration, path: Path) -> None:
    path.write_text(json.dumps(asdict(calib), indent=2))


def load_calibration(path: Path) -> TapeCalibration:
    d = json.loads(path.read_text())
    d["image_shape"] = tuple(d["image_shape"])
    return TapeCalibration(**d)
