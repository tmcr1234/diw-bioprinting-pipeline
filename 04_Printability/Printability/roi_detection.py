"""Per-photo ROI detection: logo (ruler) + identification labels.

The departmental red logo (80 x 41.7 mm) serves both as ruler and as
branding. Its position rotates between prints. White/cream paper labels
sit on the bed and must be excluded from ink analysis. This module
proposes both ROIs automatically; the optional interactive editor in
`roi_editor.py` lets the user adjust before downstream consumption.

ROI sidecar JSON schema (per foreground photo):
{
  "source": "C15-SF5.5-Black-Vp10-Fr2.HEIC",
  "image_shape": [H, W],
  "logo": {
    "bbox_xyxy": [x0, y0, x1, y1],
    "orientation": "horizontal" | "vertical",
    "long_axis_mm": 80.0,
    "short_axis_mm": 41.7,
    "mm_per_px": 0.123,
    "px_per_mm": 8.1,
    "scale_residual_pct": 1.4,
    "auto": true
  },
  "labels": [
    {"bbox_xyxy": [...], "auto": true},
    ...
  ]
}
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import cv2
import numpy as np


LOGO_LONG_MM = 80.0
LOGO_SHORT_MM = 41.7
LOGO_ASPECT = LOGO_LONG_MM / LOGO_SHORT_MM  # approx 1.918


@dataclass
class LogoROI:
    bbox_xyxy: list[int]
    orientation: str
    long_axis_mm: float
    short_axis_mm: float
    mm_per_px: float
    px_per_mm: float
    scale_residual_pct: float
    auto: bool


@dataclass
class LabelROI:
    bbox_xyxy: list[int]
    auto: bool


@dataclass
class PhotoROIs:
    source: str
    image_shape: tuple[int, int]
    logo: Optional[LogoROI]
    labels: list[LabelROI] = field(default_factory=list)


# -- Logo detection ----------------------------------------------------------

def _red_mask(rgb: np.ndarray, s_min: int = 110) -> np.ndarray:
    """HSV mask for the matte red logo material. Red wraps the hue circle
    (0/180), so we OR two ranges. `s_min` is configurable so callers can
    relax for flash photos (specular highlights desaturate the red).
    """
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    lo1 = (0, s_min, 70)
    hi1 = (10, 255, 255)
    lo2 = (170, s_min, 70)
    hi2 = (180, 255, 255)
    m1 = cv2.inRange(hsv, lo1, hi1)
    m2 = cv2.inRange(hsv, lo2, hi2)
    return ((m1 | m2) > 0).astype(np.uint8)


def _erode_mask(mask: np.ndarray, frac: float = 0.03) -> np.ndarray:
    """Erode a binary mask by a fraction of the image's shorter side."""
    h, w = mask.shape
    k = max(3, int(frac * min(h, w)))
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    return cv2.erode(mask.astype(np.uint8), kernel).astype(bool)


def _logo_candidates(red: np.ndarray, close_kernel: int,
                     aspect_tol: float, min_long_px: int) -> list[tuple]:
    """Connected components -> close stencil gaps -> convex hull ->
    minAreaRect -> aspect-filtered candidates."""
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (close_kernel, close_kernel))
    closed = cv2.morphologyEx(red, cv2.MORPH_CLOSE, kernel)
    num, labels, stats, _ = cv2.connectedComponentsWithStats(closed, connectivity=8)
    out = []
    for i in range(1, num):
        area = int(stats[i, cv2.CC_STAT_AREA])
        if area < 800:
            continue
        pts = np.column_stack(np.nonzero(labels == i))[:, ::-1]
        if len(pts) < 4:
            continue
        # Convex hull of the M-stencil pixels -> bounding rect of the
        # *logo as a whole*, not just its strokes.
        hull = cv2.convexHull(pts.astype(np.int32))
        rect = cv2.minAreaRect(hull.reshape(-1, 2).astype(np.float32))
        (_, _), (w_px, h_px), _ = rect
        long_px = max(w_px, h_px)
        short_px = min(w_px, h_px)
        if long_px < min_long_px or short_px < min_long_px / 3:
            continue
        ar = long_px / short_px
        ar_err = abs(ar - LOGO_ASPECT) / LOGO_ASPECT
        if ar_err > aspect_tol:
            continue
        # Score: prefer large hull (true logo) with low aspect-ratio error.
        hull_area = cv2.contourArea(hull)
        score = float(hull_area) * (1.0 - ar_err)
        out.append((score, rect, long_px, short_px, ar_err))
    return out


def detect_logo(rgb: np.ndarray, bed_mask: np.ndarray,
                aspect_tol: float = 0.20) -> Optional[LogoROI]:
    """Locate the red logo inside the bed ROI and derive scale.

    The logo is a stencil M-shape, so red pixels form disconnected
    strokes. Strategy:
      1. Red HSV (tight first, relaxed fallback for flash highlights).
      2. Restrict to eroded bed (no border bleed-through).
      3. Close stencil gaps, take connected components.
      4. Convex hull of each component -> true outer rectangle via
         minAreaRect.
      5. Reject candidates whose long/short ratio deviates from 80/41.7.
      6. Score by hull area x (1 - aspect_error) -> pick best.
    """
    inside_bed = _erode_mask(bed_mask, frac=0.03)
    min_long_px = int(0.10 * min(rgb.shape[:2]))  # logo > 10% of image short side

    candidates = []
    for s_min, close_k in [(140, 25), (110, 31), (80, 41)]:
        red = _red_mask(rgb, s_min=s_min) * inside_bed.astype(np.uint8)
        candidates = _logo_candidates(red, close_kernel=close_k,
                                       aspect_tol=aspect_tol,
                                       min_long_px=min_long_px)
        if candidates:
            break
    if not candidates:
        return None

    _, rect, long_px, short_px, ar_err = max(candidates, key=lambda c: c[0])
    box = cv2.boxPoints(rect).astype(np.int32)
    x0, y0 = int(box[:, 0].min()), int(box[:, 1].min())
    x1, y1 = int(box[:, 0].max()), int(box[:, 1].max())
    orientation = "horizontal" if (x1 - x0) >= (y1 - y0) else "vertical"

    mm_per_px_long = LOGO_LONG_MM / long_px
    mm_per_px_short = LOGO_SHORT_MM / short_px
    mm_per_px = 0.5 * (mm_per_px_long + mm_per_px_short)
    residual = 100.0 * abs(mm_per_px_long - mm_per_px_short) / mm_per_px

    return LogoROI(
        bbox_xyxy=[x0, y0, x1, y1],
        orientation=orientation,
        long_axis_mm=LOGO_LONG_MM,
        short_axis_mm=LOGO_SHORT_MM,
        mm_per_px=float(mm_per_px),
        px_per_mm=float(1.0 / mm_per_px),
        scale_residual_pct=float(residual),
        auto=True,
    )


# -- Label detection ---------------------------------------------------------

def detect_labels(rgb: np.ndarray, bed_mask: np.ndarray,
                  min_area_px: int = 3000) -> list[LabelROI]:
    """Locate cream/white paper labels inside the bed ROI. High-V, low-S
    achromatic patches on the dark tape. Bed is eroded first so bench
    paper just beyond the bed edge does not bleed in."""
    inside_bed = _erode_mask(bed_mask, frac=0.03)
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    S = hsv[..., 1]
    V = hsv[..., 2]
    paper = ((V > 185) & (S < 60)).astype(np.uint8) * inside_bed.astype(np.uint8)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7))
    paper = cv2.morphologyEx(paper, cv2.MORPH_OPEN, kernel)
    paper = cv2.morphologyEx(paper, cv2.MORPH_CLOSE, kernel)

    num, labels, stats, _ = cv2.connectedComponentsWithStats(paper, connectivity=8)
    out = []
    for i in range(1, num):
        if stats[i, cv2.CC_STAT_AREA] < min_area_px:
            continue
        x = int(stats[i, cv2.CC_STAT_LEFT])
        y = int(stats[i, cv2.CC_STAT_TOP])
        w = int(stats[i, cv2.CC_STAT_WIDTH])
        h = int(stats[i, cv2.CC_STAT_HEIGHT])
        out.append(LabelROI(bbox_xyxy=[x, y, x + w, y + h], auto=True))
    return out


# -- Persistence -------------------------------------------------------------

def save_rois(rois: PhotoROIs, path: Path) -> None:
    payload = {
        "source": rois.source,
        "image_shape": list(rois.image_shape),
        "logo": asdict(rois.logo) if rois.logo else None,
        "labels": [asdict(l) for l in rois.labels],
    }
    path.write_text(json.dumps(payload, indent=2))


def load_rois(path: Path) -> PhotoROIs:
    d = json.loads(path.read_text())
    logo = LogoROI(**d["logo"]) if d.get("logo") else None
    labels = [LabelROI(**l) for l in d.get("labels", [])]
    return PhotoROIs(
        source=d["source"],
        image_shape=tuple(d["image_shape"]),
        logo=logo,
        labels=labels,
    )
