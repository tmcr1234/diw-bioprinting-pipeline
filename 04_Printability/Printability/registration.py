"""Image registration utilities.

Part 1: align_pair() — match no-flash to flash by phase correlation.
Part 2: align_gcode_to_image() — affine fit of nominal raster bbox
to consensus mask foreground via rotation candidate search.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import cv2
from skimage.registration import phase_cross_correlation
from gcode_parser import RasterSpec


def align_pair(flash: np.ndarray, noflash: np.ndarray) -> tuple[np.ndarray, tuple[int, int]]:
    f_gray = flash.mean(axis=2)
    n_gray = noflash.mean(axis=2)
    shift, _, _ = phase_cross_correlation(f_gray, n_gray, upsample_factor=4)
    dy_inv, dx_inv = int(round(shift[0])), int(round(shift[1]))
    aligned = np.roll(noflash, shift=(dy_inv, dx_inv), axis=(0, 1))
    return aligned, (-dy_inv, -dx_inv)


@dataclass
class Alignment:
    affine: np.ndarray
    iou: float
    success: bool


def _bbox(mask: np.ndarray) -> tuple[int, int, int, int]:
    ys, xs = np.where(mask)
    if ys.size == 0:
        return (0, 0, 0, 0)
    return (xs.min(), ys.min(), xs.max(), ys.max())


def _rasterize_nominal(spec: RasterSpec, affine: np.ndarray, shape) -> np.ndarray:
    """Render the nominal raster footprint (bounding rectangle of all strands)
    as the alignment reference. We compare envelope-to-envelope rather than
    line-to-mask because the consensus mask is a solid foreground blob and
    line-based IoU would always be small."""
    out = np.zeros(shape, dtype=np.uint8)
    pts = []
    for s in spec.strands:
        for x_mm, y_mm in [(s.x0_mm, s.y0_mm), (s.x1_mm, s.y1_mm)]:
            p = affine @ np.array([x_mm, y_mm, 1.0])
            pts.append((int(p[0]), int(p[1])))
    if not pts:
        return out > 0
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    cv2.rectangle(out, (min(xs), min(ys)), (max(xs), max(ys)), 255, thickness=-1)
    return out > 0


def align_gcode_to_image(spec: RasterSpec, consensus_mask: np.ndarray) -> Alignment:
    bx0, by0, bx1, by1 = _bbox(consensus_mask)
    gx0, gy0, gx1, gy1 = spec.bounding_box_mm
    bbox_w_px = max(bx1 - bx0, 1)
    bbox_h_px = max(by1 - by0, 1)
    gw = max(gx1 - gx0, 1e-6)
    gh = max(gy1 - gy0, 1e-6)

    best_iou, best_affine = -1.0, None
    for rotation_deg in np.linspace(-10, 10, 21):
        for flip in (False, True):
            theta = np.deg2rad(rotation_deg)
            sx = bbox_w_px / gw
            sy = bbox_h_px / gh * (-1 if flip else 1)
            R = np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta),  np.cos(theta)]])
            S = np.diag([sx, sy])
            M = R @ S
            tx = bx0 - M[0, 0] * gx0 - M[0, 1] * gy0
            ty = by0 - M[1, 0] * gx0 - M[1, 1] * gy0
            affine = np.array([[M[0, 0], M[0, 1], tx],
                               [M[1, 0], M[1, 1], ty]])
            rendered = _rasterize_nominal(spec, affine, consensus_mask.shape)
            iou = float((rendered & consensus_mask).sum() / max((rendered | consensus_mask).sum(), 1))
            if iou > best_iou:
                best_iou, best_affine = iou, affine
    return Alignment(best_affine, best_iou, success=best_iou >= 0.7)
