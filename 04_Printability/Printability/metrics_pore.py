"""Ouyang printability index Pr = L²/(16A) per inter-strand pore."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import cv2
from gcode_parser import RasterSpec
from registration import Alignment


@dataclass
class PoreReport:
    Pr_values: list[float]
    Pr_ideal_serpentine: float
    pore_areas_px: list[float]
    pore_perimeters_px: list[float]


def ideal_pr_from_spec(spec: RasterSpec) -> float:
    """Pr of a rectangular pore of width strand_length × height strand_spacing.
    Ideal Pr = (2·(W + H))² / (16·W·H).
    """
    W = spec.strand_length_mm
    H = spec.strand_spacing_mm
    if W <= 0 or H <= 0:
        return float("nan")
    L = 2 * (W + H)
    return (L * L) / (16.0 * W * H)


def _pore_between(seg_i, seg_j, affine, mask_shape, pad: int = 10):
    p0_i = affine @ np.array([seg_i.x0_mm, seg_i.y0_mm, 1.0])
    p1_i = affine @ np.array([seg_i.x1_mm, seg_i.y1_mm, 1.0])
    p0_j = affine @ np.array([seg_j.x0_mm, seg_j.y0_mm, 1.0])
    p1_j = affine @ np.array([seg_j.x1_mm, seg_j.y1_mm, 1.0])
    pts = np.array([p0_i, p1_i, p0_j, p1_j])
    rmin = max(int(pts[:, 1].min()) - pad, 0)
    rmax = min(int(pts[:, 1].max()) + pad, mask_shape[0])
    cmin = max(int(pts[:, 0].min()) - pad, 0)
    cmax = min(int(pts[:, 0].max()) + pad, mask_shape[1])
    return rmin, rmax, cmin, cmax


def compute_pore_report(consensus_mask, spec: RasterSpec, alignment: Alignment) -> PoreReport:
    Pr_vals, areas, perims = [], [], []
    for i in range(len(spec.strands) - 1):
        rmin, rmax, cmin, cmax = _pore_between(
            spec.strands[i], spec.strands[i + 1], alignment.affine, consensus_mask.shape
        )
        roi_bg = ~consensus_mask[rmin:rmax, cmin:cmax]
        roi_bg_u8 = roi_bg.astype(np.uint8) * 255
        contours, _ = cv2.findContours(roi_bg_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            continue
        largest = max(contours, key=cv2.contourArea)
        A = float(cv2.contourArea(largest))
        L = float(cv2.arcLength(largest, True))
        if A <= 0:
            continue
        Pr = (L * L) / (16.0 * A)
        Pr_vals.append(Pr)
        areas.append(A)
        perims.append(L)
    return PoreReport(Pr_vals, ideal_pr_from_spec(spec), areas, perims)
