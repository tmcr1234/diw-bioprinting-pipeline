"""Per-strand width, continuity, centreline deviation."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from skimage.morphology import medial_axis
from gcode_parser import RasterSpec
from registration import Alignment


@dataclass
class StrandReport:
    widths_px: list[float]
    continuities: list[float]
    centreline_dev_px: list[float]


def _strand_roi(seg, affine, mask_shape, roi_half_w: int = 30):
    p0 = affine @ np.array([seg.x0_mm, seg.y0_mm, 1.0])
    p1 = affine @ np.array([seg.x1_mm, seg.y1_mm, 1.0])
    rr0, cc0 = int(p0[1]), int(p0[0])
    rr1, cc1 = int(p1[1]), int(p1[0])
    rmin = max(min(rr0, rr1) - roi_half_w, 0)
    rmax = min(max(rr0, rr1) + roi_half_w, mask_shape[0])
    cmin = max(min(cc0, cc1), 0)
    cmax = min(max(cc0, cc1), mask_shape[1])
    return rmin, rmax, cmin, cmax, (rr0, cc0), (rr1, cc1)


def compute_strand_report(consensus_mask, spec: RasterSpec, alignment: Alignment) -> StrandReport:
    widths, conts, devs = [], [], []
    for seg in spec.strands:
        rmin, rmax, cmin, cmax, p0, p1 = _strand_roi(seg, alignment.affine, consensus_mask.shape)
        roi = consensus_mask[rmin:rmax, cmin:cmax]
        if not roi.any():
            widths.append(0.0)
            conts.append(0.0)
            devs.append(float("nan"))
            continue
        skel, dist = medial_axis(roi, return_distance=True)
        skel_dists = dist[skel]
        widths.append(float(2.0 * np.median(skel_dists)))
        nominal_len = abs(p1[1] - p0[1]) + abs(p1[0] - p0[0])
        conts.append(float(skel.sum() / max(nominal_len, 1)))
        ys, _ = np.where(skel)
        if ys.size:
            devs.append(float(np.abs(ys - (skel.shape[0] / 2)).mean()))
        else:
            devs.append(float("nan"))
    return StrandReport(widths, conts, devs)
