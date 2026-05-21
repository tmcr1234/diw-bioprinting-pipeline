"""Per-turn fidelity: how cleanly was each S-turn deposited."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from gcode_parser import RasterSpec
from registration import Alignment


@dataclass
class TurnReport:
    n_clean: int
    n_total: int
    fidelity: float


def compute_turn_report(consensus_mask, spec: RasterSpec, alignment: Alignment,
                        median_strand_width_px: float) -> TurnReport:
    n_clean = 0
    for t in spec.turns:
        p = alignment.affine @ np.array([t.x_mm, t.y_mm, 1.0])
        r, c = int(p[1]), int(p[0])
        half = int(2.5 * median_strand_width_px)
        rmin = max(r - half, 0)
        rmax = min(r + half, consensus_mask.shape[0])
        cmin = max(c - half, 0)
        cmax = min(c + half, consensus_mask.shape[1])
        roi = consensus_mask[rmin:rmax, cmin:cmax]
        if not roi.any():
            continue
        col_sum = roi.any(axis=0)
        apex_width = int(col_sum.sum())
        if 0.7 * median_strand_width_px <= apex_width <= 1.5 * median_strand_width_px:
            n_clean += 1
    n_total = len(spec.turns)
    fidelity = n_clean / n_total if n_total else 1.0
    return TurnReport(n_clean, n_total, fidelity)
