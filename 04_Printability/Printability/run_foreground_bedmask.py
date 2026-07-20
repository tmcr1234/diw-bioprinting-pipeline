"""Apply the saved Black-BG calibration to a foreground print photo.

Loads one S-Test foreground (Vp10-Fr2, flash and no-flash variants),
re-segments the bed using the learned HSV bounds, and writes a side-by-
side visualization plus the binary bed mask.

This validates pose-independence: the calibration was learned on a BG
photo at one camera pose, and we apply it to a foreground at a different
pose. If the bed mask still snaps to the actual bed boundary, the
calibration-only path is viable for this dataset.
"""
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
import matplotlib.pyplot as plt

from bed_calibration import (
    bed_roi_grabcut,
    hsv_consistency_score,
    load_calibration,
    load_rgb,
)


PROJECT_ROOT = Path(
    "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE"
)
PRINTABILITY = PROJECT_ROOT / "Analises/Python/Printability"
CAL_DIR = PRINTABILITY / "results/Black-BG/_calibration"
FG_DIR = PRINTABILITY / "results/Black-BG/S-Test"
OUT_DIR = PRINTABILITY / "results/Black-BG/_bedmask"
OUT_DIR.mkdir(exist_ok=True)


SAMPLES = {
    "flash":   "C15-SF5.5-Black-Vp10-Fr2-Flash.HEIC",
    "noflash": "C15-SF5.5-Black-Vp10-Fr2.HEIC",
}


def _overlay(rgb: np.ndarray, mask: np.ndarray, colour=(0, 255, 0), alpha=0.35) -> np.ndarray:
    out = rgb.copy()
    layer = np.zeros_like(rgb)
    layer[mask] = colour
    return cv2.addWeighted(out, 1.0, layer, alpha, 0)


def main() -> None:
    for mode, fname in SAMPLES.items():
        fg_path = FG_DIR / fname
        cal_path = CAL_DIR / f"calibration_{mode}.json"
        if not fg_path.exists() or not cal_path.exists():
            print(f"[warn] missing for {mode}: fg={fg_path.exists()} cal={cal_path.exists()}")
            continue

        calib = load_calibration(cal_path)
        rgb = load_rgb(fg_path)
        bed = bed_roi_grabcut(rgb)
        area_frac = float(bed.sum() / bed.size)
        consistency = hsv_consistency_score(rgb, bed, calib)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        axes[0].imshow(rgb)
        axes[0].set_title(f"{mode}: foreground (Vp10-Fr2)")
        axes[0].axis("off")
        axes[1].imshow(_overlay(rgb, bed))
        axes[1].set_title(f"{mode}: bed ROI (GrabCut + hull) — area={area_frac:.1%}")
        axes[1].axis("off")
        fig.suptitle(
            f"Pose-independent bed segmentation — {mode}\n"
            f"BG HSV calib: H=[{calib.h_lo:.0f},{calib.h_hi:.0f}] "
            f"S=[{calib.s_lo:.0f},{calib.s_hi:.0f}] "
            f"V=[{calib.v_lo:.0f},{calib.v_hi:.0f}] | "
            f"FG consistency: H={consistency['h_within']:.0%} "
            f"S={consistency['s_within']:.0%} "
            f"V={consistency['v_within']:.0%}"
        )
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"bedmask_{mode}.png", dpi=120)
        plt.close(fig)

        # Save binary mask for downstream consumption.
        cv2.imwrite(
            str(OUT_DIR / f"bedmask_{mode}.png").replace(".png", "_binary.png"),
            (bed * 255).astype(np.uint8),
        )
        print(f"[{mode}] bed area = {area_frac:.1%} → {OUT_DIR / f'bedmask_{mode}.png'}")


if __name__ == "__main__":
    main()
