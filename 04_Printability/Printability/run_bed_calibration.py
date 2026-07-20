"""Run bed calibration on the Black-BG photos and emit a sanity-check report.

Outputs (in results/Black-BG/_calibration/):
  - calibration_flash.json
  - calibration_noflash.json
  - sanity_<mode>.png   (BG photo with learned bed mask overlay)
"""
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
import matplotlib.pyplot as plt

from bed_calibration import (
    calibrate_from_background,
    bed_mask_from_calibration,
    load_rgb,
    save_calibration,
)


PROJECT_ROOT = Path(
    "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE"
)
BG_DIR = PROJECT_ROOT / "Analises/Python/Printability/results/Black-BG"
OUT_DIR = BG_DIR / "_calibration"
OUT_DIR.mkdir(exist_ok=True)


SOURCES = {
    "flash":   ("BG-Black-Pre-Flash.jpeg", True),
    "noflash": ("BG-Black-Pre.jpeg",       False),
}


def _overlay(rgb: np.ndarray, mask: np.ndarray, colour=(0, 255, 0), alpha=0.35) -> np.ndarray:
    out = rgb.copy()
    layer = np.zeros_like(rgb)
    layer[mask] = colour
    return cv2.addWeighted(out, 1.0, layer, alpha, 0)


def main() -> None:
    results = {}
    for mode, (filename, flash) in SOURCES.items():
        bg_path = BG_DIR / filename
        if not bg_path.exists():
            print(f"[warn] missing: {bg_path}")
            continue

        calib, bootstrap_mask = calibrate_from_background(bg_path, flash=flash)
        save_calibration(calib, OUT_DIR / f"calibration_{mode}.json")

        # Re-apply learned thresholds (with closing/opening/hull) for a clean ROI.
        rgb = load_rgb(bg_path)
        refined = bed_mask_from_calibration(rgb, calib)

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        axes[0].imshow(rgb)
        axes[0].set_title(f"{mode}: raw")
        axes[0].axis("off")
        axes[1].imshow(_overlay(rgb, bootstrap_mask, colour=(255, 0, 0)))
        axes[1].set_title(f"{mode}: bootstrap (V<.45, S<.40) → LCC")
        axes[1].axis("off")
        axes[2].imshow(_overlay(rgb, refined, colour=(0, 255, 0)))
        axes[2].set_title(f"{mode}: learned HSV + hull")
        axes[2].axis("off")
        fig.suptitle(
            f"Tape calibration — {mode}\n"
            f"H:[{calib.h_lo:.1f},{calib.h_hi:.1f}]  "
            f"S:[{calib.s_lo:.1f},{calib.s_hi:.1f}]  "
            f"V:[{calib.v_lo:.1f},{calib.v_hi:.1f}]  "
            f"bed_area={calib.bed_area_fraction:.1%}"
        )
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"sanity_{mode}.png", dpi=120)
        plt.close(fig)

        results[mode] = calib
        print(f"[{mode}] H=[{calib.h_lo:.1f},{calib.h_hi:.1f}] "
              f"S=[{calib.s_lo:.1f},{calib.s_hi:.1f}] "
              f"V=[{calib.v_lo:.1f},{calib.v_hi:.1f}] "
              f"bed={calib.bed_area_fraction:.1%}")

    print(f"\nWrote → {OUT_DIR}")


if __name__ == "__main__":
    main()
