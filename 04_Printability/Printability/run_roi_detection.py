"""End-to-end ROI detection on a foreground photo.

For each foreground:
  1. Bed ROI via GrabCut + hull.
  2. Logo localization + scale calibration.
  3. Label localization.
  4. Sidecar JSON + composite visualization.
"""
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from bed_calibration import bed_roi_grabcut, load_rgb
from roi_detection import (
    PhotoROIs,
    detect_labels,
    detect_logo,
    save_rois,
)


PROJECT_ROOT = Path(
    "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE"
)
PRINTABILITY = PROJECT_ROOT / "Analises/Python/Printability"
FG_DIR = PRINTABILITY / "results/Black-BG/S-Test"
OUT_DIR = PRINTABILITY / "results/Black-BG/_rois"
OUT_DIR.mkdir(exist_ok=True)


SAMPLES = {
    "flash":   "C15-SF5.5-Black-Vp10-Fr2-Flash.HEIC",
    "noflash": "C15-SF5.5-Black-Vp10-Fr2.HEIC",
}


def _overlay_mask(rgb, mask, colour, alpha=0.30):
    out = rgb.copy()
    layer = np.zeros_like(rgb)
    layer[mask] = colour
    return cv2.addWeighted(out, 1.0, layer, alpha, 0)


def main() -> None:
    for mode, fname in SAMPLES.items():
        fg_path = FG_DIR / fname
        if not fg_path.exists():
            print(f"[warn] missing: {fg_path}")
            continue

        rgb = load_rgb(fg_path)
        bed = bed_roi_grabcut(rgb)
        logo = detect_logo(rgb, bed)
        labels = detect_labels(rgb, bed)

        rois = PhotoROIs(
            source=fname,
            image_shape=(rgb.shape[0], rgb.shape[1]),
            logo=logo,
            labels=labels,
        )
        json_path = OUT_DIR / f"rois_{mode}.json"
        save_rois(rois, json_path)

        # Visualization
        fig, ax = plt.subplots(figsize=(10, 8))
        composite = _overlay_mask(rgb, bed, (0, 255, 0), alpha=0.20)
        ax.imshow(composite)
        if logo:
            x0, y0, x1, y1 = logo.bbox_xyxy
            ax.add_patch(mpatches.Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                linewidth=2, edgecolor="red", facecolor="none",
                label=f"logo ({logo.orientation}, {logo.mm_per_px:.4f} mm/px, "
                      f"residual {logo.scale_residual_pct:.1f}%)",
            ))
        for i, lab in enumerate(labels):
            x0, y0, x1, y1 = lab.bbox_xyxy
            ax.add_patch(mpatches.Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                linewidth=2, edgecolor="yellow", facecolor="none",
                label=f"label {i+1}" if i == 0 else None,
            ))
        ax.set_title(f"ROI detection — {mode}")
        ax.axis("off")
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"rois_{mode}.png", dpi=120)
        plt.close(fig)

        if logo:
            print(f"[{mode}] logo {logo.orientation} bbox={logo.bbox_xyxy} "
                  f"scale={logo.mm_per_px:.4f} mm/px "
                  f"residual={logo.scale_residual_pct:.1f}% | "
                  f"labels found: {len(labels)}")
        else:
            print(f"[{mode}] logo NOT detected | labels: {len(labels)}")


if __name__ == "__main__":
    main()
