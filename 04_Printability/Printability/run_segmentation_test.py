"""Standalone segmentation-ensemble smoke test on a single real photo.

Runs the classical methods (no bgDiff/SAM2 — those need a background model /
checkpoint) on a bed-cropped iPhone HEIC and renders a per-method montage plus
a consensus overlay. This is a CAPABILITY test on the dark-on-dark stress case,
not a metrology run (no g-code registration, no pixel->mm scale).
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from io_heic import load_heic
from segmentation import METHODS, run_method
from consensus import compute_consensus

SRC = Path("/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/"
           "Tinta NE/Impressao 3D/Fotos Teste/Camera/2 de junho de 2026/IMG_7510.heic")
OUT = Path(__file__).parent / "results" / "_pipeline_test"
OUT.mkdir(parents=True, exist_ok=True)
BG = "Black"
# Bed crop (fractions of W,H) isolating tape+print, excluding ruler/frame/clutter.
CROP = dict(x0=0.10, x1=0.48, y0=0.26, y1=0.62)
SKIP = {"bgDiff", "SAM2"}          # need model / checkpoint
LONG_EDGE = 1100                   # downscale crop for speed


def main() -> int:
    img = load_heic(SRC)                          # RGB uint8
    H, W = img.shape[:2]
    x0, x1 = int(W*CROP["x0"]), int(W*CROP["x1"])
    y0, y1 = int(H*CROP["y0"]), int(H*CROP["y1"])
    crop = img[y0:y1, x0:x1]
    s = LONG_EDGE / max(crop.shape[:2])
    if s < 1:
        crop = cv2.resize(crop, None, fx=s, fy=s, interpolation=cv2.INTER_AREA)
    print(f"crop {crop.shape[1]}x{crop.shape[0]} px from {W}x{H} source")

    masks, times = {}, {}
    for name in METHODS:
        if name in SKIP:
            continue
        t = time.time()
        try:
            m = run_method(name, flash=crop, noflash=crop, bg=BG).astype(bool)
        except Exception as e:
            print(f"  {name:14s} FAILED: {e}")
            continue
        masks[name] = m
        times[name] = time.time() - t
        print(f"  {name:14s} fg={100*m.mean():5.1f}%  {times[name]:5.2f}s")

    consensus, agreement = compute_consensus(masks)
    print(f"\nconsensus fg={100*consensus.mean():.1f}%  agreement={agreement:.3f}  "
          f"methods={len(masks)}")

    # --- per-method montage ---
    n = len(masks) + 2
    cols = 4; rows = (n + cols - 1)//cols
    fig, ax = plt.subplots(rows, cols, figsize=(4*cols, 3.2*rows))
    ax = ax.ravel()
    ax[0].imshow(crop); ax[0].set_title("crop (bed)"); ax[0].axis("off")
    for i, (name, m) in enumerate(masks.items(), start=1):
        ax[i].imshow(m, cmap="gray"); ax[i].axis("off")
        ax[i].set_title(f"{name}  {100*m.mean():.0f}%")
    j = len(masks)+1
    ax[j].imshow(consensus, cmap="magma"); ax[j].axis("off")
    ax[j].set_title(f"CONSENSUS  agr={agreement:.2f}")
    for k in range(j+1, len(ax)): ax[k].axis("off")
    fig.tight_layout(); fig.savefig(OUT/"montage.png", dpi=110); plt.close(fig)

    # --- consensus overlay on the crop ---
    ov = crop.copy()
    ov[consensus] = (0.45*ov[consensus] + 0.55*np.array([0,255,0])).astype(np.uint8)
    plt.figure(figsize=(7,6)); plt.imshow(ov); plt.axis("off")
    plt.title(f"Consensus overlay  (fg={100*consensus.mean():.1f}%, agr={agreement:.2f})")
    plt.tight_layout(); plt.savefig(OUT/"consensus_overlay.png", dpi=130); plt.close()

    print(f"\nwrote {OUT/'montage.png'}\nwrote {OUT/'consensus_overlay.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
