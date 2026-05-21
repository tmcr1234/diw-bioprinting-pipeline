"""HEIC photo loader for the printability pipeline.

iPhone standard Photo mode emits HEIC. We decode via pillow-heif and
return a numpy uint8 RGB array. No colour-space conversion; downstream
modules handle that.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from PIL import Image
import pillow_heif

pillow_heif.register_heif_opener()


def load_heic(path: Path) -> np.ndarray:
    img = Image.open(path).convert("RGB")
    return np.asarray(img, dtype=np.uint8)


def load_pair(flash_path: Path, noflash_path: Path) -> tuple[np.ndarray, np.ndarray]:
    flash = load_heic(flash_path)
    noflash = load_heic(noflash_path)
    if flash.shape != noflash.shape:
        raise ValueError(
            f"flash/noflash shape mismatch: {flash.shape} vs {noflash.shape}"
        )
    return flash, noflash
