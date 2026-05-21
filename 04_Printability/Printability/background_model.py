"""Reference-photo-driven background model.

Loads the 4 reference HEICs from a session's _reference/ folder, computes
per-pixel μ and σ, derives a ghost mask from |post − pre|, and persists
as .npz. Method 15 (segmentation/family5_bgdiff.py) consumes this model.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from io_heic import load_heic

REF_FILES = [
    "bg_pre_flash.heic",
    "bg_pre_noflash.heic",
    "bg_post_flash.heic",
    "bg_post_noflash.heic",
]


@dataclass
class BackgroundModel:
    mu: np.ndarray
    sigma: np.ndarray
    ghost_mask: np.ndarray
    drift_score: float


def build_background_model(ref_dir: Path) -> BackgroundModel:
    refs = [load_heic(ref_dir / name).astype(np.float32) for name in REF_FILES]
    stack = np.stack(refs, axis=0)
    mu = stack.mean(axis=0)
    # Per-pixel std across only 4 samples collapses to ~0 on uniform backgrounds;
    # floor at 2.0 (≈ 8-bit camera noise floor) so Mahalanobis distances stay sane.
    sigma = np.maximum(stack.std(axis=0), 2.0)
    pre_no = refs[1]
    post_no = refs[3]
    delta = np.abs(post_no - pre_no).mean(axis=2)
    sigma_floor = float(np.median(sigma.mean(axis=2)))
    ghost_mask = delta > (5.0 * sigma_floor)
    drift_score = float(delta.mean() / 255.0)
    return BackgroundModel(mu, sigma, ghost_mask, drift_score)


def save_model(model: BackgroundModel, path: Path) -> None:
    np.savez_compressed(
        path,
        mu=model.mu,
        sigma=model.sigma,
        ghost_mask=model.ghost_mask,
        drift_score=np.array(model.drift_score),
    )


def load_model(path: Path) -> BackgroundModel:
    z = np.load(path)
    return BackgroundModel(
        z["mu"], z["sigma"], z["ghost_mask"].astype(bool), float(z["drift_score"])
    )
