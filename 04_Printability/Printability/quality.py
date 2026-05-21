"""Per-photo quality metrics. All flags advisory, never fatal."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
import cv2
from background_model import BackgroundModel


@dataclass
class QualityReport:
    noise_floor_sigma: float
    SNR: float
    var_laplacian: float
    clipped_fraction: float
    vignette_deviation: float
    reference_drift_mean: float
    pair_focus_diff: float
    flags: dict[str, bool] = field(default_factory=dict)


def _var_laplacian(img: np.ndarray) -> float:
    g = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    return float(cv2.Laplacian(g, cv2.CV_64F).var())


def _clipped_fraction(img: np.ndarray) -> float:
    return float(((img == 0) | (img == 255)).any(axis=2).mean())


def _vignette_deviation(img: np.ndarray, bg_mask: np.ndarray) -> float:
    g = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY).astype(np.float32)
    if bg_mask.sum() < 100:
        return 0.0
    ys, xs = np.where(bg_mask)
    Z = g[bg_mask]
    A = np.column_stack([np.ones_like(xs), xs, ys, xs * xs, ys * ys, xs * ys])
    coef, *_ = np.linalg.lstsq(A, Z, rcond=None)
    fitted = A @ coef
    return float(np.abs(Z - fitted).max() / 255.0)


def compute_quality(flash, noflash, foreground_mask, bg_model: BackgroundModel) -> QualityReport:
    noise_floor_sigma = float(bg_model.sigma.mean())
    fg_mean = float(noflash[foreground_mask].mean()) if foreground_mask.any() else 0.0
    bg_mean = float(noflash[~foreground_mask].mean())
    SNR = (fg_mean - bg_mean) / max(noise_floor_sigma, 1e-3)

    vl_flash = _var_laplacian(flash)
    vl_no = _var_laplacian(noflash)
    var_laplacian = vl_no

    clipped = max(_clipped_fraction(flash), _clipped_fraction(noflash))
    vignette = _vignette_deviation(noflash, ~foreground_mask)

    drift = bg_model.drift_score
    pair_focus_diff = abs(vl_flash - vl_no) / max(vl_flash, vl_no, 1e-3)

    flags = {
        "low_contrast":         SNR < 10,
        "blurry":               var_laplacian < 100,
        "clipped":              clipped > 0.01,
        "uneven_lighting":      vignette > 0.15,
        "session_drift":        drift > 0.05,
        "pair_focus_mismatch":  pair_focus_diff > 0.30,
    }
    return QualityReport(noise_floor_sigma, SNR, var_laplacian, clipped,
                         vignette, drift, pair_focus_diff, flags)
