"""Method registry: name → callable + channel preference."""
from __future__ import annotations
import numpy as np
from .family1_colour import mask_hsv, mask_lab_ab, mask_kmeans_lab, mask_kmeans_rgb
from .family2_global_thr import mask_otsu, mask_triangle, mask_li, mask_yen
from .family3_local_thr import mask_adaptive_gauss, mask_sauvola, mask_niblack
from .family4_edge_region import mask_canny, mask_felzenszwalb
from .family5_bgdiff import mask_bgdiff
from .family6_sam2 import mask_sam2

METHODS: dict[str, dict] = {
    "HSV":          {"fn": mask_hsv,            "channel": "flash",   "family": 1},
    "Lab_ab":       {"fn": mask_lab_ab,         "channel": "noflash", "family": 1},
    "kmeans_lab":   {"fn": mask_kmeans_lab,     "channel": "noflash", "family": 1},
    "kmeans_rgb":   {"fn": mask_kmeans_rgb,     "channel": "noflash", "family": 1},
    "Otsu":         {"fn": mask_otsu,           "channel": "flash",   "family": 2},
    "Triangle":     {"fn": mask_triangle,       "channel": "flash",   "family": 2},
    "Li":           {"fn": mask_li,             "channel": "flash",   "family": 2},
    "Yen":          {"fn": mask_yen,            "channel": "flash",   "family": 2},
    "adaptiveGauss":{"fn": mask_adaptive_gauss, "channel": "noflash", "family": 3},
    "Sauvola":      {"fn": mask_sauvola,        "channel": "noflash", "family": 3},
    "Niblack":      {"fn": mask_niblack,        "channel": "noflash", "family": 3},
    "Canny":        {"fn": mask_canny,          "channel": "flash",   "family": 4},
    "Felzenszwalb": {"fn": mask_felzenszwalb,   "channel": "noflash", "family": 4},
    "bgDiff":       {"fn": mask_bgdiff,         "channel": "both",    "family": 5},
    "SAM2":         {"fn": mask_sam2,           "channel": "noflash", "family": 6},
}


def list_methods(include_sam2: bool = True) -> list[str]:
    names = [n for n in METHODS if n != "SAM2"]
    if include_sam2 and "SAM2" in METHODS:
        names.append("SAM2")
    return names


def run_method(name: str, flash: np.ndarray, noflash: np.ndarray, bg: str, **kwargs) -> np.ndarray:
    spec = METHODS[name]
    if spec["channel"] == "both":
        return spec["fn"](flash=flash, noflash=noflash, bg=bg, **kwargs)
    img = flash if spec["channel"] == "flash" else noflash
    return spec["fn"](img, bg, **kwargs) if kwargs else spec["fn"](img, bg)
