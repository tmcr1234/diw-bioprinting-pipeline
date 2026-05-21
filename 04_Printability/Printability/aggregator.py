"""Append-safe CSV writer for master_printability.csv."""
from __future__ import annotations
from pathlib import Path
import pandas as pd

PER_METHOD_PR = [
    "Pr_HSV", "Pr_Lab_ab", "Pr_kmeans_lab", "Pr_kmeans_rgb",
    "Pr_Otsu", "Pr_Triangle", "Pr_Li", "Pr_Yen",
    "Pr_adaptiveGauss", "Pr_Sauvola", "Pr_Niblack",
    "Pr_Canny", "Pr_Felzenszwalb", "Pr_bgDiff", "Pr_SAM2",
]

FLAGS = [
    "low_contrast", "blurry", "clipped", "uneven_lighting",
    "session_drift", "pair_focus_mismatch",
    "alignment_failed", "low_confidence",
]

COLUMNS = [
    "ink", "bg", "shape", "Vp", "Flow", "gcode_file",
    "Pr_consensus_mean", "Pr_consensus_std",
    "Pr_consensus_min", "Pr_consensus_max",
    "Pr_ideal_serpentine",
    "strand_w_mean_px", "strand_w_std_px", "strand_w_cov",
    "n_strands_resolved", "n_strands_nominal",
    "turn_fidelity", "n_turns_clean", "n_turns_nominal",
    "agreement_score", "alignment_iou",
    *PER_METHOD_PR,
    *FLAGS,
    "overlay_path", "debug_dir",
]


def read_csv(csv_path: Path) -> pd.DataFrame:
    if not Path(csv_path).exists():
        return pd.DataFrame(columns=COLUMNS)
    return pd.read_csv(csv_path)


def write_row(csv_path: Path, row: dict) -> None:
    csv_path = Path(csv_path)
    df = read_csv(csv_path)
    key_vp = row["Vp"]
    key_flow = row["Flow"]
    if not df.empty:
        mask = (df["Vp"] == key_vp) & (df["Flow"] == key_flow)
        if mask.any():
            df = df.loc[~mask]
    df = pd.concat([df, pd.DataFrame([row], columns=COLUMNS)], ignore_index=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
