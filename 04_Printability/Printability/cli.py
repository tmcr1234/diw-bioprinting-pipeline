"""Argparse surface for the pipeline."""
from __future__ import annotations
import argparse
import os


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="DIW printability pipeline (S-Test)")
    p.add_argument("--shape", required=True)
    p.add_argument("--ink", required=True)
    p.add_argument("--bg", required=True, choices=["Black", "Blue"])
    p.add_argument("--workers", type=int, default=max(os.cpu_count() - 1, 1))
    p.add_argument("--use-sam2", dest="use_sam2", action="store_true", default=True)
    p.add_argument("--no-sam2", dest="use_sam2", action="store_false")
    p.add_argument("--debug", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--batch", action="store_true",
                   help="Iterate every (shape, ink, bg) folder with photos")
    return p.parse_args(argv)
