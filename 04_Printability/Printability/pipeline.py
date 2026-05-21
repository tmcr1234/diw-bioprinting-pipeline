"""Pipeline orchestration: discover → process → aggregate."""
from __future__ import annotations
import hashlib
import json
import multiprocessing as mp
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cli import parse_args
from io_heic import load_pair
from gcode_parser import parse_gcode, parse_filename
from registration import align_pair, align_gcode_to_image
from background_model import build_background_model, save_model, load_model
from segmentation import METHODS, run_method
from consensus import compute_consensus
from quality import compute_quality
from metrics_strand import compute_strand_report
from metrics_pore import compute_pore_report
from metrics_turn import compute_turn_report
from aggregator import COLUMNS, write_row
from cache import CacheIndex

PROJECT_ROOT = Path("/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Doutorado/Tinta NE")


@dataclass
class PairSpec:
    Vp: float
    Flow: float
    flash: Path
    noflash: Path
    gcode: Path


def discover_pairs(photo_dir: Path, shape_gcode_dir: Path) -> list[PairSpec]:
    pairs: dict[tuple[float, float], dict] = {}
    for f in photo_dir.glob("*flash.heic"):
        if f.name.endswith("noflash.heic"):
            continue
        Vp, Flow = parse_filename(f)
        pairs.setdefault((Vp, Flow), {})["flash"] = f
    for f in photo_dir.glob("*noflash.heic"):
        Vp, Flow = parse_filename(f)
        pairs.setdefault((Vp, Flow), {})["noflash"] = f
    out = []
    for (Vp, Flow), d in pairs.items():
        if "flash" not in d or "noflash" not in d:
            continue
        gcandidates = list(shape_gcode_dir.glob(f"*Vp{int(Vp)}*"))
        if not gcandidates:
            continue
        gcode = gcandidates[0]
        for g in gcandidates:
            _, gflow = parse_filename(g)
            if abs(gflow - Flow) < 1e-3:
                gcode = g
                break
        out.append(PairSpec(Vp, Flow, d["flash"], d["noflash"], gcode))
    return sorted(out, key=lambda p: (p.Vp, p.Flow))


def process_pair(pair: PairSpec, model_path: Path, sam2_client, bg: str) -> dict:
    flash, noflash = load_pair(pair.flash, pair.noflash)
    noflash_aligned, _ = align_pair(flash, noflash)
    model = load_model(model_path)

    masks = {}
    for name in METHODS:
        try:
            if name == "SAM2":
                masks[name] = run_method(name, flash=flash, noflash=noflash_aligned,
                                         bg=bg, service=sam2_client)
            elif name == "bgDiff":
                masks[name] = run_method(name, flash=flash, noflash=noflash_aligned,
                                         bg=bg, model=model)
            else:
                masks[name] = run_method(name, flash=flash, noflash=noflash_aligned, bg=bg)
        except Exception:
            masks[name] = np.zeros(flash.shape[:2], dtype=bool)

    consensus_mask, agreement = compute_consensus(masks)

    spec = parse_gcode(pair.gcode)
    alignment = align_gcode_to_image(spec, consensus_mask)
    strand = compute_strand_report(consensus_mask, spec, alignment)
    pore = compute_pore_report(consensus_mask, spec, alignment)
    turn = compute_turn_report(
        consensus_mask, spec, alignment,
        median_strand_width_px=(float(np.median(strand.widths_px)) if strand.widths_px else 1.0),
    )
    quality = compute_quality(flash, noflash_aligned, consensus_mask, model)

    row = {c: None for c in COLUMNS}
    row.update({
        "ink": pair.flash.parents[1].name,
        "bg":  pair.flash.parent.name,
        "shape": pair.flash.parents[2].name,
        "Vp": pair.Vp, "Flow": pair.Flow,
        "gcode_file": str(pair.gcode),
        "Pr_consensus_mean": float(np.mean(pore.Pr_values)) if pore.Pr_values else None,
        "Pr_consensus_std":  float(np.std(pore.Pr_values))  if pore.Pr_values else None,
        "Pr_consensus_min":  float(np.min(pore.Pr_values))  if pore.Pr_values else None,
        "Pr_consensus_max":  float(np.max(pore.Pr_values))  if pore.Pr_values else None,
        "Pr_ideal_serpentine": pore.Pr_ideal_serpentine,
        "strand_w_mean_px": float(np.mean(strand.widths_px)) if strand.widths_px else None,
        "strand_w_std_px":  float(np.std(strand.widths_px))  if strand.widths_px else None,
        "strand_w_cov": (float(np.std(strand.widths_px) / np.mean(strand.widths_px))
                         if strand.widths_px and np.mean(strand.widths_px) > 0 else None),
        "n_strands_resolved": sum(1 for c in strand.continuities if c > 0.5),
        "n_strands_nominal": len(spec.strands),
        "turn_fidelity": turn.fidelity,
        "n_turns_clean": turn.n_clean,
        "n_turns_nominal": turn.n_total,
        "agreement_score": agreement,
        "alignment_iou": alignment.iou,
        "alignment_failed": not alignment.success,
        "low_confidence": agreement < 0.7 or len(masks) < 8,
        **quality.flags,
    })
    return row


def main() -> int:
    args = parse_args()
    photo_dir = PROJECT_ROOT / "Impressao 3D" / args.shape / args.ink / args.bg
    gcode_dir = PROJECT_ROOT / "Impressao 3D" / args.shape
    out_dir = PROJECT_ROOT / "Analises/Python/Printability/results" / args.shape / args.ink / args.bg
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_dir = photo_dir / "_reference"
    model_path = ref_dir / "background_model.npz"
    if not model_path.exists() or args.force:
        if not ref_dir.exists():
            print(f"ERROR: missing reference folder {ref_dir}", file=sys.stderr)
            return 1
        save_model(build_background_model(ref_dir), model_path)

    pairs = discover_pairs(photo_dir, gcode_dir)
    if args.limit:
        pairs = pairs[:args.limit]
    if not pairs:
        print(f"No photo pairs found in {photo_dir}", file=sys.stderr)
        return 1

    csv_path = out_dir / "master_printability.csv"
    cache = CacheIndex(out_dir / ".cache_index.json")

    sam2_proc, sam2_client = None, None
    if args.use_sam2:
        from sam2_service import run_service, SAM2Client
        req_q, resp_q = mp.Queue(), mp.Queue()
        ckpt = PROJECT_ROOT / "Analises/Python/Printability/sam2_checkpoints/sam2_hiera_large.pt"
        if ckpt.exists():
            sam2_proc = mp.Process(target=run_service, args=(req_q, resp_q, ckpt))
            sam2_proc.start()
            sam2_client = SAM2Client(req_q, resp_q)
        else:
            print(f"WARN: SAM2 checkpoint missing at {ckpt}; running without SAM2", file=sys.stderr)

    try:
        if args.workers <= 1:
            # Serial path — also useful for testing
            for pair in pairs:
                key = f"Vp{int(pair.Vp)}_Flow{pair.Flow}"
                paths = {"photo_flash": pair.flash, "photo_noflash": pair.noflash,
                         "gcode": pair.gcode, "background_model": model_path}
                if not args.force and cache.is_fresh(key, paths):
                    print(f"  SKIP {key} (cached)")
                    continue
                row = process_pair(pair, model_path, sam2_client, args.bg)
                write_row(csv_path, row)
                h = hashlib.sha1(json.dumps(row, default=str, sort_keys=True).encode()).hexdigest()
                cache.update(key, paths, metrics_hash=h)
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = []
                for pair in pairs:
                    key = f"Vp{int(pair.Vp)}_Flow{pair.Flow}"
                    paths = {"photo_flash": pair.flash, "photo_noflash": pair.noflash,
                             "gcode": pair.gcode, "background_model": model_path}
                    if not args.force and cache.is_fresh(key, paths):
                        print(f"  SKIP {key} (cached)")
                        continue
                    futures.append((key, paths, pool.submit(process_pair, pair, model_path, sam2_client, args.bg)))
                for key, paths, fut in futures:
                    row = fut.result()
                    write_row(csv_path, row)
                    h = hashlib.sha1(json.dumps(row, default=str, sort_keys=True).encode()).hexdigest()
                    cache.update(key, paths, metrics_hash=h)
        cache.save()
    finally:
        if sam2_proc is not None:
            req_q.put(None)
            sam2_proc.join(timeout=5)

    print(f"OK: wrote {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
