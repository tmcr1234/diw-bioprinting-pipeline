"""SAM2 inference service for the printability pipeline.

PyTorch MPS is not fork-safe, so SAM2 cannot be loaded inside a
ProcessPoolExecutor worker. Instead, one dedicated process holds the model
and answers predict() requests via multiprocessing.Queue. Workers hold
a SAM2Client that wraps the queue.

Usage (pipeline.py):
    req_q, resp_q = mp.Queue(), mp.Queue()
    p = mp.Process(target=run_service,
                   args=(req_q, resp_q, "sam2_checkpoints/sam2_hiera_large.pt"))
    p.start()
    client = SAM2Client(req_q, resp_q)
"""
from __future__ import annotations
import uuid
import numpy as np


class SAM2Service:
    def __init__(self, checkpoint, device: str = "mps"):
        import torch
        from sam2.build_sam import build_sam2
        from sam2.automatic_mask_generator import SAM2AutomaticMaskGenerator
        self.device = device
        try:
            self.model = build_sam2("sam2_hiera_l.yaml", checkpoint, device=device)
        except Exception:
            self.device = "cpu"
            self.model = build_sam2("sam2_hiera_l.yaml", checkpoint, device="cpu")
        self.gen = SAM2AutomaticMaskGenerator(self.model, points_per_side=32)

    def predict(self, img: np.ndarray) -> np.ndarray:
        H, W = img.shape[:2]
        try:
            results = self.gen.generate(img)
        except Exception:
            return np.zeros((H, W), dtype=bool)
        cy, cx = H / 2, W / 2
        diag = (H * H + W * W) ** 0.5
        best, best_score = None, -1.0
        for r in results:
            area = r["area"]
            if not (0.05 * H * W <= area <= 0.70 * H * W):
                continue
            y, x = np.argwhere(r["segmentation"]).mean(axis=0)
            dist = ((y - cy) ** 2 + (x - cx) ** 2) ** 0.5
            if dist > 0.30 * diag:
                continue
            score = r.get("predicted_iou", 1.0) - 0.5 * (dist / diag)
            if score > best_score:
                best_score = score
                best = r["segmentation"]
        if best is None:
            return np.zeros((H, W), dtype=bool)
        return best.astype(bool)


class SAM2Client:
    def __init__(self, req_q, resp_q):
        self.req_q, self.resp_q = req_q, resp_q

    def predict(self, img: np.ndarray) -> np.ndarray:
        tag = uuid.uuid4().hex
        self.req_q.put((tag, img))
        while True:
            t, mask = self.resp_q.get()
            if t == tag:
                return mask
            self.resp_q.put((t, mask))


def run_service(req_q, resp_q, checkpoint):
    svc = SAM2Service(checkpoint)
    while True:
        item = req_q.get()
        if item is None:
            return
        tag, img = item
        try:
            mask = svc.predict(img)
        except Exception:
            mask = np.zeros(img.shape[:2], dtype=bool)
        resp_q.put((tag, mask))
