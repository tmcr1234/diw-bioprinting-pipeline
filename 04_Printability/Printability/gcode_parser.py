"""Parse PrusaSlicer G-code into a RasterSpec describing strand geometry.

Extrusion segments (G1 with positive E-delta) are grouped into strands
by constant Y. Travel moves are skipped. The filename token Vp<N> and
Fr<X>/Flow<X> populate Vp_mm_s and Flow.
"""
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class LineSegment:
    x0_mm: float
    y0_mm: float
    x1_mm: float
    y1_mm: float


@dataclass
class TurnPoint:
    x_mm: float
    y_mm: float
    in_strand_idx: int
    out_strand_idx: int


@dataclass
class RasterSpec:
    strands: list[LineSegment]
    turns: list[TurnPoint]
    strand_spacing_mm: float
    strand_length_mm: float
    bounding_box_mm: tuple[float, float, float, float]
    Vp_mm_s: float
    Flow: float
    source_gcode: Path


_VP_RE = re.compile(r"Vp(\d+(?:[.,]\d+)?)", re.IGNORECASE)
_FLOW_RE = re.compile(r"(?:Fr|Flow)(\d+(?:[.,]\d+)?)", re.IGNORECASE)


def parse_filename(path: Path) -> tuple[float, float]:
    """Extract (Vp, Flow) from a filename. Returns (0.0, 0.0) if absent."""
    name = path.name
    vp = _VP_RE.search(name)
    fl = _FLOW_RE.search(name)
    to_f = lambda m: float(m.group(1).replace(",", ".")) if m else 0.0
    return to_f(vp), to_f(fl)


def parse_gcode(path: Path) -> RasterSpec:
    Vp, Flow = parse_filename(path)
    x = y = e = 0.0
    segments: list[tuple[float, float, float, float]] = []

    for line in path.read_text().splitlines():
        if not line.startswith("G1"):
            continue
        new_x = x
        new_y = y
        new_e = e
        for tok in line.split():
            if tok.startswith("X"):
                new_x = float(tok[1:].replace(",", "."))
            elif tok.startswith("Y"):
                new_y = float(tok[1:].replace(",", "."))
            elif tok.startswith("E"):
                new_e = float(tok[1:].replace(",", "."))
        if new_e > e + 1e-9:
            segments.append((x, y, new_x, new_y))
        x, y, e = new_x, new_y, new_e

    if not segments:
        raise ValueError(f"No extrusion segments found in {path}")

    strands: list[LineSegment] = []
    turns: list[TurnPoint] = []
    current_y = segments[0][1]
    strand_start = (segments[0][0], segments[0][1])
    strand_end = (segments[0][2], segments[0][3])
    for x0, y0, x1, y1 in segments[1:]:
        if abs(y1 - current_y) < 0.05 and abs(y0 - current_y) < 0.05:
            strand_end = (x1, y1)
        else:
            strands.append(LineSegment(strand_start[0], current_y, strand_end[0], current_y))
            turns.append(TurnPoint(x0, (current_y + y1) / 2,
                                   in_strand_idx=len(strands) - 1,
                                   out_strand_idx=len(strands)))
            current_y = y1
            strand_start = (x0, y1)
            strand_end = (x1, y1)
    strands.append(LineSegment(strand_start[0], current_y, strand_end[0], current_y))

    ys = sorted(set(round(s.y0_mm, 3) for s in strands))
    spacings = [ys[i + 1] - ys[i] for i in range(len(ys) - 1)] or [0.0]
    strand_spacing = sum(spacings) / len(spacings)
    strand_length = sum(abs(s.x1_mm - s.x0_mm) for s in strands) / len(strands)

    xs = [s.x0_mm for s in strands] + [s.x1_mm for s in strands]
    ys_all = [s.y0_mm for s in strands] + [s.y1_mm for s in strands]
    bbox = (min(xs), min(ys_all), max(xs), max(ys_all))

    return RasterSpec(strands, turns, strand_spacing, strand_length,
                      bbox, Vp, Flow, path)
