"""Interactive ROI editor for foreground photos.

Loads a foreground photo and (optionally) an auto-detected ROI sidecar
JSON, then lets the user adjust on a matplotlib canvas:

  - Drag a rectangle to add a NEW ROI.
  - Click an existing rectangle to select it (red highlight).
  - Press 'd' to delete the selected ROI.
  - Press 'l' to set the next-drawn rectangle as the logo (replaces existing).
  - Press 'p' to set the next-drawn rectangle as a label.
  - Press 'o' to toggle logo orientation (horizontal <-> vertical).
  - Press 's' to save the sidecar JSON.
  - Press 'q' to quit (prompts to save unsaved edits).

Usage:
  python roi_editor.py <foreground.HEIC> [--rois path/to/rois.json]

If --rois is omitted, auto-detection runs first to seed proposals. The
edited result is saved to <photo>.rois.json next to the photo unless
--out is provided.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RectangleSelector

from bed_calibration import bed_roi_grabcut, load_rgb
from roi_detection import (
    LOGO_LONG_MM,
    LOGO_SHORT_MM,
    LabelROI,
    LogoROI,
    PhotoROIs,
    detect_labels,
    detect_logo,
    load_rois,
    save_rois,
)


class ROIEditor:
    def __init__(self, photo_path: Path, rois: PhotoROIs, out_path: Path):
        self.photo_path = photo_path
        self.rois = rois
        self.out_path = out_path
        self.dirty = False
        self.next_kind = "label"  # 'logo' or 'label'

        self.rgb = load_rgb(photo_path)
        self.fig, self.ax = plt.subplots(figsize=(12, 9))
        self.ax.imshow(self.rgb)
        self.ax.set_title(self._title_text())
        self.ax.axis("off")

        self.patches: list[tuple[mpatches.Rectangle, str, int]] = []  # (patch, kind, idx)
        self.selected: Optional[int] = None

        self._draw_existing()

        self.selector = RectangleSelector(
            self.ax,
            self._on_select,
            useblit=True,
            button=[1],
            minspanx=10, minspany=10,
            spancoords="pixels",
            interactive=True,
        )
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

        print(self.__doc__)
        plt.show()

    def _title_text(self) -> str:
        logo = self.rois.logo
        n = len(self.rois.labels)
        kind = self.next_kind.upper()
        marker = " *" if self.dirty else ""
        scale = (f"{logo.mm_per_px:.4f} mm/px" if logo else "no logo")
        ori = logo.orientation if logo else "-"
        return (f"{self.photo_path.name}  | next: {kind} | logo: {ori} {scale} | "
                f"labels: {n}{marker}\n"
                f"[L]=logo [P]=label [O]=toggle ori [D]=delete sel [S]=save [Q]=quit")

    def _refresh_title(self):
        self.ax.set_title(self._title_text())
        self.fig.canvas.draw_idle()

    def _draw_existing(self):
        if self.rois.logo:
            self._add_patch_visual(self.rois.logo.bbox_xyxy, "logo")
        for idx, lab in enumerate(self.rois.labels):
            self._add_patch_visual(lab.bbox_xyxy, "label", idx)

    def _add_patch_visual(self, bbox, kind: str, idx: int = -1):
        x0, y0, x1, y1 = bbox
        colour = "red" if kind == "logo" else "yellow"
        rect = mpatches.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                  linewidth=2, edgecolor=colour, facecolor="none")
        self.ax.add_patch(rect)
        self.patches.append((rect, kind, idx))

    def _redraw_all_patches(self):
        for rect, _, _ in self.patches:
            rect.remove()
        self.patches.clear()
        self._draw_existing()
        self.fig.canvas.draw_idle()

    def _on_select(self, eclick, erelease):
        x0, y0 = int(min(eclick.xdata, erelease.xdata)), int(min(eclick.ydata, erelease.ydata))
        x1, y1 = int(max(eclick.xdata, erelease.xdata)), int(max(eclick.ydata, erelease.ydata))
        if x1 - x0 < 10 or y1 - y0 < 10:
            return
        bbox = [x0, y0, x1, y1]
        if self.next_kind == "logo":
            long_px = float(max(x1 - x0, y1 - y0))
            short_px = float(min(x1 - x0, y1 - y0))
            mm_per_px = 0.5 * (LOGO_LONG_MM / long_px + LOGO_SHORT_MM / short_px)
            self.rois.logo = LogoROI(
                bbox_xyxy=bbox,
                orientation="horizontal" if (x1 - x0) >= (y1 - y0) else "vertical",
                long_axis_mm=LOGO_LONG_MM,
                short_axis_mm=LOGO_SHORT_MM,
                mm_per_px=mm_per_px,
                px_per_mm=1.0 / mm_per_px,
                scale_residual_pct=100.0 * abs(LOGO_LONG_MM / long_px
                                              - LOGO_SHORT_MM / short_px) / mm_per_px,
                auto=False,
            )
        else:
            self.rois.labels.append(LabelROI(bbox_xyxy=bbox, auto=False))
        self.dirty = True
        self._redraw_all_patches()
        self._refresh_title()

    def _on_click(self, event):
        if event.inaxes != self.ax or event.button != 3:  # right-click selects
            return
        for i, (rect, _, _) in enumerate(self.patches):
            x, y = rect.get_xy()
            w, h = rect.get_width(), rect.get_height()
            if x <= event.xdata <= x + w and y <= event.ydata <= y + h:
                for j, (r2, _, _) in enumerate(self.patches):
                    r2.set_edgecolor("red" if self.patches[j][1] == "logo" else "yellow")
                rect.set_edgecolor("magenta")
                self.selected = i
                self.fig.canvas.draw_idle()
                return
        self.selected = None
        self.fig.canvas.draw_idle()

    def _on_key(self, event):
        if event.key == "l":
            self.next_kind = "logo"
        elif event.key == "p":
            self.next_kind = "label"
        elif event.key == "o" and self.rois.logo:
            current = self.rois.logo.orientation
            self.rois.logo.orientation = "vertical" if current == "horizontal" else "horizontal"
            self.rois.logo.auto = False
            self.dirty = True
        elif event.key == "d" and self.selected is not None:
            _, kind, idx = self.patches[self.selected]
            if kind == "logo":
                self.rois.logo = None
            else:
                self.rois.labels.pop(idx)
            self.selected = None
            self.dirty = True
            self._redraw_all_patches()
        elif event.key == "s":
            save_rois(self.rois, self.out_path)
            self.dirty = False
            print(f"saved -> {self.out_path}")
        elif event.key == "q":
            if self.dirty:
                resp = input("unsaved edits — save before quit? [y/n] ").strip().lower()
                if resp == "y":
                    save_rois(self.rois, self.out_path)
                    print(f"saved -> {self.out_path}")
            plt.close(self.fig)
            return
        self._refresh_title()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("photo", type=Path)
    ap.add_argument("--rois", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-auto", action="store_true",
                    help="skip auto-detection seeding")
    args = ap.parse_args()

    if args.rois and args.rois.exists():
        rois = load_rois(args.rois)
        print(f"loaded existing ROIs from {args.rois}")
    elif args.no_auto:
        rgb = load_rgb(args.photo)
        rois = PhotoROIs(source=args.photo.name,
                         image_shape=(rgb.shape[0], rgb.shape[1]),
                         logo=None, labels=[])
    else:
        print(f"auto-detecting ROIs in {args.photo.name} ...")
        rgb = load_rgb(args.photo)
        bed = bed_roi_grabcut(rgb)
        rois = PhotoROIs(
            source=args.photo.name,
            image_shape=(rgb.shape[0], rgb.shape[1]),
            logo=detect_logo(rgb, bed),
            labels=detect_labels(rgb, bed),
        )
        print(f"  logo: {'detected' if rois.logo else 'NOT detected (please draw manually)'}")
        print(f"  labels auto-detected: {len(rois.labels)}")

    out_path = args.out or args.photo.with_suffix(args.photo.suffix + ".rois.json")
    ROIEditor(args.photo, rois, out_path)


if __name__ == "__main__":
    main()
