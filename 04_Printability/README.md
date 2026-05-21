# Printability Module

Image-based DIW print-quality analysis pipeline. Self-contained Python
package with its own CLI, test suite, and segmentation backends.

## When to use this

After Stage 4 of the pipeline (printing). Photograph each printed shape with
your phone (one with flash, one without) and feed the images through this
package to compute the **Ouyang Printability index** (Pr) and related
strand/pore/turn metrics.

## Quick start

```bash
cd Export/04_Printability/Printability
pip install -e .
pytest                          # 66 tests should pass
python -m diw_printability --help
```

## Per-shape benchmark suite

The expected directory structure for input images is the **per-shape suite**:

```
photos/
├── Cubes/         <ink_name>/<image>.HEIC
├── Area-Cube/
├── S-Test/
├── Tower/
├── Scaffold/
├── Bridge/
└── Misc/
```

Two iPhone photos per print: one with flash, one without. The pipeline
chooses the better one automatically via the `consensus` module.

## Segmentation methods (14)

Six families:

1. `family1_colour` — HSV/Lab thresholding
2. `family2_global_thr` — Otsu, Yen, Li
3. `family3_local_thr` — Sauvola, Niblack
4. `family4_edge_region` — Canny + watershed
5. `family5_bgdiff` — background subtraction (needs paired flash/no-flash)
6. `family6_sam2` — SAM2 (placeholder — install separately)

## Outputs

- Per-image segmentation overlays
- Per-shape aggregate metrics (mean Pr, strand width CV, pore-area CV, etc.)
- CSV summary across the whole suite

## Operational SOP

See `Export/03_SOPs/Printing_Parameters_SOP_v4_EN.pdf` for the bench-side
printing protocol that produces the inputs to this pipeline.

---

## Status

This is the most recently developed module in the Export. It is functional
(22 commits, 66 passing tests at time of export) but expect more rapid
iteration than the rheology layer. Check the package's own `CHANGELOG` if it
exists; otherwise consult the maintainer.
