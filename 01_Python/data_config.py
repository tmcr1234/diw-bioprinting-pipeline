"""
data_config.py
==============

Project-local data configuration for the shared Python layer — the Python
equivalent of ``inks_local.m``.

THE PROBLEM THIS SOLVES
-----------------------
``Export/`` is a single shared clone symlinked into every project. Until
now the Python scripts hard-coded their input folders and ink lists as
module-level literals (``INPUT_FOLDER = "./Reologia/Viscosity"``), so a
project whose folders were named differently, or that carried a different
set of inks, had to FORK the script. A fork that changes nothing but a
folder name and a list is the exact mechanism by which two copies of an
analysis silently drift apart. MATLAB solved this with ``inks_local.m``;
this is the missing Python half.

HOW IT WORKS
------------
Each project keeps a ``data_config_local.py`` in its own root (next to
``inks_local.m``), and the shared scripts call::

    from data_config import load
    cfg = load()
    folder = cfg.folder("flow")

``load()`` searches the current working directory first — scripts are run
from the project root, per the pipeline rule — then falls back to built-in
defaults that reproduce today's hard-coded literals exactly. A project that
does nothing therefore behaves exactly as before.

Copy ``data_config_local.template.py`` to your project root as
``data_config_local.py`` to start.

Author: T.M.C. Rodrigues - 2026-08-26
"""

from __future__ import annotations

import fnmatch
import importlib.util
import sys
from dataclasses import dataclass, field
from pathlib import Path

# The literals the shared scripts used before this module existed. Keeping
# them here — and ONLY here — means "no local config" is a defined state,
# not an accident.
_DEFAULT_FOLDERS = {
    "flow":      "./Reologia/Viscosity",
    "amplitude": "./Reologia/Strain Sweep - 3D Bioprinting",
    "frequency": "./Reologia/Frequency Sweep - 3D Bioprinting",
    "recovery":  "./Reologia/Recovery - 3D Bioprinting",
    "force":     "./Reologia/Flow with force",
}

_DEFAULT_RESULTS = "./Analises/Python/Results"

_CONFIG_FILENAME = "data_config_local.py"


@dataclass
class DataConfig:
    """Resolved data configuration for one project."""
    folders: dict = field(default_factory=lambda: dict(_DEFAULT_FOLDERS))
    results_dir: str = _DEFAULT_RESULTS
    sample_include: list | None = None
    exclude_patterns: list = field(default_factory=list)
    geometry_filter: str | None = None
    source: str = "built-in defaults"

    def folder(self, kind: str) -> str:
        """Folder for one test type. Raises on an unknown kind rather than
        returning something plausible-looking."""
        if kind not in self.folders:
            raise KeyError(
                f"unknown data folder kind {kind!r}. Known: "
                f"{sorted(self.folders)}. Add it to {_CONFIG_FILENAME}."
            )
        return self.folders[kind]

    def describe(self) -> str:
        """One block, for a report header — so a result carries its own
        provenance instead of relying on someone remembering the run."""
        lines = [f"Data configuration: {self.source}"]
        for k in sorted(self.folders):
            lines.append(f"  {k:<10} -> {self.folders[k]}")
        lines.append(f"  results    -> {self.results_dir}")
        if self.sample_include:
            lines.append(f"  include    -> {self.sample_include}")
        if self.exclude_patterns:
            lines.append(f"  exclude    -> {self.exclude_patterns}")
        lines.append(f"  geometry   -> {self.geometry_filter or '(no filter — see A13)'}")
        return "\n".join(lines)


def load(start: str | Path | None = None) -> DataConfig:
    """
    Load ``data_config_local.py`` from the project root, or fall back to the
    built-in defaults.

    Searches `start` (default: cwd) and then its parents, so running from a
    subfolder still finds the project's config.
    """
    here = Path(start) if start else Path.cwd()
    for candidate in [here, *here.parents]:
        cfg_path = candidate / _CONFIG_FILENAME
        if cfg_path.is_file():
            return _load_from(cfg_path)
    return DataConfig()


def _load_from(cfg_path: Path) -> DataConfig:
    spec = importlib.util.spec_from_file_location("data_config_local", cfg_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {cfg_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["data_config_local"] = mod
    spec.loader.exec_module(mod)

    folders = dict(_DEFAULT_FOLDERS)
    folders.update(getattr(mod, "FOLDERS", {}) or {})

    return DataConfig(
        folders=folders,
        results_dir=getattr(mod, "RESULTS_DIR", _DEFAULT_RESULTS),
        sample_include=getattr(mod, "SAMPLE_INCLUDE", None),
        exclude_patterns=list(getattr(mod, "EXCLUDE_PATTERNS", []) or []),
        geometry_filter=getattr(mod, "GEOMETRY_FILTER", None),
        source=str(cfg_path),
    )


# ----------------------------------------------------------------------
# A13 — manifest-based filtering
# ----------------------------------------------------------------------

class UnclassifiedFileError(RuntimeError):
    """Raised when a file matches neither an include rule nor an explicit
    exclude rule. Failing here is the entire point: the alternative is a
    PP50 parallel-plate run quietly entering a CP50-only analysis and
    producing a self-consistent, physically wrong answer."""


def filter_files(paths, cfg: DataConfig, strict: bool = True) -> list[Path]:
    """
    Apply the project's exclude patterns and geometry filter to a list of
    candidate files.

    Rules, in order:
      1. Anything matching ``EXCLUDE_PATTERNS`` is dropped, silently — it
         was explicitly declared unwanted (NEG / Amb blanks, aborted runs).
      2. If ``GEOMETRY_FILTER`` is set, a file whose name names THAT
         geometry is kept, and a file naming a DIFFERENT known geometry is
         dropped with a warning.
      3. A file that names no geometry at all is unclassified. With
         ``strict=True`` (the default) this RAISES. Mixing measuring systems
         is not a formatting inconvenience; it changes the physics.

    Set ``strict=False`` only for exploratory listing, never for a run whose
    numbers will be reported.
    """
    known_geometries = ("CP50", "PP50", "CP25", "PP25", "CC27", "CP60")
    kept, dropped, unclassified = [], [], []

    for p in (Path(x) for x in paths):
        name = p.name
        if any(fnmatch.fnmatch(name, pat) or pat.lower() in name.lower()
               for pat in cfg.exclude_patterns):
            dropped.append((p, "excluded by pattern"))
            continue

        if cfg.geometry_filter:
            up = name.upper()
            wanted = cfg.geometry_filter.upper()
            if wanted in up:
                kept.append(p)
            elif any(g in up for g in known_geometries):
                other = next(g for g in known_geometries if g in up)
                dropped.append((p, f"geometry {other} != {cfg.geometry_filter}"))
            else:
                unclassified.append(p)
        else:
            kept.append(p)

    if cfg.geometry_filter is None and paths:
        print(
            "  [!] No GEOMETRY_FILTER set. Files from different measuring "
            "systems will be mixed without complaint, which gives a "
            "self-consistent but physically wrong answer. Set "
            f"GEOMETRY_FILTER in {_CONFIG_FILENAME}."
        )

    if unclassified:
        msg = (
            f"{len(unclassified)} file(s) name no known measuring system, so "
            f"they can be neither included nor excluded under "
            f"GEOMETRY_FILTER={cfg.geometry_filter!r}:\n  "
            + "\n  ".join(str(p.name) for p in unclassified[:10])
            + f"\nEither rename them to carry the geometry, or add them to "
              f"EXCLUDE_PATTERNS in {_CONFIG_FILENAME}."
        )
        if strict:
            raise UnclassifiedFileError(msg)
        print(f"  [!] {msg}")
        kept.extend(unclassified)

    for p, why in dropped:
        print(f"  [skip] {p.name}  ({why})")

    return kept
