"""
data_config_local.py  (TEMPLATE — copy to your PROJECT ROOT, next to inks_local.m)
=================================================================================

Project-local data configuration for the shared Python scripts under
``Export/01_Python/``. This is the Python counterpart of ``inks_local.m``:
the shared code stays generic, the project-specific names live here.

Nothing in this file is required. Anything you leave out falls back to the
built-in default in ``Export/01_Python/data_config.py``, which reproduces
the folder names the scripts used before this mechanism existed.

RULE: do NOT edit folder names back into the shared scripts. That is a
fork, and forks drift.
"""

# ----------------------------------------------------------------------
# Where the raw rheology lives, per test type. Relative to the project root
# (scripts are always run from the project root, never from 01_Python/).
# ----------------------------------------------------------------------
FOLDERS = {
    "flow":      "./Reologia/Viscosity",
    "amplitude": "./Reologia/Strain Sweep - 3D Bioprinting",
    "frequency": "./Reologia/Frequency Sweep - 3D Bioprinting",
    "recovery":  "./Reologia/Recovery - 3D Bioprinting",
    "force":     "./Reologia/Flow with force",
}

# Where reports and CSVs are written.
RESULTS_DIR = "./Analises/Python/Results"

# ----------------------------------------------------------------------
# Restrict to a subset of samples. None = every raw file in the folder.
# Entries match the SHORT filename prefix (the text before " - ").
# ----------------------------------------------------------------------
SAMPLE_INCLUDE = None
# e.g. SAMPLE_INCLUDE = ["C10", "C15", "C20", "C25", "gel"]

# ----------------------------------------------------------------------
# Files to drop before anything reads them. Substring or glob. Use this for
# blanks, aborted runs, and anything the lab notebook says to ignore — so
# the reason is recorded once instead of being re-derived every run.
# ----------------------------------------------------------------------
EXCLUDE_PATTERNS = [
    # "NEG",      # negative / blank runs
    # "Amb",      # ambient-only runs
    # "*_ABORT*",
]

# ----------------------------------------------------------------------
# Measuring system this analysis is valid for. Files naming a DIFFERENT
# known geometry are dropped with a warning; files naming NO geometry raise,
# because silently mixing a PP50 parallel-plate run into a CP50 analysis
# gives an answer that looks self-consistent and is physically wrong.
# Set it to None only if your filenames genuinely never carry the geometry
# — and then accept that the guard cannot help you.
# ----------------------------------------------------------------------
GEOMETRY_FILTER = None
# e.g. GEOMETRY_FILTER = "CP50"
