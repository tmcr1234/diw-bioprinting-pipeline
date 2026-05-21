"""mtime-keyed cache index for the pipeline."""
from __future__ import annotations
import json
from pathlib import Path


class CacheIndex:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.entries: dict[str, dict] = {}
        if self.path.exists():
            self.entries = json.loads(self.path.read_text())

    def is_fresh(self, key: str, paths: dict[str, Path]) -> bool:
        if key not in self.entries:
            return False
        e = self.entries[key]
        for name, p in paths.items():
            mtime_key = f"{name}_mtime"
            if mtime_key not in e:
                return False
            if not Path(p).exists():
                return False
            if abs(Path(p).stat().st_mtime - e[mtime_key]) > 1e-3:
                return False
        return True

    def update(self, key: str, paths: dict[str, Path], metrics_hash: str) -> None:
        rec = {f"{n}_mtime": Path(p).stat().st_mtime for n, p in paths.items() if Path(p).exists()}
        rec["metrics_hash"] = metrics_hash
        self.entries[key] = rec

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.entries, indent=2, sort_keys=True))
