import time
from pathlib import Path
from cache import CacheIndex


def test_fresh_returns_false_for_unknown_key(tmp_path):
    idx = CacheIndex(tmp_path / "cache.json")
    assert not idx.is_fresh("k", {"a": tmp_path / "missing"})


def test_fresh_true_after_update(tmp_path):
    p = tmp_path / "file.txt"
    p.write_text("x")
    idx = CacheIndex(tmp_path / "cache.json")
    idx.update("k", {"a": p}, metrics_hash="abc")
    assert idx.is_fresh("k", {"a": p})


def test_fresh_false_after_file_touch(tmp_path):
    p = tmp_path / "file.txt"
    p.write_text("x")
    idx = CacheIndex(tmp_path / "cache.json")
    idx.update("k", {"a": p}, metrics_hash="abc")
    time.sleep(0.05)
    p.write_text("y")
    assert not idx.is_fresh("k", {"a": p})


def test_save_and_reload(tmp_path):
    p = tmp_path / "file.txt"
    p.write_text("x")
    idx = CacheIndex(tmp_path / "cache.json")
    idx.update("k", {"a": p}, metrics_hash="abc")
    idx.save()
    idx2 = CacheIndex(tmp_path / "cache.json")
    assert idx2.is_fresh("k", {"a": p})
