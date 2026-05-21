import numpy as np
from registration import align_pair


def test_align_pair_recovers_known_shift():
    rng = np.random.default_rng(42)
    base = rng.integers(0, 255, (200, 200, 3), dtype=np.uint8)
    shifted = np.roll(base, shift=(5, -3), axis=(0, 1))
    aligned, (dy, dx) = align_pair(base, shifted)
    assert abs(dy - 5) <= 1
    assert abs(dx + 3) <= 1
    diff = np.abs(aligned[20:180, 20:180].astype(int) - base[20:180, 20:180].astype(int))
    assert diff.mean() < 5
