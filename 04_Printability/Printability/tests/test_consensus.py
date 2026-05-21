import numpy as np
from consensus import compute_consensus, is_degenerate


def test_consensus_unanimous():
    mask = np.zeros((50, 50), dtype=bool)
    mask[10:40, 10:40] = True
    masks = {"a": mask, "b": mask, "c": mask}
    consensus, agree = compute_consensus(masks)
    np.testing.assert_array_equal(consensus, mask)
    assert agree == 1.0


def test_consensus_majority_wins():
    big = np.zeros((50, 50), dtype=bool)
    big[10:40, 10:40] = True
    small = np.zeros((50, 50), dtype=bool)
    small[20:30, 20:30] = True
    masks = {"a": big, "b": big, "c": small}
    consensus, _ = compute_consensus(masks)
    assert consensus[15, 15] == True


def test_is_degenerate_flags_empty():
    assert is_degenerate(np.zeros((10, 10), dtype=bool))


def test_is_degenerate_flags_all_foreground():
    assert is_degenerate(np.ones((10, 10), dtype=bool))


def test_consensus_drops_degenerate_masks():
    good = np.zeros((50, 50), dtype=bool)
    good[10:40, 10:40] = True
    empty = np.zeros((50, 50), dtype=bool)
    masks = {"a": good, "b": good, "c": empty}
    consensus, _ = compute_consensus(masks)
    np.testing.assert_array_equal(consensus, good)
