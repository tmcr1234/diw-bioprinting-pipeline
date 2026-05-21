import pytest
from tests.fixtures.synthetic import make_synthetic
from segmentation.family3_local_thr import mask_adaptive_gauss, mask_sauvola, mask_niblack

METHODS = [mask_adaptive_gauss, mask_sauvola, mask_niblack]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("bg", ["Black", "Blue"])
def test_local_threshold_produces_plausible_mask(method, bg):
    """Wiring check only. Local thresholds are tuned for textured photos and
    perform unevenly on uniform synthetic backgrounds; their real validation
    comes from the multi-method consensus on real S-Test photos. Here we just
    confirm the method runs, produces a binary mask of the correct shape, and
    is non-degenerate (some but not all pixels foreground)."""
    img, gt = make_synthetic(bg)
    mask = method(img, bg)
    assert mask.dtype == bool
    assert mask.shape == gt.shape
    fg = mask.sum()
    assert 0 < fg < mask.size, f"{method.__name__} on {bg}: degenerate mask, fg={fg}"
