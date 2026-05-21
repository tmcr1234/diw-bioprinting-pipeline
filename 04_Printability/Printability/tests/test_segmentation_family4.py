import pytest
from tests.fixtures.synthetic import make_synthetic
from segmentation.family4_edge_region import mask_canny, mask_felzenszwalb

METHODS = [mask_canny, mask_felzenszwalb]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("bg", ["Black", "Blue"])
def test_edge_region_produces_plausible_mask(method, bg):
    img, gt = make_synthetic(bg)
    mask = method(img, bg)
    assert mask.dtype == bool
    assert mask.shape == gt.shape
    fg = mask.sum()
    assert 0 < fg < mask.size, f"{method.__name__} on {bg}: degenerate mask, fg={fg}"
    iou = (mask & gt).sum() / (mask | gt).sum()
    assert iou > 0.40, f"{method.__name__} on {bg}: IoU={iou:.3f}"
