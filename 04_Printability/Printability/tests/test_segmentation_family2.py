import pytest
from tests.fixtures.synthetic import make_synthetic
from segmentation.family2_global_thr import mask_otsu, mask_triangle, mask_li, mask_yen

METHODS = [mask_otsu, mask_triangle, mask_li, mask_yen]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("bg", ["Black", "Blue"])
def test_global_threshold_recovers_square(method, bg):
    img, gt = make_synthetic(bg)
    mask = method(img, bg)
    iou = (mask & gt).sum() / (mask | gt).sum()
    assert iou > 0.80, f"{method.__name__} on {bg}: IoU={iou:.3f}"
