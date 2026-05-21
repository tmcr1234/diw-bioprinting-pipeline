import pytest
from tests.fixtures.synthetic import make_synthetic
from segmentation.family1_colour import mask_hsv, mask_lab_ab, mask_kmeans_lab, mask_kmeans_rgb

METHODS = [mask_hsv, mask_lab_ab, mask_kmeans_lab, mask_kmeans_rgb]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("bg", ["Black", "Blue"])
def test_method_recovers_synthetic_square(method, bg):
    img, gt = make_synthetic(bg)
    mask = method(img, bg)
    assert mask.dtype == bool
    assert mask.shape == gt.shape
    iou = (mask & gt).sum() / (mask | gt).sum()
    assert iou > 0.85, f"{method.__name__} on {bg}: IoU={iou:.3f}"
