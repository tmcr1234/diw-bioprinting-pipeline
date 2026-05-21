import numpy as np
from segmentation.family6_sam2 import mask_sam2


def test_mask_sam2_returns_empty_mask_when_service_none():
    img = np.zeros((100, 100, 3), dtype=np.uint8)
    mask = mask_sam2(img, "Black", service=None)
    assert mask.shape == (100, 100)
    assert mask.dtype == bool
    assert mask.sum() == 0
