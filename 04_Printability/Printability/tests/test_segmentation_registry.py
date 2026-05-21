import numpy as np
from segmentation import METHODS, list_methods, run_method
from tests.fixtures.synthetic import make_synthetic


def test_methods_registered():
    names = list_methods(include_sam2=False)
    expected_min = {"HSV", "Lab_ab", "kmeans_lab", "kmeans_rgb",
                    "Otsu", "Triangle", "Li", "Yen",
                    "adaptiveGauss", "Sauvola", "Niblack",
                    "Canny", "Felzenszwalb"}
    assert expected_min.issubset(set(names))


def test_run_method_dispatches_with_correct_channel():
    img, gt = make_synthetic("Black")
    mask = run_method("HSV", flash=img, noflash=img, bg="Black")
    assert mask.shape == gt.shape
    assert mask.dtype == bool
