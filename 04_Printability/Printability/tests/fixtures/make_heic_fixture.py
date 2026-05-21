"""One-shot script to write tests/fixtures/tiny.heic from a synthetic numpy array.
Run once; the resulting .heic is committed."""
import numpy as np
from PIL import Image
import pillow_heif
pillow_heif.register_heif_opener()


def main():
    rng = np.random.default_rng(0)
    img = (rng.uniform(0, 255, size=(64, 64, 3))).astype(np.uint8)
    Image.fromarray(img).save("tests/fixtures/tiny.heic", format="HEIF")


if __name__ == "__main__":
    main()
