import mahotas as mh
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import pandas as pd  # noqa


def test_rgb2gray():
    img2 = npia.Image('PtnanoCenk0.5et_171218_AM.0_00002_1.png', 3)
    img2.image = np.array([[[255, 231, 97], [45, 23, 18]]])
    test_output = img2.rgb2gray()
    assert np.isclose(test_output[0, 0], 222.9), 'bad grayscale'
    assert np.isclose(test_output[0, 1], 29.008)
    return

def test_apply_threshold():
    img2 = npia.Image('PtnanoCenk0.5et_171218_AM.0_00002_1.png', 3)
    img2.image = np.array([[[300, 231, 97], [45, 23, 18]]])
    img2.apply_threshold()
    test_output2 = img2.image
    assert np.isclose(test_output2[0, 0, 0], 255), 'bad thresholding'
    assert np.isclose(test_output2[0, 1, 2], 12)
    return