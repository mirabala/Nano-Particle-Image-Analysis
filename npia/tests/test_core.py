from os.path import abspath, dirname, join

import mahotas as mh  # noqa
import matplotlib.pyplot as plt  # noqa
import numpy as np

from ..core import Image


def test_rgb2gray():
    default = abspath(join(dirname(__file__),
                           'PtnanoCenk0.5et_171218_AM.0_00002_1.png'))
    img2 = Image(default, 3)
    img2.image = np.array([[[255, 231, 97], [45, 23, 18]]])
    test_output = img2.rgb2gray()
    assert np.isclose(test_output[0, 0], 222.9), 'bad grayscale'
    assert np.isclose(test_output[0, 1], 29.008)
    return


def test_apply_threshold():
    default = abspath(join(dirname(__file__),
                           'PtnanoCenk0.5et_171218_AM.0_00002_1.png'))
    img2 = Image(default, 3)
    img2.image = np.array([[[300, 231, 97], [45, 23, 18]]])
    img2.apply_threshold()
    test_output2 = img2.image
    assert np.isclose(test_output2[0, 0, 0], 255), 'bad thresholding'
    assert np.isclose(test_output2[0, 1, 2], 0)
    return
