# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from __future__ import absolute_import, division

from os.path import isfile

from six import string_types
from contextlib import contextmanager
import numpy as np
from scipy._lib._util import check_random_state


@contextmanager
def maybe_open(file_handle, mode="r"):
    """Context manager to open if passed a filename and seek/set on fileobject

    Parameter
    --------
    file_handle: str or IO
        open as IO object if str type

    """
    is_file = False
    if isinstance(file_handle, string_types):
        file_handle = open(file_handle, mode)
        is_file = True
    else:
        curpos = file_handle.tell()
        file_handle.seek(0)

    try:
        yield file_handle
    finally:
        if is_file:
            file_handle.close()
        else:
            file_handle.seek(curpos)


def check_file_exists(fname):
    """check if file exists or raise error"""
    if not isfile(fname):
        raise IOError("File does not exist: {}".format(fname))


class TowerSample(object):
    def __init__(self, distribution, start, stop, ngrid=1000, random_state=None):
        """Draw sample from an arbitrary distribution using tower sampling.

        Parameters
        ----------
        distribution : callable
            distribution to sample. It should accept only 1 parameter
        start, stop : float
            define a range over which the distribution should be sampled.
        ngrid : int, optional
            number of elements in the tower
        random_state : np.random.RandomState, optional
            random state to use. If ``None`` use np.random singleton

        """
        self._linspace = np.linspace(start, stop, ngrid)
        self._tower = distribution(self._linspace).cumsum()
        self._random_state = check_random_state(random_state)

    def _value(self):
        random = self._random_state.uniform(high=self._tower[-1])
        idx = np.where(random < self._tower)[0][0]
        return self._linspace[idx]

    def draw(self, size=1):
        """
        size : int, optional
            number of samples to draw
        """
        if size == 1:
            return self._value()
        else:
            return np.array([self._value() for _ in range(size)])


def closest_image(p1, p2, box):
    """find closest image in periodic box"""
    box = np.asarray(box)
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)

    box2 = .5 * box
    res = p2.copy()
    diff = p1 - p2
    for i in range(3):
        if diff[i] < -box2[i]:
            res[i] -= box[i]
        elif diff[i] > box2[i]:
            res[i] += box[i]
    return res


def norm(v):
    """faster norm calculation"""
    return np.sqrt(np.dot(v, v))


def triangle_inequality(x, y, z):
    """evaluate triangle inequality for 3 points"""
    return x < (y + z) and y < (z + x) and z < (y + x)


def spherical_to_cartesion(r, theta, phi):
    """following physics convention"""
    s_theta = np.sin(theta)
    return np.array(
        [r * s_theta * np.cos(phi), r * s_theta * np.sin(phi), r * np.cos(theta)]
    )


def law_cosines(a, b, c):
    """calculate an angle cosine gamma in a triangle based on law of cosines"""
    return (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
