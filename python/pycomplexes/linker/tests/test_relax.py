from pycomplexes.linker import relax

import numpy as np


def test_norm():
    for _ in range(10):
        v = np.random.random(size=3)
        np.testing.assert_almost_equal(relax.norm(v), np.linalg.norm(v))


def test_norm2():
    for _ in range(10):
        v = np.random.random(size=3)
        np.testing.assert_almost_equal(relax.norm2(v), np.linalg.norm(v) ** 2)


def test_harmonic():
    np.testing.assert_almost_equal(relax.harmonic(0, 0, 1), 0)


def test_cross():
    for _ in range(10):
        a = np.random.random(size=3)
        b = np.random.random(size=3)
        np.testing.assert_almost_equal(relax.cross(a, b), np.cross(a, b))


def test_potential():
    xyz = np.random.random(size=(10, 3))
    np.testing.assert_almost_equal(
        relax.potential(xyz),
        relax.torsion_potential(xyz)
        + relax.bond_potential(xyz)
        + relax.angle_potential(xyz),
    )
