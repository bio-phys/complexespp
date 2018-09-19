import numpy as np

from scipy.spatial.distance import cdist

from pycomplexes.linker.collision_detector import CollisionDetector

import pytest

SEED = 508842746
BOX = np.asarray([100, 100, 100], dtype=np.float32)
cell_size = 7


def test_collision():
    rng = np.random.RandomState(SEED)
    coords = rng.uniform(0, 100, size=(1000, 3))
    detector = CollisionDetector(coords, cell_size, BOX)
    assert detector.collision(coords[0])


def test_no_collision():
    rng = np.random.RandomState(SEED)
    coords = rng.uniform(0, 10, size=(1000, 3))
    detector = CollisionDetector(coords, cell_size, BOX)
    assert not detector.collision([50, 50, 50])


def test_pbc_collision():
    rng = np.random.RandomState(SEED)
    coords = rng.uniform(0, 10, size=(1000, 3))
    detector = CollisionDetector(coords, cell_size, BOX)
    assert detector.collision([99, 99, 99])
