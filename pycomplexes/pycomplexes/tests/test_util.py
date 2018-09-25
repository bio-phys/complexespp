import numpy as np

from scipy import stats

from pycomplexes import util

SEED = 508842746


class TestTowerSample(object):
    def test_size_1(self):
        sample = util.TowerSample(stats.norm.pdf, -5, 5, random_state=SEED).draw()
        assert isinstance(sample, float)

    def test_size_n(self):
        size = 10
        sample = util.TowerSample(stats.norm.pdf, -5, 5, random_state=SEED).draw(
            size=size
        )
        assert sample.size == size

    def test_correct_sampling(self):
        sample = util.TowerSample(stats.norm.pdf, -5, 5, random_state=SEED).draw(
            size=10000
        )
        assert np.all(sample < 5)
        assert np.all(sample > -5)
        np.testing.assert_almost_equal(sample.mean(), 0.0, decimal=1)
        np.testing.assert_almost_equal(sample.var(), 1.0, decimal=1)

    def test_kstest(self):
        for _ in range(10):
            rvs = util.TowerSample(
                stats.norm.pdf, -5, 5, ngrid=5000, random_state=SEED
            ).draw(size=10000)
            ks = stats.kstest(rvs, stats.norm.cdf)
            assert ks.pvalue > .5


def test_norm():
    for _ in range(100):
        x = np.random.normal(size=3)
        np.testing.assert_almost_equal(np.linalg.norm(x), util.norm(x))


def test_closest_image():
    box = [10, 10, 10]
    p1 = [1, 1, 1]
    p2 = [9, 1, 1]
    ci = util.closest_image(p1, p2, box)
    np.testing.assert_equal(ci, [-1, 1, 1])


def test_triangle_inequality():
    assert util.triangle_inequality(2, 3, 3)
    assert not util.triangle_inequality(2, 3, 5)


def test_spherical_to_cartesion():
    np.testing.assert_almost_equal(
        util.spherical_to_cartesion(1, np.pi / 2, 0), [1, 0, 0]
    )


def test_law_cosines():
    a, b, c = 2, 3, 3
    np.testing.assert_almost_equal(
        np.arccos(util.law_cosines(a, b, c)), 1.2309594173407747
    )
    np.testing.assert_almost_equal(
        np.arccos(util.law_cosines(b, c, a)), 0.6796738189082439
    )
    np.testing.assert_almost_equal(
        np.arccos(util.law_cosines(c, a, b)), 1.2309594173407747
    )
