import numpy as np

from scipy import stats
from scipy.special import erf
from scipy.spatial.distance import cdist

from pycomplexes.linker import gauss
from pycomplexes.util import norm

import pytest

SEED = 50884274


def test_gauss_pdf():
    b = 3.81
    x = np.array(
        [
            0.,
            0.55555556,
            1.11111111,
            1.66666667,
            2.22222222,
            2.77777778,
            3.33333333,
            3.88888889,
            4.44444444,
            5.,
        ]
    )
    expected = np.array(
        [
            0.,
            0.0224105,
            0.0814626,
            0.1562736,
            0.2222318,
            0.2605958,
            0.2642221,
            0.2375745,
            0.1923176,
            0.1415329,
        ]
    )
    np.testing.assert_almost_equal(gauss.gauss_pdf(x, b), expected)


def test_get_bead_coord():
    ps = np.array([0, 0, 0])
    pe = np.array([10, 0, 0])
    p3 = np.array([1, 0, 2])

    # add pi / 4 offset due do internal coordinate system offset used in conversion.
    point = gauss.get_bead_coord(ps, pe, norm(ps - p3), norm(pe - p3), np.pi / 4)
    np.testing.assert_almost_equal(point, p3)


@pytest.fixture
def gauss_linker():
    return gauss.GaussianChainDistribution(3.81, random_state=SEED, ngrid=1000)


def gauss_cdf(r, bond_length):
    rms = bond_length ** 2
    a = 3 / (2 * rms)
    fac = 4 * np.pi * (a / np.pi) ** 1.5
    return fac * (
        (np.sqrt(np.pi) * erf(np.sqrt(a) * r)) / (4 * np.power(a, 1.5))
        - (r * np.exp(-a * r * r)) / (2 * a)
    )


def test_gaussian_chain_sphere(gauss_linker):
    rvs = [gauss_linker.sphere() for _ in range(10000)]
    ks = stats.kstest(rvs, lambda x: gauss_cdf(x, gauss_linker._bond_length))
    # 20 % of the time we expect the difference between the CDF estimated from
    # rvs and the analytic cdf is as large as the observed one this time. This
    # mean the NULL-hypothesis is confirmed and we really sample a normal
    # distribution. The value here is small because we can't achieve a better
    # sampling for this function with the discretization.
    assert ks.pvalue > .2
    # For detail see
    # https://sites.google.com/a/ucsc.edu/krumholz/teaching-and-courses/ast119_w15/class-10


def cleaved_gauss_cdf(x, start, stop, bond_length):
    scale = 1 / (gauss_cdf(stop, bond_length) - gauss_cdf(start, bond_length))
    return scale * (gauss_cdf(x, bond_length) - gauss_cdf(start, bond_length))


@pytest.mark.parametrize("d12, d13", ((8.519, 5.672), (5.672, 8.519)))
def test_gaussian_chain_circle_on_sphere(gauss_linker, d12, d13):
    nEl = 5
    N = 10000

    gl = gauss_linker
    rvs = [gl.circle_on_sphere(d12, d13, nEl - 1) for _ in range(N)]

    start = d12 - d13 if d12 > d13 else d13 - d12
    stop = d12 + d13
    ks = stats.kstest(
        rvs,
        lambda x: cleaved_gauss_cdf(x, start, stop, np.sqrt(nEl - 1) * gl._bond_length),
    )
    assert ks.pvalue > .2


def test_grow_gauss_chain():
    start = [0, 0, 0]
    stop = [10, 10, 10]
    box = [20, 20, 20]
    N = 50
    bond_length = 5
    xyz = gauss.grow_gauss_chain(
        start, stop, box, N, random_state=SEED, bond_length=bond_length
    )

    np.testing.assert_equal(xyz[0], start)
    np.testing.assert_equal(xyz[-1], stop)
    np.testing.assert_equal(xyz.shape, (N, 3))
    assert not np.any(cdist(xyz, xyz) + np.eye(N) * 10 < bond_length)
