import numpy as np

from scipy.spatial.distance import cdist

from .. import util
from .collision_detector import CollisionDetector


def gauss_pdf(r, b):
    """calculate PDF to find a bead at distance r in  gaussian chain
    Parameters
    ----------
    r : float
        distance
    b : float
        bond length

    Returns
    -------
    float : probability
    """
    b2 = b * b
    return (
        4
        * np.pi
        * r ** 2
        * (3 / (2 * np.pi * b2)) ** (3 / 2)
        * np.exp(-3 * r ** 2 / 2 / b2)
    )


class GaussianChainDistribution(object):
    """
    Helper class to draw correct distributed random numbers to build gaussian chain
    """

    def __init__(self, bond_length, start=0, ngrid=400, random_state=None):
        self._bond_length = bond_length
        self._ngrid = ngrid

        self._random_state = util.check_random_state(random_state)
        # stop being 3 times the bond_length is a rule of thumb for acceptable results.
        self._tower = util.TowerSample(
            lambda r: gauss_pdf(r, bond_length),
            start,
            3 * bond_length,
            ngrid=ngrid,
            random_state=self._random_state,
        )

    def sphere(self):
        """draw distance to next sphere"""
        return self._tower.draw()

    def circle_on_sphere(self, sphere_radius, distance_to_sphere, nbeads):
        """Draw distance distance from bead 0 to bead N-1 so that it lies on sphere around bead N
        """
        start = (
            distance_to_sphere - sphere_radius
            if distance_to_sphere > sphere_radius
            else sphere_radius - distance_to_sphere
        )
        stop = distance_to_sphere + sphere_radius
        chain_length = np.sqrt(nbeads) * self._bond_length

        return util.TowerSample(
            lambda r: gauss_pdf(r, chain_length),
            start=start,
            stop=stop,
            ngrid=self._ngrid,
            random_state=self._random_state,
        ).draw()

    def draw(self, start, stop, nbeads, itmax=100):
        """
        select new bead coordinaets based on a start and end point and chain length
        """
        distance = util.norm(start - stop)
        d_to_start = self.sphere()
        d_to_end = self.circle_on_sphere(d_to_start, distance, nbeads)

        i = 0
        while (
            not (util.triangle_inequality(distance, d_to_start, d_to_end)) and i < itmax
        ):
            d_to_start = self.sphere()
            d_to_end = self.circle_on_sphere(d_to_start, distance, nbeads)
            i += 1

        if i == itmax:
            raise RuntimeError("couldn't find points to fulfill triangle inequality")

        phi = self._random_state.uniform(0, 2 * np.pi)
        return get_bead_coord(start, stop, d_to_start, d_to_end, phi)


def get_bead_coord(p_start, p_end, dist_to_start, dist_to_end, phi):
    # build a coordinate system with z-axis along start-end vector
    az = (p_end - p_start) / util.norm(p_end - p_start)
    ax = np.ones(3)
    if az[0] != 0:
        ax[0] = -(az[1] + az[2]) / az[0]
    elif az[1] != 0:
        ax[1] = -(az[0] + az[2]) / az[1]
    elif az[2] != 0:
        ax[2] = -(az[0] + az[1]) / az[2]
    ax /= util.norm(ax)
    ay = np.cross(az, ax)
    a = np.array([ax, ay, az]).T

    # calculate angle between z-axis and vector using law of cosines
    d = util.norm(p_end - p_start)
    cos_theta = util.law_cosines(dist_to_start, d, dist_to_end)
    # catch edge cases for numerical stability
    if cos_theta > 1:
        theta = 0
    elif cos_theta < -1:
        theta = np.pi
    else:
        theta = np.arccos(cos_theta)

    # convert point into lab coordinate system
    r = util.spherical_to_cartesion(dist_to_start, theta, phi)
    return np.dot(a, r) + p_start


def grow_gauss_chain(
    start,
    stop,
    box,
    nbeads,
    bond_length=3.81,
    check_collision=True,
    other_coords=None,
    random_state=None,
    ntrials=100,
):
    """grow a linker of length N between the point start and stop. Here start and stop are being fixed.

    Parameters
    ----------
    start, stop : array_like (3, )
         start and end point
    box : array_like (3, )
         periodic box information
    nbeads : int
         chain length including start and stop points
    bond_length : float, optional
         length of a single bond
    check_collision : bool, optional
         check if chain grows into itself.
    other_coords : array_like (3, ), optional
         other coordinates to check against clashes Only works if
         check_collision is ``True``.
    random_state : int, optional
         int or numpy random state.
    ntrials : int, optional
         number of retries for growing a single bead

    Returns
    -------
    xyz : array_like (N, 3)
        bead positions including start and stop

    """
    random_state = util.check_random_state(random_state)
    # random growing direction
    grow_from_stop = random_state.uniform(0, 1) < .5
    if grow_from_stop:
        start, stop = stop, start
    # account for periodic boundary conditions
    stop = util.closest_image(start, stop, box)

    # setup detector if required
    col_detector = None
    if other_coords is not None:
        col_detector = CollisionDetector(other_coords, bond_length, box)

    chain_pdf = GaussianChainDistribution(
        bond_length,
        start=bond_length if check_collision else 0,
        random_state=random_state,
    )
    xyz = np.ones((nbeads, 3)) * -1000
    xyz[0] = start
    xyz[-1] = stop

    bl = bond_length
    for i in range(nbeads - 2):
        nbeads_left = nbeads - i - 2
        # greatly increases changes of successful growth if we retry failed attempts
        new_pos = chain_pdf.draw(start, stop, nbeads_left)
        if check_collision:
            for _ in range(ntrials):
                try:
                    new_pos = chain_pdf.draw(start, stop, nbeads_left)
                except RuntimeError:
                    continue

                self_collision = np.any(cdist([new_pos], xyz) < bl)
                other_collision = False
                if col_detector is not None:
                    other_collision = col_detector.collision(new_pos)
                if not self_collision and not other_collision:
                    break
        start = new_pos
        xyz[i + 1] = new_pos

    if grow_from_stop:
        xyz = xyz[::-1]

    if np.any(xyz == -1000):
        raise RuntimeError("Grown into something")

    return xyz
