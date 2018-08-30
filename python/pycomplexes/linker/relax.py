import numba

import numpy as np
from scipy import constants
from MDAnalysis.lib import transformations

from .collision_detector import CollisionDetector

UNITS_ENERGY = constants.k * 300


@numba.jit(nopython=True)
def norm2(x):
    return np.sum(x * x)


@numba.jit(nopython=True)
def norm(x):
    return np.sqrt(norm2(x))


@numba.jit(nopython=True)
def pseudo_angle(theta):
    gamma = 0.1  # mol  / kcal
    ka = 106.4  # kcal / mol / rad^2
    kb = 26.3  # kcal / mol / rad^2
    ea = 4.3  # kcal / mol
    tha = 1.60  # rad
    thb = 2.27  # rad
    # conversion factor goes from kcal/mol to joule and then to kT
    unit_conversion = 6.9477e-21 / UNITS_ENERGY
    return (
        -np.log(
            np.exp(-gamma * (ka * (theta - tha) * (theta - tha) + ea))
            + np.exp(-gamma * kb * (theta - thb) * (theta - thb))
        )
        / gamma
        * unit_conversion
    )


@numba.jit(nopython=True, fastmath=True)
def angle_potential(xyz):
    u = 0
    for i in range(len(xyz) - 2):
        d1 = xyz[i] - xyz[i + 1]
        d2 = xyz[i + 1] - xyz[i + 2]
        angle = np.arccos(np.dot(d1, d2) / np.sqrt(norm2(d2) * norm2(d1)))
        u += pseudo_angle(angle)
    return u


@numba.jit(nopython=True)
def harmonic(x, x0, k):
    return .5 * k * (x - x0) ** 2


@numba.jit(nopython=True)
def bond_potential(xyz):
    bond = 3.81
    k = 634.0572472420426
    u = 0
    for i in range(len(xyz) - 1):
        u += harmonic(norm(xyz[i] - xyz[i + 1]), bond, k)
    return u


@numba.jit(nopython=True)
def pseudo_torsion(phi):
    # Here we use values for ALA-ALA
    delta = [287.354830, 271.691192, 180.488748, 108.041256]
    V = [0.936472, 2.307767, 0.131743, 0.613133]
    s = 0
    for i in range(4):
        s += V[i] * (1 + np.cos((i + 1) * phi - delta[i]))
    return s


@numba.jit(nopython=True)
def cross(v1, v2):
    out = np.empty(3)
    out[0] = v1[1] * v2[2] - v1[2] * v2[1]
    out[1] = v1[2] * v2[0] - v1[0] * v2[2]
    out[2] = v1[0] * v2[1] - v1[1] * v2[0]
    return out


@numba.jit(nopython=True)
def torsion_potential(xyz):
    u = 0
    for i in range(2, len(xyz) - 1):
        d1 = xyz[i - 2] - xyz[i - 1]
        d2 = xyz[i - 1] - xyz[i]
        d3 = xyz[i] - xyz[i + 1]

        c12 = cross(d1, d2)
        c23 = cross(d2, d3)
        c123 = cross(c12, c23)

        phi = -np.arctan2(np.dot(c123, d2) / norm(d2), np.dot(c12, c23))
        u += pseudo_torsion(phi)
    return u


@numba.jit(nopython=True)
def potential(xyz):
    return torsion_potential(xyz) + bond_potential(xyz) + angle_potential(xyz)


@numba.jit(nopython=True)
def translation(xyz, idx, width):
    xyz[idx] = xyz[idx] + (np.random.random(3) - .5) * width
    return xyz


@numba.jit
def rotation(xyz, idx, max_angle):
    angle = np.random.uniform(0, max_angle)
    vec = xyz[idx + 1] - xyz[idx - 1]
    vec /= norm(vec)
    x = xyz[idx] - xyz[idx-1]
    mat = transformations.rotation_matrix(angle, vec)[:3, :3]
    xyz[idx] = np.dot(mat, x) + xyz[idx-1]
    return xyz


def relax(xyz, temp=300, nsweep=100, target=.3, other_coords=None, box=None):
    beta = 300 / temp
    tmp = xyz.copy()
    E = potential(xyz)
    N = len(xyz)
    energy = np.empty(nsweep)
    acc = 0
    nstep = N - 2
    # use 3.81 the standard bond width
    width = 3.81 * .2

    # setup detector if required
    col_detector = None
    if other_coords is not None:
        col_detector = CollisionDetector(other_coords, 3.81, box)

    for i in range(nsweep):
        # do a sweep
        for j in range(nstep):
            idx = np.random.randint(N-2) + 1
            if np.random.uniform(1) < 0.5:
                tmp = translation(tmp, idx, width)
            else:
                tmp = rotation(tmp, idx, 2 * np.pi)

            # if a collision occurs this means the trials is rejected
            if col_detector is not None and col_detector.collision(tmp[idx]):
                tmp = xyz.copy()
                continue

            Enew = potential(tmp)
            dE = Enew - E
            if dE < 0 or np.random.random() < np.exp(-beta * dE):
                xyz = tmp.copy()
                E = Enew
                acc += 1
            else:
                tmp = xyz.copy()

        energy[i] = E

        # update step width
        rate = acc / ((i + 1) * nstep)
        if rate > target:
            width *= 1.1
        else:
            width *= .9

    return energy, xyz
