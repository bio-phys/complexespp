from MDAnalysis.lib.pkdtree import PeriodicKDTree as PKDTree

import numpy as np


class CollisionDetector(object):
    """Given a set of N particles this class allows to query if any new particle
    will overlap. If particles overlap this is called a collision

    Example
    -------
    >>> cd = CollisionDetector(coords, 3.81, [100, 100, 100])
    >>> cd.collision([50, 50, 50])

    """

    def __init__(self, coords, diameter, box):
        """
        Parameters
        ----------
        coords : array_like (N, 3)
            coordinates of N particles
        diameter : float
            diameter of particles
        box : array_like (3)
            box dimensions
        """
        box = np.asarray(box, dtype=np.float32)
        coords = np.asarray(coords, dtype=np.float32)
        self._kdt = PKDTree(box)
        self._kdt.set_coords(coords)
        self._diameter = diameter

    def collision(self, pos):
        """check is there is a collision

        Parameters
        ----------
        pos : array_like (3)
             coordinates of new particles

        Returns
        -------
        collision : bool
             ``True`` if there is a collision,
        """
        pos = np.asarray(pos, dtype=np.float32)
        self._kdt.search(pos, self._diameter)
        indices = self._kdt.get_indices()
        return len(indices) != 0
