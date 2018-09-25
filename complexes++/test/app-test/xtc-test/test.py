#!/usr/bin/env python
# Copyright (c) 2018 the complexes++ development team and contributors
# (see the file AUTHORS for the full list of names)
#
# This file is part of complexes++.
#
# complexes++ is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# complexes++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with complexes++.  If not, see <https://www.gnu.org/licenses/>

import sys

import MDAnalysis as mda
import numpy as np
from six.moves import range

pdb = mda.Universe('traj.pdb')
xtc = mda.Universe('traj_reference.pdb', 'traj.xtc')
trr = mda.Universe('traj_reference.pdb', 'traj.trr')

p = pdb.select_atoms('all')
x = xtc.select_atoms('all')
t = trr.select_atoms('all')

res = np.bool(True)
for i in range(len(xtc.trajectory)):
    xtc.trajectory[i]
    pdb.trajectory[i]
    trr.trajectory[i]
    res = np.logical_and(res, np.allclose(x.positions, p.positions,
                                          rtol=0, atol=1e-3))
    res = np.logical_and(res, np.allclose(t.positions, p.positions,
                                          rtol=0, atol=1e-3))

print("all ok? = {}".format(res))
sys.exit(not res)
