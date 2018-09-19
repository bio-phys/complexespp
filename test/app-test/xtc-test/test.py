#!/usr/bin/env python
# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------

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
