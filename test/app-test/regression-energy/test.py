#!/usr/bin/env python
# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd

energy = pd.read_csv('energy-rerun.stat').energy
ref = np.loadtxt('reference-energy.dat')

res = True
try:
    np.testing.assert_almost_equal(ref, energy, decimal=4)
except:
    res = False

print("all ok? = {}".format(res))
sys.exit(not res)
