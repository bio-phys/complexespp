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
import os
from glob import glob

folders = glob('Residue-*')


res = True

for f in folders:
    cpp = pd.read_csv(os.path.join(f, 'energy.stat')).energy
    ref = np.load(os.path.join(f, 'ref.npz'))['arr_0']
    try:
        np.testing.assert_almost_equal(ref, cpp, decimal=5)
    except:
        res = False

print("all ok? = {}".format(res))
sys.exit(not res)
