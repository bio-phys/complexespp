#!/usr/bin/env python
# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------

import numpy as np
import pandas as pd
import sys
from glob import glob

all_ok=True
for i in range(0,3):
    org = pd.read_csv('test' + str(i) + '/27-remc.stat')
    other = pd.read_csv('test' + str(i) + '/27-seq.stat')
    all_ok &= np.allclose(org.energy, other.energy, rtol=0, atol=.05)

print("allclose = {}".format(all_ok))
sys.exit(not all_ok)
