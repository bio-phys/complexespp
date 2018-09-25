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
for i in range(3):
    remc = pd.read_csv('test{}/27-remc.stat'.format(i))
    hrex = pd.read_csv('test{}/27-hrex.stat'.format(i))
    all_ok &= np.allclose(remc.energy, hrex.energy, rtol=0, atol=.05)

with open('test1/complexes-remc.log') as fh:
    remc_nlines = len(fh.readlines())
with open('test1/complexes-hrex.log') as fh:
    hrex_nlines = len(fh.readlines())

# used different verbosity settings
all_ok &= remc_nlines != hrex_nlines


print("allclose = {}".format(all_ok))
sys.exit(not all_ok)
