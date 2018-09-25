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
