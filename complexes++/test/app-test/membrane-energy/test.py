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
