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
import MDAnalysis as mda


def main():
    data = pd.read_csv('run-energy.stat')

    # check we write every frame
    if np.alltrue(data['# frame'] != np.arange(len(data))):
        print("number of rows in run-enery.stat incorrect")
        sys.exit(1)

    data_restart = pd.read_csv('run-energy-restart.stat')
    if data_restart['# frame'].iloc[0] != 49:
        print("wrong initial sweep number of restart file")
        sys.exit(1)

    if data_restart['# frame'].iloc[-1] != 99:
        print("simulated wrong number of sweeps in restart")
        sys.exit(1)

    # check if atom positions in restart overlap at least one frame
    top = 'run-test_reference.pdb'
    u = mda.Universe(top, 'run-test.trr')
    u.trajectory[-1]
    u_restart = mda.Universe(top, 'run-test-restart.trr')
    np.testing.assert_almost_equal(u.atoms.positions,
                                   u_restart.atoms.positions)


if __name__ == '__main__':
    main()
    print("all OK")
