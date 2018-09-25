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
