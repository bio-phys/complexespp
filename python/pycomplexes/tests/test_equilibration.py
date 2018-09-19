# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
"""
module to test the equilibration module...

"""
from __future__ import absolute_import, print_function
from six.moves import cStringIO as StringIO

import yaml
import MDAnalysis as mda
import collections

from os.path import join as pjoin, dirname
from numpy.testing import assert_equal
import pytest

import pycomplexes.equilibration as equ
from pycomplexes.testing import datafiles

CUR_DIR = dirname(__file__)
DATA_DIR = pjoin(CUR_DIR, "data")


def test_update_coordinates(datafiles):
    """
    testing function update_coordinates

    """
    start_cplx = datafiles["lj-sys.cplx"]
    with open(start_cplx) as f:
        system = yaml.load(f)
    test_reference = datafiles["test_traj_reference.pdb"]
    test_traj = datafiles["test_traj.xtc"]
    eq_u = mda.Universe(test_reference, test_traj)
    eq_u.trajectory[-1]
    system = equ.update_coordinates(system, eq_u)
    updated_cplx = datafiles["lj-sys_updated.cplx"]
    with open(updated_cplx, "r") as f:
        updated_reference = yaml.load(f)
    assert_equal(system, updated_reference)


# Following the preparation of the different possible cases for argument_parsing:
Arguments = collections.namedtuple(
    "Arguments", "equilibrate_cplx trajectory reference_pdb cplx complexes_config frame"
)
args_console = Arguments(
    equilibrate_cplx="equilibrate_cplx",
    trajectory="trajectory",
    reference_pdb="reference_pdb",
    cplx="cplx",
    complexes_config=None,
    frame=None,
)
expected_console = "equilibrate_cplx", "trajectory", "reference_pdb", "cplx", None
# In the following reference_pdb is given as a test. Should not be processed...
# Using the same test config, as for testing visualize
args_com_config = Arguments(
    equilibrate_cplx=None,
    trajectory=None,
    reference_pdb="any",
    cplx="cplx",
    complexes_config=pjoin(DATA_DIR, "test_vis_com_config.conf"),
    frame=None,
)
expected_com_config = "cplx.cplx", "xtc.xtc", "xtc_reference.pdb", "cplx", None


@pytest.mark.parametrize(
    "test_input, expected",
    [(args_console, expected_console), (args_com_config, expected_com_config)],
)
def test_argument_processing(test_input, expected):
    assert equ.argument_processing(test_input) == expected
