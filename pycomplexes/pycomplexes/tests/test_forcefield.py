# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from os.path import join as pjoin, dirname, basename
from glob import glob
import six
from numpy import array

import pytest
from numpy.testing import assert_array_equal

from pycomplexes import forcefield
from pycomplexes.testing import data, datafiles


def test_scale_interaction():
    # yapf: disable
    ff = {'energies':
           {'ALA':
             {'ALA': -2}
           },
          'other': 'test'}
    # yapf: enable
    scale = 2
    offset = 2
    new_ff = forcefield.scale_interaction(ff, scale, offset)
    assert new_ff["other"] == "test"
    assert new_ff["energies"]["ALA"]["ALA"] == -8


def test_read_forcefield(datafiles):
    from_fname = forcefield.read_forcefield(datafiles["testff"])
    from_name = forcefield.read_forcefield(forcefield.get_forcefield_file("KimHummer"))
    assert from_name.keys() == from_fname.keys()


def test_write_forcefield(tmpdir):
    # yapf: disable
    ff = {'energies':
           {'ALA':
             {'ALA': -2}
           },
          'other': 'test'}
    # yapf: enable

    # see if I can write if file doesn't exists
    forcefield.write_forcefield(ff, pjoin(tmpdir.dirname, "the-answer-is-42"))

    # raise error if file exists
    with pytest.raises(IOError):
        forcefield.write_forcefield(ff, pjoin(tmpdir.dirname, "the-answer-is-42"))


def test_known_forcefields():
    ff = forcefield.known_forcefields()
    assert isinstance(ff, list)


def test_get_forcefield_file(datafiles):
    from_fname = forcefield.get_forcefield_file(datafiles["testff"])
    assert from_fname == datafiles["testff"]

    from_name = forcefield.get_forcefield_file("KimHummer")
    assert basename(from_name) == "KimHummer"
