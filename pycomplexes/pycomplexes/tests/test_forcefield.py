# Copyright (c) 2018 the pycomplexes development team and contributors
# (see the file AUTHORS for the full list of names)
#
# This file is part of pycomplexes.
#
# pycomplexes is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pycomplexes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pycomplexes.  If not, see <https://www.gnu.org/licenses/>
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
