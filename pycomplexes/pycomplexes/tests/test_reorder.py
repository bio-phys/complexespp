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
from os.path import join as pjoin, dirname
import MDAnalysis as mda
from numpy.testing import assert_array_equal
from pycomplexes import reorder
from collections import namedtuple
from pycomplexes.testing import data


def test_create_selection_strings(data):
    """Testing the create_selection_string function with a simple
    example-pdb-file.
    """
    pdbfile = "reorder_test.pdb"
    fname = data[pdbfile]
    u = mda.Universe(fname)
    test_string_array = reorder.create_selection_strings(u)
    expected_string_array = [
        "segid A and resid 1",
        "segid A and resid 2",
        "segid A and resid 3",
        "segid A and resid 4",
        "segid A and resid 5",
        "segid A and resid 6",
        "segid A and resid 7",
        "segid A and resid 10",
        "segid B and resid 8",
        "segid B and resid 9",
        "segid C and resid 11",
        "segid C and resid 12",
    ]
    assert_array_equal(test_string_array, expected_string_array)


def test_main(tmpdir, data):
    """Testing the main function of reorder.py with a simple
    example-pdb-file. Especially the successful writing of
    a pdb is tested (in tmpdir).
    """
    pdbfile = "reorder_test.pdb"
    fname = data[pdbfile]
    Args = namedtuple("Args", "pdbfile, outputfile")
    output_filename = pjoin(tmpdir.dirname, "test.pdb")
    args = Args(fname, output_filename)
    reorder.Reorder.main(args)  # Testint the class method
    with open(pjoin(tmpdir.dirname, "test.pdb"), "r") as outputfile:
        test_lines = [l for l in outputfile.readlines() if "ATOM" in l]
    expected_lines = (
        "ATOM      1  CA  LYS A   1      -7.998   1.609  -1.490  1.00  0.00      A    C\n",
        "ATOM      2  CA  THR A   2      -4.754   3.630  -1.008  1.00  0.00      A    C\n",
        "ATOM      3  CA  TRP A   3      -1.368   3.031  -2.738  1.00  0.00      A    C\n",
        "ATOM      4  CA  ASN A   4       1.663   2.174  -0.525  1.00  0.00      A    C\n",
        "ATOM      5  CA  PRO A   5       4.950   3.274  -2.280  1.00  0.00      A    C\n",
        "ATOM      6  CA  ALA A   6       6.957   1.731   0.626  1.00  0.00      A    C\n",
        "ATOM      7  CA  THR A   7       5.931  -1.892  -0.279  1.00  0.00      A    C\n",
        "ATOM      8  CA  TRP A  10      -2.632  -1.198  -2.639  1.00  0.00      A    C\n",
        "ATOM      9  CA  GLY B   8       4.569  -1.407  -3.861  1.00  0.00      B    C\n",
        "ATOM     10  CA  LYS B   9       0.918  -2.504  -3.180  1.00  0.00      B    C\n",
        "ATOM     11  CA  THR C  11      -4.250  -1.170   0.860  1.00  0.00      C    C\n",
        "ATOM     12  CA  GLU C  12      -7.968  -1.794   1.715  1.00  0.00      C    C\n",
    )
    assert_array_equal(test_lines, expected_lines)
