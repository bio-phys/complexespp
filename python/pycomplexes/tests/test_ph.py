# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
import MDAnalysis as mda
import numpy as np
import yaml
import pytest
import os
from os.path import dirname, join as pjoin
from numpy.testing import assert_equal
from numpy.testing import assert_array_equal

from collections import namedtuple
from pycomplexes import ph
from pycomplexes.testing import data, datafiles

top_without_copies = """
box: [300.0, 300.0, 300.0]
topology:
    b0-0:
        coordinate-file: {0}/chmp3_model.pdb
        connections: []
        domains:
            Chm:
                type : "rigid"
                selection: protein
    b0-1:
        coordinate-file: {0}/chmp3_model.pdb
        connections: []
        domains:
            Chm:
                type : "rigid"
                selection: protein
"""

top_with_copies = """
box: [300.0, 300.0, 300.0]
topology:
    b0:
        coordinate-file: {}/chmp3_model.pdb
        copies: 2
        domains:
            Chm:
                type : "rigid"
                selection: protein
        connections: []
"""


def test_degree_of_dissociation():
    phs = [0., 6., 12.]
    test_list = [ph.degree_of_dissociation(x, 6.) for x in phs]
    expected_list = [0.0, 0.5, 1.0]
    # We assert approximately. typical float comparison problem...
    # since the result of this function is rounded in net_charge anyways,
    # this approximate test should be very sufficient
    assert pytest.approx(test_list, abs=1e-5) == expected_list


input_tuple = (6., -1)
expected_list = [
    -0.0,
    -0.0,
    -0.0,
    -0.0,
    -0.01,
    -0.03,
    -0.09,
    -0.24,
    -0.5,
    -0.76,
    -0.91,
    -0.97,
    -0.99,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
    -1.0,
]
input_tuple2 = (7., 1)
expected_list2 = [
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    0.99,
    0.97,
    0.91,
    0.76,
    0.5,
    0.24,
    0.09,
    0.03,
    0.01,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]


@pytest.mark.parametrize(
    "test_input,expected",
    [(input_tuple, expected_list), (input_tuple2, expected_list2)],
)
def test_net_charge(test_input, expected):
    phs = np.arange(2, 12, .5)
    test_list = [ph.net_charge(x, *test_input) for x in phs]
    assert test_list == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (0., [-0.0, 1.0, 1.0, 1.0, -0.0, -0.0, 0]),
        (6.5, [-1.0, 1.0, 0.5, 1.0, -0.99, -0.0, 0]),
        (15., [-1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0]),
    ],
)
def test_change_charges_in_domain(test_input, expected):
    some_domain_dict = {
        "beads": ["ASP", "LYS", "HIS", "ARG", "GLU", "TYR", "ALA"],
        "charges": [-1, 1, 1, 1, -1, -1, 0],
        "name": "Spongebob",
    }
    ph.change_charges_in_domain(some_domain_dict, test_input, ph.CHARGED_SIDECHAINS)
    # NOTE: charges are changed in place
    assert some_domain_dict["charges"] == expected


def test_change_charges(data):
    expected_cplx = """
box: [100.0, 100.0, 100.0]
definitions:
  some:
    definitions
topologies:
- domains:
    0:
      beads: [ARG]
      charges: [0.0]
    1:
      beads: [LYS]
      charges: [0.0]
  name: rick
  ndomains: 2
- domains:
    2:
      beads: [GLU]
      charges: [-1.0]
    3:
      beads: [HIS]
      charges: [0.0]
  name: trick
  ndomains: 2
- domains:
    4:
      beads: [TYR]
      charges: [-1.0]
    5:
      beads: [ASP]
      charges: [-1.0]
  name: track
  ndomains: 2
- domains:
    6:
      beads: [ALA]
      charges: [0.0]
    7:
      beads: [ALA]
      charges: [0.0]
  name: donald
  ndomains: 2
   """
    expected_cplx = yaml.load(expected_cplx)
    with open(data["test_ph.cplx"]) as f:
        test_cplx = yaml.load(f)
    test_cplx = ph.change_charges(test_cplx, 15.)
    assert_equal(test_cplx, expected_cplx)
