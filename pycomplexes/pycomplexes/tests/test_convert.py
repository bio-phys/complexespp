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
import MDAnalysis as mda
import numpy as np
import yaml
import pytest
import os
from os.path import dirname, join as pjoin
from numpy.testing import assert_equal
from numpy.testing import assert_array_equal

from collections import namedtuple
from pycomplexes import convert
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


@pytest.mark.parametrize("top", [top_with_copies, top_without_copies])
def test_parse_top_file(top, data):
    top = yaml.load(top.format(data.folder))

    cplx = convert.convert(top, random=False)

    with open(data["test_convert.cplx"], "r") as f:
        expected_cplx = yaml.load(f)

    assert_equal(cplx, expected_cplx)


@pytest.fixture
def topology():
    top_yaml = """
box: [100, 100, 100]
topology:
    Binding:
        coordinate-file: foobar.pdb
        domains:
            foo:
                type: rigid
                selection: 'all'
            bar:
                type: rigid
                selection: 'all'
            foobar:
                type: gaussian
                selection: 'all'
    """
    return yaml.load(top_yaml)


def test_chain_ids(datafiles):
    # - Objective: Get a list from segid and resid from an atom group
    #
    # - Expected Result: A list of segid and resid formated as follow ["Segid"
    # - "Resid", ...]
    #
    # - Implementation: Call the function with an atom group formed by 3
    # - residus

    # Creation of an atom group with 3 residues (note residu 1 to 8 are not
    # described on the model of the chmp3
    u = mda.Universe(datafiles["chmp3_model.pdb"])
    ag = u.select_atoms("resid 9-11")
    IdsInput = namedtuple("IdsInput", ["segids", "resids"])
    # definition of an input list ['segid resid'] for the chmp3 and residues
    # 9-132, residues 1 to 8 are not described on the model
    test_input_list = np.arange(9, 12, 1)
    segment = "A"
    ids_list = [IdsInput(segment, n) for n in test_input_list]
    # the output of chain_ids is formated in a list with the format ['"segid""
    # ""resid"']
    ids_list_expected = ["{} {}".format(segment, elem.resids) for elem in ids_list]
    ids = convert.chain_ids(ag)
    assert_array_equal(ids, ids_list_expected)


def test_normalize_resnames():
    test_input = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "CYM",
        "HIP",
        "HID",
        "HIE",
        "GLH",
        "ASH",
        "LYN",
        "LJ",
    ]
    expected_resnames = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "CYS",
        "HIS",
        "HIS",
        "HIS",
        "GLU",
        "ASP",
        "LYS",
        "LJ",
    ]
    assert_equal(convert.normalize_resnames(test_input), expected_resnames)


def test_charges():
    # - Objective: Get the charge for the 20 common amino acids + amber + dummy (residues)
    #
    # - Expected Result: All amino acids should have the charge 0 except for
    #                    the amino acids 'HIS', 'ASP', 'GLU', 'LYS', 'ARG' that
    #                    should have respectively the charges 0.5, -1, -1, 1
    #                    and 1
    #
    # - Implementation: Call the function with the 20 common residues + amber + dummy
    test_input = np.array(
        [
            ["ALA", 0],
            ["ARG", 1],
            ["ASN", 0],
            ["ASP", -1],
            ["CYS", 0],
            ["GLN", 0],
            ["GLU", -1],
            ["GLY", 0],
            ["HIS", .5],
            ["ILE", 0],
            ["LEU", 0],
            ["LYS", 1],
            ["MET", 0],
            ["PHE", 0],
            ["PRO", 0],
            ["SER", 0],
            ["THR", 0],
            ["TRP", 0],
            ["TYR", 0],
            ["VAL", 0],
            ["CYM", -1],
            ["HIP", 1],
            ["HID", 0],
            ["HIE", 0],
            ["GLH", 0],
            ["ASH", 0],
            ["LYN", 0],
            ["LJ", 0],
        ],
        dtype=object,
    )  # included a dummy here, just for clarity.
    # testing the neutral non canonic names is in a way redundant
    MDAResidueMock = namedtuple("MDAResidueMock", "resname")
    charges_out = convert.charges([MDAResidueMock(n) for n in test_input[:, 0]])
    assert_array_equal(charges_out, test_input[:, 1])


def test_ProteinDomain(datafiles):
    u = mda.Universe(datafiles["chmp3_model.pdb"])
    node = {"selection": "all", "type": "rigid", "mc-moves": {"foo": "bar"}}
    dom = convert.ProteinDomain(node, "r", u, 0)
    assert_equal(dom.moves, {"foo": "bar"})


def test_ProteinTopology(datafiles):
    pdb_file = datafiles["chmp3_model.pdb"]
    topology = yaml.load(
        """
    coordinate-file: '{}'
    domains:
        rigid1:
            type: rigid
            selection: 'bynum 1-132'
        link1:
            type: gaussian
            selection: 'bynum 133-148'
            start_connection: [rigid1, 'bynum 132']
            end_connection: [rigid2, 'bynum 149']
        rigid2:
            type: rigid
            selection: 'bynum 149-164'
        link2:
            type: gaussian
            selection: 'bynum 165-198'
            start_connection: [rigid2, 'bynum 164']
            end_connection: [rigid3, 'bynum 199']
        rigid3:
            type: rigid
            selection: 'bynum 199-214'
    """.format(
            pdb_file
        )
    )
    ptop = convert.ProteinTopology(topology, "test", False, [100, 100, 100])
    assert len(ptop.connections) == 6


def test_1000_copies(datafiles):
    top_dict = {
        "box": [10, 10, 10],
        "topology": {
            "LJ": {
                "domains": {"Lys": {"selection": "bynum 1", "type": "rigid"}},
                "connections": [],
                "coordinate-file": datafiles["chmp3_model.pdb"],
                "copies": 1000,
            }
        },
    }

    cplx = convert.convert(top_dict, random=True)
    assert len(cplx["topologies"]) == 1000
