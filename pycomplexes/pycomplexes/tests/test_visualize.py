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
"""
module to test the visualize module...

"""
from __future__ import absolute_import, print_function
from six.moves import cStringIO as StringIO

import yaml
import six
import collections

from os.path import join as pjoin, dirname
from numpy.testing import assert_array_equal
import pytest

import pycomplexes.visualize as vis
from pycomplexes import convert


DATA_DIR = pjoin(dirname(__file__), "data")


def test_vmd_sel_chainid():
    """
    tests the getchainid function
    """
    test_chainids_array = ["A 1", "B 2", "L 12"]
    test_output_string = vis.vmd_sel_chainid(test_chainids_array[1])
    expected_output_string = "chain {} and resid {}".format("B", 2)
    assert expected_output_string == test_output_string


def test_vmd_selection_string():
    """
    tests the vmd_selection_string function
    """
    test_topology = {
        "topologies": {0: {"domains": {0: {"chain-ids": ["A 1", "B 2", "L 12"]}}}}
    }
    chainids = test_topology["topologies"][0]["domains"][0]["chain-ids"]
    test_output_string = vis.vmd_selection_string(chainids)
    expected_output_string = (
        "chain {} and resid {}".format("A", 1)
        + " or "
        + "chain {} and resid {}".format("B", 2)
        + " or "
        + "chain {} and resid {}".format("L", 12)
    )
    assert expected_output_string == test_output_string


def test_vmd_selection_string_trace():
    """
    tests the vmd_selection_string function
    """
    test_topology = {
        "topologies": {
            0: {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "gaussian"}
                }
            }
        }
    }
    chainids = test_topology["topologies"][0]["domains"][0]["chain-ids"]
    test_output_string = vis.vmd_selection_string_trace(chainids)
    expected_output_string = (
        "chain {} and resid {}".format("A", 0)
        + " or "
        + "chain {} and resid {}".format("A", 1)
        + " or "
        + "chain {} and resid {}".format("B", 2)
        + " or "
        + "chain {} and resid {}".format("L", 12)
        + " or "
        + "chain {} and resid {}".format("L", 13)
    )
    assert expected_output_string == test_output_string


def test_vmd_rep_gen():
    """
    testing the vmd_rep_gen function.
    """
    # creating test topology
    #############################
    test_topology = {
        "topologies": [
            {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "rigid"}
                }
            }
        ]
    }
    test_topology["topologies"][0]["domains"].update(
        {1: {"chain-ids": ["B 0", "B 1", "C 2", "M 12"], "type": "rigid"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {2: {"chain-ids": ["C 0", "C 1", "D 2", "N 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["ndomains"] = 3
    test_topology["topologies"][0]["full-move"] = False

    test_topology["topologies"][0]["domains"][0]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][1]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][2]["nbeads"] = 3
    test_topology["definitions"] = convert.DEFAULT_DEFINITIONS
    ##############################
    # assigning certain test-values
    coloring = "random"
    ###############################
    test_output_string = vis.vmd_rep_gen(test_topology, coloring)
    expected_output_string = (
        "mol addrep 0 \n"
        "mol modselect 0 0 chain A and resid 0 or chain A and resid 1 or chain B "
        "and resid 2 or chain L and resid 12 \nmol modcolor 0 0 ColorID 0 \n"
        "mol modstyle 0 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \n"
        "mol modselect 1 0 chain B and resid 0 or chain B and resid 1 or chain C and resid 2 "
        "or chain M and resid 12 \nmol modcolor 1 0 ColorID 1 \n"
        "mol modstyle 1 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \n"
        "mol modselect 2 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 "
        "or chain N and resid 12 or chain N and resid 13 \n"
        "mol modcolor 2 0 ColorID 2 \nmol modstyle 2 0 Trace 0.300000 12.000000 \n\n"
        "mol addrep 0 \n"
        "mol modselect 3 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 or "
        "chain N and resid 12 \nmol modcolor 3 0 ColorID 2 \n"
        "mol modstyle 3 0 VDW 1.000000 12.000000 \n\n"
        "#hiding automatically generated all-rep\n# delrep syntax:  #rep #mol \nmol delrep 4 0 "
    )
    assert expected_output_string == test_output_string


def test_domain_rep():
    """
    testing the domain_rep_function
    """
    # creating test topology
    #############################
    test_topology = {
        "topologies": [
            {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "rigid"}
                }
            }
        ]
    }
    test_topology["topologies"][0]["domains"].update(
        {1: {"chain-ids": ["B 0", "B 1", "C 2", "M 12"], "type": "rigid"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {2: {"chain-ids": ["C 0", "C 1", "D 2", "N 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["ndomains"] = 3
    test_topology["topologies"][0]["full-move"] = False

    test_topology["topologies"][0]["domains"][0]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][1]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][2]["nbeads"] = 3
    test_topology["definitions"] = convert.DEFAULT_DEFINITIONS
    ##############################
    # assigning certain test-values and creating test-string in necessary loop!
    coloring = "random"
    out = StringIO()
    bead_counter = 1
    rep_counter = 0
    for j, top in enumerate(test_topology["topologies"]):
        domains = top["domains"]
        for i, dom in six.iteritems(domains):
            chainids = dom["chain-ids"]
            rep_number = rep_counter
            first_bead = bead_counter
            last_bead = first_bead + dom["nbeads"]
            domain_number = i
            bead_counter = last_bead + 1
            dom_rep, rep_counter = vis.vmd_domain_rep(
                dom["type"], chainids, coloring, rep_counter, rep_number, i, dom
            )
            out.write(dom_rep)
    test_output_string = out.getvalue()
    out.close()
    expected_output_string = (
        "mol addrep 0 \n"
        "mol modselect 0 0 chain A and resid 0 or chain A and resid 1 or chain B and resid 2 or "
        "chain L and resid 12 \n"
        "mol modcolor 0 0 ColorID 0 \n"
        "mol modstyle 0 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \n"
        "mol modselect 1 0 chain B and resid 0 or chain B and resid 1 or chain C and resid 2 or "
        "chain M and resid 12 \n"
        "mol modcolor 1 0 ColorID 1 \n"
        "mol modstyle 1 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \nmol modselect 2 0 chain C and resid 0 or chain C and resid 1 or "
        "chain D and resid 2 or chain N and resid 12 or chain N and resid 13 \n"
        "mol modcolor 2 0 ColorID 2 \nmol modstyle 2 0 Trace 0.300000 12.000000 \n\n"
        "mol addrep 0 \n"
        "mol modselect 3 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 or "
        "chain N and resid 12 \nmol modcolor 3 0 ColorID 2 \n"
        "mol modstyle 3 0 VDW 1.000000 12.000000 \n\n"
    )
    assert expected_output_string == test_output_string


def test_vmd_visualize_script():
    """
    testing the generate vvmd_visualize_script function
    """

    #############################
    test_topology = {
        "topologies": [
            {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "rigid"}
                }
            }
        ]
    }
    test_topology["topologies"][0]["domains"].update(
        {1: {"chain-ids": ["B 0", "B 1", "C 2", "M 12"], "type": "rigid"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {2: {"chain-ids": ["C 0", "C 1", "D 2", "N 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["ndomains"] = 3
    test_topology["topologies"][0]["full-move"] = False

    test_topology["topologies"][0]["domains"][0]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][1]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][2]["nbeads"] = 3
    test_topology["definitions"] = convert.DEFAULT_DEFINITIONS
    ##############################
    # assigning certain test-values
    pdb_filename = "test.pdb"
    trajectory_filename = "test.xtc"
    cplx_filename = "test.cplx"
    coloring = "random"
    step = 1
    ###############################
    test_output_string = vis.vmd_visualize_script(
        pdb_filename, trajectory_filename, test_topology, cplx_filename, coloring, step
    )
    expected_output_string = (
        "# file generated: 14:18:46.890530 \n"
        "# topology-input: test.cplx \n# pdb-input: test.pdb \n"
        "# trajectory-input: test.xtc \n"
        " #################################### \n"
        "# Please be sure, that vmd operates in the directory, that your files are located...\n"
        " #################################### \n"
        " #################################### \n"
        "# opening file \n #################################### \n\n"
        " mol new test.pdb \n# Using plugin xtc for coordinates from file test.xtc\n"
        "mol addfile test.xtc type xtc first 0 last -1 step 1 waitfor -1 0\n"
        "animate style Loop\n#adding representation to molecule 0 \n"
        " #################################### \n"
        "mol addrep 0 \nmol modselect 0 0 chain A and resid 0 or chain A and resid 1 or "
        "chain B and resid 2 or chain L and resid 12 \n"
        "mol modcolor 0 0 ColorID 0 \nmol modstyle 0 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \n"
        "mol modselect 1 0 chain B and resid 0 or chain B and resid 1 or chain C and resid 2 or "
        "chain M and resid 12 \nmol modcolor 1 0 ColorID 1 \n"
        "mol modstyle 1 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \nmol modselect 2 0 chain C and resid 0 or chain C and resid 1 or "
        "chain D and resid 2 or chain N and resid 12 or chain N and resid 13 \n"
        "mol modcolor 2 0 ColorID 2 \nmol modstyle 2 0 Trace 0.300000 12.000000 \n\n"
        "mol addrep 0 \n"
        "mol modselect 3 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 or "
        "chain N and resid 12 \nmol modcolor 3 0 ColorID 2 \n"
        "mol modstyle 3 0 VDW 1.000000 12.000000 \n\n#hiding automatically generated all-rep\n"
        "# delrep syntax:  #rep #mol \nmol delrep 4 0 "
        "\npbc box"
    )
    assert_array_equal(
        expected_output_string.split("\n")[-34:], test_output_string.split("\n")[-34:]
    )


def test_visualize_unified():
    """
    Testing the vmd_visualize_vdw_particle-function
    """
    #############################
    test_topology = {
        "topologies": [
            {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "rigid"}
                }
            }
        ]
    }
    test_topology["topologies"][0]["domains"].update(
        {1: {"chain-ids": ["B 0", "B 1", "C 2", "M 12"], "type": "rigid"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {2: {"chain-ids": ["C 0", "C 1", "D 2", "N 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["ndomains"] = 3
    test_topology["topologies"][0]["full-move"] = False

    test_topology["topologies"][0]["domains"][0]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][1]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][2]["nbeads"] = 3
    ##############################
    # assigning certain test-values
    pdb_filename = "test.pdb"
    trajectory_filename = "test.xtc"
    cplx_filename = "test.cplx"
    step = 1
    ###############################
    test_output_string = vis.vmd_visualize_unified(
        pdb_filename, trajectory_filename, cplx_filename, step
    )
    expected_output_string = (
        "# file generated: 11:18:48.460495 \n"
        "# topology-input: test.cplx \n"
        "# pdb-input: test.pdb \n"
        "# trajectory-input: test.xtc \n "
        "#################################### \n"
        "# Please be sure, that vmd operates in the directory, that your files are located...\n "
        "#################################### \n "
        "#################################### \n"
        "# opening file \n #################################### \n\n"
        " mol new test.pdb \n"
        "# Using plugin xtc for coordinates from file test.xtc\n"
        "mol addfile test.xtc type xtc first 0 last -1 step 1 waitfor -1 0\n"
        "animate style Loop\n"
        "#adding representation to molecule 0 \n "
        "#################################### \n"
        "mol addrep 0 \n"
        "mol modselect 0 0 all \n"
        "mol modcolor 0 0 ColorID 1 \nmol modstyle 0 0 VDW 1.0 12.000000 \n\n"
        "#hiding automatically generated all-rep\n"
        "# delrep syntax:  #rep #mol \n"
        "mol delrep 1 0 "
        "\npbc box"
    )
    assert_array_equal(
        expected_output_string.split("\n")[1:], test_output_string.split("\n")[1:]
    )


def test_vmd_gaussian_rep():
    """
    testing vmd_gaussian_rep with random and with domain-specific coloring
    """
    # creating test topology
    #############################
    test_topology = {
        "topologies": [
            {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "rigid"}
                }
            }
        ]
    }
    test_topology["topologies"][0]["domains"].update(
        {1: {"chain-ids": ["B 0", "B 1", "C 2", "M 12"], "type": "rigid"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {2: {"chain-ids": ["C 0", "C 1", "D 2", "N 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {3: {"chain-ids": ["D 0", "D 1", "E 2", "O 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["ndomains"] = 4
    test_topology["topologies"][0]["full-move"] = False

    test_topology["topologies"][0]["domains"][0]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][1]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][2]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][3]["nbeads"] = 3
    ##############################
    # assigning certain test-values and creating test-string in necessary loop!
    coloring = "random"
    out = StringIO()
    bead_counter = 1
    rep_counter = 0
    top = test_topology["topologies"][0]
    domains = top["domains"]
    i = 2
    dom = domains[i]
    chainids = dom["chain-ids"]
    rep_number = rep_counter
    first_bead = bead_counter
    last_bead = first_bead + dom["nbeads"]
    bead_counter = last_bead + 1
    rigid_rep_output_string, rep_counter = vis.vmd_gaussian_rep(
        chainids, coloring, rep_counter, rep_number, i
    )
    out.write(rigid_rep_output_string)
    # Trying the domain-coloring scheme, as well the correct increasing of the rep_number is tested
    rep_number = rep_counter
    coloring = "domain"
    i = 3
    rigid_rep_output_string, rep_number = vis.vmd_gaussian_rep(
        chainids, coloring, rep_counter, rep_number, i
    )
    out.write(rigid_rep_output_string)
    test_output_string = out.getvalue()
    out.close()
    expected_output_string = (
        "mol addrep 0 \n"
        "mol modselect 0 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 "
        "or chain N and resid 12 or chain N and resid 13 \n"
        "mol modcolor 0 0 ColorID 2 \n"
        "mol modstyle 0 0 Trace 0.300000 12.000000 \n\n"
        "mol addrep 0 \n"
        "mol modselect 1 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 or "
        "chain N and resid 12 \n"
        "mol modcolor 1 0 ColorID 2 \n"
        "mol modstyle 1 0 VDW 1.000000 12.000000 \n\n"
        "mol addrep 0 \n"
        "mol modselect 2 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 or "
        "chain N and resid 12 or chain N and resid 13 \n"
        "mol modcolor 2 0 ColorID 0 \n"
        "mol modstyle 2 0 Trace 0.300000 12.000000 \n\n"
        "mol addrep 0 \n"
        "mol modselect 3 0 chain C and resid 0 or chain C and resid 1 or chain D and resid 2 or "
        "chain N and resid 12 \n"
        "mol modcolor 3 0 ColorID 0 \n"
        "mol modstyle 3 0 VDW 1.000000 12.000000 \n\n"
    )
    assert expected_output_string == test_output_string


def test_vmd_rigid_rep():
    """
    testing vmd_rigid_rep with random and with domain-specific coloring.
    """
    # creating test topology
    #############################
    test_topology = {
        "topologies": [
            {
                "domains": {
                    0: {"chain-ids": ["A 0", "A 1", "B 2", "L 12"], "type": "rigid"}
                }
            }
        ]
    }
    test_topology["topologies"][0]["domains"].update(
        {1: {"chain-ids": ["B 0", "B 1", "C 2", "M 12"], "type": "rigid"}}
    )
    test_topology["topologies"][0]["domains"].update(
        {2: {"chain-ids": ["C 0", "C 1", "D 2", "N 12"], "type": "gaussian"}}
    )
    test_topology["topologies"][0]["ndomains"] = 3
    test_topology["topologies"][0]["full-move"] = False

    test_topology["topologies"][0]["domains"][0]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][1]["nbeads"] = 3
    test_topology["topologies"][0]["domains"][2]["nbeads"] = 3
    ##############################
    # assigning certain test-values and creating test-string in necessary loop!
    coloring = "random"
    out = StringIO()
    bead_counter = 1
    rep_counter = 0
    top = test_topology["topologies"][0]
    domains = top["domains"]
    i = 0
    dom = domains[i]
    chainids = dom["chain-ids"]
    rep_number = rep_counter
    first_bead = bead_counter
    last_bead = first_bead + dom["nbeads"]
    bead_counter = last_bead + 1
    rigid_rep_output_string, rep_counter = vis.vmd_rigid_rep(
        chainids, coloring, rep_counter, rep_number, i
    )
    out.write(rigid_rep_output_string)
    # Trying the domain-coloring scheme, as well the correct increasing of the rep_number is tested
    rep_number = rep_counter
    coloring = "domain"
    i = 1
    rigid_rep_output_string, rep_number = vis.vmd_rigid_rep(
        chainids, coloring, rep_counter, rep_number, i
    )
    out.write(rigid_rep_output_string)
    test_output_string = out.getvalue()
    out.close()
    expected_output_string = (
        "mol addrep 0 \n"
        "mol modselect 0 0 chain A and resid 0 or chain A and resid 1 or chain B and resid 2 or "
        "chain L and resid 12 \n"
        "mol modcolor 0 0 ColorID 0 \n"
        "mol modstyle 0 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
        "mol addrep 0 \n"
        "mol modselect 1 0 chain A and resid 0 or chain A and resid 1 or chain B and resid 2 "
        "or chain L and resid 12 \nmol modcolor 1 0 ColorID 10 \n"
        "mol modstyle 1 0 QuickSurf 1.000000 0.5 1.0 1.0 \n\n"
    )
    assert expected_output_string == test_output_string


# Following the preparation of the different possible cases for argument_parsing:
Arguments = collections.namedtuple(
    "Arguments",
    "cplx xtc vmd_top coloring run_vmd complexes_config vis_config output unify step",
)
args_console = Arguments(
    cplx="cplx",
    xtc="xtc",
    vmd_top="vmd_top",
    coloring="domain",
    run_vmd=False,
    complexes_config=None,
    vis_config=None,
    output=None,
    unify=None,
    step=2,
)
expected_console = "cplx", "vmd_top", "xtc", "domain", "cplx_output.vmd", False, None, 2

args_com_config = Arguments(
    cplx=None,
    xtc=None,
    vmd_top=None,
    coloring="domain",
    run_vmd=True,
    complexes_config=pjoin(DATA_DIR, "test_vis_com_config.conf"),
    vis_config=None,
    output=None,
    unify=None,
    step=1,
)
expected_com_config = (
    "cplx.cplx",
    "xtc_reference.pdb",
    "xtc.xtc",
    "domain",
    "cplx_output.vmd",
    True,
    None,
    1,
)

args_vis_config = Arguments(
    cplx=None,
    xtc=None,
    vmd_top=None,
    coloring="domain",
    run_vmd=None,
    complexes_config=None,
    vis_config=pjoin(DATA_DIR, "test_vis_vis_config.conf"),
    output=None,
    unify=None,
    step=1,
)
expected_vis_config = (
    "cplx.cplx",
    "vmd_top.pdb",
    "xtc.xtc",
    "random",
    "test.vmd",
    True,
    True,
    1,
)
# coloring key will be be overwritten by the vis-config, is set to test, that console input is ignored
# for the vis_conf-case


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (args_console, expected_console),
        (args_com_config, expected_com_config),
        (args_vis_config, expected_vis_config),
    ],
)
def test_argument_processing(test_input, expected):
    """testing the argument_processing function
    With different input arguments, to evaluate all possible input cases:
    -console input
    -vis-conf
    -com-conf

    see test_input above
    """
    assert vis.argument_processing(test_input) == expected


def test_write_to_file(tmpdir):
    """
        testing the write_to_file function
    """
    outputfilename = pjoin(tmpdir.dirname, "test_output.vmd")

    test_string = "This is only test-content"
    vis.write_to_file(test_string, outputfilename)
    with open(outputfilename, "r") as f:
        read_string = yaml.load(f)
    assert read_string == test_string


def test_check_files_exists(tmpdir):
    files = ["test.xtc", "test.pdb"]
    with tmpdir.as_cwd():
        # touch files
        for f in files:
            with open(f, "w") as fh:
                pass
        vis.check_files_exist(*files)
        with pytest.raises(RuntimeError):
            vis.check_files_exist("some", "other")


def test_check_file_size(tmpdir):
    limit = 10 / 1024  # x MB in units of GB
    with tmpdir.as_cwd():
        with open("small", "w") as fh:
            pass
        # https://stackoverflow.com/a/6497779/2207958
        with open("large", "w") as fh:
            fh.truncate(limit * 1024**3 + 1024)
        vis.check_file_size(["small", ], limit=limit)
        with pytest.warns(UserWarning):
            vis.check_file_size(["large", ], limit=limit)


# def test_vmd_runner():
#     ############
