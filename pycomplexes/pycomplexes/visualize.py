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
"""Module to create file, that includes instructions for vmd how to cope with
cplx input and to automatically display included groups. Loads pdb topology
also, as well as trajectory. yaml has to be installed!!

Notes
-----
* Concerning use with vmd:
- If not directly opened with argument -r the command "play" has to be used to
        open the output-file. (This equals the "load visualization state..."
        button in context menu "file")

- For "playing" output files vmd has to be newly started. At least there cannot
        have been a molecule loaded before... (This is due to reusing the mol
        index 0, which seems to be only reseted after restarting vmd...)

* Further notices:

- is currently only implemented for one topology, respectively one molecule, but
might easily be implemented for more

"""
from __future__ import absolute_import, print_function, division
from six.moves import cStringIO as StringIO

import yaml
import datetime
import os
import numpy as np
import platform
import six
import warnings

from .scripts import _ScriptMeta

ADDREP_FORMAT = (
    "mol addrep 0 \n"
    "mol modselect {rep_number} {mol_number} {resid_string} \n"
    "mol modcolor {rep_number} {mol_number} ColorID {colorid} \n"
    "mol modstyle {rep_number} {mol_number} {rep_type} {rep_param} \n\n"
)

HEADER_FORMAT = (
    "# file generated: {generation_time} \n"
    "# topology-input: {topology_inputfile} \n"
    "# pdb-input: {pdb_inputfile} \n"
    "# trajectory-input: {trajectory_inputfile} \n"
    " #################################### \n"
    "# Please be sure, that vmd operates in the directory, that your files are located...\n"
    " #################################### \n"
)

LOADTRAJ_FORMAT = (
    "# Using plugin xtc for coordinates from file {traj_file}\n"
    "mol addfile {traj_file} type {traj_type} first 0 last -1 step {step} waitfor -1 0\n"
    "animate style Loop\n"
)

LOAD_VMD_TOP_FORMAT = (
    " #################################### \n"
    "# opening file \n"
    " #################################### \n\n "
    "mol new {vmd_top} \n"
)

DEL_REP_FORMAT = (
    "#hiding automatically generated all-rep\n"
    "# delrep syntax:  #rep #mol \n"
    "mol delrep {rep_counter} 0 "
)

SHOW_SIM_BOX = "\npbc box"

FUNCTIONS = """
##################################
# Add our own drawing primitives #
##################################
proc vmd_draw_plane {mol z} {
    if {$mol=="top"} then {set mol [molinfo top]}
    set cell [lindex [pbc get -molid $mol] 0]

    set A [list 0 0 $z]
    set B [list [lindex $cell 0] 0 $z]
    set C [list [lindex $cell 0] [lindex $cell 1] $z]
    set D [list 0 [lindex $cell 1] $z]
    set E [list [expr .5 * [lindex $cell 0]] [lindex $cell 1] $z]

    graphics $mol triangle $A $B $E
    graphics $mol triangle $A $D $E
    graphics $mol triangle $C $B $E
}

proc vmd_draw_tube {mol x y r} {
    if {$mol=="top"} then {set mol [molinfo top]}
    set cell [lindex [pbc get -molid $mol] 0]
    set z [lindex $cell 2]
    set top [list $x $y 0]
    set bottom [list $x $y $z]
    graphics $mol cylinder $top $bottom radius $r resolution 15 filled no
}

"""


def vmd_sel_chainid(chainid):
    """gets the tuple of chain and resid for each bead in a certain chain-ids-list
    returns a selection string for a certain BEAD in a format like resid ___ in
    chain ___

    Parameters:
    ----------
    chainid : string-array
        Formatted string like 'A 9'

    Returns:
    ----------
    string
        String that can be understood by vmd as a selction of a specific residue
        in a specific chain

    """
    return "chain {} and resid {}".format(*chainid.split())


def vmd_selection_string(chainids):
    """Creates a VERY long selection string for a certain WHOLE "DOMAIN". Format
    like: resid ___ and chain ___ or resid ...

    Parameter:
    ----------
    chainids : string-array
        Is a list of Strings. Get it from a yaml dictionary.

    Returns:
    ----------
    string
        Long concatenation of single residue selection strings, concatenated
        with or (logical or for vmd)

    """
    if len(chainids) == 1:
        return vmd_sel_chainid(chainids[0])
    else:
        return (
            vmd_sel_chainid(chainids[0])
            + " or "
            + " or ".join("{}".format(vmd_sel_chainid(el)) for el in chainids[1:])
        )


def vmd_selection_string_trace(chainids):
    """creates a VERY long selection string for a certain WHOLE "DOMAIN". Format
    like: resid ___ and chain ___ or resid ... adds one more bead at the
    beginning and at the end for connection between rigid domains and gaussians
    domains

    Parameter:
    ----------
    chainids : string-array
         A list of Strings. Get it from a yaml dictionary.

    Returns:
    --------
        long concatenation of single residue selection strings, concatenated
        with or (logical or for vmd). Before the first and after the last
        residue is one additional residue selected, for linking domains with
        traces.

    """
    middle_sel = vmd_selection_string(chainids)
    end_sel = "chain {} and resid {}".format(
        chainids[-1].split()[0], int(chainids[-1].split()[1]) + 1
    )
    # catching the case: resid_start= (-1): vmd has no problem with
    # residue-selection 0, but it errors on selection of negative residue
    # numbers, The following if statement prevents the making of a selection
    # string including negative residue numbers.
    if int(chainids[0].split()[1]) - 1 >= 0:
        start_sel = "chain {} and resid {}".format(
            chainids[0].split()[0], int(chainids[0].split()[1]) - 1
        )
        return start_sel + " or " + middle_sel + " or " + end_sel
    else:
        return middle_sel + " or " + end_sel


def get_type_translater(definitions):
    """
    Parameter
    ---------
    definitions : dict
        simulation definition section of a topology

    Returns
    -------
    translator : dict
         dictionary to convert domain names to movement types
    """
    translator = {}
    for dom, params in six.iteritems(definitions["domains"]):
        translator[dom] = params["move"]
        # So now this is tied to the expectations that the cplx is produced by
        # convert. If one deviates from these then the whole thing falls apart.
        # This is potentially very brittle
        if dom == "gaussian":
            translator[dom] = dom
        if params["move"] == "membrane":
            translator[dom] = "{}-membrane".format(params["defaults"]["type"])
    return translator


def vmd_rep_gen(topology, coloring):
    """Goes through the given topology and creates representations from that. The
    domains will be in a certain coloring scheme (specified by COLORING
    input-param). Different domain types will be shown in different rep_types by
    vmd.

    Parameter:
    ----------
    topology : dictionary
        Dictionary from the cplx read by yaml
    coloring : string
        random or domain. Parsed as input argument for choosing color-scheme

    Returns:
    --------
    string
        The string containing all the interpretable stuff for vmd to create the
        representations (also HEADER)

    """
    out = StringIO()
    bead_counter = 1
    rep_counter = 0
    name2type = get_type_translater(topology["definitions"])
    for j, top in enumerate(topology["topologies"]):
        domains = top["domains"]
        for i, dom in six.iteritems(domains):
            chainids = dom["chain-ids"]
            rep_number = rep_counter
            first_bead = bead_counter
            last_bead = first_bead + dom["nbeads"]
            bead_counter = last_bead + 1
            dom_rep, rep_counter = vmd_domain_rep(
                name2type[dom["type"]],
                chainids,
                coloring,
                rep_counter,
                rep_number,
                i,
                dom,
            )
            out.write(dom_rep)
    # The Following has to be the LAST included instruction in the file !!
    # Because otherwise vmd will automatically generate an additional
    # representation including all atoms, displaying them as lines.
    out.write(DEL_REP_FORMAT.format(rep_counter=rep_counter))
    concatenated_string = out.getvalue()
    out.close()
    return concatenated_string


def vmd_domain_rep(type_of_dom, chainids, coloring, rep_counter, rep_number, i, dom):
    """
    function that writes different rep-styles for different domain-types.
    so far types: rigid and gaussian are included
    is meant to make it easier to include more coloring schemes.

    Parameters:
    -----------
    type_of_dom : string
        Value behind the type-key in sub-dictionary:domain
    chainids : string_array
        List of Strings. Get it from a yaml dictionary. is passed on
    coloring : string
        Chooser for the coloring scheme: random, or domain
    rep_counter : Integer
        These are counters/indices for outer loop
    rep_number : Integer
        These are counters/indices for outer loop
    i: Integer
        These are counters/indices for outer loop
    dom : dict
        complete domain description

    Returns:
    --------
    vmd_rep_string : string
    updated counters : int
    """
    if type_of_dom == "rigid":
        return vmd_rigid_rep(chainids, coloring, rep_counter, rep_number, i)
    elif type_of_dom == "gaussian":
        return vmd_gaussian_rep(chainids, coloring, rep_counter, rep_number, i)
    elif type_of_dom == "flat-membrane":
        return vmd_flat_membrane_rep(
            chainids, coloring, rep_counter, rep_number, i, dom
        )
    elif type_of_dom == "tube-membrane":
        return vmd_tube_membrane_rep(
            chainids, coloring, rep_counter, rep_number, i, dom
        )
    else:
        raise RuntimeError("Unkown domain type: {}".format(type_of_dom))


def vmd_flat_membrane_rep(chainids, coloring, rep_counter, rep_number, i, dom):
    rep = "draw plane {}\n".format(dom["zaxis"])
    return rep, rep_counter + 1


def vmd_tube_membrane_rep(chainids, coloring, rep_counter, rep_number, i, dom):
    rep = "draw tube {} {} {}\n".format(dom["x"], dom["y"], dom["radius"])
    return rep, rep_counter + 1


def check_files_exist(*files):
    all_isfile = np.array([os.path.isfile(f) for f in files])
    if not np.all(all_isfile):
        raise RuntimeError(
            "Couldn't find files: {}".format(np.array(files)[~all_isfile])
        )


def check_file_size(files, limit=10, message=""):
    for fname in files:
        size = os.stat(fname).st_size / 1024 ** 3
        if size > limit:
            warnings.warn(
                "{}\n{} is {} GB big.".format(message, fname, size), UserWarning
            )


def vmd_rigid_rep(chainids, coloring, rep_counter, rep_number, i):
    """
    creates rep-string for rigid domains
    is meant to make it easier to include more coloring schemes.

    Parameters:
    -----------
    type_of_dom : string
        Value behind the type-key in sub-dictionary:domain
    chainids : string_array
        List of Strings. Get it from a yaml dictionary. is passed on
    coloring : string
        Chooser for the coloring scheme: random, or domain
    rep_counter : Integer
        These are counters/indices for outer loop
    rep_number : Integer
        These are counters/indices for outer loop
    i: Integer
        These are counters/indices for outer loop

    Returns:
    --------
    vmd_rep_string : string
    updated counters : int
    """
    if coloring == "domain":
        colorid = 10
    elif coloring == "random":
        colorid = i % 32
    else:
        raise ValueError(
            "Coloring method '{}' not supported for gaussian domain".format(coloring)
        )
    rep_type = "QuickSurf"
    rep_param = "1.000000 0.5 1.0 1.0"
    rep = ADDREP_FORMAT.format(
        colorid=colorid,
        rep_number=rep_number,
        mol_number=0,
        resid_string=vmd_selection_string(chainids),
        rep_type=rep_type,
        rep_param=rep_param,
    )

    return rep, rep_counter + 1


def vmd_gaussian_rep(chainids, coloring, rep_counter, rep_number, i):
    """
    creates rep-string for gaussian domains
    is meant to make it easier to include more coloring schemes.

    Parameters:
    -----------
    chainids : string_array
        List of Strings. Get it from a yaml dictionary. is passed on
    coloring : string
        Chooser for the coloring scheme: random, or domain
    rep_counter : Integer
        These are counters/indices for outer loop
    rep_number : Integer
        These are counters/indices for outer loop
    i: Integer
        These are counters/indices for outer loop

    Returns:
    --------
    rep : string
    updated counters : int
    """
    out = StringIO()
    if coloring == "domain":
        colorid = 0
    elif coloring == "random":
        colorid = i % 32
    else:
        raise ValueError(
            "Coloring method '{}' not supported for gaussian domain".format(coloring)
        )
    rep_type = "Trace"
    rep_param = "0.300000 12.000000"
    out.write(
        ADDREP_FORMAT.format(
            colorid=colorid,
            rep_number=rep_number,
            mol_number=0,
            resid_string=vmd_selection_string_trace(chainids),
            rep_type=rep_type,
            rep_param=rep_param,
        )
    )
    rep_counter += 1
    rep_type = "VDW"
    rep_param = "1.000000 12.000000"
    rep_number = rep_counter
    out.write(
        ADDREP_FORMAT.format(
            colorid=colorid,
            rep_number=rep_number,
            mol_number=0,
            resid_string=vmd_selection_string(chainids),
            rep_type=rep_type,
            rep_param=rep_param,
        )
    )
    rep = out.getvalue()
    return rep, rep_counter + 1


def vmd_visualize_script(
    pdb_filename, trajectory_filename, topology, cplx_filename, coloring, step
):
    """Function that creates the string in a syntax, that vmd can interpret. This
    string can be saved to a text/data file.

    Includes a comment-header. The very last printed lines prevent vmd from
    automatically creating an additional Representation.

    Parameters:
    -----------
    pdb_filename: String
        Containing name of pdb-file (parsed argument)
    trajectory_filename: String
        Containing name of trajectory file (parsed argument)
    topology: dictionary
        yaml-dictionary from cplx
    cplx_filename: String
        Containing name of cplx-file (parsed argument)
    coloring: String
        (parsed argument) for choosing coloring scheme

    Returns:
    --------
    String
        the complete string to be written to data-file

    """

    trajectory_type = os.path.splitext(trajectory_filename)[-1][1:]
    out = StringIO()
    out.write(
        HEADER_FORMAT.format(
            generation_time=datetime.datetime.time(datetime.datetime.now()),
            topology_inputfile=cplx_filename,
            pdb_inputfile=pdb_filename,
            trajectory_inputfile=trajectory_filename,
        )
    )
    out.write(FUNCTIONS)
    out.write(LOAD_VMD_TOP_FORMAT.format(vmd_top=pdb_filename))
    out.write(
        LOADTRAJ_FORMAT.format(
            traj_file=trajectory_filename, traj_type=trajectory_type, step=step
        )
    )
    out.write(
        "#adding representation to molecule 0 \n #################################### \n"
    )
    out.write(vmd_rep_gen(topology, coloring))
    out.write(SHOW_SIM_BOX)
    rep = out.getvalue()
    out.close()
    return rep


def vmd_visualize_unified(pdb_filename, trajectory_filename, cplx_filename, step):
    """Function that creates a string in a syntax that can be interpreted by vmd.
    This string can be saved to a text/data file.

    Includes a comment-header. The very last printed lines prevent vmd from
    automatically creating an additional Representation.
    This Function creates a very basic representation of the given topology:
    All simulation beads are displayed as vdw,particles of defined radius
    (default = 1.0)

    Parameters
    ----------
    pdb_filename: String
        Containing name of pdb-file (parsed argument)
    trajectory_filename: String
        Containing name of trajectory file (parsed argument)
    cplx_filename: String
        Containing name of cplx-file (parsed argument)

    Returns
    -------
    String
        the complete string to be written to data-file

    """
    vdw_radius = 1.0
    trajectory_type = os.path.splitext(trajectory_filename)[-1][1:]
    out = StringIO()
    out.write(
        HEADER_FORMAT.format(
            generation_time=datetime.datetime.time(datetime.datetime.now()),
            topology_inputfile=cplx_filename,
            pdb_inputfile=pdb_filename,
            trajectory_inputfile=trajectory_filename,
        )
    )
    out.write(LOAD_VMD_TOP_FORMAT.format(vmd_top=pdb_filename))
    out.write(
        LOADTRAJ_FORMAT.format(
            traj_file=trajectory_filename, traj_type=trajectory_type, step=step
        )
    )
    out.write(
        "#adding representation to molecule 0 \n #################################### \n"
    )
    colorid = 1
    rep_type = "VDW"
    rep_param = "{} 12.000000".format(vdw_radius)
    rep_number = 0
    out.write(
        ADDREP_FORMAT.format(
            colorid=colorid,
            rep_number=rep_number,
            mol_number=0,
            resid_string="all",
            rep_type=rep_type,
            rep_param=rep_param,
        )
    )
    rep_counter = 1
    out.write(DEL_REP_FORMAT.format(rep_counter=rep_counter))
    out.write(SHOW_SIM_BOX)
    rep = out.getvalue()
    out.close()
    return rep


def write_to_file(string, outputfilename):
    with open(outputfilename, "w") as outputfile:
        outputfile.write(string)


def vmd_runner(outputfilename):
    """Luxurious tool that executes vmd for you. Depends on your OS how far it gets.

    - In case of linux vmd has to be included in your .bashrc with alias="vmd"
    - In case of Windows the vmd-executable should be in the directory behind
      variable: vmd_directory
    - Has not been tested on mac-OS yet. But there is a certain chance, that it
      already works like this...

    Parameters:
    -----------
    outputfilename: str
        filename given to open ("play") file directly

    """
    print("You decided to run vmd directly.")
    if platform.system() == "Linux" or platform.system() == "Darwin":
        os.system("vmd -e " + outputfilename)
    else:
        print(
            "Starting vmd interactively isn't supported for platform '{}'. "
            "Please start vmd manually.".format(platform.system())
        )


def argument_processing(args):
    """Function the parsed arguments. Keeps main() clean. Mainly distinguishes
     between the direct parsing of arguments, the vis-config and the the extraction
     from the complexes-config.

    Parameters
    ----------
    args : namedtuple
        parsed arguments (from argparse)

    Returns
    -------
    processed arguments

    """
    # check whether complexes_config is chosen to provide cplx, xtc and vmd_top
    if args.complexes_config:
        print("Extracting Information from complexes-config: " + args.complexes_config)
        with open(args.complexes_config) as f:
            com_config = yaml.load(f)
        cplx = com_config["structure"]
        xtc = com_config["output"]["file"]
        vmd_top = os.path.splitext(com_config["output"]["file"])[0] + "_reference.pdb"
        if args.output:
            output_filename = args.output
        else:
            output_filename = os.path.splitext(cplx)[0] + "_output.vmd"
        # arguments that are not given in the complexes-config are simply passed through
        return (
            cplx,
            vmd_top,
            xtc,
            args.coloring,
            output_filename,
            args.run_vmd,
            args.unify,
            args.step,
        )
    if args.vis_config:  # seperate vis-conf given
        print("Extracting Information from visualization-config: " + args.vis_config)
        print("Note that any additional console argument will be ignored")
        with open(args.vis_config) as f:
            vis_config = yaml.load(f)
        cplx = vis_config["cplx"]
        xtc = vis_config["xtc"]
        vmd_top = vis_config["vmd_top"]
        if "coloring" in vis_config:
            coloring = vis_config["coloring"]
        else:
            coloring = "random"
        output_filename = vis_config["output_filename"]
        run_vmd = vis_config["run_vmd"]
        if "unify" in vis_config:
            unify = vis_config["unify"]
        else:
            unify = None
        if "step" in vis_config:
            step = vis_config["unify"]
        else:
            step = args.step
        return cplx, vmd_top, xtc, coloring, output_filename, run_vmd, unify, step
    else:  # No config given, direct argument input
        print("Reading console arguments.")
        if args.output:
            output_filename = args.output
        else:
            output_filename = os.path.splitext(args.cplx)[0] + "_output.vmd"
        return (
            args.cplx,
            args.vmd_top,
            args.xtc,
            args.coloring,
            output_filename,
            args.run_vmd,
            args.unify,
            args.step,
        )


class Visualize(six.with_metaclass(_ScriptMeta)):
    description = (
        "create vmd file from cplx-, pdb- and trajectoryfile. "
        "There is 3 options how to povide these files:\n"
        "direct argument-input,\n"
        "extract from complexes-simulation-config (-cc)\n"
        "use a special config for visualize (-vc)\n"
        "The latter provides the possibility to give any argument in there.\n"
        "Other optional inputs: change outputfilename (-o),"
        "run vmd directly(-r), change coloring scheme(-c), choose unified representation (-u)."
    )

    @staticmethod
    def parser(p):
        files = p.add_argument_group(
            "files",
            "files to be passed to visualize. Are only passed if neither "
            "-cc, nor -vc is given",
        )
        files.add_argument(
            "cplx",
            type=str,
            nargs="?",
            help="the name of the file that contains the topology (.cplx)",
        )
        files.add_argument(
            "vmd_top",
            type=str,
            nargs="?",
            help="the name of the file, taht contains the topology for vmd (f.e. pdb-file)",
        )
        files.add_argument(
            "xtc",
            type=str,
            nargs="?",
            help="the name of the file that contains the trajectory (.xtc)",
        )
        # forbid -c and unify together:
        representation = p.add_mutually_exclusive_group()
        representation.add_argument(
            "-c",
            "--coloring",
            type=str,
            default="random",
            help="coloring scheme can be chosen. Possible choices: random, domain \n"
            "- random (default):   every domain has a different color \n"
            "- domain:             gaussian domains have one color, rigid domains another. \n"
            " Flag can be overwritten with the alternative argument --unify",
        )
        representation.add_argument(
            "-u",
            "--unify",
            action="store_true",
            help="overwrites the coloring scheme and creates a single representation, where all "
            "beads are represented as vdw-particles with radius set to 1.0.",
        )
        p.add_argument(
            "-o",
            "--output",
            type=str,
            default=None,
            help="output_filename can be specified. Add an appropriate filename extension. "
            "Per default this would be 'CPLX'_output.vmd",
        )
        p.add_argument(
            "-r",
            "--run_vmd",
            help="If used: script directly starts vmd to open the newly created output-file."
            "Currently only fully available on Linux-Systems."
            "On windows at least starts vmd in script directory.",
            action="store_true",
        )
        p.add_argument(
            "-s",
            "--step",
            type=int,
            default=1,
            help="Every n-th frame should be loaded.",
        )
        # forbid to give two configs
        config = p.add_mutually_exclusive_group()
        config.add_argument(
            "-cc",
            "--complexes_config",
            type=str,
            default=None,
            help="Use the input-file of the complexes-simulation to extract the names of the:\n"
            "-cplx-file\n"
            "-vmd-top\n"
            "-xtc-file\n"
            "from there. The rest is still given by optional arguments",
        )
        config.add_argument(
            "-vc",
            "--vis_config",
            type=str,
            default=None,
            help="Use an extra config file only created for visualize to extract info from there. "
            "This config can contain all parsed arguments. It has to be in a yaml-dict format."
            "All additionally given Arguments (in the console) will be ignored!",
        )

    @staticmethod
    def main(args):
        """Main program. Calls functions to process and parse arguments. Opens topology
            file and calls function to create output string.

            Parameters:
            -----------
            args : namedtuple
                parsed arguments (from argparse)
            """

        cplx, vmd_top, xtc, coloring, output_filename, run_vmd, unify, step = argument_processing(
            args
        )
        check_files_exist(cplx, vmd_top, xtc)
        check_file_size(
            [xtc, ],
            message="Files may to big to load completely into VMD. Consider using '-s' option to skip frames.",
        )

        with open(cplx) as f:
            topology = yaml.load(f)
        print("Given topology: " + cplx)
        # check whether coloring was overwritten by unify-argument:
        if unify:
            print("You chose unified representation")
            output_string = vmd_visualize_unified(vmd_top, xtc, cplx, step)
        else:
            if coloring == "domain":
                print("You chose domain-type-specific coloring")
            #     output_string = vmd_visualize_script(vmd_top, xtc, topology, cplx, coloring)
            # else:
            output_string = vmd_visualize_script(
                vmd_top, xtc, topology, cplx, coloring, step
            )
        write_to_file(output_string, output_filename)
        print("Your output will be in file: " + output_filename)
        if run_vmd:
            vmd_runner(output_filename)


# # produces a generally okay output file, but not very user friendly in vmd
# probably better if the outputstring producing function was smarter,
# so that it recognizes series of residues in the same chain...
