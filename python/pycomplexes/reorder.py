# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from __future__ import absolute_import, print_function

import MDAnalysis as mda
import numpy as np
import os

import six
from .scripts import _ScriptMeta


def create_selection_strings(atoms):
    """Functions that creates the selection string for a given pdb-topology.

    The returned selection string is ordered in a way, that the resids of a
    segments are ordered together by resid, so that no segments are mixed in
    the reordered pdb-file.

    Parameters
    ----------
    atoms : mda.Universe or mda.Atomgroup
        This is the given MDAnalysis-Universe

    Returns
    -------
    selections_strings : list of str
       This is the selection string to be given to the MDAnalysis
       atom-selector (u.atoms.select_atoms(...))

    """
    # always use atomgroup regardless of input type (universe or atomgroup)
    ag = atoms.atoms

    segids_unique = [
        "{}".format(seg) for seg in np.sort(np.unique([s.segid for s in ag.segments]))
    ]
    selection_strings = []
    for segid in segids_unique:
        resids = np.sort(ag.select_atoms("segid {}".format(segid)).resids)
        for resid in resids:
            selection_strings.append("segid {} and resid {}".format(segid, resid))
    return selection_strings


# generate a parser for this class
class Reorder(six.with_metaclass(_ScriptMeta)):
    description = "reorder pdb files, so that linkers are really linked after visualization (with visualize.py)"

    @staticmethod
    def parser(p):
        # print("adding to parser")
        p.add_argument("pdbfile", type=str, help="pdb file to be reordered")
        p.add_argument(
            "--outputfile",
            "-o",
            type=str,
            default=None,
            help="optional change of outputfilename. Has to have an appropriate filename extension for "
            "MDA-topologies (usually .pdb, but there is a lot of other possibilities. Can be looked"
            " up at MDA-homepage.) If nothing given, default: INPUTPDB_reordered.pdb",
        )

    @staticmethod
    def main(args):
        """Main-function that processes the parsed arguments and writes a re-
        ordered pdb-file using MDAnalysis.

        Notes
        -----
            Outputfilename has to have an appropriate filename extension for MDA-
            topologies (usually .pdb, but there is a lot of other possibilities.
            Can be looked up at MDA-homepage.)

        Parameters
        ----------
        args : argparse-tuple
            parsed arguments (pdbfile, outputfile)

        """
        u = mda.Universe(args.pdbfile)
        print("Given pdb to be reordered: " + args.pdbfile)
        selection_strings = create_selection_strings(u)
        if args.outputfile:
            output_filename = args.outputfile
        else:
            output_filename = os.path.splitext(args.pdbfile)[0] + "_reordered.pdb"
        print("Your output will be in file: " + output_filename)
        core = u.atoms.select_atoms(*selection_strings)
        core.write(output_filename)
        print("Writing successful.")
