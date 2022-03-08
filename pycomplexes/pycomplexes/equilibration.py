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
from __future__ import absolute_import, print_function
import MDAnalysis as mda
import yaml

import six
from .scripts import _ScriptMeta

# TODO: be a bit more specific about filter MDA warnings
import warnings

warnings.filterwarnings("ignore")


def update_coordinates(system, universe):
    """update the coordinates in system, with universe-coordinates

    Parameters
    ----------
    system : dict
        dict loaded from yaml-file (.cplx)
    universe : object(MDA-universe)
        MDA universe. Should nly the desired frame (in most cases probably last)

    Returns
    -------
    system : dict
        yaml-dict with updated coordinates

    """
    index_bead_low = 0
    index_bead_high = 0
    # renaming the universe (maybe bad habit, but should improve readability)
    u = universe
    # processing each topology
    for cur_topology in system["topologies"]:
        # processing the universe with equilibrated coordinates replacing now
        # the non equilibrated coordinates with the equilibrated ones
        for key in sorted(cur_topology["domains"]):
            # update the index of the last bead of current top
            index_bead_high += cur_topology["domains"][key]["nbeads"]
            # replace the coordinates
            cur_topology["domains"][key]["coordinates"] = u.atoms.positions[
                index_bead_low:index_bead_high
            ].tolist()
            # update the index of the last bead of the current top
            index_bead_low += cur_topology["domains"][key]["nbeads"]
    return system


def argument_processing(args):
    """Function to distinguish, whether input should be extracted from
    complexes config, or from console input

    Parameters
    ----------
    args : Namedtuple
        parsed arguments

    Returns
    -------
    processed arguments

    """
    if args.complexes_config:
        with open(args.complexes_config, "r") as f:
            config = yaml.safe_load(f)
        equilibrate_cplx = config["structure"]
        trajectory = config["output"]["file"]
        reference_pdb = trajectory[:-4] + "_reference.pdb"
        return equilibrate_cplx, trajectory, reference_pdb, args.cplx, args.frame
    else:
        return (
            args.equilibrate_cplx,
            args.trajectory,
            args.reference_pdb,
            args.cplx,
            args.frame,
        )


# generate a parser for this class
class Equilibration(six.with_metaclass(_ScriptMeta)):
    description = "reuse coordinates from equilibration run"

    @staticmethod
    def parser(p):
        files = p.add_argument_group(
            "files",
            "files to be passed to equilibration. Are only passed if -cc is not given",
        )
        files.add_argument(
            "equilibrate_cplx",
            type=str,
            nargs="?",
            help="cplx of equilibration simulation",
        )
        files.add_argument(
            "trajectory", type=str, nargs="?", help="equilibration trajectory"
        )
        files.add_argument(
            "reference_pdb", type=str, nargs="?", help="PDB of equilibration simulation"
        )
        # This is still positional
        p.add_argument("cplx", type=str, help="cplx file to write")
        p.add_argument(
            "-cc",
            "--complexes_config",
            type=str,
            help="complexes-config, from where information should be extracted.",
        )
        p.add_argument(
            "--frame", type=int, default=-1, help="frame to use, defaults to last"
        )

    @staticmethod
    def main(args):
        """Parse a generated cplx and equilibrate it

        Parameter
        ---------
        args : namedtuple
            parsed arguments (from argparse)

        Return
        ------
        file:
            in which the cplx topologies are equilibrated

        """
        equilibrate_cplx, trajectory, reference_pdb, cplx, frame = argument_processing(
            args
        )

        with open(equilibrate_cplx) as f:
            system = yaml.safe_load(f)

        # update the current system-coordinates to coordinates at specified frame (parsed argument)
        eq_u = mda.Universe(reference_pdb, trajectory)
        eq_u.trajectory[frame]

        system = update_coordinates(system, eq_u)

        # writing the updated cplx
        with open(cplx, "w") as f:
            yaml.dump(system, f)
