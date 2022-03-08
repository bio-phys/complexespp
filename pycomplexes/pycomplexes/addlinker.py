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

from os.path import splitext
from collections import namedtuple

import MDAnalysis as mda
import numpy as np
import six
import yaml
from tqdm import tqdm
import warnings

from .scripts import _ScriptMeta
from . import util
from .linker import gauss


def get_type_ags(cplx, sim, typename):
    tops = cplx["topologies"]
    ags = []
    total_beads = 0
    for t in tops:
        for dom_id, d in six.iteritems(t["domains"]):
            nbeads = d["nbeads"]
            if d["type"] == typename:
                a = sim.atoms[total_beads : total_beads + nbeads]
                ags.append(a)
            # update counter of total beads seen so far
            total_beads += nbeads
    return ags


class Link(object):
    def __init__(self, start, stop, ag, name):
        self.start = start
        self.stop = stop
        self.ag = ag
        self.name = name


def get_linkers(cplx, sim):
    """get gaussian linkers as atom groups and the start and end positions placed
    on a rigid domain

    """
    all_linkers = []
    total_beads = 0
    for top in cplx["topologies"]:
        connections = top['connections']
        linkers = []
        for dom in six.itervalues(top["domains"]):
            # running update of corret atom selection for current domain
            nbeads = dom['nbeads']
            ag = sim.atoms[total_beads: total_beads + nbeads]
            total_beads += nbeads

            # create Link object if gaussian domain
            if dom['type'] != 'gaussian':
                continue
            c = connections[dom['meta-data']['connection_id']]
            start = sim.select_atoms('segid {} and resid {}'.format(*c['domain-a'][1].split()))[0]
            stop = sim.select_atoms('segid {} and resid {}'.format(*c['domain-b'][1].split()))[0]
            linkers.append(Link(start, stop, ag, dom['name']))

        all_linkers.extend(linkers)

    return all_linkers


def addlinker(cplx, sim, out, start=None, stop=None, step=None):
    """
    add linker domains for a given trajectory
    """
    try:
        from .linker import relax
    except Exception as e:
        relax = None
        warnings.warn("Couldn't import relax algorithm for linkers: {}".format(str(e)))

    linkers = get_linkers(cplx, sim)
    rigid = sum(get_type_ags(cplx, sim, "rigid"))
    box = sim.dimensions[:3]

    with mda.Writer(out, n_atoms=sim.atoms.n_atoms) as w:
        for ts in tqdm(sim.trajectory[start:stop:step]):
            other_coords = rigid.positions
            for link in linkers:
                xyz = gauss.grow_gauss_chain(
                    start=link.start.position,
                    stop=link.stop.position,
                    nbeads=link.ag.n_atoms + 2,
                    other_coords=other_coords,
                    box=box,
                )
                # first parameter is the energy relaxation I don't need
                if relax is not None:
                    _, xyz = relax.relax(xyz, other_coords=other_coords, box=box)
                # ignore start and stop they are not part of this atomgroup
                link.ag.positions = xyz[1:-1]

                # add to list of stuff I don't want to overlap with
                other_coords = np.concatenate((other_coords, xyz))

            w.write(sim.atoms)


# generate a parser for this class
class Addlinker(six.with_metaclass(_ScriptMeta)):
    description = "Generate Linker configurations"

    @staticmethod
    def parser(p):
        p.add_argument("config", type=str, help="config file of simulation")
        p.add_argument("out", type=str, help="name of new trajectory")
        p.add_argument("--start", type=int, default=None)
        p.add_argument("--stop", type=int, default=None)
        p.add_argument("--step", type=int, default=None)

    @staticmethod
    def main(args):
        util.check_file_exists(args.config)

        with open(args.config) as f:
            config = yaml.safe_load(f)

        cplx = config["structure"]
        traj = config["output"]["file"]
        top = "{}_reference.pdb".format(splitext(traj)[0])

        util.check_file_exists(cplx)
        util.check_file_exists(traj)

        with open(cplx) as fh:
            cplx = yaml.safe_load(fh)

        sim = mda.Universe(top, traj)
        addlinker(cplx, sim, args.out, args.start, args.stop, args.step)
