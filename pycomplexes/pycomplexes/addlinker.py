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


class LinkEnd(object):
    def __init__(self, name, bead):
        self.name = name
        self.bead = bead


class LinkerConnection(object):
    def __init__(self, linker_name, start_bead, stop_bead, start, stop):
        self.linker = (linker_name, start_bead, stop_bead)
        self.start = start
        self.stop = stop

    def __repr__(self):
        return "LinkerCon: {}".format(self.linker[0])


def get_linker_connections(top):
    """parse the connections into a more suitable data structure
    """
    linker_connections = []
    # collect all gaussian domains
    connections = []
    for con in six.itervalues(top["connections"]):
        if con["type"] == "gaussian":
            connections.append(con)
    # now we look for the domains connected to those
    for con in connections:

        start_name = con["domain-a"][0]
        end_name = con["domain-b"][0]
        for c in six.itervalues(top["connections"]):
            if c == con:
                continue
            if c["domain-a"][0] == start_name or c["domain-b"][0] == start_name:
                start = c
            elif c["domain-a"][0] == end_name or c["domain-b"][0] == end_name:
                end = c

        if start["domain-a"][0] == start_name:
            linker, linker_start = start["domain-b"]
        else:
            linker, linker_start = start["domain-a"]

        if end["domain-a"][0] == end_name:
            linker_end = end["domain-b"][1]
        else:
            linker_end = end["domain-a"][1]

        # yapf: disable
        lc = LinkerConnection(linker, linker_start, linker_end,
                              LinkEnd(start_name, con['domain-a'][1]),
                              LinkEnd(end_name, con['domain-b'][1]))
        # yapf: enable
        linker_connections.append(lc)

    return linker_connections


class Link(object):
    def __init__(self, start, stop, ag, name, lc=None):
        self.start = start
        self.stop = stop
        self.ag = ag
        self.name = name
        self._lc = lc

    @classmethod
    def empty(cls, link_connection):
        name = link_connection.linker[0]
        return cls(start=None, stop=None, ag=None, name=name, lc=link_connection)


def get_linkers(cplx, sim):
    """get gaussian linkers as atom groups and the start and end positions placed
    on a rigid domain

    """
    all_linkers = []
    total_beads = 0
    for top in cplx["topologies"]:
        # get some raw linker connections parsed into a more suitable data
        # structure
        linker_connections = get_linker_connections(top)
        linkers = [Link.empty(lc) for lc in linker_connections]

        # go through domains and update linker structures when we see a
        # suitable domain. This avoid reiterating over the domains. Can be
        # avoided in the future if we store domain id's in the linker
        # connections
        for dom_id, dom in six.iteritems(top["domains"]):
            dom_name = dom["name"]
            nbeads = dom["nbeads"]
            ag = sim.atoms[total_beads : total_beads + nbeads]

            # Do I need information of this domain
            is_part_of_linker = False
            found_links = []
            for link in linkers:
                lc = link._lc
                possible_names = [lc.linker[0], lc.start.name, lc.stop.name]
                if dom_name in possible_names:
                    is_part_of_linker = True
                    found_links.append(link)

            # update Link entry if yes
            if is_part_of_linker:
                for link in found_links:
                    if dom_name == link._lc.start.name:
                        bead = link._lc.start.bead
                        atom = ag.select_atoms(
                            "segid {} and resid {}".format(*bead.split())
                        )
                        link.start = atom[0]
                    elif dom_name == link._lc.stop.name:
                        bead = link._lc.stop.bead
                        atom = ag.select_atoms(
                            "segid {} and resid {}".format(*bead.split())
                        )
                        link.stop = atom[0]
                    elif dom_name == link._lc.linker[0]:
                        link.ag = ag

            total_beads += dom["nbeads"]

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
            config = yaml.load(f)

        cplx = config["structure"]
        traj = config["output"]["file"]
        top = "{}_reference.pdb".format(splitext(traj)[0])

        util.check_file_exists(cplx)
        util.check_file_exists(traj)

        with open(cplx) as fh:
            cplx = yaml.load(fh)

        sim = mda.Universe(top, traj)
        addlinker(cplx, sim, args.out, args.start, args.stop, args.step)
