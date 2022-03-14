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
from six.moves import range, zip

from functools import lru_cache

import copy
from os.path import isfile

import MDAnalysis as mda
import numpy as np
import six
import yaml
import abc

from .scripts import _ScriptMeta
from .forcefield import (
    get_forcefield_file,
    known_forcefields,
    read_forcefield,
    ListForcefieldsAction,
)


AMBER_ALIASES = {
    "CYM": "CYS",
    "HIP": "HIS",
    "HID": "HIS",
    "HIE": "HIS",
    "GLH": "GLU",
    "ASH": "ASP",
    "LYN": "LYS",
}

CHARGE_ASSIGNMENT = {
    "ASP": -1,
    "GLU": -1,
    "CYM": -1,
    "LYS": 1,
    "ARG": 1,
    "HIP": 1,
    "HIS": .5,
}  # assigning zero charge is the default for resnames not specified here


_DEFAULT_DEFINITIONS = """
  domains:
    gaussian:
      defaults: {rotation: 0, translation: 0}
      move: rigid
    rigid:
      defaults: {rotation: 3.1415, translation: 3}
      move: rigid
    flat:
      move: membrane
      defaults: { type: flat, z0: 15, ePSI0: -30 }
    tube:
      move: membrane
      defaults: { type: tube, z0: 15, ePSI0: -30 }
  pair-interaction:
  - domain-type-pair: [default, default]
    function: LJH
  - domain-type-pair: [rigid, gaussian]
    function: None
  # Add membrane interactions here
  - domain-type-pair: [rigid, flat]
    function: Membrane
  - domain-type-pair: [gaussian, flat]
    function: None
  - domain-type-pair: [flat, flat]
    function: None
  - domain-type-pair: [rigid, tube]
    function: Membrane
  - domain-type-pair: [gaussian, tube]
    function: None
  - domain-type-pair: [tube, tube]
    function: None
"""

DEFAULT_DEFINITIONS = yaml.safe_load(_DEFAULT_DEFINITIONS)


def charges(residues, charge_assignment=CHARGE_ASSIGNMENT):
    """
    Set charges according to the nature of the residue

    Parameter
    ---------
    residues : mda.Residues
        group od residues
    charge_assignment : dict
        lookup resname -> charge

    Return
    ------
    charges : np.ndarray
        charges for residues
    """
    charges = np.zeros(len(residues))
    for i, res in enumerate(residues):
        charges[i] = charge_assignment.get(res.resname, 0)
        # all not specifically declared bead-types will be neutral
    return charges


def chain_ids(ag):
    """
    Return a list of segment id and resids from an atom group

    Parameter
    ---------
    ag : mda.AtomGroup

    Return
    ------
    ids : list
        a list of segment id and resids

    """
    return ["{} {}".format(seg, res) for seg, res in zip(ag.segids, ag.resids)]


def normalize_resnames(resnames, aliases=AMBER_ALIASES):
    """
    Convert resnames from AMBER to KimHummer names. This allows to use
    more resnames to define charges and use the same LJ parameters.

    Parameters
    ----------
    resnames : list
        residue names from topology
    aliases : dict
        translates non-canonical resnames to canonical

    Returns
    -------
    canonified : list
        translated resnames

    """
    return [aliases.get(r, r) for r in resnames]


class AbstractDomain(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def to_cplx(self):
        pass


class ProteinDomain(AbstractDomain):
    """
    Protein Domain internal representation
    """

    def __init__(self, node, name, universe, id, randomize=False, box=None):
        self.type = node["type"]
        self.name = name
        self.id = id
        self.moves = node.get("mc-moves", dict())
        self.node = node
        ag = universe.select_atoms(node["selection"])
        if randomize:
            if box is None:
                raise RuntimeError("Need box info for random placement")
            point = [np.random.uniform(0, h) for h in box]
            ag.translate(-ag.centroid()).translate(point)
        self._metadata = {'remark': 'created with pycomplexes'}
        self._populate(ag)

    def _populate(self, ag):
        self.n_beads = ag.n_atoms
        self.beads = normalize_resnames(ag.resnames.tolist())
        self.charges = charges(ag.residues).tolist()
        self.positions = ag.positions.copy().tolist()
        self.chain_ids = chain_ids(ag)

    def to_cplx(self):
        return {
            "beads": self.beads,
            "chain-ids": self.chain_ids,
            "charges": self.charges,
            "coordinates": self.positions,
            "name": self.name,
            "type": self.type,
            "nbeads": self.n_beads,
            "meta-data": self._metadata,
            "mc-moves": self.moves,
        }


class MembraneDomain(AbstractDomain):
    """
    Membrane Domain internal representation
    """

    def __init__(self, node, name, id):
        node = copy.deepcopy(node)
        self.type = node["type"]
        node.pop("type")
        self.args = node
        self.id = id
        self.name = name

    def to_cplx(self):
        base = {"name": self.name, "type": self.type, "chain-ids": ["Z 1"], "nbeads": 1}
        base.update(self.args)
        return base


class AbstractTopology(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def to_cplx(self):
        pass


@lru_cache(maxsize=32)
def universe_cache(filename):
    return mda.Universe(filename)


def atom2beadid(atom):
    return "{} {}".format(atom.segid, atom.resid)


class ProteinTopology(AbstractTopology):
    """
    protein topology internal representation
    """

    def __init__(self, node, name, randomize=False, box=None):
        u = universe_cache(node["coordinate-file"])
        self.name = name
        self.domains = []
        self.move = node.get("move", False)

        for i, name in enumerate(sorted(six.iterkeys(node["domains"]))):
            n = node["domains"][name]
            self.domains.append(ProteinDomain(n, name, u, i, randomize, box))

        # guess connections for gaussian domains
        guessed_connections = []
        flag_linkers_plus_random_placement=False
        for dom in self.domains:
            n = dom.node
            if dom.type == "gaussian":
                if randomize:
                    flag_linkers_plus_random_placement=True
                start_atom = u.select_atoms(n["start_connection"][1])
                if start_atom.n_atoms != 1:
                    raise RuntimeError(
                        "Start connection for gaussian domain '{}' selects more then one atom. Selection: {}".format(
                            dom.name, n["start_connection"][1]
                        )
                    )
                start_bead = atom2beadid(start_atom[0])

                end_atom = u.select_atoms(n["end_connection"][1])
                if end_atom.n_atoms != 1:
                    raise RuntimeError(
                        "End connection for gaussian domain '{}' selects more then one atom. Selection: {}".format(
                            dom.name, n["end_connection"][1]
                        )
                    )
                end_bead = atom2beadid(end_atom[0])

                con3 = {
                    "domain-a": [n["start_connection"][0], start_bead],
                    "domain-b": [n["end_connection"][0], end_bead],
                    "type": "gaussian",
                    "params": {"N": dom.n_beads, "bond-length": 3.81},
                }
                guessed_connections.append(con3)
                dom._metadata['connection_id'] = len(guessed_connections) - 1
        if flag_linkers_plus_random_placement:
            print("WARNING: You used the randomized starting position in combination with Gaussian linker domains!")
            print("The random placement does not respect the distance restraints, " +
                  "which are implied by the Gaussian linker domains")
            print("We strongly recommend not to combine random initial placement with Gaussian Linker domains!")

        # update connections
        connections = node.get("connections", dict())
        offset = len(connections)
        for i, c in enumerate(guessed_connections):
            connections[i + offset] = c
        self.connections = connections

    def to_cplx(self):
        cplx = {
            "connections": self.connections,
            "name": self.name,
            "ndomains": len(self.domains),
            "full-move": self.move,
            "translation": 20,
            "rotation": .5,
            "domains": {d.id: d.to_cplx() for d in self.domains},
        }
        return cplx


class MembraneTopology(AbstractTopology):
    """
    membrane topology internal representation
    """

    def __init__(self, node, name):
        self.name = name
        self.domains = [MembraneDomain(node["membrane"], name, 0)]

    def to_cplx(self):
        return {
            "connections": [],
            "name": self.name,
            "ndomains": 1,
            "full-move": False,
            "domains": {self.domains[0].id: self.domains[0].to_cplx()},
        }


def expand(topology):
    """
    expand certain keywords like 'copies' in the topology.

    Parameters
    ----------
    topology : dict
       Topology in TOP file format

    Returns
    -------
    dict
       Topology in TOP file format
    """
    new_tops = {}

    # expand 'copies' keyword
    for name, t in six.iteritems(topology):
        if "copies" in t:
            for i in range(t["copies"]):
                new_tops["{}-{}".format(name, i)] = copy.deepcopy(t)
        else:
            new_tops[name] = t
    return new_tops


def parse(node, randomize=False, box=None):
    """
    Parameters
    ----------
    node : dict
        topology in TOP file format
    randomize : bool, optional
        random placement of topologies
    box : array-like, optional
        box information, needed for randomize==True

    Returns
    -------
    list of topology objects
    """
    universe_cache.cache_clear()
    t = []
    for name in sorted(six.iterkeys(node)):
        topology = node[name]
        if ("type" not in topology) or (topology["type"] == "protein"):
            t.append(ProteinTopology(topology, name, randomize, box))
        elif topology["type"] == "membrane":
            t.append(MembraneTopology(topology, name))
    return t


def normalize(tops):
    """
    Normalize topology.

    Normalization tasks:
      - ensure each domain has a unique id

    Parameters
    ----------
    tops : list
       list of topology objects

    Returns
    -------
    tops : list
       list of topology objects
    """
    i = 0
    for t in tops:
        for d in t.domains:
            d.id = i
            i += 1
    return tops


def to_cplx(tops):
    """
    convert list of topology to the CPLX file format
    """
    return [t.to_cplx() for t in tops]


def convert(
    top,
    forcefield="KimHummer",
    alpha=1,
    debye_length=10,
    dielectric_constant=80,
    random=False,
):
    """
    convert a topology to a cplx

    Parameters
    ----------
    top : dict
        topology description in TOP file format
    random : bool, optional
        decide if topologies should be placed randomly in box
    forcefield : str, optional
        forcefield to use
    alpha : float, optional
        softcore scaling

    Returns
    -------
    dict in CPLX format
    """
    box = top["box"]
    topol = expand(top["topology"])
    topol = parse(topol, randomize=random, box=box)
    topol = normalize(topol)
    cplx = to_cplx(topol)

    ff = read_forcefield(get_forcefield_file(forcefield))
    ff["alpha"] = alpha
    ff["debye-length"] = debye_length
    ff["dielectric-constant"] = dielectric_constant

    return {
        "box": box,
        "definitions": DEFAULT_DEFINITIONS,
        "topologies": cplx,
        "forcefield": ff,
    }


# generate a parser for this class
class Convert(six.with_metaclass(_ScriptMeta)):
    description = "Convert Topology files to cplx"

    @staticmethod
    def parser(p):
        p.add_argument("top", type=str, help="topology file")
        p.add_argument("cplx", type=str, help="cplx file to write")
        p.add_argument(
            "--forcefield",
            type=str,
            default="KimHummer",
            help="forcefield to use. Can be a included name or a relative file-path",
        )
        p.add_argument(
            "--alpha", type=float, default=1, help="alpha value for softcore potential"
        )
        p.add_argument("--debye_length", type=float, default=10)
        p.add_argument("--dielectric_constant", type=float, default=80)
        p.add_argument(
            "--random", action="store_true", help="random placement of the topologies in the box"
        )
        p.add_argument(
            "--list_forcefields",
            action=ListForcefieldsAction,
            help="show known forcefields",
            nargs=0,
        )

    @staticmethod
    def main(args):
        if not isfile(args.top):
            raise IOError("File does not exist: {}".format(args.top))

        with open(args.top) as f:
            top = yaml.safe_load(f)

        cplx = convert(
            top,
            random=args.random,
            forcefield=args.forcefield,
            alpha=args.alpha,
            dielectric_constant=args.dielectric_constant,
            debye_length=args.debye_length,
        )

        with open(args.cplx, "w") as f:
            yaml.dump(cplx, f)
