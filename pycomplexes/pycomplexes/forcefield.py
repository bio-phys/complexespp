# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from __future__ import absolute_import, print_function
import six
import yaml
from collections import defaultdict
from copy import deepcopy
from pkg_resources import resource_filename, resource_exists, resource_listdir
import argparse
import os

from .scripts import _ScriptMeta


def get_forcefield_file(name):
    if resource_exists(__name__, "forcefields/{}".format(name)):
        return resource_filename(__name__, "forcefields/{}".format(name))
    else:
        return name


def known_forcefields():
    return resource_listdir(__name__, "forcefields")


def read_forcefield(ff_file):
    with open(ff_file) as fh:
        return yaml.load(fh)


def scale_interaction(old_forcefield, scale, offset):
    # from Kim, Y. C., & Hummer, G. (2008).: Coarse-grained Models for
    # Simulations of Multiprotein Complexes: Application to Ubiquitin Binding
    # equation 5: \epsilon_{ij} = \lambda * (e_{ij} - e_0)
    energies = defaultdict(dict)
    # yapf: disable
    for a in old_forcefield['energies']:
        for b in old_forcefield['energies'][a]:
            energies[a][b] = scale * (old_forcefield['energies'][a][b] - offset)
    # yapf: enable
    new_ff = deepcopy(old_forcefield)
    new_ff["energies"] = energies
    return new_ff


def write_forcefield(ff, fname):
    if os.path.exists(fname):
        raise IOError("file '{}' already exists.".format(fname))
    with open(fname, "w") as fh:
        yaml.dump(ff, fh)


class ListForcefieldsAction(argparse.Action):
    """
    Action to list known forcefields with argparse

        p.add_argument(
            '--list_forcefields',
            action=ListForcefieldsAction,
            help='show known forcefields',
            nargs=0)

    """

    def __call__(self, parser, namespace, values, option_string=None):
        print("Known Forcefields:")
        for ff in known_forcefields():
            print("\t{}".format(ff))
        parser.exit()


class Forcefield(six.with_metaclass(_ScriptMeta)):
    description = """Quickly change forcefield parameters. Provide parameters
                     lambda (scale) and e_0 (offset) and scale and existing
                     forcefield (e_ij) to receive new forcefield (E_ij)
                     with :\n
                     E_ij = lambda * (e_ij - e_0)"""

    @staticmethod
    def parser(p):
        p.add_argument("forcefield", type=str, help="name")
        p.add_argument("new_forcefield", type=str, help="new name")
        p.add_argument(
            "--list_forcefields",
            action=ListForcefieldsAction,
            help="show known forcefields",
            nargs=0,
        )
        p.add_argument(
            "--scale",
            type=float,
            default=1.0,
            help="scale interaction energies, default 1.0",
        )
        p.add_argument(
            "--offset",
            type=float,
            default=0.0,
            help="offset interaction energies, default 0.0",
        )

    @staticmethod
    def main(args):
        forcefield = read_forcefield(get_forcefield_file(args.forcefield))
        forcefield = scale_interaction(forcefield, args.scale, args.offset)
        write_forcefield(forcefield, args.new_forcefield)
