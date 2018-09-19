# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from __future__ import absolute_import, print_function
from six.moves import range, zip

import numpy as np
import yaml
import six
from os.path import isfile, splitext
import warnings

from .scripts import _ScriptMeta


# -1 corresponds to an acidic sidechain
# +1 corresponds to an basic sidechain
# from propka.cfg
CHARGED_SIDECHAINS = {
    "ASP": {"charge": -1, "pK": 3.8},
    "HIS": {"charge": 1, "pK": 6.5},
    "LYS": {"charge": 1, "pK": 10.5},
    "GLU": {"charge": -1, "pK": 4.5},
    "ARG": {"charge": 1, "pK": 12.5},
    "TYR": {"charge": -1, "pK": 10.0},
}


def degree_of_dissociation(pH, pK):
    """From Henderson-Hasselbalch Eqn.

    Parameters
    ----------
    pH : float
        pH value
    pK : float
        pK of titrable group

    Returns
    -------
    degr_of_diss : float

    """
    return 1. / (10 ** (pK - pH) + 1)


def net_charge(pH, pK, charge=-1):
    """Returns net_charge of a titrable group at a certain pH.
    Rounded to 2 decimals.

    Parameters
    ----------
    pH : float
        pH value
    pK : float
        pK of titrable group
    charge : either 1 or -1
        specifies whether it is acidic or basic

    Returns
    -------
    net_charge : float

    """
    if charge == -1:
        net_charge = -degree_of_dissociation(pH, pK)
    elif charge == 1:
        net_charge = -degree_of_dissociation(pH, pK) + 1
    else:
        warnings.warn(
            "only -1, or 1 are allowed for charge. "
            "This variable specifies whether the titrable group is acidic (-1) or basic (1)."
        )
        net_charge = 0
    return float(np.round(net_charge, 2))


def change_charges_in_domain(domain, ph, charged_sidechains):
    """
    update charges in domains according to Henderson-Hasselbalch
    """
    for i, aa in enumerate(domain["beads"]):
        if aa in charged_sidechains:
            domain["charges"][i] = net_charge(ph, **charged_sidechains[aa])


def change_charges(cplx, ph=7., charged_sidechains=CHARGED_SIDECHAINS):
    """
    update charges in cplx according to Henderson-Hasselbalch
    """
    for top in cplx["topologies"]:
        for domain in six.itervalues(top["domains"]):
            change_charges_in_domain(domain, ph, charged_sidechains)
    return cplx


# generate a parser for this class
class PH(six.with_metaclass(_ScriptMeta)):
    description = (
        "Change charges in cplx according to external pH using Henderson-Hasselbalch."
    )

    @staticmethod
    def parser(p):
        p.add_argument("cplx", type=str, help="cplx file to read")
        p.add_argument("ph", type=float, help="enternal pH")
        p.add_argument(
            "-pk",
            "--pKs",
            type=str,
            default=None,
            help="Give a file with custom pK values. Expect dict with {<resname>: {'charge': <>, 'pK': <>}} in the yaml format",
        )
        p.add_argument(
            "-o",
            "--output",
            type=str,
            default=None,
            help="specify name of output cplx. Default append chosen pH to name",
        )

    @staticmethod
    def main(args):
        if not isfile(args.cplx):
            raise IOError("File does not exist: {}".format(args.cplx))

        with open(args.cplx) as f:
            cplx = yaml.load(f)
        if args.pKs:
            if not isfile(args.pKs):
                raise IOError("File does not exist: {}".format(args.pKs))
            with open(args.pKs, "r") as f:
                charged_sidechains = yaml.load(f)
            if not isinstance(charged_sidechains, dict):
                raise RuntimeError(
                    "Cannot parse content of file: >{}< as yaml-dictionary.".format(
                        args.pKs
                    )
                )
        else:
            charged_sidechains = CHARGED_SIDECHAINS

        cplx = change_charges(cplx, args.ph, charged_sidechains=charged_sidechains)
        out_fname = args.output or "{}_pH{:2.2f}.cplx".format(
            splitext(args.cplx)[0], args.ph
        )
        with open(out_fname, "w") as f:
            yaml.dump(cplx, f)
