# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from __future__ import absolute_import, print_function

from string import Template
import os
import six
import yaml
import numpy as np

from .scripts import _ScriptMeta

DEFAULT_CONFIG = Template(
    """
structure: $cplx
montecarlo:
    algorithm: nvt
    algorithm-params:
        accept-func: metropolis
        temperatur: 300
    seed: $seed
    initial-position:
        strategy: none
    short-range-cutoff:
        enable: True
        container: sparse
        radius: 30
output:
    log: ${name}.log
    file: ${name}.trr
    freq: 1
    nstructures: $nstructures
    stat-file: ${name}.stat
    restart-freq: -1
"""
)


def generate_config(cplx, nstructures=100):
    name = os.path.basename(cplx)
    return yaml.load(
        DEFAULT_CONFIG.substitute(
            cplx=cplx,
            name=name,
            nstructures=nstructures,
            seed=np.random.randint(int(1e6)),
        )
    )


class Generate_Config(six.with_metaclass(_ScriptMeta)):
    description = "Generate a default configuration for complexes++"

    @staticmethod
    def parser(p):
        p.add_argument("config", type=str, help="name of output configuration file")
        p.add_argument("--cplx", default="structure.cplx", type=str, help="cplx file")
        p.add_argument(
            "--nstructures", type=int, help="number of frames in output", default=1000
        )

    @staticmethod
    def main(args):

        c = generate_config(args.cplx, args.nstructures)

        with open(args.config, "w") as f:
            yaml.dump(c, f, default_flow_style=False)
