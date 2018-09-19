# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
from six import StringIO
import numpy as np
import pandas as pd

import pytest

from pycomplexes import demux


# yapf: disable
@pytest.mark.parametrize('att, sol', (('2-3=x', (2, 3, True)),
                                      ('1-2=', (1, 2, False)),
                                      ('4-5=o', (4, 5, False))))
def test_parse_attemp(att, sol):
    x, y, suc = demux.parse_attempt(att)
    assert x == sol[0]
    assert y == sol[1]
    assert suc == sol[2]
# yapf: enable


def test_parse_remc_conf():
    log = """[LOG] [Replica_Conf]
  number_replica: 6
  replica_directories: id = 0 : simulation-temp-298, id = 1 : simulation-temp-318, id = 2 : simulation-temp-338, id = 3 : simulation-temp-358, id = 4 : simulation-temp-378, id = 5 : simulation-temp-398
  config_file_regex: config
  exchange_rate: 5
  statistics_rate: 1000
  total_mc_sweep: 50"""
    log = StringIO(log)
    nrep, ex_rate = demux.parse_remc_conf(log)
    assert nrep == 6
    assert ex_rate == 5


@pytest.mark.parametrize(
    "attempts, sol",
    [
        ("0-1=x 1-2= 2-3=x 3-4= 4-5=x", [1, 0, 3, 2, 5, 4]),
        ("0-1= 1-2=x 2-3= 3-4=x 4-5=", [0, 2, 1, 4, 3, 5]),
        ("0-1=x 1-2= 2-3=o 3-4= 4-5=x", [1, 0, 2, 3, 5, 4]),
    ],
)
def test_update_sims(attempts, sol):
    sims = np.arange(6)
    sims = demux.update_sims(attempts.split(), sims)
    np.testing.assert_equal(sims, sol)


def test_ExchangeParser():
    log = """  exchange_rate: 5
  statistics_rate: 1000
  total_mc_sweep: 50
[LOG] [Replica_Attempt]
  sweep: 5
  delta_energy: -0.664201 (kT)
  attempted_exchange: 0-1
  exchanges: 0-1=x, 1-2=, 2-3=x, 3-4=, 4-5=x
  probabilities: 0-1=0.958818, 1-2=, 2-3=1, 3-4=, 4-5=1
[LOG] [Replica_Attempt]
  sweep: 10
  exchanges: 0-1=, 1-2=x, 2-3=, 3-4=x, 4-5=
  probabilities: 0-1=, 1-2=0.959594, 2-3=, 3-4=1, 4-5=
[LOG] [Replica_Attempt]
  sweep: 15
  delta_energy: -1.38796 (kT)
  attempted_exchange: 0-1
  exchanges: 0-1=x, 1-2=, 2-3=x, 3-4=, 4-5=x
  probabilities: 0-1=0.915872, 1-2=, 2-3=0.945745, 3-4=, 4-5=1
[LOG] [Replica_Attempt]
  sweep: 20
  exchanges: 0-1=, 1-2=x, 2-3=, 3-4=x, 4-5=
  probabilities: 0-1=, 1-2=1, 2-3=, 3-4=1, 4-5=
[LOG] [Replica_Attempt]
  sweep: 25
  delta_energy: 1.68188 (kT)
  attempted_exchange: 0-1
  exchanges: 0-1=x, 1-2=, 2-3=x, 3-4=, 4-5=x
  probabilities: 0-1=1, 1-2=, 2-3=1, 3-4=, 4-5=0.967101"""
    log = StringIO(log)
    ep = demux.ExchangeParser.from_log(log)
    np.testing.assert_equal(ep.sweeps, [5, 10, 15, 20, 25])
    # yapf: disable
    np.testing.assert_equal(ep.exchanges, ['0-1=x 1-2= 2-3=x 3-4= 4-5=x',
                                           '0-1= 1-2=x 2-3= 3-4=x 4-5=',
                                           '0-1=x 1-2= 2-3=x 3-4= 4-5=x',
                                           '0-1= 1-2=x 2-3= 3-4=x 4-5=',
                                           '0-1=x 1-2= 2-3=x 3-4= 4-5=x'])
    # yapf: enable


def test_ExchangeParser():
    exchanges = [
        "0-1=x 1-2= 2-3=x 3-4= 4-5=x",
        "0-1= 1-2=x 2-3= 3-4=x 4-5=",
        "0-1=x 1-2= 2-3=x 3-4= 4-5=x",
    ]
    nreplicas = 6
    df = demux.untangle(exchanges, nreplicas)
    # yapf: disable
    res = [[0., 1., 2., 3., 4., 5.],
           [1., 0., 3., 2., 5., 4.],
           [1., 3., 0., 5., 2., 4.],
           [3., 1., 5., 0., 4., 2.]]
    # yapf: enable
    np.testing.assert_equal(df.values, res)


def test_structure_cont():
    df = pd.read_csv(
        StringIO(
            """sim-12,sim-1,sim-10,sim-11,frame,time
0.0,1.0,2.0,3.0,0,0.0
1.0,0.0,3.0,2.0,3,3000.0
1.0,3.0,0.0,2.0,6,6000.0
"""
        )
    )
    cont = demux.structure_cont(df)

    ref = pd.read_csv(
        StringIO(
            """
repl-12,repl-1,repl-10,repl-11,frame,time
2,3,0,1,0,0.0
1,0,3,2,3,3000.0
0,1,3,2,6,6000.0
    """
        )
    )

    np.testing.assert_equal(cont.values, ref.values)
    assert np.all(cont.columns == ref.columns)
