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
from six.moves import range

import numpy as np
import re
import six
from six import string_types
import yaml
import os
from os.path import isfile
import pandas as pd

from .scripts import _ScriptMeta
from .util import maybe_open


def parse_remc_conf(log):
    """parse replica configuration data

    Returns
    -------
    n_replicas, exchange_rate
    """
    with maybe_open(log) as fh:
        for line in fh:
            if "number_replica" in line:
                n_replicas = int(line.split()[-1])
            elif "exchange_rate" in line:
                exchange_rate = int(line.split()[-1])
                break
    return n_replicas, exchange_rate


def parse_attempt(att):
    """parse a single attempt"""
    sims, res = att.split("=")
    idx, idy = np.int64(sims.split("-"))
    return idx, idy, res == "x"


def update_sims(attemps, sims):
    """update sim order based on attempts"""
    sims = sims.copy()
    for att in attemps:
        idx, idy, suc = parse_attempt(att)
        if suc:
            sims[idx], sims[idy] = sims[idy], sims[idx]
    return sims


class ExchangeParser(object):
    """

    The only goal of this class is to have a simple parser
    that filters the exchange events from a log file

    The filtering is solved using a finite state machine.

    B -> S -> E -> P <- N <- (enter here)
    ^              |
    |--------------|

    B: Beginning
    S: Sweep
    E: Exchanges
    P: Probability
    N: None (entry point helper)

    Members
    -------
    from_log:
    sweeps:
    exchanges:
    """

    # enum trick works in python 2.7 and up (change to Enum when dropping 2.7)
    class _states(object):
        beginning = 1
        sweep = 2
        exchange = 3
        probability = 4
        none = 5

    def __init__(self):
        self._state = self._states.none
        # yapf: disable
        # Format: [regex, fromself._state, toself._state]
        self._transitions = [('\s*total_mc_sweep', self._states.none, self._states.probability),
                             (".*Replica_Attempt", self._states.probability, self._states.beginning),
                             ("\s+sweep", self._states.beginning, self._states.sweep),
                             ("\s+exchanges", self._states.sweep, self._states.exchange),
                             ("\s+probabilities", self._states.exchange, self._states.probability)]
        # yapf: enable
        self._transitions = [
            (re.compile(regex), f, t) for (regex, f, t) in self._transitions
        ]
        self._sweeps = []
        self._exchanges = []

    def advance(self, line):
        for regex, from_state, to_state in self._transitions:
            if regex.match(line):
                if not self._state == from_state:
                    raise RuntimeError(
                        "Illegal log file record expected line to match regex pattern: {}".format(
                            regex.pattern
                        )
                    )
                self._state = to_state
                # update diverse counters
                if self._state == self._states.sweep:
                    self._sweeps.append(int(line.split()[-1]))
                if self._state == self._states.exchange:
                    self._exchanges.append(line.replace(",", "").split(":")[-1].strip())

    @property
    def sweeps(self):
        return np.asarray(self._sweeps)

    @property
    def exchanges(self):
        return self._exchanges

    @classmethod
    def from_log(cls, log):
        ep = cls()
        with maybe_open(log) as fh:
            for i, line in enumerate(fh):
                try:
                    ep.advance(line)
                except RuntimeError as e:
                    raise type(e)(
                        "Error on line: {}\n Original Error Message: {}".format(
                            i, str(e)
                        )
                    )

        return ep


def log_path_and_dt(config_file):
    with open(config_file) as fh:
        config = yaml.load(fh)
    log = config["output"]["log"]
    folder = os.path.split(config_file)[0]
    log = os.path.join(folder, log)
    dt = config["montecarlo"]["algorithm-params"].get("dt", 1)
    return log, dt


def untangle(exchanges, nreplicas):
    """parse the exchanges into a dataframe
    """
    sims = np.arange(nreplicas)
    res = np.zeros((len(exchanges) + 1, nreplicas))
    res[0] = sims

    for i, ex in enumerate(exchanges):
        j = i + 1  # offset parameter
        attemps = ex.split()
        sims = update_sims(attemps, sims)
        res[j] = sims

    return pd.DataFrame(res, columns=["sim-{}".format(i) for i in range(nreplicas)])


def demux(config):
    """
    demux a REMC simulation

    Parameters
    ----------
    config : str
        filename of config to use

    Returns
    -------
    dataframe with demux information
    """
    log, dt = log_path_and_dt(config)
    nreplicas, exchange_rate = parse_remc_conf(log)
    exchanges = ExchangeParser.from_log(log)

    df = untangle(exchanges.exchanges, nreplicas)
    df["frame"] = (np.arange(len(df))) * exchange_rate
    df["time"] = df["frame"] * dt * 1e3
    return df


def structure_cont(df):
    """Change dataframe from thermodynamically continuous to structurally
    continuous.

    Parameters
    ----------
    df : pandas data frame
        Columns contain the ids of the thermodynamic configuration and rows the
        structure id. Also called thermodynamically continuous.

    Returns
    -------
    pandas data frame
        Columns contain the structure id and rows the id of the current
        thermodynamic configuration. Also called structurally continuous.

    """

    def update_name(n):
        if "sim" not in n:
            return n
        else:
            return "repl-" + n.split("-")[-1]

    # We have to do this for NOT sorting like: sim-1,sim-10,sim-11,sim-2...
    keys = df.columns.difference(["frame", "time"])
    sorted_keys = [
        "sim-{}".format(idx)
        for idx in np.sort([int(x.replace("sim-", "")) for x in keys])
    ]
    # obtain thermodynamically sorted array
    values = df[sorted_keys].values
    # argsort each row to obtain the configuration index for each replica
    df[sorted_keys] = np.argsort(values, axis=1)
    df.rename(columns={k: update_name(k) for k in df.columns}, inplace=True)
    return df


def write_xvg(df, file_handle):
    """
    Write dataframe to file_handle in xvg format
    """
    for i in range(len(df)):
        row = df.iloc[i]
        file_handle.write(
            "{:.2f} {}\n".format(
                row["time"], " ".join([str(el) for el in row.values[:-2].astype(int)])
            )
        )


class Demux(six.with_metaclass(_ScriptMeta)):
    description = """Demux replica exchange simulation to analyze replica jumps or generate input
    for 'gmx trjcat' to obtain continuous trajectories. By default output
    suitable for 'gmx trjcat' will be produced. In this mode the columns will
    contain the id of the replica in each folder/simulation at each exchange
    step, this is also called thermodynamically continuous. Using the switch
    `--structure_cont` the columns will contain in which folder/simulation a
    given replica is at each exchange step. The last option is optimal to
    analyze how a replica diffuses through, e.g. temperature space. A different
    name for that option is structurally continuous."""

    @staticmethod
    def parser(p):
        p.add_argument("config", type=str, help="config file")
        p.add_argument("out", type=str, help="xvg to write")
        p.add_argument(
            "--format",
            type=str,
            help="format of output either xvg(default) or csv",
            default="xvg",
        )
        p.add_argument(
            "--structure_cont",
            action="store_true",
            help="Sort by replica, (structure continuous)",
        )

    @staticmethod
    def main(args):
        if not isfile(args.config):
            raise IOError("File does not exist: {}".format(args.log))

        exchanges = demux(args.config)

        if args.structure_cont:
            exchanges = structure_cont(exchanges)

        if args.format == "csv":
            exchanges.to_csv(args.out, index=False)
        elif args.format == "xvg":
            with open(args.out, "w") as fh:
                write_xvg(exchanges, fh)
        else:
            print("Unknown format type please choose 'xvg' or 'csv'.")
