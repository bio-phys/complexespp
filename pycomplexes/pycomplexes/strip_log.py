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
from io import StringIO
import six
import yaml

from .util import maybe_open
from .scripts import _ScriptMeta


def process_n_hits(n_hits, log_fname, keyword="PROGRAM"):
    flag = False
    if n_hits == 0:  # This should never happen
        msg = "Could not find keyword \"{}\" in file {}".format(keyword, log_fname)
        msg = "\n".join([msg, "Not a complexes log file?"])
        print(msg)
        raise RuntimeError(msg)
    elif n_hits > 1:  # This should never happen
        msg = ["Warning!\nFound keyword \"{}\" {} times in file {}".format(keyword, n_hits, log_fname),
               "Apparently you are processing an appended log file...",
               "Will extract last run within this log file. This will increase the runtime of demux. ",
               "You should consider clearing this log file by means of the `pycomplexes strip_log` tool..."]
        print("\n".join(msg))
        # raise RuntimeError(msg)
        flag = True
    return flag


def skim_log_file(log_fname, keyword="PROGRAM"):
    """skim through log file to check for multiple simulation starts

    Returns
    -------
    multiple_hits_flag
    lidxs
    """
    # pattern = re.compile(keyword) regex is actually slower
    with maybe_open(log_fname) as fh:
        lidxs = [lidx for (lidx, line) in enumerate(fh) if keyword in line]
        # lidxs = [lidx for (lidx, line) in enumerate(fh) if pattern.search(line)]  # actually slower
    n_hits = len(lidxs)
    multiple_hits_flag = process_n_hits(n_hits, log_fname, keyword=keyword)
    return multiple_hits_flag, lidxs


def strip_log(log, line_ndx):
    stripped = StringIO()
    with open(log) as f:
        [stripped.write(line) for line in f.readlines()[line_ndx:]]
    return stripped


# generate a parser for this class
class Strip_Log(six.with_metaclass(_ScriptMeta)):
    description = "Strip last simulation run from log file."

    @staticmethod
    def parser(p):
        p.add_argument("log", type=str, help="log file of simulation")
        p.add_argument("out", type=str, help="name of stripped log")
        p.add_argument("bak", type=str, help="name of backed up log")

    @staticmethod
    def main(args):
        util.check_file_exists(args.log)

