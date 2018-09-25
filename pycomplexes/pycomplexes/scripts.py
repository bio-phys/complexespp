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
import argparse


_parser = argparse.ArgumentParser()
_subparsers = _parser.add_subparsers(help="Commands ")


# this is the meta class
# IF NEW SUBPARSER IS TO BE ADDED, INCLUDE IN __init__
class _ScriptMeta(type):
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        super(_ScriptMeta, cls).__init__(name, bases, classdict)

        subparser = _subparsers.add_parser(name.lower(), help=classdict["description"])

        cls.parser(subparser)
        subparser.set_defaults(func=cls.main)


def main():
    args = _parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
