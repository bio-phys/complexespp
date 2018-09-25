# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------
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
