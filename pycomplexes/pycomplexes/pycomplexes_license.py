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

import six

from .scripts import _ScriptMeta

class License(six.with_metaclass(_ScriptMeta)):
    description = "Print License Information"

    @staticmethod
    def parser(p):
        pass

    @staticmethod
    def main(args):
        print(" Copyright (c) 2018 the pycomplexes development team and contributors")
        print(" (see the file AUTHORS for the full list of names)")
        print("")
        print(" pycomplexes is free software: you can redistribute it and/or modify")
        print(" it under the terms of the GNU General Public License as published by")
        print(" the Free Software Foundation, either version 3 of the License, or")
        print(" (at your option) any later version.")
        print("")
        print(" pycomplexes is distributed in the hope that it will be useful,")
        print(" but WITHOUT ANY WARRANTY; without even the implied warranty of")
        print(" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the")
        print(" GNU General Public License for more details.")
        print("")
        print(" You should have received a copy of the GNU General Public License")
        print(" along with pycomplexes.  If not, see <https://www.gnu.org/licenses/>")
