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

from setuptools import setup, find_packages

setup(name='pycomplexes',
      version='0.1dev',
      description='python helper library for complexes++',
      author='Max Linke',
      packages=find_packages(),
      install_requires=['numpy', 'pyyaml', 'six', 'pandas', 'tqdm', 'numba',
                        'MDAnalysis>=0.17.0', 'backports.functools_lru_cache'],
      entry_points={'console_scripts':
                    ['pycomplexes=pycomplexes.scripts:main']},
      package_data={'pycomplexes': ['forcefields/*']}
      )
