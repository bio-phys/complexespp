# -------------------------------------------------------------------------
# Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# The code comes without warranty of any kind
# Please refer to Kim and Hummer J.Mol.Biol. 2008
# -------------------------------------------------------------------------

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
