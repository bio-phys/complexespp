#!/bin/bash

set -e

cwd=$(pwd)

# install python library
cd python
python setup.py install

# build complexes
cd $cwd/scripts/build_scripts
./build_dynamic.sh

# add complexes-pp and complexes to path
export PATH=$HOME/.local/bin:/root/complexes-pp/bin:$PATH

cd $cwd/tutorials
for nb in `find . -name '*ipynb' -not -name '*checkpoint*'`
do
echo $nb
jupyter nbconvert --execute --ExecutePreprocessor.timeout=1000 --ExecutePreprocessor.kernel_name=python3 $nb
done

