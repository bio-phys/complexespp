#!/bin/bash

set -e

cwd=`pwd`

# test python library
cd pycomplexes
python -m pytest

cd $cwd
