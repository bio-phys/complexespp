#!/bin/bash

set -e

cwd=`pwd`

# test python library
cd python
python -m pytest

cd $cwd
