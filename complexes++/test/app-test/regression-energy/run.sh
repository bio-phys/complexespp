#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++
    $complexes --backup=false --config=rerun.conf --rerun=True --nb-threads=1
}

function testit {
    python test.py
}

function clean {
    rm -f *data *pdb *xtc
}

run
testit
res=$?
#clean
exit $res
