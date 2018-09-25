#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++
    $complexes --backup=false --config=pdb.yaml
    $complexes --backup=false --config=xtc.yaml
    $complexes --backup=false --config=trr.yaml
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
clean
exit $res
