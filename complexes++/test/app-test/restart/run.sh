#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++
    $complexes --backup=false --config=conf.yaml
    $complexes --restart=restart.nc --config=conf_restart.yaml
}

function testit {
    python test.py
}

function clean {
    rm -f *data *pdb *nc *trr
}

clean
run
testit
res=$?
clean
exit $res
