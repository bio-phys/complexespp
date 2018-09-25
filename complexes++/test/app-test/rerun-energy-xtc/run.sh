#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++
    $complexes --backup=true --config=conf.yaml
    cp run-energy.stat backup_run-energy.stat
    $complexes --config=conf2.yaml --rerun=True
}

function testit {
    python test.py
}

function clean {
    rm -f *data *pdb backup*
}

run
testit
res=$?
clean
exit $res
