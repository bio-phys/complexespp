#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../../src/complexes++
    for folder in Residue-*
    do
        cd $folder
        $complexes -c config.yaml --rerun=True --nb-thread=1 --backup=false
        cd ..
    done
}

function testit {
    python test.py
}


run
testit
res=$?
exit $res
