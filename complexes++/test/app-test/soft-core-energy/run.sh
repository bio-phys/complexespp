#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++
    $complexes -c test.config --rerun=True --nb-threads=1 --backup=false
}

function testit {
    python test.py
}


run
testit
res=$?
exit $res
