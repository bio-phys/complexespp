#!/bin/bash

# stop on any error
set -e

function run {
    cp -r ./test ./test0
    cp -r ./test ./test1
    cp -r ./test ./test2
    cp -r ./test ./test3
    complexes=../../../src/complexes++
    $complexes --backup=false --config=27.config --multidir ./test0/ ./test1/ ./test2/ ./test3/
}

function testit {
    python test.py
}

function clean {
    rm -f ./test*/*.xtc ./test*/*.pdb ./test*/*.log
    rm -r ./test0 ./test1 ./test2 ./test3
}

run
testit
res=$?
clean
exit $res
