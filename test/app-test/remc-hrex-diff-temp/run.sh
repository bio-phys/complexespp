#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++
    $complexes --backup=false --config=27-remc.config --multidir ./test0/ ./test1/ ./test2/ ./test3/ --replex 3 --replex-accept remc --replex-verbosity=all > remc.log.out 2>&1
    $complexes --backup=false --config=27-hrex.config --multidir ./test0/ ./test1/ ./test2/ ./test3/ --replex 3 --replex-accept hrex > hrex.log.out 2>&1
}

function testit {
    python test.py
}

function clean {
    rm -f ./test*/*.xtc ./test*/*.pdb ./test*/*.log *.log.out *.data
}

run
testit
res=$?
clean
exit $res
