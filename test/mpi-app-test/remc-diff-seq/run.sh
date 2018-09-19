#!/bin/bash

# stop on any error
set -e

function run {
    run_as_root_docker=""
    if [[ $(mpirun --version 2>&1 | grep "Open MPI" | wc -l) == "1" ]] && [[ "$USER" == "" ]] ; then
        run_as_root_docker="--allow-run-as-root"
    fi
    mpirun -np 3 $run_as_root_docker ../../../src/complexes++_mpi --backup=false --config=27-remc.config --multidir ./test0/ ./test1/ ./test2/ ./test3/ --replex 3 --replex-accept remc > remc.log.out 2>&1
    cat remc.log.out
    ../../../src/complexes++ --backup=false --config=27-seq.config --multidir ./test0/ ./test1/ ./test2/ ./test3/ --replex 3 --replex-accept remc > seq.log.out 2>&1
    cat seq.log.out
}

function testit {
    python test.py
}

function clean {
    rm -f ./test*/*.xtc ./test*/*.pdb ./test*/*.log *.log.out
}

run
testit
res=$?
clean
exit $res
