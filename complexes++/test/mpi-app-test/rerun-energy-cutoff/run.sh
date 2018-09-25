#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++_mpi
    run_as_root_docker=""
    if [[ $(mpirun --version 2>&1 | grep "Open MPI" | wc -l) == "1" ]]  && [[ "$USER" == "" ]] ; then
        run_as_root_docker="--allow-run-as-root"
    fi
    mpirun -np 1 $run_as_root_docker $complexes --backup=true --config=conf.yaml
    cp run-energy.stat backup_run-energy.stat
    mpirun -np 1 $run_as_root_docker $complexes --config=conf2.yaml --rerun=True
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
