#!/bin/bash

# stop on any error
set -e

function run {
    complexes=../../../src/complexes++_mpi
    run_as_root_docker=""
    if [[ $(mpirun --version 2>&1 | grep "Open MPI" | wc -l) == "1" ]]  && [[ "$USER" == "" ]] ; then
        run_as_root_docker="--allow-run-as-root"
    fi
    mpirun -np 1 $run_as_root_docker $complexes --backup=false --config=conf.yaml
    mpirun -np 1 $run_as_root_docker $complexes --restart=restart.nc --config=conf_restart.yaml
}

function testit {
    python test.py
}

function clean {
    rm -f *data *pdb *nc *.trr
}

clean
run
testit
res=$?
clean
exit $res
