#!/bin/bash

function run {
    complexes=../../../src/complexes++_mpi
    
    run_as_root_docker=""
    if [[ $(mpirun --version 2>&1 | grep "Open MPI" | wc -l) == "1" ]]  && [[ "$USER" == "" ]] ; then
        run_as_root_docker="--allow-run-as-root"
    fi
    output=$(mpirun -np 1 $run_as_root_docker valgrind --tool=memcheck --leak-check=full  --show-leak-kinds=all --errors-for-leak-kinds=definite,indirect --error-exitcode=1 \
             $complexes --backup=False --config=setup-valgrind.yaml 2&>1 )

    echo $output
    if [[ $(echo $output | grep "operator new" | wc -l) != "0" ]] ; then
        echo "valgrind contains possible C++ leaks"
        return 1
    else
        echo "valgrind looks fine"
        return 0
    fi
}

function clean {
    rm -f *out
}

run
res=$?
clean
exit $res
