#!/bin/bash

function run {
    complexes=../../../src/complexes++
    valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all \
             --suppressions=libomp.supp --gen-suppressions=all \
             --errors-for-leak-kinds=definite,indirect --error-exitcode=1 \
             $complexes --backup=False --config=setup-valgrind.yaml
}

function clean {
    rm -f *out
}

run
res=$?
clean
exit $res
