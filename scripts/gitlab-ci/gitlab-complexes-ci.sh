#!/bin/bash

set -ex

function build_and_test_in_mode() {
    cd $cwd/scripts/build_scripts
    ./build_dynamic.sh -t $1

    # run complexes tests
    cd /tmp/complexes-pp/build/complexes-pp_dynamic
    if [[ $RUN_TEST == "true" ]]; then
        ctest -V
    fi
}

echo $USE_DEVTOOLS
if [[ $USE_DEVTOOLS == "true" ]]; then
    scl enable devtoolset-4 bash
fi

export NPROC=1
USER=root
cwd=`pwd`
build_and_test_in_mode $BUILD_TYPE
