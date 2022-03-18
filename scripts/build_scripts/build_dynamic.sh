#!/bin/bash

#
# Script to build of Complexes++ in dynamic mode (using shared libs, if available)
# 2016  Klaus Reuter, khr@rzg.mpg.de
#

# quit on any error
set -e

# --- option parsing

help() { echo "Usage: $0 [-f] [-t] " 1>&2;
         echo "Use -f to force build of deps"
         echo "Use -t to build tests"
         exit 1; }

force="FALSE"
test="OFF"
while getopts ":fht" o; do
    case "${o}" in
        f)
            force="TRUE"
            ;;
        h)
            help
            ;;
        t)
            test="ON"
            ;;
        *)
            ;;
    esac
done
shift $((OPTIND-1))

# --- include configuration file
source build.cfg

# --- check for the very basic prerequisites
type cmake
type $CC
type $CXX


BASE=`pwd`/../../complexes++
# readlink does not work on MAC so we turn towards the following Python hack
#BASE=`readlink -f $BASE`
BASE=`python -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' $BASE`
BUILD=$BUILD/complexes-pp_dynamic

# --- remove old build artifacts if we force building

if [ "${force}" == "TRUE" ]; then
    rm -rf $PREFIX/bin
    rm -rf $BUILD
fi

mkdir -p $PREFIX
mkdir -p $BUILD

# --- system-wide installations (--> loaded environment modules at MPCDF) override local dependency installation
if [ x"$FMT_ROOT" == x"" ]; then
  export FMT_ROOT=$PREFIX/dep
fi
if [ x"$YAML_CPP_ROOT" == x"" ]; then
  export YAML_CPP_ROOT=$PREFIX/dep
fi
if [ x"$PKG_CONFIG_PATH" == x"" ]; then
  export PKG_CONFIG_PATH=$PREFIX/dep/lib/pkgconfig
fi

BUILDTYPE=Release
# check that there was a extra option and if it spells debug when all lower case
if [[ -n $1 && $(echo $1 | tr '[:upper:]' '[:lower:]') == "debug" ]]
then
    BUILDTYPE=Debug
fi


cd $BUILD
cmake -DCMAKE_BUILD_TYPE=${BUILDTYPE} \
      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DTESTS=${test} \
      $BASE
echo "make -j $NPROC"
VERBOSE=1 make -j $NPROC
make install
