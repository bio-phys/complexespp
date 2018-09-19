#!/bin/bash

#
# Script to build of Complexes++ in pure static mode (standalone executable).
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


BASE=`pwd`/../..
# readlink does not work on MAC so we turn towards the following Python hack
#BASE=`readlink -f $BASE`
BASE=`python -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' $BASE`
BUILD=$BUILD/complexes-pp_static

# --- remove old build artifacts if we force building

if [ "${force}" == "TRUE" ]; then
    rm -rf $PREFIX/bin
    rm -rf $BUILD
fi

mkdir -p $PREFIX
mkdir -p $BUILD

export PKG_CONFIG_PATH=$PREFIX/dep/lib/pkgconfig
export BOOST_ROOT=$PREFIX/dep
export FMT_ROOT=$PREFIX/dep
export YAML_CPP_ROOT=$PREFIX/dep

cd $BUILD
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXE_LINKER_FLAGS="-static" \
      -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DBUILD_SHARED_LIBS=OFF \
      -DTESTS=${test} \
      $BASE
echo "make -j $NPROC"
VERBOSE=1 make -j $NPROC
make install
