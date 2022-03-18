#!/bin/bash

#
# Script to build the dependencies of Complexes+ as static libraries
# 2016  Klaus Reuter, khr@rzg.mpg.de
#

if [ x"$USER" == x"gitlab-runner" ]; then
  echo "GITLAB-CI of the Complexes++ dependencies is currently unsupported."
  exit 1
fi

# quit on any error
set -e

# --- option parsing

help() { echo "Usage: $0 [-f] " 1>&2;
         echo "Use -f to force build of deps"
         exit 1; }

force="FALSE"
while getopts ":fh" o; do
    case "${o}" in
        f)
            force="TRUE"
            ;;
        h)
            help
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

BASE=$(pwd)
PREFIX=${PREFIX}/dep

# --- remove old build artifacts if we force building

if [ "${force}" == "TRUE" ]; then
    rm -rf $TARGZ
    rm -rf $PREFIX
    rm -rf $BUILD
fi

mkdir -p $PREFIX
mkdir -p $BUILD
mkdir -p $TARGZ

# determine the number of CPU cores used for the builds
if hash nproc &>/dev/null; then
  NPROC=$(nproc)
else
  NPROC=1
fi

# --- download, if necessary
function download {
    base_url=$2
    name=$3
    base_dir=$1

    dest=${base_dir}/${name}
    # there are cases where I want to rename the file otherwise assume filename is part of url
    if [ ! -z "$4" -a "$4" == "rename" ]; then
        url=${base_url}
    else
        url=${base_url}/${name}
    fi

    if [ ! -e ${dest} ]; then
        if hash curl &>/dev/null; then
            # use -L to follow redirects
            curl -o ${dest} -L ${url}
        else
            wget -O ${dest} ${url}
        fi
    fi
}

download ${TARGZ} https://github.com/jbeder/yaml-cpp/archive/ ${YAML_CPP}.tar.gz

# --- Get submodules

git submodule update --init --recursive

# --- Check MD5 Sums

cd $TARGZ
if hash sha512sum &>/dev/null; then
  sha512sum --quiet --check ${BASE}/deps.md5
fi
cd $BASE

export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig
mkdir -p $PKG_CONFIG_PATH


# --- build yaml-cpp
YAML_CPP_OK=$PREFIX/.yaml_cpp_ok
if [ ! -e $YAML_CPP_OK ]; then
  cd $BUILD
  if [ ! -d $YAML_CPP ]; then
    mkdir ${YAML_CPP}
    tar -xzf ${TARGZ}/${YAML_CPP}.tar.gz -C ${YAML_CPP}
  fi
  cd $YAML_CPP/yaml-cpp-${YAML_CPP}
  mkdir -p build && cd build
  cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DBUILD_SHARED_LIBS=OFF ..
  make -j $NPROC
  make install
  # cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
  #       -DBUILD_SHARED_LIBS=ON ..
  # make -j `nproc`
  # make install
  touch $YAML_CPP_OK
fi
cd $BASE

echo "OK!"
