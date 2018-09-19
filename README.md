<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Complexes](#complexes)
- [complexes++](#complexes)
    - [Installation](#installation)
        - [Dependencies](#dependencies)
        - [Advanced installation](#advanced-installation)
    - [Help](#help)
- [pycomplexes](#pycomplexes)
    - [Installation](#installation-1)
    - [Help](#help-1)

<!-- markdown-toc end -->
# Complexes 

Complexes consists of two separate programs. *complexes++* a Monte-Carlo engine
and a helper toolbox pycomplexes.

To find detailed description of theory and implemented Monte-Carlo algorithms
please have a look into the manual. Tutorials how to use complexes++ and
pycomplexes can be found in the tutorials folder.

# complexes++

complexes++ is a Monte-Carlo engine for coarse grained protein models. For an
explanation of the implemented algorithms please have a look into the manual.
Usage examples are given in the tutorials directory in the form of jupyter
notebooks.

## Installation

In case you cannot install the [dependencies](#dependencies) of complexes yourself
we also provide scripts that download and build the dependencies for you. Those
require that you also build complexes with our scripts.

    cd scripts/build_scripts
    ./build_deps.sh
    ./build_dynamic.sh

Using these scripts complexes++ will be installed at `~/complexes-pp/bin`. You
can change the installation directory in `build.cfg`. To find the executable in
bash by name you have to update your `PATH` variable.

    export PATH=$HOME/complexes-pp/bin:$PATH

To update to a new version of complexes with our scripts run it can happen that the
dependencies have to be rebuild. To rebuild complexes++ with all dependencies run
the following:

    cd scripts/build_scripts
    ./build_deps.sh -f
    ./build_dynamic.sh -f


### Dependencies

Complexes++ has the following external dependencies.

- boost >= 1.55
- fmt
- yaml-cpp

### Advanced installation

Assuming all [build dependencies](#dependencies) are installed on your system you
can build and install complexes++ using cmake directly

    mkdir build && cd build
    cmake ..
    make install

The default installation path is `/usr`. To change the installation directory to
`/your/path` change the second line to

    cmake -DCMAKE_INSTALL_PREFIX:PATH=/your/path ..

This will install complexes in `/your/path/bin`.

## Help

To find help for the complexes++ program options use

    complexes++ --help


# pycomplexes

Pycomplexes is a python toolbox to help setup and visualize simulations. It can
be used as a standalone application and as a python library for scripting.

## Installation

To install pycomplexes in your main python environment you can use.

    cd python
    pip install .

If you get warnings about file-permissions you can install pycomplexes in your
home directory by appending the `--user` flag to the last command. Please keep
in mind that adding this flag means you also have to update your `PATH` variable
to include `~/.local/bin` to ensure the pycomplexes executable is found.

    export PATH=$HOME/.local/bin:$PATH

## Help

To find help about the tools included in pycomplexes run

    pycomplexes --help

For help for an individual tool run

    pycomplexes <command> --help
