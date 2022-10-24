[![complexespp CI](https://github.com/bio-phys/complexespp/actions/workflows/complexes-ci.yml/badge.svg)](https://github.com/bio-phys/complexespp/actions/workflows/complexes-ci.yml)

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Installation](#installation)
- [Dependencies](#dependencies)
- [Advanced installation](#advanced-installation)
- [Help](#help)

<!-- markdown-toc end -->


complexes++ is a Monte-Carlo engine for coarse grained protein models. For an
explanation of the implemented algorithms please have a look into the [docs](../docs/).


You can find examples of using complexes++ in the pycomplexes included [tutorials](../pycomplexes/tutorials).

# Installation

To install all build requirements you can run the following commands on ubuntu/debian. 

	git submodule update --init --recursive # to get inastemp

To build complexes use

	mkdir build && cd build
	cmake ..
	make -j
	sudo make install

The default installation path is `/usr`. To change the installation directory to
`/your/path` change the second line to

    cmake -DCMAKE_INSTALL_PREFIX:PATH=/your/path ..

This will install complexes in `/your/path/bin`.

# Help

To find help for the complexes++ program options use

    complexes++ --help

# Development

Complexes++ is written in C++. Do develop

## run tests

We use a mixture of unit and integration tests. You can find them all in the tests folder. 
Assuming you already build complexes you can prepare the test environment and run tests 
with the following commands.

	conda create -n complexestest -f conda-linux-64.lock
	conda activate complexestest
	cd build
	make tests


This will create a conda environment with the python packages we need to our integration tests.
If you want to add new run time dependencies for the tests please add them to the pyproject.toml 
file and create a new lock file using [conda-lock](https://github.com/conda-incubator/conda-lock).
