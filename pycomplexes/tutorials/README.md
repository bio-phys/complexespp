# Complexes Tutorials

This directory contains several use-case examples for complexes++ and
pycomplexes. The examples are written as jupyter notebooks to interleave bash
commands and an explanation of what is being done.


# Install Jupyter 

To install the jupyter program we recommend the anaconda distribution. Please
note that some examples use pycomplexes as a python library. Therefore please
ensure that pycomplexes is installed in your conda environment.


# Setup environment for notebook

The notebooks assume that you have installed `complexes++` and `pycomplexes` and that
both commands can be found in the PATH. To check if either executable can be found
use `!which complexes` and `!which pycomplexes` in a notebook. If the output is empty
the program can not be found. 

The quickest way to get jupyter and pycomplexes running is to use the conda environment
we ship in `conda-linux-64.lock`, using `conda create --file conda-linux-64.lock -n complexes`.

Complexes++ has to be installed separately. Please refer to the complexes++ docs for details.
If you cannot find complexes++ check your PATH variable and ensure that the binary is in a 
directly of PATH.
