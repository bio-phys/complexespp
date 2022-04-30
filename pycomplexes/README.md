# pycomplexes 

This python package contains the **pycomplexes** toolbox. This toolbox can be
used to prepare input files for complexes++ from topology and structure files.
For a full list of commands run.

    pycomplexes --help

For help for an individual command run:

    pycomplexes <command> --help

## Installation

To install and run pycomplexes we recommend to use conda environments. To create 
a environment called pycomplexes and install it run:

    conda create -n pycomplexes -f conda-linux-64.lock
    conda activate pycomplexes
    python setup.py install

To make full use of pycomplexes we assume you installed [complexes++](../complexes++).

## Tutorials

pycomplexes includes a few jupyte notebooks to demonstrate how to use pycomplexes and complexes++ together.
You can find those in the tutorials folder.

## Help

To find help about the tools included in pycomplexes run

    pycomplexes --help

For help for an individual tool run

    pycomplexes <command> --help

## Development Notes

The unit test are based on the py.test framework. To run them call:

    python -m pytest -vv

This ensures that are running against the current version in the checkout. To
add coverage information as well add `--cov=pycomplexes`.

If you want to add new dependencies please add them to the pyproject.toml 
file and create a new lock file using [conda-lock](https://github.com/conda-incubator/conda-lock).
