# pycomplexes 

This python package contains the **pycomplexes** toolbox. This toolbox can be
used to prepare input files for complexes++ from topology and structure files.
For a full list of commands run.

    pycomplexes --help

For help for an individual command run:

    pycomplexes <command> --help

## Installation

To install pycomplexes in your main python environment you can use.

    pip install .

If you get warnings about file-permissions you can install pycomplexes in your
home directory by appending the `--user` flag to the last command. Please keep
in mind that adding this flag means you also have to update your `PATH` variable
to include `~/.local/bin` to ensure the pycomplexes executable is found.

    export PATH=$HOME/.local/bin:$PATH

## Development Notes

The unit test are based on the py.test framework. To run them call:

    python -m pytest -vv

This ensures that are running against the current version in the checkout. To
add coverage information as well add `--cov=pycomplexes`.

## Help

To find help about the tools included in pycomplexes run

    pycomplexes --help

For help for an individual tool run

    pycomplexes <command> --help
