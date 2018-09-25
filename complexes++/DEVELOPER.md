<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Out Of Source Builds](#out-of-source-builds)
- [Release/Debug builds](#releasedebug-builds)
- [Source Code Formatting](#source-code-formatting)
- [Compilers](#compilers)
- [Speed Up Recompilation](#speed-up-recompilation)
- [Compilation Flags](#compilation-flags)
- [Tests](#tests)
    - [Running](#running)
    - [Writing](#writing)
    - [Leak Checking](#leak-checking)

<!-- markdown-toc end -->

# Out Of Source Builds

For development builds we recommend out of source builds.

    mkdir build && cd build
    cmake ..
    make
    
# Release/Debug builds

By default cmake will generate a release build. This sets some optimization
flags and other stuff for use. But sometimes it can be useful to run Complexes++
in a debugger with debug symbols. For this you can specify the build type with

    cmake -DCMAKE_BUILD_TYPE=Release/Debug ..

During development it is a good idea to use debug builds. Using out of source
builds you can permanently set up a release and a debug build. Just create two
folders with the config-options you want and then build Complexes++ in the
respective folder.

# Source Code Formatting

To keep a consistent code style please use the following git precommit hook

https://github.com/andrewseidl/githook-clang-format

# Compilers

To change the compiler run the cmake step with

    CXX=<compiler> cmake ..

For development we recommend to use clang, since it provides much better error
messages. Before any code is merged it is tested with a variety of different
compilers.

# Speed Up Recompilation

CMake is configured to automatically make use of ccache if available. This can
help to speed up rebuilds after changing a branch to less then half a minute.

# Compilation Flags

Complexes++ has several 

1. TESTS (default on)
   Build testsuite

2. USE_LOGOUTPUT (default off)
   Enable `LOG_OUTPUT(...)` macros for more verbose logging.

3. USE_TIMINGOUTPUT (default off)
   Enable `TIMEZONE` macros to measure algorithm times.
   
4. USE_DEBUGINFO (default off)
  attach debug symbols in release builds.

# Tests

## Running 

We have setup several tests for check the individual classes of the program and
others to test Complexes++ as a whole. All tests can be run in the build
directory with.

    make test

## Writing

If you want to test classes and API interfaces of the Code you should use the
GoogleTest-Framework we have set up in the folder `test/unit-test`. Tests that
run a short Complexes++ simulation and verify the results should be placed in
individual directories in `test/app-test`. We also have some statistical tests.
They are not run by the CI server. However they can be run individually between
releases to ensure everything works as expected.

## Leak Checking

Valgrind is a powerful tool for performance analysis and memory leak checking.
To run the memory checking tool use:

    valgrind --tool=memcheck --leak-check=full complexes
    
Complexes++ uses OpenMP which can lead to many false positives. Have a look into
suppression files to remove any leaks reported about OpenMP.
