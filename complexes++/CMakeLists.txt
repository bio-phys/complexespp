cmake_minimum_required (VERSION 3.18)
set(PROJECT_NAME_STR complexes++)
project(${PROJECT_NAME_STR} CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# options
option(TESTS "Enable tests." ON)
option(USE_TIMINGOUTPUT "Enable output time of important sections." OFF)
option(USE_DEBUGINFO "Attach source code to binary to use debuger/profiler even in release." OFF)

# check if the compiler supports c++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# test if a mpi file can be compiled
OPTION( COMPLEXES_USE_MPI "Set to OFF to disable MPI" OFF )
include(testmpi)
TestMpi()
if( MPI_WORKS )
  set(COMPLEXES_USE_MPI ON)
  message("-- Build mpi executable. Compiler warnings are turned off.")
endif()

# https://stackoverflow.com/a/24305849/2207958
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# look for external dependencies that we require
find_package(yaml-cpp REQUIRED)
find_package(OpenMP)

# --- set up compiler flags, depending on the compiler (g++, icpc) and the build type ---
# To set the compiler, set the environment variable `CXX`
# to the compiler's executable name before invoking CMAKE.

# debug info
if( USE_DEBUGINFO )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# -- Warnings --
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith -Wcast-qual -Wconversion")

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion -Wold-style-cast")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report -diag-disable 1682")
else() # GCC
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast")
endif()

# --- OpenMP ---
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "Openmp is disabled")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# --- Disable Error only for non MPI build (needed on Fedora 25) ---
if ( NOT COMPLEXES_USE_MPI )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

# --- compiler specific build-dependent flags ---
if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math") # Seems to have a bug: -funroll-loops
  if(CMAKE_EXE_LINKER_FLAGS MATCHES "-static")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2")
  else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  endif()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ip")
  if(CMAKE_EXE_LINKER_FLAGS MATCHES "-static")
    # generate a binary for distribution that runs on virtually any x86_64 Linux platform
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2 -axAVX,SSE4.2 -qopenmp-link-static")
  else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -xHost")
  endif()
endif()



# add external libraries that live in the repository
add_subdirectory(ext/xdrfile-1.1.4)
add_subdirectory(ext/fmt-8.1.1 EXCLUDE_FROM_ALL)

include_directories(SYSTEM
  # FMT should not need this according to docs
  ext/fmt-8.1.1/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${XDRFILE_INCLUDE_DIR}
  ${yaml-cpp_INCLUDE_DIR}
  ext/CLI11
  )

# inastemp
set(INASTEMP_JUST_LIB TRUE)
add_subdirectory(ext/inastemp)
INCLUDE_DIRECTORIES(
         ${INASTEMP_BINARY_DIR}/Src
         ${INASTEMP_SOURCE_DIR}/Src
         ${INASTEMP_INCLUDE_DIR}
    )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INASTEMP_CXX_FLAGS}")

# set rpath, following the hints at https://cmake.org/Wiki/CMake_RPATH_handling
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# add actual source code CMakeLists.txt
subdirs(src)

if (TESTS)
  add_subdirectory(ext/gtest-1.7.0)
  # there is a warning in gtest out of our controll.
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set_target_properties(gtest PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error -diag-disable 2304 -diag-disable 1682")
    set_target_properties(gtest_main PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error -diag-disable 2304 -diag-disable 1682")
  else()
    set_target_properties(gtest PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")
  endif()
  enable_testing()
  subdirs(test)
endif()

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    message( STATUS "Doxygen has been found, use \"make doc\" to generate the documentation." )
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../manual/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(
        doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
else()
    message( STATUS "Doxygen has not been found. The creation of the doc is not possible." )
endif(DOXYGEN_FOUND)
