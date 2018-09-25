# Try to find yaml-cpp
# https://code.google.com/p/yaml-cpp/
#
# yaml-cpp_FOUND - system has yaml-cpp
# yaml-cpp_INCLUDE_DIR - yaml-cpp include directory
# yaml-cpp_LIBRARY - link symbols for yaml-cpp

include(LibFindMacros)

libfind_pkg_check_modules(yaml-cpp_PKGCONF yaml-cpp)

find_path(yaml-cpp_INCLUDE_DIR yaml-cpp/yaml.h
  #NAMES yaml-cpp
  HINTS $ENV{YAML_CPP_ROOT} $ENV{YAML_CPP_HOME}
  PATH_SUFFIXES include
  PATHS ${yaml-cpp_PKGCONF_INCLUDE_DIRS}
)

find_library(yaml-cpp_LIBRARY
  NAMES yaml-cpp
  HINTS $ENV{YAML_CPP_ROOT} $ENV{YAML_CPP_HOME}
  PATH_SUFFIXES lib64 lib
  PATHS ${yaml-cpp_PKGCONF_LIBRARY_DIRS}
)

set(yaml-cpp_PROCESS_INCLUDES yaml-cpp_INCLUDE_DIR)
set(yaml-cpp_PROCESS_LIBS yaml-cpp_LIBRARY)
libfind_process(yaml-cpp)
