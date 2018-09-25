# Try to find c++format
# https://fmt.github.io/
#
# fmt_FOUND - system has fmt
# fmt_INCLUDE_DIR - fmt include directory
# fmt_LIBRARY - link symbols for fmt

include(LibFindMacros)

libfind_pkg_check_modules(fmt_PKGCONF fmt)

find_path(fmt_INCLUDE_DIR fmt/format.h
  HINTS $ENV{FMT_ROOT} $ENV{FMT_HOME}
  PATH_SUFFIXES include include/fmt
  PATHS ${fmt_PKGCONF_INCLUDE_DIRS}
)

#set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
find_library(fmt_LIBRARY
  NAMES fmt
  HINTS $ENV{FMT_ROOT} $ENV{FMT_HOME}
  PATH_SUFFIXES lib64 lib
  PATHS ${fmt_PKGCONF_LIBRARY_DIRS}
)

set(fmt_PROCESS_INCLUDES fmt_INCLUDE_DIR)
set(fmt_PROCESS_LIBS fmt_LIBRARY)
libfind_process(fmt)
