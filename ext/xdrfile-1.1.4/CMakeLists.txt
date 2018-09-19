cmake_minimum_required (VERSION 2.8.12)

include(ExternalProject)

ExternalProject_Add(xdrfile
  PREFIX 1.1.4
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_COMMAND CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/configure --prefix=<INSTALL_DIR>
  BUILD_COMMAND ${MAKE}
  )

ExternalProject_Get_Property(xdrfile install_dir)
set(XDRFILE_INSTALL_DIR ${install_dir} CACHE INTERNAL "install")
set(XDRFILE_INCLUDE_DIR ${XDRFILE_INSTALL_DIR}/include CACHE INTERNAL "include")
set(XDRFILE_LIB_DIR ${XDRFILE_INSTALL_DIR}/lib/libxdrfile.a CACHE INTERNAL "lib")