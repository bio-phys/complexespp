set(SOURCES
  "domains/abstractdomain.cpp"
  "domains/beads.cpp"
  "domains/connection.cpp"
  "domains/rigiddomain.cpp"
  "energy/energy.cpp"
  "io/cplx.cpp"
  "io/ffparameters.cpp"
  "io/gro.cpp"
  "io/io.cpp"
  "io/pdb.cpp"
  "io/reader.cpp"
  "io/restart.cpp"
  "io/trajectoryfile.cpp"
  "io/xdr.cpp"
  "mc/accept_func.cpp"
  "mc/mc.cpp"
  "mc/npt.cpp"
  "mc/nvt.cpp"
  "mc/simulation.cpp"
  "pairkernels/pairkernelbuilder.cpp"
  "parallelization/taskslimiter.cpp"
  "setup/application.cpp"
  "setup/config.cpp"
  "util/array.cpp"
  "util/file.cpp"
  "util/linalg.cpp"
  "util/log.cpp"
  "util/random.cpp"
  "util/string.cpp"
  "util/util.cpp"
  "util/demangle.cpp"
  )

set(LIBRARIES
  ${XDRFILE_LIB_DIR}
  ${yaml-cpp_LIBRARY}
  ${PROJECT_BINARY_DIR}/prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}yaml-cpp${CMAKE_STATIC_LIBRARY_SUFFIX}
  fmt::fmt-header-only
  ${complexes++_LIBRARIES}
  CACHE STRING "link time dependencies")

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.cpp.in
  ${CMAKE_BINARY_DIR}/src/version.cpp
)
list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/version.cpp" version.h)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/complexesconfig.h.in
  ${CMAKE_BINARY_DIR}/src/complexesconfig.h
)
list(APPEND SOURCES complexesconfig.h)

include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_CURRENT_LIST_DIR})

# build files as a library so that I can reuse the object files
add_library(${PROJECT_NAME_STR}_objects OBJECT ${SOURCES})
add_dependencies(${PROJECT_NAME_STR}_objects xdrfile yaml_cpp_project)

add_executable (${PROJECT_NAME_STR} $<TARGET_OBJECTS:${PROJECT_NAME_STR}_objects> "main.cpp")
set_target_properties(${PROJECT_NAME_STR} PROPERTIES
	IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}yaml-cpp${CMAKE_STATIC_LIBRARY_SUFFIX}
)

target_link_libraries(${PROJECT_NAME_STR} PRIVATE ${LIBRARIES})
  
install (TARGETS ${PROJECT_NAME_STR}
  RUNTIME DESTINATION bin)

if(COMPLEXES_USE_MPI)
    add_executable (${PROJECT_NAME_STR}_mpi $<TARGET_OBJECTS:${PROJECT_NAME_STR}_objects> "mainmpi.cpp")
    target_link_libraries(${PROJECT_NAME_STR}_mpi ${LIBRARIES})
    install (TARGETS ${PROJECT_NAME_STR}_mpi RUNTIME DESTINATION bin)
endif()

