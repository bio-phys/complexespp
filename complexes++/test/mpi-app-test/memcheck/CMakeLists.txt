file(COPY "complex.cplx" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "setup-valgrind.yaml" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "run.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(mpi-app-memcheck run.sh)
