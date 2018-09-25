file(COPY "test" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

file(COPY "test.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "run.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(mpi-app-multidir run.sh)
