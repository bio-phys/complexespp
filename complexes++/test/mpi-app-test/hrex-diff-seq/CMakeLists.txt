file(COPY "test0" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY "test1" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY "test2" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY "test3" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

file(COPY "test.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "test-log.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "run.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(mpi-app-hrex-diff-seq run.sh)
