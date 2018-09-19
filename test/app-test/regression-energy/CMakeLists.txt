file(COPY "system.cplx" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "rerun.conf" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "structures.pdb" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "reference-energy.dat" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "run.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "test.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(app-energy-regression run.sh)
