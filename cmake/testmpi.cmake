###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
# This goes with the getCpuInfos.cpp
# This will create one CMAKE value per output option from the cpp file.
# For example the output of the CPP file can be:
# SSE3=TRUE;AVX=FALSE
# Then it will create:
# CPU_INFO_SSE3 = TRUE
# CPU_INFO_AVX = FALSE
#
# The binary should return 0 on success.
###########################################################################################
macro(TestMpi)
# The original CPP file
set(TestMpiFile "${PROJECT_SOURCE_DIR}/cmake/testmpi.cpp")

# Fatal error if the file does not exist
if(NOT EXISTS ${TestMpiFile})
	message(FATAL_ERROR "The TestMpiFile does not exist (${TestMpiFile})")
endif()


# Simply try and run
try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
      ${CMAKE_CURRENT_BINARY_DIR} ${TestMpiFile}
      COMPILE_OUTPUT_VARIABLE comp
      RUN_OUTPUT_VARIABLE run)

# If it has successfuly compiled an run
if(COMPILE_RESULT_VAR)
    if(RUN_RESULT_VAR EQUAL 0)
        message(STATUS "MPI is turned ON : test file can be compiled and executed")
    else()
        message(STATUS "MPI is turned ON : test file can be compiled, but execution was not possible")
    endif()
	set( MPI_WORKS TRUE )
else()
	message(STATUS "MPI is turned OFF : test file cannot be compiled, error is :")
	message(STATUS "${comp}")
	set( MPI_WORKS FALSE )
endif()

endmacro(TestMpi)
