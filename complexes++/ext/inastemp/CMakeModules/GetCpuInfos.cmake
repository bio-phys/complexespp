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
macro(GetCpuInfos)
# The original CPP file
set(GetCpuInfosFile "${PROJECT_SOURCE_DIR}/CMakeModules/getCpuInfos.cpp")

# Fatal error if the file does not exist
if(NOT EXISTS ${GetCpuInfosFile})
	message(FATAL_ERROR "The GetCpuInfosFile does not exist (${GetCpuInfosFile})")
endif()

OPTION( INASTEMP_ISDE_CPU  "Set to ON to run the CPU detection over sde64" OFF )

# Compile and execute the file
if(INASTEMP_ISDE_CPU)
    SET( INASTEMP_ISDE_CPU_ARGS "-knl" CACHE STRING "Arguments for sde64"  )

    if($ENV{VERBOSE})
    	message(STATUS "GetCpuInfosFile -- use intel SDE")
    	message(STATUS "GetCpuInfosFile -- INASTEMP_ISDE_CPU_ARGS = ${INASTEMP_ISDE_CPU_ARGS}")
    endif()
    
    get_filename_component(
		GetCpuInfosFileExec ${GetCpuInfosFile}
		NAME_WE
	)

    set(GetCpuInfosFileExec "${CMAKE_CURRENT_BINARY_DIR}/${GetCpuInfosFileExec}")

    try_compile(COMPILE_RESULT_VAR
            ${CMAKE_CURRENT_BINARY_DIR}/GetCpuInfos
            ${GetCpuInfosFile}
            OUTPUT_VARIABLE comp
            COPY_FILE ${GetCpuInfosFileExec})

    if(COMPILE_RESULT_VAR)
        if(NOT EXISTS ${GetCpuInfosFileExec})
	        message(FATAL_ERROR "The GetCpuInfosFile compiled file does not exist (${GetCpuInfosFileExec})")
        endif()

        exec_program("sde64 ${INASTEMP_ISDE_CPU_ARGS} -- ${GetCpuInfosFileExec}" ${CMAKE_CURRENT_BINARY_DIR}
                OUTPUT_VARIABLE run
                RETURN_VALUE RUN_RESULT_VAR)
    endif() 
else()
    # Simply try and run
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
          ${CMAKE_CURRENT_BINARY_DIR} ${GetCpuInfosFile}
          COMPILE_OUTPUT_VARIABLE comp
          RUN_OUTPUT_VARIABLE run)
endif()

# If it has successfuly compiled an run
if(COMPILE_RESULT_VAR AND (RUN_RESULT_VAR EQUAL 0) )
	set( CPU_OPTIONS ${run} )
	# For each value
	foreach(optionNode ${run})
		# Get name and value
		string(REPLACE "=" ";" optionNameAndValue ${optionNode})
		list(LENGTH optionNameAndValue optionLength)
		# If we get both
		if(optionLength EQUAL 2)
			list(GET optionNameAndValue 0 optionName)
			list(GET optionNameAndValue 1 optionValue)
			# create cmake variable
			set(CPU_INFO_${optionName} ${optionValue})
		else()
			message(WARNING "GetCpuInfosFile wrong format for ${optionNode}.")
		endif()
	endforeach()
	# output the sentence from the binrary
    if($ENV{VERBOSE})
    	message(STATUS "GetCpuInfosFile -- results : ${CPU_OPTIONS}")
    endif()
else()
	message(WARNING "GetCpuInfosFile -- did not return correctly.")
	message(WARNING "GetCpuInfosFile -- compilation output : ${comp}.")
	message(WARNING "GetCpuInfosFile -- execution output : ${run}.")
endif()

endmacro(GetCpuInfos)
