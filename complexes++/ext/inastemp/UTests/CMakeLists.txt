###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
project(UTESTS_INASTEMP C CXX)

ADD_DEFINITIONS(${INASTEMP_CXX_FLAGS})

OPTION( INASTEMP_BUILD_UTESTS  "Set to ON to build the unit tests" ON )

if(INASTEMP_BUILD_UTESTS)
    OPTION( INASTEMP_ISDE_UTESTS  "Set to ON to run the unit test over sde64" OFF )

    if($ENV{VERBOSE})
        MESSAGE(STATUS "UTests -- INASTEMP_CXX_FLAGS   : ${INASTEMP_CXX_FLAGS}")
        MESSAGE(STATUS "UTests -- INASTEMP_ISDE_UTESTS : ${INASTEMP_ISDE_UTESTS}")
    endif()

    if(INASTEMP_ISDE_UTESTS)
        SET( INASTEMP_ISDE_UTESTS_ARGS "-knl" CACHE STRING "Arguments for sde64"  )
    endif()

    # Find all code files
    file(	
	    GLOB_RECURSE
	    source_tests_files
	    ./*.cpp
	    )


    # Adding the project sources dir as an include dir
    INCLUDE_DIRECTORIES(
         ${INASTEMP_BINARY_DIR}/Src    
         ${INASTEMP_SOURCE_DIR}/Src    
         ${UTESTS_INASTEMP_SOURCE_DIR}
    )


    # Add execs - 1 cpp = 1 exec
    foreach(exec ${source_tests_files}) 
	    get_filename_component(
		    execname ${exec}
		    NAME_WE
	    )

        file(STRINGS "${exec}" generictest REGEX "\@TYPE\@")

	    if(generictest)
            foreach(TYPE ${ALL_USED_TYPES})
                set(typexec "${execname}-${TYPE}")
                configure_file( ${exec} "${typexec}.cpp" )

                if($ENV{VERBOSE})
                        MESSAGE(STATUS "UTests -- Add ${typexec}")
                endif()

	            add_executable(
			            ${typexec}
			            "${typexec}.cpp"
			            )

	            target_link_libraries(
		            ${typexec}
		            ${INASTEMP_LIBRARIES}
		            )

                if(INASTEMP_ISDE_UTESTS)
                    add_test(${typexec} sde64 ${INASTEMP_ISDE_UTESTS_ARGS} -- ${CMAKE_CURRENT_BINARY_DIR}/${typexec})
                else()
                    add_test(${typexec} ${CMAKE_CURRENT_BINARY_DIR}/${typexec})
                endif()
            endforeach()
        else()
            if($ENV{VERBOSE})
                MESSAGE(STATUS "UTests -- Add ${execname}")
            endif()

            add_executable(
				    ${execname}
				    ${exec}
			    )

            target_link_libraries(
	            ${execname}
	            ${INASTEMP_LIBRARIES}
	            )
            if(INASTEMP_ISDE_UTESTS)
                add_test(${typexec} sde64 ${INASTEMP_ISDE_UTESTS_ARGS} -- ${CMAKE_CURRENT_BINARY_DIR}/${execname})
            else()
                add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
            endif()
        endif()
    endforeach(exec)
endif()
