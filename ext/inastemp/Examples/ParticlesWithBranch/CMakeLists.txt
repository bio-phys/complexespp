###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
project(EXAMPLES_PARTICLES_INASTEMP C CXX)

ADD_DEFINITIONS(${INASTEMP_CXX_FLAGS})

OPTION( INASTEMP_BUILD_PATTERNS  "Set to ON to build pattern examples" ON )


if($ENV{VERBOSE})
    MESSAGE(STATUS "Examples -- PARTICLES")
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
     ./
)


# Add execs - 1 cpp = 1 exec
foreach(exec ${source_tests_files}) 
    get_filename_component(
	    execname ${exec}
	    NAME_WE
    )

    set(execname "particles-${execname}")

    if($ENV{VERBOSE})
        MESSAGE(STATUS "Examples -- PARTICLES ${execname}")
    endif()

    add_executable(
		    ${execname}
		    ${exec}
	    )

    target_link_libraries(
        ${execname}
        ${INASTEMP_LIBRARIES}
        )
endforeach(exec)

