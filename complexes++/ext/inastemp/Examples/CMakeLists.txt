###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
project(EXAMPLES_INASTEMP C CXX)

OPTION( INASTEMP_BUILD_EXAMPLES  "Set to ON to build examples" ON )

if($ENV{VERBOSE})
    MESSAGE(STATUS "Examples -- INASTEMP_BUILD_EXAMPLES = ${INASTEMP_BUILD_EXAMPLES}")
endif()

if( INASTEMP_BUILD_EXAMPLES )
    file(GLOB all_examples ./*)

    foreach(example_name ${all_examples})
        if(IS_DIRECTORY ${example_name})
          if($ENV{VERBOSE})
              message(STATUS "Examples -- ${example_name}" )
          endif()
          add_subdirectory(${example_name})
        endif()
    endforeach()
endif()
