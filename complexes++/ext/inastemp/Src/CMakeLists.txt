###########################################################################
# Inastemp - Berenger Bramas MPCDF - 2016
# Under MIT Licence, please you must read the LICENCE file.
###########################################################################
project(LIB_INASTEMP C CXX)

ADD_DEFINITIONS(${INASTEMP_CXX_FLAGS})

# Adding the entire project dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_CURRENT_SOURCE_DIR}
     ${CMAKE_CURRENT_BINARY_DIR}
)

file(
    GLOB_RECURSE 	 	 	 	
 	hpp_in_dir
	./*.hpp ./*.h
)
FOREACH(hpp_file ${hpp_in_dir})	
    file(RELATIVE_PATH relative_hpp_file ${CMAKE_CURRENT_SOURCE_DIR} ${hpp_file})
    if(INASTEMP_DO_INSTALL)
        INSTALL( FILES ${hpp_file} DESTINATION include/${relative_hpp_file} )
    endif()

    if($ENV{VERBOSE})
        MESSAGE(STATUS "Src -- Add ${relative_hpp_file}")
    endif()
ENDFOREACH()

if(INASTEMP_DO_INSTALL)
    INSTALL( FILES "${INASTEMP_BINARY_DIR}/Src/InastempConfig.h" DESTINATION include/Src )
endif()

