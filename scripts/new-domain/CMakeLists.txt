
add_custom_target(add-domain
                   COMMAND ./add-domain.sh "${CMAKE_SOURCE_DIR}/scripts/new-domain"
                   WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/scripts/new-domain")
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES add-domain)
