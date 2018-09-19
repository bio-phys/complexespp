# enable GTest
include_directories(SYSTEM
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)

# build all test files
file(GLOB TEST_SOURCES "*.cpp")
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
# ${Project_Name_Str}_objects is defined in src/CMakeLists.txt
add_executable(${PROJECT_TEST_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:${PROJECT_NAME_STR}_objects>)
target_link_libraries(${PROJECT_TEST_NAME} gtest gtest_main ${LIBRARIES})
add_test(unit-test ${PROJECT_TEST_NAME})
file(COPY "data" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
