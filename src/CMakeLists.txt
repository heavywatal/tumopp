execute_process(
  COMMAND git describe --always --dirty
  OUTPUT_VARIABLE GIT_DESCRIBE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(GIT_DESCRIBE)
  set(PROJECT_VERSION ${GIT_DESCRIBE})
endif()
message(STATUS "PROJECT_VERSION: ${PROJECT_VERSION}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp @ONLY
)

if(sfmt_FOUND)
  set(WTL_MT_HPP sfmt.hpp)
else()
  set(WTL_MT_HPP wtl/random.hpp)
endif()
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/random.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/random.hpp @ONLY
)

# waiting for CMake 3.12
get_property(wtl_INCLUDE_DIRS TARGET wtl PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(sfmt_INCLUDE_DIRS TARGET sfmt PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "wtl_INCLUDE_DIRS: ${wtl_INCLUDE_DIRS}")

file(GLOB srcfiles *.cpp)
add_library(obj OBJECT ${srcfiles})
target_include_directories(obj BEFORE
  PRIVATE ${Boost_INCLUDE_DIRS} ${wtl_INCLUDE_DIRS} ${sfmt_INCLUDE_DIRS}
)
set_property(TARGET obj PROPERTY POSITION_INDEPENDENT_CODE TRUE)
