execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git describe --always --abbrev=0 --dirty
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git show -s --format=%ci
  OUTPUT_VARIABLE GIT_COMMIT_TIME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "commit ${GIT_COMMIT_HASH} [${GIT_BRANCH}]")
message(STATUS "Date:  ${GIT_COMMIT_TIME}")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp @ONLY
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
