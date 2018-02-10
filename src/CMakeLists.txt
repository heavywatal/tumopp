# waiting for CMake 3.12
get_property(wtl_INCLUDE_DIRS TARGET wtl PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(sfmt_INCLUDE_DIRS TARGET sfmt PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "wtl_INCLUDE_DIRS: ${wtl_INCLUDE_DIRS}")

file(GLOB srcfiles *.cpp)
add_library(obj OBJECT ${srcfiles})
target_include_directories(obj BEFORE
  PRIVATE ${Boost_INCLUDE_DIRS} ${wtl_INCLUDE_DIRS} ${sfmt_INCLUDE_DIRS}
)
