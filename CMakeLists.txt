cmake_minimum_required(VERSION 3.1)
project(tumopp CXX)

include(GNUInstallDirs)
option(PROFILE "Use gperftools" FALSE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_FLAGS "-O3 -march=native -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(~/local/include)
include_directories(/usr/local/include)
link_directories(~/local/lib)
link_directories(/usr/local/lib)

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams system)
message("Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message("Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

aux_source_directory(src srcfiles)
message("Source files: ${srcfiles}")

add_executable(${PROJECT_NAME} main.cpp)
add_library(lib${PROJECT_NAME} SHARED ${srcfiles})

target_link_libraries(${PROJECT_NAME} lib${PROJECT_NAME})
target_link_libraries(lib${PROJECT_NAME} ${Boost_LIBRARIES} sfmt)
set_target_properties(lib${PROJECT_NAME} PROPERTIES PREFIX "")

if(PROFILE)
  target_link_libraries(${PROJECT_NAME} profiler)
endif()

install(TARGETS ${PROJECT_NAME} lib${PROJECT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp")
