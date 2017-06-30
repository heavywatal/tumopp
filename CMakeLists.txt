cmake_minimum_required(VERSION 3.1)
project(tumopp CXX)
option(PROFILE "Use gperftools" FALSE)

include_directories(~/local/include)
include_directories(/usr/local/include)
link_directories(~/local/lib)
link_directories(/usr/local/lib)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams system)
message("Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message("Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "~/local" CACHE PATH "install prefix" FORCE)
endif()
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

if(APPLE)
    set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
    set(CMAKE_SHARED_LINKER_FLAGS "-install_name ${CMAKE_INSTALL_PREFIX}/lib/lib${PROJECT_NAME}.dylib")
else()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-soname,${CMAKE_INSTALL_PREFIX}/lib/lib${PROJECT_NAME}.so")
endif()

aux_source_directory(src srcfiles)
message("Source files: ${srcfiles}")

add_executable(${PROJECT_NAME} main.cpp)
add_library(lib${PROJECT_NAME} SHARED ${srcfiles})

target_link_libraries(${PROJECT_NAME} lib${PROJECT_NAME})
target_link_libraries(lib${PROJECT_NAME} ${Boost_LIBRARIES} sfmt)
set_target_properties(lib${PROJECT_NAME} PROPERTIES
  NO_SONAME TRUE
  PREFIX ""
)

if(PROFILE)
  target_link_libraries(${PROJECT_NAME} profiler)
endif()

install(TARGETS ${PROJECT_NAME} lib${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)

install(DIRECTORY src/
  DESTINATION include/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp")
