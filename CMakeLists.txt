cmake_minimum_required(VERSION 3.1)
project(tumopp CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "~/local" CACHE PATH "install prefix" FORCE)
endif()
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
set(CMAKE_BUILD_TYPE release)

include_directories(~/local/include)
include_directories(/usr/local/include)
link_directories(~/local/lib)

add_executable(a.out main.cpp)
target_link_libraries(a.out ${PROJECT_NAME})

add_subdirectory(src)

#install(TARGETS a.out RUNTIME DESTINATION bin)

install(DIRECTORY src/ DESTINATION include/${PROJECT_NAME}
        FILES_MATCHING PATTERN "*.hpp")
