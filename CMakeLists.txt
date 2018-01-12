cmake_minimum_required(VERSION 3.1)
project(tumopp CXX)

include(GNUInstallDirs)
option(PROFILE "Use gperftools" FALSE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g)
else()
  add_definitions(-DNDEBUG)
endif()
add_compile_options(-O3 -march=native -Wall -Wextra -Wpedantic)

find_package(wtl REQUIRED)
link_libraries(wtl)
find_package(sfmt REQUIRED)
link_libraries(sfmt)

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams system)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
message("Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message("Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

aux_source_directory(src srcfiles)
message("Source files: ${srcfiles}")

set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_library(sharedlib SHARED ${srcfiles})
target_link_libraries(sharedlib ${Boost_LIBRARIES} wtl sfmt)
set_target_properties(sharedlib PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} sharedlib)
if(PROFILE)
  target_link_libraries(${PROJECT_NAME} profiler)
endif()

install(TARGETS ${PROJECT_NAME} sharedlib
  EXPORT ${PROJECT_NAME}-config
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp")
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}
)
