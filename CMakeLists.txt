cmake_minimum_required(VERSION 3.1)
project(tumopp CXX)

include(GNUInstallDirs)
option(PROFILE "Use gperftools" FALSE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams system)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")

find_package(wtl)
find_package(sfmt)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)
add_compile_options(-march=native -Wall -Wextra -Wpedantic)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)

add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:obj>)
target_include_directories(${PROJECT_NAME} BEFORE
  INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_link_libraries(${PROJECT_NAME} PRIVATE wtl sfmt ${Boost_LIBRARIES})

add_executable(${PROJECT_NAME}-exe main.cpp)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME})
if(PROFILE)
  target_link_libraries(${PROJECT_NAME}-exe PRIVATE profiler)
endif()
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)

install(TARGETS ${PROJECT_NAME}-exe ${PROJECT_NAME}
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

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
