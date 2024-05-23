cmake_minimum_required(VERSION 3.15)
project(tumopp
  VERSION 0.8.6
  LANGUAGES CXX)

cmake_policy(SET CMP0076 NEW)
set(CMAKE_VERBOSE_MAKEFILE ON)
include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
cmake_print_variables(CMAKE_INSTALL_PREFIX)
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
cmake_print_variables(CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")

set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

function(import_env variable)
  if(DEFINED ENV{${variable}})
    set(${variable} $ENV{${variable}} PARENT_SCOPE)
  endif()
endfunction()

function(find_package_or_fetch package version repo)
  find_package(${package} ${version} QUIET)  # 3.24
  if(${package}_FOUND)
    cmake_print_variables(${package}_DIR)
    return()
  endif()
  string(TOUPPER "${package}" upackage)
  import_env("FETCHCONTENT_SOURCE_DIR_${upackage}")
  if(IS_DIRECTORY "${FETCHCONTENT_SOURCE_DIR_${upackage}}")
    set(tag "HEAD")
  else()
    set(tag "v${version}")
  endif()
  cmake_print_variables(FETCHCONTENT_SOURCE_DIR_${upackage} tag)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_Declare(
    ${package}
    GIT_REPOSITORY https://github.com/${repo}.git
    GIT_TAG ${tag}
    GIT_SHALLOW ON
    # FIND_PACKAGE_ARGS  # 3.24
  )
  # FetchContent_MakeAvailable(${package})  # 3.28
  if(NOT ${package}_POPULATED)
    FetchContent_Populate(${package})
    add_subdirectory(${${package}_SOURCE_DIR} ${${package}_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
  cmake_print_variables(${package}_SOURCE_DIR ${package}_BINARY_DIR)
endfunction()

find_package_or_fetch(wtl 0.9.0 heavywatal/cxxwtl)
find_package_or_fetch(sfmt 0.7.5 heavywatal/sfmt-class)
find_package_or_fetch(clippson 0.8.6 heavywatal/clippson)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
add_library(${PROJECT_NAME})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
add_subdirectory(src)

add_library(common INTERFACE)
target_compile_options(common INTERFACE
  -Wall -Wextra -pedantic
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>:-march=native>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>:-march=armv8.3-a+sha3>
)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PROJECT_NAME}
  PUBLIC wtl::sfmt
  PRIVATE common wtl::wtl clippson::clippson
)

option(BUILD_EXECUTABLE "Build executable file" ON)
if(BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME}-exe src/main.cpp)
  target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME} common wtl::wtl)
  find_package(ZLIB)
  if(ZLIB_FOUND)
    target_compile_definitions(${PROJECT_NAME}-exe PRIVATE ZLIB_FOUND)
    target_link_libraries(${PROJECT_NAME}-exe PRIVATE ZLIB::ZLIB)
  endif()
  set_target_properties(${PROJECT_NAME}-exe PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
  )
  install(TARGETS ${PROJECT_NAME}-exe
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()

set(VERSION_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${VERSION_CONFIG} COMPATIBILITY AnyNewerVersion
)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp")
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
  NAMESPACE ${PROJECT_NAME}::
)
install(FILES ${VERSION_CONFIG}
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

include(CTest)
if(BUILD_TESTING AND ${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  add_subdirectory(test)
endif()
