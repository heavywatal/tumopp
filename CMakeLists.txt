cmake_minimum_required(VERSION 3.28)
project(tumopp
  VERSION 0.9.3
  LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
cmake_print_variables(CMAKE_INSTALL_PREFIX)

add_compile_definitions(
  _LIBCPP_REMOVE_TRANSITIVE_INCLUDES=1
  FMT_USE_LOCALE=0
  FMT_OPTIMIZE_SIZE=1
)
add_compile_options(
  -Wall -Wextra -Wpedantic
  -Wconversion -Wsign-conversion -Wshadow
  $<IF:$<BOOL:$<CONFIG>>,,-Werror>
  $<IF:$<BOOL:$<CONFIG>>,,-O2>
  $<IF:$<BOOL:$<CONFIG>>,,-g>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>:-march=native>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>:-march=armv8.4-a>
)

macro(Pkg_MakeAvailable)
  cmake_parse_arguments(arg "" "" "" ${ARGN})
  FetchContent_MakeAvailable(${arg_UNPARSED_ARGUMENTS})
  foreach(pkg IN LISTS arg_UNPARSED_ARGUMENTS)
    cmake_print_variables(${pkg}_FOUND ${pkg}_DIR)
    if(NOT ${pkg}_FOUND)
      cmake_print_variables(${pkg}_SOURCE_DIR ${pkg}_BINARY_DIR)
    endif()
  endforeach()
endmacro()

FetchContent_Declare(fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/12.2.0.tar.gz
  URL_HASH SHA256=8b852bb5aa6e7d8564f9e81394055395dd1d1936d38dfd3a17792a02bebd7af0
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 12.2.0
)

FetchContent_Declare(pcglite
  URL https://github.com/heavywatal/pcglite/archive/refs/tags/v0.2.2.tar.gz
  URL_HASH SHA256=f082db23389a653630a9596704bdea9e587486a2cb3a100e106529a4de7b0193
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.2.2
)

FetchContent_Declare(clippson
  URL https://github.com/heavywatal/clippson/archive/refs/tags/v0.8.12.tar.gz
  URL_HASH SHA256=aea4fa10a078a41e31ad2980e883c65025f35510465e8ef7a72a0535c0c4f555
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.8.12
)

FetchContent_Declare(wtl
  URL https://github.com/heavywatal/cxxwtl/archive/refs/tags/v0.11.3.tar.gz
  URL_HASH SHA256=b8df5b253ddb7d8772a432deef7fa5d92eed05e16720c9a1eb14e630a47aea70
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.11.3 COMPONENTS zlib
)

Pkg_MakeAvailable(fmt pcglite clippson wtl)

add_library(${PROJECT_NAME} STATIC)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
add_subdirectory(src)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PROJECT_NAME}
  PUBLIC pcglite::pcglite
  PRIVATE fmt::fmt wtl::wtl wtl::zlib clippson::clippson
)

add_executable(${PROJECT_NAME}-exe src/main.cpp)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}-exe
  EXPORT ${PROJECT_NAME}-config
)

set(config_destination ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${config_destination}
  NAMESPACE ${PROJECT_NAME}::
)

set(version_file ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${version_file}
  COMPATIBILITY AnyNewerVersion
)
install(FILES ${version_file}
  DESTINATION ${config_destination}
)

include(CTest)
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  add_subdirectory(test)
endif()
