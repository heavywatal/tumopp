cmake_minimum_required(VERSION 3.28)
project(tumopp
  VERSION 0.9.2
  LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
cmake_print_variables(CMAKE_INSTALL_PREFIX)

add_compile_definitions("_LIBCPP_REMOVE_TRANSITIVE_INCLUDES")
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
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 407c905e45ad75fc29bf0f9bb7c5c2fd3475976f
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 12.1.0
)

FetchContent_Declare(pcglite
  GIT_REPOSITORY https://github.com/heavywatal/pcglite.git
  GIT_TAG 4eff9c9088e8b503c63c5147eb2e0a994b07a57a
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.2.2
)

FetchContent_Declare(clippson
  GIT_REPOSITORY https://github.com/heavywatal/clippson.git
  GIT_TAG df72e9dc1ceac38a46db452e8102e807301784c5
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.8.12
)

FetchContent_Declare(wtl
  GIT_REPOSITORY https://github.com/heavywatal/cxxwtl.git
  GIT_TAG f66a31048c0d934c29882082099988d9423b6bf3
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.11.2 COMPONENTS zlib
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
target_compile_definitions(${PROJECT_NAME}
  PRIVATE FMT_USE_LOCALE=0
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
