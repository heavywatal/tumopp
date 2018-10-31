cmake_minimum_required(VERSION 3.4)
project(tumopp
  VERSION 0.6.5
  LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)

function(find_package_or_git package version repo)
  find_package(Git)
  find_package(${package} ${version} QUIET)
  if(${package}_FOUND)
    return()
  endif()
  set(subdir external/${package}-v${version})
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/.git)
    execute_process(COMMAND
      ${GIT_EXECUTABLE} clone --recursive --depth=1 --branch=v${version} https://github.com/${repo}.git ${subdir}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  endif()
  message(STATUS "add_subdirectory(${subdir})")
  add_subdirectory(${subdir} EXCLUDE_FROM_ALL)
endfunction()

find_package(ZLIB REQUIRED)
find_package(sfmt 0.7.1)
find_package_or_git(wtl 0.7.4 heavywatal/cxxwtl)
find_package_or_git(clippson 0.7.3 heavywatal/clippson)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-march=native -Wall -Wextra -pedantic)

if(NOT CMAKE_BUILD_TYPE)
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
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

add_library(${PROJECT_NAME} SHARED)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
add_subdirectory(src)

set_target_properties(${PROJECT_NAME} PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PROJECT_NAME} PRIVATE wtl::wtl clippson::clippson ZLIB::ZLIB)
if(sfmt_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE SFMT_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE wtl::sfmt)
endif()

add_executable(${PROJECT_NAME}-exe main.cpp)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)

set(VERSION_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${VERSION_CONFIG} COMPATIBILITY AnyNewerVersion
)

install(TARGETS ${PROJECT_NAME}-exe ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
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
