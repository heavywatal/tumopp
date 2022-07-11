cmake_minimum_required(VERSION 3.13)
project(tumopp
  VERSION 0.8.4
  LANGUAGES CXX)

cmake_policy(SET CMP0076 NEW)
set(CMAKE_VERBOSE_MAKEFILE ON)
include(GNUInstallDirs)
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")
add_compile_options(-march=native -Wall -Wextra -pedantic)

set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

function(find_package_or_git package version repo)
  find_package(${package} ${version} QUIET)
  if(${package}_FOUND)
    return()
  endif()
  set(subdir ${PROJECT_BINARY_DIR}/${package}-v${version})
  if(NOT EXISTS ${subdir})
    find_package(Git)
    execute_process(COMMAND
      ${GIT_EXECUTABLE} clone --recursive --depth=1 --branch=v${version} https://github.com/${repo}.git ${subdir}
    )
  endif()
  message(STATUS "add_subdirectory(${subdir})")
  add_subdirectory(${subdir} EXCLUDE_FROM_ALL)
endfunction()

find_package_or_git(wtl 0.8.5 heavywatal/cxxwtl)
find_package_or_git(clippson 0.8.3 heavywatal/clippson)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
add_library(${PROJECT_NAME})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
add_subdirectory(src)

target_compile_features(${PROJECT_NAME}
  PUBLIC cxx_std_11
  PRIVATE cxx_std_14
)
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
target_link_libraries(${PROJECT_NAME} PRIVATE wtl::wtl clippson::clippson)
find_package(sfmt 0.7.2)
if(sfmt_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE SFMT_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE wtl::sfmt)
endif()

option(BUILD_EXECUTABLE "Build executable file" ON)
if(BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME}-exe src/main.cpp)
  target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME} wtl::wtl)
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
