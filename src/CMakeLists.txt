execute_process(
  COMMAND git describe --always --dirty
  OUTPUT_VARIABLE GIT_DESCRIBE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(GIT_DESCRIBE)
  set(PROJECT_VERSION ${GIT_DESCRIBE})
endif()
message(STATUS "PROJECT_VERSION: ${PROJECT_VERSION}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp @ONLY
)

# Be patient until 3.13 is popularized
target_sources(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/coord.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tissue.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp
)
