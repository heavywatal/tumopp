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
  version.cpp.in
  version.cpp @ONLY
)

target_sources(${PROJECT_NAME} PRIVATE
  cell.cpp
  coord.cpp
  simulation.cpp
  tissue.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
