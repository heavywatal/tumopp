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
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp @ONLY
)

file(GLOB srcfiles *.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${srcfiles})
set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
