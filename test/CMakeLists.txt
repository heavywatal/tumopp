aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
foreach(src IN LISTS source_files)
  cmake_path(GET src STEM name_we)
  add_executable(test-${name_we} ${src})
  set_target_properties(test-${name_we} PROPERTIES CXX_EXTENSIONS OFF)
  target_link_libraries(test-${name_we} PRIVATE ${PROJECT_NAME} fmt::fmt wtl::wtl)
  add_test(NAME ${name_we} COMMAND $<TARGET_FILE:test-${name_we}>)
endforeach()

find_program(BASH_EXE bash)
if(BASH_EXE)
  add_test(NAME cli.sh
    COMMAND ${BASH_EXE} ${CMAKE_CURRENT_SOURCE_DIR}/cli.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
endif(BASH_EXE)
