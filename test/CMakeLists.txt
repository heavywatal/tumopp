link_libraries(${PROJECT_NAME}::${PROJECT_NAME} common wtl::wtl)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
foreach(src IN LISTS source_files)
  get_filename_component(name_we ${src} NAME_WE)
  add_executable(test-${name_we} ${src})
  set_target_properties(test-${name_we} PROPERTIES CXX_EXTENSIONS OFF)
  add_test(NAME ${name_we}
    COMMAND $<TARGET_FILE:test-${name_we}>
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
endforeach()
