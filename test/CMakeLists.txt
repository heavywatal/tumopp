include_directories(${CMAKE_SOURCE_DIR}/src)
link_libraries(sharedlib)

macro(add_test_ filename)
  get_filename_component(name_we ${filename} NAME_WE)
  add_executable(test-${name_we} ${filename})
  add_test(NAME ${name_we} COMMAND $<TARGET_FILE:test-${name_we}>)
endmacro()

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
foreach(src IN LISTS source_files)
  add_test_(${src})
endforeach()
