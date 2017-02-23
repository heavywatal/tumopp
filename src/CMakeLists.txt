include_directories(~/local/include)
add_library(${PROJECT_NAME} SHARED
  simulation.cpp
  tissue.cpp
  cell.cpp
)

find_package(Boost REQUIRED COMPONENTS program_options filesystem system iostreams)
message(${Boost_INCLUDE_DIRS})
message(${Boost_LIBRARY_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} sfmt)

set_target_properties(${PROJECT_NAME} PROPERTIES
  INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

install(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION lib)
