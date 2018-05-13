
add_executable       ( test_StructParser test_StructParser.cpp )
target_link_libraries( test_StructParser )

add_executable       ( test_CommandParser test_CommandParser.cpp )
target_link_libraries( test_CommandParser )

add_custom_command(
     TARGET ${PROJECT_NAME} test_StructParser test_CommandParser
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)



