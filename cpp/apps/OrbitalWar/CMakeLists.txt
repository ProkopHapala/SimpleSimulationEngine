
add_executable( orbitEditor 
	orbitEditor.cpp 
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:Splines>
)
target_link_libraries( orbitEditor ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( test_OptContinuousThrust
	test_OptContinuousThrust
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:DynamicOpt>
)
target_link_libraries( test_OptContinuousThrust  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${PROJECT_NAME} orbitEditor
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)
