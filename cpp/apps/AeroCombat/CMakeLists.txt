
set( NAME "AeroCombat" )

add_executable( ${NAME} 
	main.cpp FieldPatch.cpp AeroCraftWorld.cpp AeroCraftGUI.cpp #AeroCraft.cpp
	#$<TARGET_OBJECTS:${NAME}_objects> 
	$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:AeroSurf>
	$<TARGET_OBJECTS:AeroCraft>
	$<TARGET_OBJECTS:SDL2OGL> 
)

target_link_libraries( ${NAME} ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${NAME}
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
