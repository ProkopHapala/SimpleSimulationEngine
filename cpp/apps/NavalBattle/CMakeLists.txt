
set( NAME "NavalBattle" )

add_executable( ${NAME}_single 
	${NAME}_single.cpp Battleship.cpp
#	$<TARGET_OBJECTS:Shooter>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_single ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_custom_command(
    TARGET ${NAME}_single
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
