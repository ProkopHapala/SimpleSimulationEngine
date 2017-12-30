
set( NAME "AeroCombat" )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/engines
        ${COMMON_SRCS}/SDL2OGL
        ${OpenGL_INCLUDE_DIR}
        ${GLU_INCLUDE_DIR}
        ${SDL2_INCLUDE_DIR}
)

add_executable( ${NAME}
	main.cpp #FieldPatch.cpp  #AeroCraft.cpp AeroCraftWorld.cpp AeroCraftGUI.cpp
	#$<TARGET_OBJECTS:${NAME}_objects>
	#$<TARGET_OBJECTS:Body>
	#$<TARGET_OBJECTS:AeroSurf>
	#$<TARGET_OBJECTS:AeroCraft>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME} ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable( AeroCraft_editor
	editor.cpp
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:AeroSurf>
	$<TARGET_OBJECTS:AeroCraft>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( AeroCraft_editor ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${NAME} AeroCraft_editor
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
