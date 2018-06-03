
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

set( NAME "AeroCombat" )

add_executable( ${NAME}_main
	${NAME}_main.cpp 
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable( AeroCraft_editor
	AeroCraft_editor.cpp
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:AeroSurf>
	$<TARGET_OBJECTS:AeroCraft>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( AeroCraft_editor ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
endif()
add_executable( AeroCombatOGL3
	AeroCombatOGL3.cpp 
	$<TARGET_OBJECTS:TerrainGrid2D> 
	$<TARGET_OBJECTS:Noise> 
	$<TARGET_OBJECTS:SDL2OGL3>
)
target_link_libraries( AeroCombatOGL3 ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_custom_command(
    TARGET ${NAME}_main AeroCraft_editor AeroCombatOGL3
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
