
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

set( NAME "SwordPlay" )

add_executable( ${NAME}_main
	${NAME}_main.cpp 
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_custom_command(
    TARGET ${NAME}_main
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
