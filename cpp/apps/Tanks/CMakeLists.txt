
find_package( OpenGL REQUIRED )
find_package( GLU REQUIRED    )
find_package( SDL2 REQUIRED   )

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

set( NAME "Tanks" )

add_executable( ${NAME}_single 
	${NAME}_single.cpp Tank.cpp
#	$<TARGET_OBJECTS:Shooter>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_single ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
