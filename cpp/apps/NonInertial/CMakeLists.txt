
find_package( OpenGL REQUIRED )
find_package( GLU REQUIRED    )
find_package( SDL2 REQUIRED   )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
        ${OpenGL_INCLUDE_DIR}
        ${GLU_INCLUDE_DIR}
        ${SDL2_INCLUDE_DIR}
)

set( NAME "NonInert" )

add_executable( ${NAME}_seats 
	${NAME}_seats.cpp ${NAME}World.cpp Warrior2D.cpp Projectile2D.cpp 
	$<TARGET_OBJECTS:Body2D>
	$<TARGET_OBJECTS:NBodyWorld2D>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_seats ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
