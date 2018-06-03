
set( NAME "NonInert" )

add_executable( ${NAME}_main
	${NAME}_main.cpp ${NAME}World.cpp Warrior2D.cpp Projectile2D.cpp 
	$<TARGET_OBJECTS:Body2D>
	$<TARGET_OBJECTS:NBodyWorld2D>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
