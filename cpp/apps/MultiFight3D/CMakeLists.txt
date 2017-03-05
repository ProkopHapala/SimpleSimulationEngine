
set( NAME "MultiFight3D" )

add_executable( ${NAME}_single 
	${NAME}_single.cpp ${NAME}World.cpp Warrior3D.cpp Projectile3D.cpp 
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_single ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
