
set( NAME "MultiFight3D" )

add_executable( ${NAME}_main 
	${NAME}_main.cpp ${NAME}World.cpp #Warrior3D.cpp Projectile3D.cpp 
	$<TARGET_OBJECTS:Warrior3D>
	$<TARGET_OBJECTS:Projectile3D>
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
