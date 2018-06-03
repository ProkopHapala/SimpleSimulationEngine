
set( NAME "SailWar" )

#add_library( ${NAME}_objects OBJECT Frigate2D Yacht2D Projectile Gun GameWorld GameScreen )

add_executable( ${NAME}_main

	${NAME}_main.cpp Frigate2D.cpp Yacht2D.cpp Projectile.cpp Gun.cpp GameWorld.cpp GameScreen.cpp

	$<TARGET_OBJECTS:Convex2d>
	$<TARGET_OBJECTS:Voronoi>
 
	$<TARGET_OBJECTS:Collisions> 
	$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:Body2D> 
	$<TARGET_OBJECTS:AeroSurf2D> 

	$<TARGET_OBJECTS:SDL2OGL> 
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


# === small test

add_executable       ( test_buoyancy test_buoyancy.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Convex2d> )
target_link_libraries( test_buoyancy ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SailPolar test_SailPolar.cpp Yacht2D.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Body2D> $<TARGET_OBJECTS:AeroSurf2D>  )
target_link_libraries( test_SailPolar ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

# === misc

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_command(
    TARGET ${NAME}_main test_buoyancy test_SailPolar
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)



