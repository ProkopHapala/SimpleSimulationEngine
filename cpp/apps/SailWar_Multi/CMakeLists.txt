
find_package( SDL2_net REQUIRED )

# === server

add_executable( SailWar_server

	SailWar_server.cpp Frigate2D.cpp Yacht2D.cpp Projectile.cpp Gun.cpp SailWarWorld.cpp

	$<TARGET_OBJECTS:Convex2d>
	$<TARGET_OBJECTS:Voronoi>
 
	$<TARGET_OBJECTS:Collisions> 
	$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:Body2D> 
	$<TARGET_OBJECTS:AeroSurf2D> 

	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( SailWar_server ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${SDL2_NET_LIBRARY} )

# === client

add_executable( SailWar_client

	SailWar_client.cpp Frigate2D.cpp Yacht2D.cpp Projectile.cpp Gun.cpp SailWarWorld.cpp

	$<TARGET_OBJECTS:Convex2d>
	$<TARGET_OBJECTS:Voronoi>
 
	$<TARGET_OBJECTS:Collisions> 
	$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:Body2D> 
	$<TARGET_OBJECTS:AeroSurf2D> 

	$<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( SailWar_client ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${SDL2_NET_LIBRARY} )

# === misc

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
add_custom_command(
    TARGET SailWar_client SailWar_server
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)



