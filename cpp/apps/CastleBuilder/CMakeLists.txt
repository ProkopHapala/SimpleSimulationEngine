
set( NAME "CastleBuilder" )

add_executable( ${NAME}_main
	${NAME}_main.cpp CastleWorld.cpp
	$<TARGET_OBJECTS:Noise> 
	$<TARGET_OBJECTS:Body2D> 
	$<TARGET_OBJECTS:SDL2OGL> 
	$<TARGET_OBJECTS:TerrainSimplex> 
	$<TARGET_OBJECTS:TiledView> 
)

target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
