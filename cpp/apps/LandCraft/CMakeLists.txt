
set( NAME "LandCraft" )

add_executable( ${NAME}_main 
	${NAME}_main.cpp #LCWorld.cpp #LTFaction.cpp LTUnit.cpp LTUnitType.cpp LTcommon.h
	#$<TARGET_OBJECTS:Body2D> 
	$<TARGET_OBJECTS:SDL2OGL> 
	$<TARGET_OBJECTS:TerrainGrid2D> 
	$<TARGET_OBJECTS:Noise> 
	#$<TARGET_OBJECTS:TerrainCubic> 
	#$<TARGET_OBJECTS:TiledView> 
)

target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

file(GLOB ${NAME}_HEADERS_ *.h )
add_custom_target( ${NAME}_headers SOURCES ${${NAME}_HEADERS_} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

add_custom_command(
    TARGET ${NAME}_main
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
