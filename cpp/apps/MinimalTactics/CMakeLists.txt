
set( NAME "MinimalTactics" )

add_executable( ${NAME}_main
	${NAME}_main.cpp TacWorld.cpp Faction.cpp Unit.cpp UnitType.cpp MinimalTacticsCommon.h
	$<TARGET_OBJECTS:Body2D> 
	$<TARGET_OBJECTS:SDL2OGL> 
	$<TARGET_OBJECTS:TerrainCubic> 
	$<TARGET_OBJECTS:TiledView> 
)

target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file(GLOB math_HEADERS_ *.h )
#add_custom_target(math_headers SOURCES ${math_HEADERS_} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

add_custom_command(
    TARGET ${NAME}_main
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
