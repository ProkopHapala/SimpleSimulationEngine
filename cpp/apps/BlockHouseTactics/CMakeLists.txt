
set( NAME "BlockHouseTactics" )

add_executable( ${NAME}_main
	${NAME}_main.cpp BlockHouseWorld.cpp 
	$<TARGET_OBJECTS:SDL2OGL> 
	$<TARGET_OBJECTS:SoftBody> 
)

target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
