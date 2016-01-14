
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

set( NAME "AeroCraft" )


add_library( ${NAME}_objects OBJECT AeroCraft FieldPatch GameWorld GameScreen )

add_executable( ${NAME} 
	main.cpp 
	$<TARGET_OBJECTS:${NAME}_objects> 

	$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:AeroSurf> 
	$<TARGET_OBJECTS:SDL2OGL> 
)

target_link_libraries( ${NAME} ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
