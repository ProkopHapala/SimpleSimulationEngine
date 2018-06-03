
#find_package( SDL2_image REQUIRED )
#find_package( SDL2_ttf   REQUIRED )

set( NAME "MolecularEditor" )

#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )

include_directories( ${COMMON_SRCS}/dynamics/Molecular )

add_executable( ${NAME}_main
	${NAME}_main.cpp 
	#$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)

target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${NAME}_main
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resource
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
