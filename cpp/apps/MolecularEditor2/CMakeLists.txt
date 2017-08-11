
#find_package( SDL2_image REQUIRED )
#find_package( SDL2_ttf   REQUIRED )

set( NAME "MolecularEditor2" )

include_directories( ${COMMON_SRCS}/dynamics/Molecular )

add_executable( ${NAME} 
	molecularEditor2.cpp 
	#$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)

target_link_libraries( ${NAME} ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${NAME}
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
