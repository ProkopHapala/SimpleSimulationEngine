
find_package( OpenGL    REQUIRED )
find_package( GLU       REQUIRED )
find_package( SDL2      REQUIRED )
#find_package( SDL2_image REQUIRED )
#find_package( SDL2_ttf   REQUIRED )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
        ${OpenGL_INCLUDE_DIR}
        ${GLU_INCLUDE_DIR}
        ${SDL2_INCLUDE_DIR}
)

set( NAME "MolecularEditor" )

#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )

add_executable( ${NAME} 
	molecularGUI.cpp 
	#$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)

target_link_libraries( ${NAME} ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${NAME}
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resource
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
