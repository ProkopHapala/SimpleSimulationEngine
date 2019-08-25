
set( NAME "MolecularEditor2" )

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable( ${NAME}_main
	${NAME}_main.cpp 
	#$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( ConfSearch 
	ConfSearch.cpp 
	#$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( ConfSearch ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#add_executable( test_reset 	test_reset.cpp 	$<TARGET_OBJECTS:DynamicOpt> 	$<TARGET_OBJECTS:SDL2OGL>   )
#target_link_libraries( test_reset ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( test_MMFF  test_MMFF.cpp  $<TARGET_OBJECTS:DynamicOpt> )
target_link_libraries( test_MMFF )


add_custom_command(
    TARGET ${NAME}_main ConfSearch test_MMFF
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
