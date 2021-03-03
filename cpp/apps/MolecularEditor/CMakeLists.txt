
set( NAME "MolecularEditor" )

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable( ${NAME}_main
	${NAME}_main.cpp 
	#$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( ${NAME}_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( ${NAME}_old_main
	${NAME}_old_main.cpp 
	#$<TARGET_OBJECTS:Body> 
	$<TARGET_OBJECTS:MolecularEngine>  
	$<TARGET_OBJECTS:DynamicOpt>
	$<TARGET_OBJECTS:SDL2OGL>  
)
target_link_libraries( ${NAME}_old_main ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


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

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
endif()
add_executable       ( test_CLCFGO_OGL3 test_CLCFGO_OGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_CLCFGO_OGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )
#target_link_libraries( test_CLCFGO_OGL3 ${OpenGL_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )


add_custom_command(
    TARGET ${NAME}_main ${NAME}_old_main ConfSearch test_MMFF test_CLCFGO_OGL3
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources     ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
)
