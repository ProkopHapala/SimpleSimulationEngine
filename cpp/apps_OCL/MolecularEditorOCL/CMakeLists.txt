
set( NAME "MolecularEditorOCL" )

#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )

include_directories( ${COMMON_SRCS}/dynamics/Molecular )

add_executable( ${NAME}_main
    ${NAME}_main.cpp 
    $<TARGET_OBJECTS:MolecularEngine>
    $<TARGET_OBJECTS:DynamicOpt>
    $<TARGET_OBJECTS:SDL2OGL>
)

target_link_libraries( ${NAME}_main ${OpenCL_LIBRARY} ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

SET( MOLECULAR_EDITOR_SRCS       ${MY_SRC_DIR}/apps/MolecularEditor2 )
#add_custom_target( ${NAME}_src_dep SOURCES ${MOLECULAR_EDITOR_SRCS} )

add_custom_command(
    TARGET ${NAME}_main
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources  ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/inputs ${CMAKE_CURRENT_SOURCE_DIR}/inputs
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/data   ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/data   ${CMAKE_CURRENT_SOURCE_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/doc    ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MOLECULAR_EDITOR_SRCS}/doc    ${CMAKE_CURRENT_SOURCE_DIR}/doc
)
