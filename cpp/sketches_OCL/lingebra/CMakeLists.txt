
#ADD_CUSTOM_TARGET( shaders ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )

add_executable       ( test_mmul test_mmul.cpp )
target_link_libraries( test_mmul ${OpenCL_LIBRARY} )

#add_executable       ( test_mmul test_mmul.cpp $<TARGET_OBJECTS:SDL2OGL3> )
#target_link_libraries( test_mmul ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

#add_custom_command(
#    TARGET ${PROJECT_NAME} test_DiffractShader test_SphereShader test_OrbitalRayMarch test_Texture test_SSAO test_ScreenOGL3 test_MeshOGL3
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
#)

# instead do softling manually
#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#ADD_CUSTOM_TARGET(link_target ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )


