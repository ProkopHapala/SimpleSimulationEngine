
#ADD_CUSTOM_TARGET( shaders ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )

add_executable       ( test_DiffractShader test_DiffractShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_DiffractShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereShader test_SphereShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_SphereShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_ShaderDepth test_ShaderDepth.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_ShaderDepth ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MeshOGL3 test_MeshOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_MeshOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_OrbitalRayMarch test_OrbitalRayMarch.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_OrbitalRayMarch ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Texture test_Texture.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_Texture ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SSAO test_SSAO.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_SSAO ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

#add_executable       ( test_ScreenOGL3 test_ScreenOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Body> )
add_executable       ( test_ScreenOGL3 test_ScreenOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_ScreenOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} ${_LIBRARY} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_DiffractShader test_SphereShader test_OrbitalRayMarch test_Texture test_SSAO test_ScreenOGL3 test_MeshOGL3
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)

# instead do softling manually
#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#ADD_CUSTOM_TARGET(link_target ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )

