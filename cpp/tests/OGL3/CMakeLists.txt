
add_executable       ( test_DiffractShader test_DiffractShader.cpp )
target_link_libraries( test_DiffractShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereShader test_SphereShader.cpp )
target_link_libraries( test_SphereShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Texture test_Texture.cpp )
target_link_libraries( test_Texture ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SSAO test_SSAO.cpp )
target_link_libraries( test_SSAO ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

# instead do softling manually
#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )


