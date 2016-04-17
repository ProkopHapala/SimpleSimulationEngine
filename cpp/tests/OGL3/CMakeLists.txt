
add_executable       ( test_DiffractShader test_DiffractShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_DiffractShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereShader test_SphereShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_SphereShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Texture test_Texture.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_Texture ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SSAO test_SSAO.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_SSAO ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

#add_executable       ( test_ScreenOGL3 test_ScreenOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Body> )
add_executable       ( test_ScreenOGL3 test_ScreenOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_ScreenOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} ${_LIBRARY} )

# instead do softling manually
#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )


