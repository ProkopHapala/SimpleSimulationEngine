
project(sketches_SDL_OGL3)

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
endif()

#ADD_CUSTOM_TARGET( shaders ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )

add_executable       ( test_DiffractShader test_DiffractShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_DiffractShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereShader test_SphereShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_SphereShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_ShaderDepth test_ShaderDepth.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_ShaderDepth ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_LandScape test_LandScape.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_LandScape ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_DrawOGL3 test_DrawOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_DrawOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Instances test_Instances.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_Instances ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Atoms test_Atoms.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_Atoms ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Vegetation test_Vegetation.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_Vegetation ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Horizont test_Horizont.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_Horizont ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_VAOs test_VAOs.cpp $<TARGET_OBJECTS:SDL2OGL3>  )
target_link_libraries( test_VAOs ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_PatchesOGL3 test_PatchesOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_PatchesOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_StencilTextures test_StencilTextures.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_StencilTextures ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_VolumetricTexture test_VolumetricTexture.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_VolumetricTexture ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RenderToTexture test_RenderToTexture.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_RenderToTexture ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Sprites test_Sprites.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_Sprites ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_AntiAliasing test_AntiAliasing.cpp $<TARGET_OBJECTS:SDL2OGL3>  )
target_link_libraries( test_AntiAliasing ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MeshOGL3 test_MeshOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_MeshOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_GeometryShader test_GeometryShader.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_GeometryShader ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Tubes test_Tubes.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_Tubes ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_OrbitalRayMarch test_OrbitalRayMarch.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_OrbitalRayMarch ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Texture test_Texture.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_Texture ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SSAO test_SSAO.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_SSAO ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

#add_executable       ( test_ScreenOGL3 test_ScreenOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> $<TARGET_OBJECTS:Body> )
add_executable       ( test_ScreenOGL3 test_ScreenOGL3.cpp $<TARGET_OBJECTS:SDL2OGL3> )
target_link_libraries( test_ScreenOGL3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY} ${SDL2_LIBRARY} )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} ${_COMMON_SHADERS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_DiffractShader test_SphereShader test_OrbitalRayMarch test_Texture test_SSAO test_ScreenOGL3 test_MeshOGL3 test_Instances
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)

# instead do softling manually
#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#ADD_CUSTOM_TARGET(link_target ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )


