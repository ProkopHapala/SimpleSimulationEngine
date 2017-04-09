
add_executable       ( test_clConvolve2D test_clConvolve2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_clConvolve2D ${OpenCL_LIBRARY} ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_clSpaceFrame2D test_clSpaceFrame2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_clSpaceFrame2D ${OpenCL_LIBRARY} ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_clInterpolate2D test_clInterpolate2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_clInterpolate2D ${OpenCL_LIBRARY} ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_clConvolve2D test_clSpaceFrame2D
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/cl ${CMAKE_CURRENT_BINARY_DIR}/cl
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/py ${CMAKE_CURRENT_BINARY_DIR}/py
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)
