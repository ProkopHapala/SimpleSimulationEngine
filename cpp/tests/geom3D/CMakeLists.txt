
add_executable       ( test_QuatRotSampling test_QuatRotSampling.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_QuatRotSampling ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

