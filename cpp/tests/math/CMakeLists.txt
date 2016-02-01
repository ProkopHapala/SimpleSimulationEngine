
add_executable       ( test_TresholdFunc test_TresholdFunc.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_TresholdFunc ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )