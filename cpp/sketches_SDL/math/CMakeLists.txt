
add_executable       ( test_TresholdFunc test_TresholdFunc.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_TresholdFunc ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_FastMath test_FastMath.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_FastMath ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Convex2d test_Convex2d.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Convex2d> )
target_link_libraries( test_Convex2d ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_opt2d test_opt2d.cpp )
target_link_libraries( test_opt2d ${SDL2_LIBRARY} )


