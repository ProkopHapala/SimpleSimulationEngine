
find_package( OpenGL REQUIRED )
find_package( GLU    REQUIRED )
find_package( SDL2   REQUIRED )

add_executable       ( test_AppSDL2OGL test_AppSDL2OGL.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_AppSDL2OGL ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_HashMap2D test_HashMap2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_HashMap2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_HashMap2D_uniformity test_HashMap2D_uniformity.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_HashMap2D_uniformity ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )
