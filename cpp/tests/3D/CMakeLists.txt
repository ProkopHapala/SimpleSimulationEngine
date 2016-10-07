
add_executable       ( test_QuatRotSampling test_QuatRotSampling.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_QuatRotSampling ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Solids test_Solids.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Solids ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_GUI test_GUI.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_GUI ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MousePicking test_MousePicking.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_MousePicking ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Raytracing test_Raytracing.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Raytracing ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_BlockBuilder test_BlockBuilder.cpp $<TARGET_OBJECTS:SDL2OGL>  )
target_link_libraries( test_BlockBuilder ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TrussBuilder   test_TrussBuilder.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TrussBuilder> $<TARGET_OBJECTS:SoftBody> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_TrussBuilder  ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


add_executable       ( test_SphereTree test_SphereTree.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:SoftBody> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_SphereTree  ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

