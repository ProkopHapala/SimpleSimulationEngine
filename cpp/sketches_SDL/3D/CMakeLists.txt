
add_executable       ( test_QuatRotSampling test_QuatRotSampling.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_QuatRotSampling ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_ConfDynamics test_ConfDynamics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> )
target_link_libraries( test_ConfDynamics ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Solids test_Solids.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Solids ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Patches test_Patches.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Patches ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_GUI test_GUI.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_GUI ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MousePicking test_MousePicking.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_MousePicking ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Raytracing test_Raytracing.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Raytracing ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Mesh test_Mesh.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Mesh> )
target_link_libraries( test_Mesh ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RigidBody test_RigidBody.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Collisions> )
target_link_libraries( test_RigidBody ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_Collision test_Collision.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Collisions> )
target_link_libraries( test_Collision ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_BlockBuilder test_BlockBuilder.cpp $<TARGET_OBJECTS:SDL2OGL>  )
target_link_libraries( test_BlockBuilder ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TrussBuilder   test_TrussBuilder.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TrussBuilder> $<TARGET_OBJECTS:SoftBody> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_TrussBuilder  ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereTree test_SphereTree.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:SoftBody> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_SphereTree  ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SoftMolecularDynamics.cpp  test_SoftMolecularDynamics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_SoftMolecularDynamics.cpp  ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_GUI test_Solids test_Mesh
    #COMMAND ln -sf ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)
