
project(sketches_SDL_3D)

add_executable       ( test_QuatRotSampling test_QuatRotSampling.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_QuatRotSampling ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Projection test_Projection.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Projection ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Solids test_Solids.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Solids ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SceneGraph test_SceneGraph.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_SceneGraph ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereSampling test_SphereSampling.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_SphereSampling ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Patches test_Patches.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Patches ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_GUI test_GUI.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_GUI ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_EditorGizmo test_EditorGizmo.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_EditorGizmo ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MusculeEditor test_MusculeEditor.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_MusculeEditor ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Camera test_Camera.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Camera ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MousePicking test_MousePicking.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_MousePicking ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereTree test_SphereTree.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_SphereTree ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

#add_executable       ( test_ShotHit test_ShotHit.cpp $<TARGET_OBJECTS:SDL2OGL> )
#target_link_libraries( test_ShotHit ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Raytracing test_Raytracing.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Raytracing ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Radiosity test_Radiosity.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Radiosity ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Scatterer test_Scatterer.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Scatterer ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RayScattererMMC test_RayScattererMMC.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_RayScattererMMC ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Elasticity test_Elasticity.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Truss> $<TARGET_OBJECTS:SoftBody> $<TARGET_OBJECTS:DynamicOpt> )
target_link_libraries( test_Elasticity ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY}  )

add_executable       ( test_Electromagnetic test_Electromagnetic.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Electromagnetic ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY}  )

add_executable       ( test_MultipoleAccel test_MultipoleAccel.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_MultipoleAccel ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY}  )

add_executable       ( test_Mesh test_Mesh.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Mesh> )
target_link_libraries( test_Mesh ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RigidBody test_RigidBody.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Collisions> )
target_link_libraries( test_RigidBody ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_VortexLattice test_VortexLattice.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_VortexLattice ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Collision test_Collision.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Collisions> )
target_link_libraries( test_Collision ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_CompressiveParticles test_CompressiveParticles.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_CompressiveParticles ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_AABBTree test_AABBTree.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_AABBTree ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_BlockBuilder test_BlockBuilder.cpp $<TARGET_OBJECTS:SDL2OGL>  )
target_link_libraries( test_BlockBuilder ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TrussBuilder   test_TrussBuilder.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TrussBuilder> $<TARGET_OBJECTS:SoftBody> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_TrussBuilder  ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereGaussSeidel test_SphereGaussSeidel.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_SphereGaussSeidel ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )
add_executable       ( test_Stick test_Stick.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>   )
target_link_libraries( test_Stick ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_GUI test_Solids test_Mesh
    #COMMAND ln -sf ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
