
add_executable       ( test_AppSDL2OGL test_AppSDL2OGL.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_AppSDL2OGL ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_HashMap2D test_HashMap2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_HashMap2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SphereTree2D test_SphereTree2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_SphereTree2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TileTree2D test_TileTree2D.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_TileTree2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )
#target_compile_options( test_TileTree2D PRIVATE "-O2")

add_executable       ( test_HashMap2D_uniformity test_HashMap2D_uniformity.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_HashMap2D_uniformity ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_HashMap2D_3 test_HashMap2D_3.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_HashMap2D_3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_NBodyColHashMap test_NBodyColHashMap.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Body2D> )
target_link_libraries( test_NBodyColHashMap ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_NBodyWorld test_NBodyWorld.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:NBodyWorld2D> $<TARGET_OBJECTS:Body2D> )
target_link_libraries( test_NBodyWorld ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_PolyLine test_PolyLine.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_PolyLine ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SimplexGrid test_SimplexGrid.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:SimplexGrid> )
target_link_libraries( test_SimplexGrid ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Fluid2D test_Fluid2D.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Fluid2D> )
target_link_libraries( test_Fluid2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TerrainHydraulics test_TerrainHydraulics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TerrainHydraulics> $<TARGET_OBJECTS:Noise> )
target_link_libraries( test_TerrainHydraulics ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TerrainCubic test_TerrainCubic.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TerrainCubic> $<TARGET_OBJECTS:TiledView> )
target_link_libraries( test_TerrainCubic ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_TerrainRBF test_TerrainRBF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TerrainRBF> $<TARGET_OBJECTS:TiledView> )
target_link_libraries( test_TerrainRBF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_CityGen test_CityGen.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TerrainRBF> $<TARGET_OBJECTS:TiledView> )
target_link_libraries( test_CityGen ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


add_executable       ( test_Tris2Rect test_Tris2Rect.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:TiledView> )
target_link_libraries( test_Tris2Rect ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

include_directories( ${COMMON_SRCS}/dynamics/Molecular )
add_executable       ( test_GlobOpt2D test_GlobOpt2D.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:MoleculeWorld2D> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Body2D> )
target_link_libraries( test_GlobOpt2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Plotting2D test_Plotting2D.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:MoleculeWorld2D> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Body2D> )
target_link_libraries( test_Plotting2D ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_Plotting2D
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)
