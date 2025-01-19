
include_directories( ${MY_SRC_DIR}/libs/KosmoSuite/cpp ${COMMON_SRCS}/Orbital )

add_executable( orbitEditor
	orbitEditor.cpp
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:Splines>
)
target_link_libraries( orbitEditor ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( spaceTactics
	spaceTactics.cpp
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:Splines>
)
target_link_libraries( spaceTactics ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( constructionBlockApp
    constructionBlockApp.cpp
	$<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:MeshBuilder2>
)
target_link_libraries( constructionBlockApp ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( SolarSystemMap
	SolarSystemMap.cpp
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:Splines>
    #$<TARGET_OBJECTS:MeshBuilder2>
)
target_link_libraries( SolarSystemMap ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


add_executable( test_OptContinuousThrust
	test_OptContinuousThrust.cpp
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:DynamicOpt>
)
target_link_libraries( test_OptContinuousThrust  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( test_SpaceFlightODE
    test_SpaceFlightODE.cpp
    $<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( test_SpaceFlightODE  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


# ====== WITH_LUA    START
if( WITH_LUA )
# === Lua # https://www.gamedev.net/forums/topic/681869-linking-lua-via-cmake/
find_package( Lua52   REQUIRED )
add_definitions(-DLUA)
include_directories( ${LUA_INCLUDE_DIR}  ${COMMON_SDL_SRCS}/Lua )

add_executable( spaceCraftEditor
    spaceCraftEditor.cpp
    $<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:Truss>
    $<TARGET_OBJECTS:SoftBody>
    $<TARGET_OBJECTS:DynamicOpt>
    $<TARGET_OBJECTS:MeshBuilder2>
    $<TARGET_OBJECTS:TrussDynamics_d>
    $<TARGET_OBJECTS:TrussDynamics_f>
)
target_link_libraries( spaceCraftEditor  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} )

add_executable( spaceCraftEditorNew
    SpaceCraftEditorNew.cpp
    $<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:Truss>
    $<TARGET_OBJECTS:SoftBody>
    $<TARGET_OBJECTS:DynamicOpt>
    $<TARGET_OBJECTS:MeshBuilder2>
    $<TARGET_OBJECTS:TrussDynamics_d>
    $<TARGET_OBJECTS:TrussDynamics_f>
)
target_link_libraries( spaceCraftEditorNew  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} )

add_executable( spaceCraftDynamics
    spaceCraftDynamics.cpp
    $<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:Truss>
    $<TARGET_OBJECTS:SoftBody>
    $<TARGET_OBJECTS:DynamicOpt>
    $<TARGET_OBJECTS:MeshBuilder2>
    $<TARGET_OBJECTS:TrussDynamics_d>
    $<TARGET_OBJECTS:TrussDynamics_f>
)
target_link_libraries( spaceCraftDynamics  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} )

if(  WITH_OPENCL )
    include_directories( ${COMMON_SRCS}/OCL )
    add_executable( spaceCraftDynamicsOCL
        spaceCraftDynamicsOCL.cpp
        $<TARGET_OBJECTS:SDL2OGL>
        $<TARGET_OBJECTS:Truss>
        $<TARGET_OBJECTS:SoftBody>
        $<TARGET_OBJECTS:DynamicOpt>
        #$<TARGET_OBJECTS:TrussDynamics_f>
    )
    target_link_libraries( spaceCraftDynamicsOCL  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} ${OpenCL_LIBRARY} )
endif()   # WITH_OPENCL

endif()   # WITH_LUA
# ====== WITH_LUA   END


file(GLOB OrbitalWar_HEADERS_ *.h )
add_custom_target( OrbitalWar_HEADERS SOURCES ${OrbitalWar_HEADERS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} orbitEditor test_OptContinuousThrust test_SpaceFlightODE
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/py ${CMAKE_CURRENT_BINARY_DIR}/py
)
