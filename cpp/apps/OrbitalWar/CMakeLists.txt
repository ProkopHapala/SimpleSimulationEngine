
# === Lua # https://www.gamedev.net/forums/topic/681869-linking-lua-via-cmake/
find_package( Lua52   REQUIRED )
set(INCLUDE_DIRS ${LUA_INCLUDE_DIR} ${MY_SRC_DIR}/libs/KosmoSuite/cpp ${COMMON_SDL_SRCS}/Lua )
add_definitions(-DLUA)
include_directories( ${INCLUDE_DIRS} )
#include_directories( ${MY_SRC_DIR}/libs/KosmoSuite/cpp  ${COMMON_SDL_SRCS}/Lua )

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


add_executable( SolarSystemMap
	SolarSystemMap.cpp
	$<TARGET_OBJECTS:Body>
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:Splines>
)
target_link_libraries( SolarSystemMap ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )


add_executable( test_OptContinuousThrust
	test_OptContinuousThrust
	$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:DynamicOpt>
)
target_link_libraries( test_OptContinuousThrust  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( test_SpaceFlightODE
    test_SpaceFlightODE.cpp
    $<TARGET_OBJECTS:SDL2OGL>
)
target_link_libraries( test_SpaceFlightODE  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable( spaceCraftEditor
    spaceCraftEditor.cpp
    $<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:Truss>
    $<TARGET_OBJECTS:SoftBody>
    $<TARGET_OBJECTS:DynamicOpt>
)
target_link_libraries( spaceCraftEditor  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} )

add_executable( spaceCraftDynamics
    spaceCraftDynamics.cpp
    $<TARGET_OBJECTS:SDL2OGL>
    $<TARGET_OBJECTS:Truss>
    $<TARGET_OBJECTS:SoftBody>
    $<TARGET_OBJECTS:DynamicOpt>
)
target_link_libraries( spaceCraftDynamics  ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} )

file(GLOB OrbitalWar_HEADERS_ *.h )
add_custom_target( OrbitalWar_HEADERS SOURCES ${OrbitalWar_HEADERS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} orbitEditor test_OptContinuousThrust test_SpaceFlightODE spaceCraftEditor
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/py ${CMAKE_CURRENT_BINARY_DIR}/py
)
