
# === Lua # https://www.gamedev.net/forums/topic/681869-linking-lua-via-cmake/
find_package( Lua52   REQUIRED )
set(INCLUDE_DIRS ${LUA_INCLUDE_DIR} ${MY_SRC_DIR}/libs/KosmoSuite/cpp ${COMMON_SDL_SRCS}/Lua )
add_definitions(-DLUA)
include_directories( ${INCLUDE_DIRS} )

add_executable       ( cad2d cad2d.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( cad2d ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${LUA_LIBRARY} )

file(GLOB cad_HEADERS_ *.h )
add_custom_target( cad_HEADERS SOURCES ${cad_HEADERS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} cad2d
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/py   ${CMAKE_CURRENT_BINARY_DIR}/py
)
