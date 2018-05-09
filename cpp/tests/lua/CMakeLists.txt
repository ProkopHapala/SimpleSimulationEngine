
# https://www.gamedev.net/forums/topic/681869-linking-lua-via-cmake/
find_package( Lua52   REQUIRED )
#FIND_PACKAGE(Lua 5.3 REQUIRED)
set(INCLUDE_DIRS ${LUA_INCLUDE_DIR})
add_definitions(-DLUA)
include_directories(${INCLUDE_DIRS})

include_directories(
#        ${COMMON_SRCS}/math
        ${COMMON_SDL_SRCS}/Lua
)

# put this after add_executable
#target_link_libraries(PROJECT_NAME ${LUA_LIBRARY})

add_executable       ( test_Lua_getVar test_Lua_getVar.cpp )
target_link_libraries( test_Lua_getVar ${LUA_LIBRARY} )

add_executable       ( test_Lua_callback test_Lua_callback.cpp )
target_link_libraries( test_Lua_callback ${LUA_LIBRARY} )

add_executable       ( test_Lua_macros test_Lua_macros.cpp )
target_link_libraries( test_Lua_macros ${LUA_LIBRARY} )

#add_library(test_Lua SHARED test_Lua.cpp )

add_custom_command(
     TARGET ${PROJECT_NAME} test_Lua_getVar test_Lua_callback test_Lua_macros
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)



