#ifndef LuaHelpers_H
#define LuaHelpers_H

#include <string>
#include <vector>
#include <iostream>

#ifdef __cplusplus
# include <lua5.2/lua.hpp>
#else
# include <lua5.2/lua.h>
# include <lua5.2/lualib.h>
# include <lua5.2/lauxlib.h>
#endif


// https://stackoverflow.com/questions/41387796/access-nested-tables-from-lua-to-c-to-get-values
// https://stackoverflow.com/questions/25940366/passing-array-to-c-as-argument-in-the-stack

void lua_getVec3(lua_State *L, int idx, Vec3d& vec){
    // universal helper function to get Vec3 function argument from Lua to C++ function
    luaL_checktype(L, idx, LUA_TTABLE);
    lua_rawgeti(L, idx, 1); vec.x = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_rawgeti(L, idx, 2); vec.y = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_rawgeti(L, idx, 3); vec.z = lua_tonumber(L, -1); lua_pop(L, 1);
    //lua_pop(L, 3);
}

void lua_getMat3(lua_State *L, int idx, Mat3d& mat){
    // universal helper function to get Vec3 function argument from Lua to C++ function

    luaL_checktype(L, idx, LUA_TTABLE);
    lua_pushinteger(L, 1); lua_gettable(L, idx); lua_getVec3(L, -1, mat.a ); lua_pop(L, 1);
    lua_pushinteger(L, 2); lua_gettable(L, idx); lua_getVec3(L, -1, mat.b ); lua_pop(L, 1);
    lua_pushinteger(L, 3); lua_gettable(L, idx); lua_getVec3(L, -1, mat.c ); lua_pop(L, 1);

    //luaL_checktype(L, idx, LUA_TTABLE);   printf("DEBUG 0\n");
    //lua_pushinteger(L, 1); lua_gettable(L, idx); printf("DEBUG 0.1\n"); lua_getVec3(L, -1, mat.a ); printf("DEBUG 1\n"); lua_pop(L, 1);
    //lua_pushinteger(L, 2); lua_gettable(L, idx); printf("DEBUG 1.1\n"); lua_getVec3(L, -1, mat.b ); printf("DEBUG 2\n"); lua_pop(L, 1);
    //lua_pushinteger(L, 3); lua_gettable(L, idx); printf("DEBUG 2.1\n"); lua_getVec3(L, -1, mat.c ); printf("DEBUG 3\n"); lua_pop(L, 1);
    //lua_rawgeti(L, idx, 1); lua_gettable(L, idx); printf("DEBUG 0.1\n"); lua_getVec3(L, -1, mat.a ); printf("DEBUG 1\n"); lua_pop(L, 1);
    //lua_rawgeti(L, idx, 2); lua_gettable(L, idx); printf("DEBUG 1.1\n"); lua_getVec3(L, -1, mat.b ); printf("DEBUG 2\n"); lua_pop(L, 1);
    //lua_rawgeti(L, idx, 3); lua_gettable(L, idx); printf("DEBUG 2.1\n"); lua_getVec3(L, -1, mat.c ); printf("DEBUG 3\n"); lua_pop(L, 1);
}

#endif

