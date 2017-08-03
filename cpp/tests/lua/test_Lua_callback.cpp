// see tutorial
// https://csl.name/post/lua-and-cpp/

#include <stdio.h>

#ifdef __cplusplus
# include <lua5.2/lua.hpp>
#else
# include <lua5.2/lua.h>
# include <lua5.2/lualib.h>
# include <lua5.2/lauxlib.h>
#endif

#include "Vec3.h"
#include "Mat3.h"
#include "LuaHelpers.h"

//lua_State *state = NULL;

extern "C"
int howdy(lua_State* state){
  // The number of function arguments will be on top of the stack.
  int args = lua_gettop(state);
  printf("howdy() was called with %d arguments:\n", args);
  for ( int n=1; n<=args; ++n) {
    printf("  argument %d: '%s'\n", n, lua_tostring(state, n));
  }
  // Push the return value on top of the stack. NOTE: We haven't popped the
  // input arguments to our function. To be honest, I haven't checked if we
  // must, but at least in stack machines like the JVM, the stack will be
  // cleaned between each function call.
  lua_pushnumber(state, 123);
  // Let Lua know how many return values we've passed
  return 1;
}

void doSomePhysics(int n, Vec3d pos, Vec3d vel){
    //printf( "%i (%g,%g,%g) (%g,%g,%g) \n", n, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z );
    Vec3d G = {0.0,0.0,-9.81};
    double dt = 0.03;
    double restitution = -1.0;
    for(int i=0; i<n; i++){
        vel.add_mul(G,dt);
        pos.add_mul(vel,dt);
        printf("%i pos=(%g,%g,%g) (%g,%g,%g)\n", i, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z );
        if( (pos.z<0)&&(vel.z<0) ) vel.z*=restitution;
    }
};

void doSomePhysics2(int n, Vec3d pos, Mat3d mat){
    printf( " -------------- \n" );
    printf( " n=%i pos=(%g,%g,%g) \n", n, pos.x, pos.y, pos.z );
    printf("mat = \n");
    mat.print();
};

extern "C"
int l_doSomePhysics(lua_State* L){
    // lua interface for doSomePhysics
    Vec3d pos,vel;
    int n = lua_tointeger(L, 1);
    lua_getVec3(L, 2, pos );
    lua_getVec3(L, 3, vel );
    doSomePhysics(n,pos,vel);
    //lua_pushnumber(L, 123);
    return 3;
}

    extern "C"
    int l_doSomePhysics2(lua_State* L){
        // lua interface for doSomePhysics
        Vec3d pos;
        Mat3d mat;
        int n = lua_tointeger(L, 1);
        lua_getVec3(L, 2, pos );
        lua_getMat3(L, 3, mat );
        doSomePhysics2(n,pos,mat);
        //lua_pushnumber(L, 123);
        return 3;
    }


void print_error(lua_State* state) {
  // The error message is on top of the stack.
  // Fetch it, print it and then pop it off the stack.
  const char* message = lua_tostring(state, -1);
  puts(message);
  lua_pop(state, 1);
}

void execute(lua_State* state, const char* filename){
    // Make standard libraries available in the Lua object
    luaL_openlibs(state);
    int result;
    // Load the program; this supports both source code and bytecode files.
    result = luaL_loadfile(state, filename);
    if ( result != LUA_OK ) { print_error(state); return; }
    // Finally, execute the program by calling into it.
    // Change the arguments if you're not running vanilla Lua code.
    result = lua_pcall(state, 0, LUA_MULTRET, 0);
    if ( result != LUA_OK ) { print_error(state); return; }
}

int main(int argc, char** argv){
    lua_State *state = luaL_newstate();
    lua_register(state, "doSomePhysics",  l_doSomePhysics);
    lua_register(state, "doSomePhysics2", l_doSomePhysics2);

    //execute( state, "data/hello.lua" );
    //execute( state, "data/callback.lua" );
    execute( state, "data/doSomePhysics.lua" );
}
