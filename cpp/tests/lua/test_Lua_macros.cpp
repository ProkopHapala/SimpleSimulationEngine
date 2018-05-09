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

// DUMP_MACRO
// https://stackoverflow.com/questions/985403/seeing-expanded-c-macros
#define CALL(macro, arguments) macro arguments
#define STR(...) STR_(__VA_ARGS__)
#define STR_(...) # __VA_ARGS__
#define DUMP_MACRO(macro, ...) \
    do { \
        puts ( \
            "'" \
            # macro STR(DUMP_MACRO_ARGS_ ## __VA_ARGS__) \
            "'\n expands to \n'" \
            STR(CALL(macro, DUMP_MACRO_ARGS_ ## __VA_ARGS__)) \
            "'\n" \
        ); \
    } while (0)
/* helpers for DUMP_MACRO, add more if required */
#define DUMP_MACRO_ARGS_
#define DUMP_MACRO_ARGS_0 ()
#define DUMP_MACRO_ARGS_1 (<1>)
#define DUMP_MACRO_ARGS_2 (<1>, <2>)
#define DUMP_MACRO_ARGS_3 (<1>, <2>, <3>)


// References:
// https://codecraft.co/2014/11/25/variadic-macros-tricks/
// http://ptspts.blogspot.cz/2013/11/how-to-apply-macro-to-all-arguments-of.html
// https://stackoverflow.com/questions/824639/variadic-recursive-preprocessor-macros-is-it-possible
// http://saadahmad.ca/cc-preprocessor-metaprogramming-lists-and-for_each/
// https://stackoverflow.com/questions/1872220/is-it-possible-to-iterate-over-arguments-in-variadic-macros

#define _GET_NTH_ARG(_1,_2,_3,_4,_5,N,...) N
#define MAKE_FUNC_1(fname,TR,T1)             TR l_##fname(T1 x1)            {return (TR)(x1); }
#define MAKE_FUNC_2(fname,TR,T1,T2)          TR l_##fname(T1 x1,T2 x2)      {return (TR)(x1+x2); }
#define MAKE_FUNC_3(fname,TR,T1,T2,T3)       TR l_##fname(T1 x1,T2 x2,T3 x3){return (TR)(x1+x2+x3); }
#define MAKE_FUNC_4(fname,TR,T1,T2,T3,T4)    TR l_##fname(T1 x1,T2 x2,T3 x3,T4 x4){return (TR)(x1+x2+x3+x4); }
#define MAKE_FUNC_5(fname,TR,T1,T2,T3,T4,T5) TR l_##fname(T1 x1,T2 x2,T3 x3,T4 x4,T5 x5){return (TR)(x1+x2+x3+x4+x5); }

#define MAKE_FUNC(fname,TR,...) _GET_NTH_ARG(__VA_ARGS__, \
        MAKE_FUNC_5(fname,TR,__VA_ARGS__), \
        MAKE_FUNC_4(fname,TR,__VA_ARGS__), \
        MAKE_FUNC_3(fname,TR,__VA_ARGS__), \
        MAKE_FUNC_2(fname,TR,__VA_ARGS__), \
        MAKE_FUNC_1(fname,TR,__VA_ARGS__), \
    )

MAKE_FUNC(func1,double,float)
MAKE_FUNC(func2,double,float,int)
MAKE_FUNC(func3,double,float,int,long)
MAKE_FUNC(func4,double,float,int,long,float)
MAKE_FUNC(func5,double,float,int,long,float,double)
//MAKE_FUNC_1(func1,double,float)
//MAKE_FUNC_2(func2,double,float,int)
//MAKE_FUNC_3(func3,double,float,int,long)

#define LUA_GET_int(i)     Lua::getInt   (L,i)
#define LUA_GET_float(i)   (float)Lua::getDouble(L,i)
#define LUA_GET_double(i)  Lua::getDouble(L,i)
#define LUA_GET_string(i)  Lua::getString(L,i)

#define LUA_PUSH_int(a)    lua_pushnumber(L,a)
#define LUA_PUSH_float(a)  lua_pushnumber(L,a)
#define LUA_PUSH_double(a) lua_pushnumber(L,a)
#define LUA_PUSH_float(a)  lua_pushstring(L,a)

//#define LUA_ARGS_2(T1,T2)  LUA_GET_##T1(1), LUA_GET_##T2(2)

#define LUA_ARG_(T)
#define LUA_GET_(T)
#define LUA_ARG_1(T,...)   LUA_GET_##T(1) LUA_ARG_2(__VA_ARGS__)
#define LUA_ARG_2(T,...)  ,LUA_GET_##T(2) LUA_ARG_3(__VA_ARGS__)
#define LUA_ARG_3(T,...)  ,LUA_GET_##T(3) LUA_ARG_4(__VA_ARGS__)
#define LUA_ARG_4(T,...)  ,LUA_GET_##T(4) LUA_ARG_5(__VA_ARGS__)
#define LUA_ARG_5(T,...)  ,LUA_GET_##T(5) LUA_ARG_6(__VA_ARGS__)
#define LUA_ARG_6(T,...)  ,LUA_GET_##T(6) LUA_ARG_7(__VA_ARGS__)
#define LUA_ARG_7(T,...)  ,LUA_GET_##T(7) LUA_ARG_8(__VA_ARGS__)
#define LUA_ARG_8(T,...)  ,LUA_GET_##T(8) LUA_ARG_ (__VA_ARGS__)

//#define  MAKE_LUA_FUNC_2(TR,fname,T1,T2) int l_##fname(lua_State * L){ LUA_PUSH_##TR( fname( LUA_GET_##T1(1), LUA_GET_##T2(2) ) ); return 1; }

#define  MAKE_LUA_FUNC(TR,fname,...) int l_##fname(lua_State * L){ LUA_PUSH_##TR( fname( LUA_ARG_1(__VA_ARGS__) ) ); return 1; }

double add(float x, float y){
    return x+y;
};

//MAKE_LUA_FUNC_2( double, add, float, float )
//MAKE_LUA_FUNC  ( double, add, float, float )

int main(int argc, char** argv){
    lua_State *state = luaL_newstate();
    //lua_register(state, "doSomePhysics",  l_doSomePhysics);
    //lua_register(state, "doSomePhysics2", l_doSomePhysics2);
    //execute( state, "data/doSomePhysics.lua" );

    //DUMP_MACRO( MAKE_LUA_FUNC_2( double, add, float, float ) );
    DUMP_MACRO( MAKE_LUA_FUNC( double, add, float, float ) );

    /*
    DUMP_MACRO( MAKE_FUNC(func1,double,float) );
    DUMP_MACRO( MAKE_FUNC(func2,double,float,int) );
    DUMP_MACRO( MAKE_FUNC(func3,double,float,int,long));
    DUMP_MACRO( MAKE_FUNC(func4,double,float,int,long,float) );
    DUMP_MACRO( MAKE_FUNC(func5,double,float,int,long,float,double) );
    */

    printf("result %f \n", l_func1(15.100) );
    printf("result %f \n", l_func2(15.100,3) );
    printf("result %f \n", l_func3(15.100,3,6) );
    printf("result %f \n", l_func4(15.100,3,6,9.6) );
    printf("result %f \n", l_func5(15.100,3,6,9.6,11.2) );
}
