#ifndef LuaClass_h
#define LuaClass_h

#include <lua5.2/lua.hpp>

class LuaClass{ public:
    lua_State* L;

    void getError(          int i,  const char * s ){ printf( "LuaERROR [%i] : %s\n", s     ); }
    void clean   () { int n = lua_gettop(L); lua_pop(L, n); }

    void dumpStack(){
        int n = lua_gettop( L);  printf("LuaStack[%i] = {",n);
        for (int i=1; i<=n; i++){
            int t = lua_type( L, i);
            switch (t) {
                case LUA_TSTRING:  printf(" %i:'%s';", i, lua_tostring ( L, i)); break;
                case LUA_TBOOLEAN: printf(" %i:%s;",   i, lua_toboolean( L, i) ? "true" : "false"); break;
                case LUA_TNUMBER:  printf(" %i:%g;",   i, lua_tonumber ( L, i)); break;
                default:           printf(" %i:%s;",   i, lua_typename ( L, t)); break;
            }
        }
        printf("}\n");
    }

    bool getBool(int i){ return (bool)lua_toboolean(L, i); }
    bool getBool()     { return getBool(-1); }

    double getDouble(int i ) {
        if(!lua_isnumber(L, i)) { getError( i, "Not a double"); }
        return lua_tonumber(L, i);
    }
    double getDouble() { return getDouble(-1); };

    long getInt(int i) {
        long li=0;
        if(lua_isnumber(L, i)) {
            double f = lua_tonumber(L, i);
            long   li = floor(f);
            if( li != f ) getError( i, "Not an integer");
        }else{
            getError( i, "Not a number");
        }
        return li;
    }
    long getInt() { return getInt(-1); };

    const char * getString(int i) {
        const char * s = NULL;
        if(lua_isstring(L, i)) {  return lua_tostring(L, i);       }
        else                   {  getError( i, "Not a string"); return NULL; }
    }
    const char * getString(){ return getString(-1); };

    double getDoubleField( const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        //if (!lua_isnumber(L, -1)) getErro( -1, "field not a number" );
        //int result = (int)lua_tonumber(L, -1);
        double f= getDouble(-1);
        lua_pop(L, 1);        // remove number
        return f;
    }

    long getIntField( const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        long li= getDouble(-1);
        lua_pop(L, 1);        // remove number
        return li;
    }

    const char * getStringField( const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        //if (!lua_isnumber(L, -1)) getErro( -1, "field not a number" );
        //int result = (int)lua_tonumber(L, -1);
        const char * s= getString(-1);
        lua_pop(L, 1);        // remove number
        return s;
    }

    void getVec3( int idx, Vec3d& vec){
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); vec.x = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); vec.y = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 3); vec.z = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_pop(L, 3);
    }
    void getVec3( Vec3d& vec){  getVec3( -1, vec); }

    void getMat3( int idx, Mat3d& mat){
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_pushinteger(L, 1); lua_gettable(L, idx); getVec3(-1, mat.a ); lua_pop(L, 1);
        lua_pushinteger(L, 2); lua_gettable(L, idx); getVec3(-1, mat.b ); lua_pop(L, 1);
        lua_pushinteger(L, 3); lua_gettable(L, idx); getVec3(-1, mat.c ); lua_pop(L, 1);
    }
    void getMat3( Mat3d& mat){ getMat3( -1, mat); };

    void getDoubleVector( std::vector<double>& vec ) {
        luaL_checktype(L, -1, LUA_TTABLE);
        lua_pushnil(L);
        while(lua_next(L, -2)) { vec.push_back( lua_tonumber(L, -1) ); }
        clean();
    }

    bool getVar(const char * nameStr) {
        /*
        char wstr[64];
        const char * ch0 = strpbrk(nameStr, "." );
        //nameStr = NULL;

        lua_getglobal(L, wstr);
        if(lua_isnil(L, -1)) { getError(-1, "no such global"); return false; }

        bool next = true;
        while (next){
            const char * ch1 = strpbrk(nameStr, "." );
            if( ch1 == NULL ){ next=false; ch1=ch0+strlen(ch0); };
            strncpy(wstr,ch0,ch1-ch0);
        };
        */
        int level = 0;
        int i0 = 0;
        char wstr[64];
        for(int i=0; i<256; i++) {
            char ch = nameStr[i];
            if( (ch=='.') || (ch=='\n') ){
                strncpy(wstr, nameStr+i0, i-i0 );
                if(level == 0)       { lua_getglobal(L,     wstr); }
                else                 { lua_getfield (L, -1, wstr); }
                if(lua_isnil(L, -1)) { printf( "LuaERROR: %s of %s not found \n", wstr, nameStr ); return false; }
                else                 { level++; }
            }
        }
        return true; // success
        /*
        for( int i = 0; i < variableName.size(); i++) { // walk down the hierarchy
            if(variableName.at(i) == '.') {
                if(level == 0) { lua_getglobal(L,     var.c_str()); }
                else           { lua_getfield (L, -1, var.c_str()); }
                if(lua_isnil(L, -1)) { printError(variableName, var + " is not defined"); return false; }
                else                 { var = ""; level++; }
            } else { var += variableName.at(i);  }
        }
        if(level == 0) { lua_getglobal(L,     var.c_str()); }
        else           { lua_getfield (L, -1, var.c_str()); }
        if(lua_isnil(L, -1)) { printError(variableName, var + " is not defined"); return false; }
        return true; // success
        */
    }

    virtual void init(){ L = luaL_newstate(); };
};

#endif

