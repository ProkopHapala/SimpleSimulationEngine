
#ifndef  EditSpaceCraft_h
#define  EditSpaceCraft_h

#include <string>
#include <unordered_map>

#include "LuaHelpers.h"
#include "SpaceCraft.h"
//#include "LuaClass.h"

//using namespace SpaceCrafting;

namespace SpaceCrafting{

char str[4096];
int fontTex = 0;
lua_State  * theLua=0;
SpaceCraft * theSpaceCraft=0;
std::unordered_map<std::string,Material*> materials;
std::unordered_map<std::string,Commodity*> comodities;

int l_Material (lua_State * L){
    Material *mat = new Material();
    Lua::dumpStack(L);
    strcpy( mat->name,  Lua::getStringField(L, "name"    ) );
    mat->density      = Lua::getDoubleField(L, "density" );
    mat->Spull        = Lua::getDoubleField(L, "Spull"   );
    mat->Spush        = Lua::getDoubleField(L, "Spush"   );
    mat->Kpull        = Lua::getDoubleField(L, "Kpull"   );
    mat->Kpush        = Lua::getDoubleField(L, "Kpush"   );
    mat->reflectivity = Lua::getDoubleField(L, "reflectivity" );
    mat->Tmelt        = Lua::getDoubleField(L, "Tmelt"   );
    auto it = materials.find( mat->name );
    if( it == materials.end() ){ materials.insert({mat->name,mat}); }else{ printf( "Material `%s` replaced\n", mat->name ); delete it->second; it->second = mat;  }
    printf( "Material %s dens %g Strength (%g,%g) Stiffness (%g,%g) Refl %g Tmelt %g \n", mat->name, mat->density, mat->Kpull,mat->Kpush,   mat->Spull,mat->Spush,  mat->reflectivity, mat->Tmelt );
    return 0;
};

int l_Node    (lua_State * L){
    Vec3d pos;
    Lua::getVec3(L, 1, pos );

    int id =  theSpaceCraft->nodes.size();
    theSpaceCraft->nodes.push_back( Node(pos) );
    printf( "Node (%g,%g,%g)  ->  %i\n",  pos.x, pos.y, pos.z, id );
    lua_pushnumber(L, id);
    return 1;
};

int l_Rope    (lua_State * L){
    Rope r;
    r.p0    = Lua::getInt   (L,1);
    r.p1    = Lua::getInt   (L,2);
    r.thick = Lua::getDouble(L,3);
    const char * matn = Lua::getString(L,4);
    auto it = materials.find(matn);
    if( it != materials.end() ){r.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    r.id   = theSpaceCraft->ropes.size();
    printf( "Rope (%i,%i) %g %s ->  %i\n", r.p0, r.p1, r.thick, matn, r.id );
    theSpaceCraft->ropes.push_back( r );
    lua_pushnumber(L, r.id);
    return 1;
};

int l_Girder  (lua_State * L){
    //Lua::dumpStack(L);
    Girder g;
    g.p0   = Lua::getInt (L,1);
    g.p1   = Lua::getInt (L,2);
             Lua::getVec3(L,3, g.up );
    g.nseg = Lua::getInt (L,4);
    g.mseg = Lua::getInt (L,5);
             Lua::getVec2(L,6, g.wh );
    const char * matn = Lua::getString(L,7);
    auto it  = materials.find(matn);
    if ( it != materials.end() ){g.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    g.id     = theSpaceCraft->girders.size();
    printf( "Girder (%i,%i) (%g,%g,%g) (%i,%i), (%g,%g) %s ->  %i\n", g.p0, g.p1, g.up.x, g.up.y, g.up.z, g.nseg, g.mseg, g.wh.x, g.wh.y, matn, g.id );
    theSpaceCraft->girders.push_back( g );
    lua_pushnumber(L, g.id);
    return 1;
};

int l_Ring    (lua_State * L){
    //Lua::dumpStack(L);
    Ring o;
    o.p0   = Lua::getInt (L,1); printf("DEBUG 1 \n");
    o.p1   = Lua::getInt (L,2); printf("DEBUG 2 \n");
             Lua::getVec3(L,3, o.up ); printf("DEBUG 3 \n");
    o.nseg = Lua::getInt (L,4); printf("DEBUG 4 \n");
    o.R    = Lua::getDouble(L,5); printf("DEBUG 5 \n");
             Lua::getVec2(L,6, o.wh ); printf("DEBUG 6 \n");
    const char * matn = Lua::getString(L,7);
    auto it  = materials.find(matn);
    if ( it != materials.end() ){o.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id     = theSpaceCraft->rings.size();
    printf( "Ring (%i,%i) (%g,%g,%g) (%i), %g, (%g,%g) %s ->  %i\n", o.p0, o.p1, o.up.x, o.up.y, o.up.z, o.nseg, o.R, o.wh.x, o.wh.y, matn, o.id );
    theSpaceCraft->rings.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Radiator (lua_State * L){
    Radiator rd;
    rd.g1          = Lua::getInt   (L,1);
    rd.g1span.x    = Lua::getDouble(L,2);
    rd.g1span.y    = Lua::getDouble(L,3);
    rd.g2          = Lua::getInt   (L,4);
    rd.g2span.x    = Lua::getDouble(L,5);
    rd.g2span.y    = Lua::getDouble(L,6);
    rd.temperature = Lua::getDouble(L,7);
    //const char * matn = Lua::getString(L,4);
    //auto it = materials.find(matn);
    //if( it != materials.end() ){r.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    rd.id   = theSpaceCraft->radiators.size();
    printf( "Radiator %i(%.2f,%.2f) %i(%.2f,%.2f) %f -> %i\n", rd.g1, rd.g1span.x, rd.g1span.y,   rd.g2, rd.g2span.x, rd.g2span.y,  rd.temperature, rd.id );
    theSpaceCraft->radiators.push_back( rd );
    lua_pushnumber(L, rd.id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Tank (lua_State * L){
    Tank o;
    Lua::getVec3(L, 1, o.pose.pos   );
    Lua::getVec3(L, 2, o.pose.rot.c );
    o.pose.rot.c.getSomeOrtho( o.pose.rot.b, o.pose.rot.a );
    Lua::getVec3(L, 3, o.span );
    const char * matn = Lua::getString(L,4);
    o.id   = theSpaceCraft->tanks.size();
    printf( "Tank pos (%f,%f,%f) dir (%f,%f,%f) l %f r %f -> %i\n",  o.pose.pos.x,o.pose.pos.x,o.pose.pos.x,  o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,   o.span.b, o.span.c, o.id );
    theSpaceCraft->tanks.push_back( o );

    lua_pushnumber(L, o.id);
    return 1;
};


int l_Gun     (lua_State * L){ return 0; };
int l_Truster (lua_State * L){ return 0; };



//int l_Tank    (lua_State * L){ return 0; };
//int l_Radiator(lua_State * L){ return 0; };

void initSpaceCraftingLua(){
    theLua = luaL_newstate();
    lua_State  * L = theLua;
    luaL_openlibs(L);
    lua_register(L, "Material", l_Material );
    lua_register(L, "Node",     l_Node     );
    lua_register(L, "Rope",     l_Rope     );
    lua_register(L, "Girder",   l_Girder   );
    lua_register(L, "Ring",     l_Ring     );
    lua_register(L, "Gun",      l_Gun      );
    lua_register(L, "Truster",  l_Truster  );
    lua_register(L, "Tank",     l_Tank     );
    lua_register(L, "Radiator", l_Radiator );

}

}; // namespace SpaceCrafting

#endif
