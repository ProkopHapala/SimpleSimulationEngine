
#ifndef  EditSpaceCraft_h
#define  EditSpaceCraft_h

#include <string>
#include <unordered_map>

#include "LuaHelpers.h"
#include "SpaceCraft.h"
//#include "LuaClass.h"

#include "TriangleRayTracer.h"
#include "Radiosity.h"

//using namespace SpaceCrafting;

namespace SpaceCrafting{

char str[4096];
int fontTex = 0;
lua_State  * theLua=0;
SpaceCraft * theSpaceCraft=0;
std::unordered_map<std::string,Material*>  materials;
std::unordered_map<std::string,Commodity*> comodities;

Radiosity radiositySolver;


constexpr int NBIT_kind = 8;
constexpr int NBIT_id   = sizeof(int)*8-NBIT_kind;
constexpr int MASK_id   = (1<<NBIT_id)-1;

inline void i2idKind( int i, int& kind, int& id ){ kind=(i>>NBIT_id);  id=(i&&MASK_id); }
inline int  idKind2i( int kind, int id ){ return (kind<<NBIT_id)|(MASK_id&&id); }

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
    Rope o;
    o.p0    = Lua::getInt   (L,1);
    o.p1    = Lua::getInt   (L,2);
    o.thick = Lua::getDouble(L,3);
    const char * matn = Lua::getString(L,4);
    auto it = materials.find(matn);
    if( it != materials.end() ){o.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id   = theSpaceCraft->ropes.size();
    printf( "Rope (%i,%i) %g %s ->  %i\n", o.p0, o.p1, o.thick, matn, o.id );
    theSpaceCraft->ropes.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

int l_Girder  (lua_State * L){
    //Lua::dumpStack(L);
    Girder o;
    o.p0   = Lua::getInt (L,1);
    o.p1   = Lua::getInt (L,2);
             Lua::getVec3(L,3, o.up );
    o.nseg = Lua::getInt (L,4);
    o.mseg = Lua::getInt (L,5);
             Lua::getVec2(L,6, o.wh );
    const char * matn = Lua::getString(L,7);
    auto it  = materials.find(matn);
    if ( it != materials.end() ){o.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id     = theSpaceCraft->girders.size();
    printf( "Girder (%i,%i) (%g,%g,%g) (%i,%i), (%g,%g) %s ->  %i\n", o.p0, o.p1, o.up.x, o.up.y, o.up.z, o.nseg, o.mseg, o.wh.x, o.wh.y, matn, o.id );
    theSpaceCraft->girders.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

int l_Ring    (lua_State * L){
    //Lua::dumpStack(L);
    Ring o;
    //o.p0   = Lua::getInt (L,1);
    //o.p1   = Lua::getInt (L,2);
    // Lua::getVec3(L,3, o.up );
    Lua::getVec3(L,1, o.pose.pos   ); //print(o.pose.pos);   printf("pos\n");
    Lua::getVec3(L,2, o.pose.rot.c ); //print(o.pose.rot.c); printf("dir\n"); // ax
    Lua::getVec3(L,3, o.pose.rot.b ); //print(o.pose.rot.b); printf("up\n");// up
    o.pose.rot.fromDirUp(o.pose.rot.c,o.pose.rot.b);
    o.nseg = Lua::getInt (L,4);
    o.R    = Lua::getDouble(L,5);
             Lua::getVec2(L,6, o.wh );
    const char * matn = Lua::getString(L,7);
    auto it  = materials.find(matn);
    if ( it != materials.end() ){o.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id     = theSpaceCraft->rings.size();
    printf( "Ring (%i,%i) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%i), %g, (%g,%g) %s ->  %i\n", o.nseg,
        o.pose.pos.x, o.pose.pos.y, o.pose.pos.z,
        o.pose.rot.a.x, o.pose.rot.a.y, o.pose.rot.a.z,
        o.pose.rot.b.x, o.pose.rot.b.y, o.pose.rot.b.z,
        o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,
        o.R, o.wh.x, o.wh.y, matn, o.id );
    //printf( "Ring (%i,%i) (%g,%g,%g) (%i), %g, (%g,%g) %s ->  %i\n", o.p0, o.p1, o.up.x, o.up.y, o.up.z, o.nseg, o.R, o.wh.x, o.wh.y, matn, o.id );
    theSpaceCraft->rings.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Radiator (lua_State * L){
    Radiator o;
    o.g1          = Lua::getInt   (L,1);
    o.g1span.x    = Lua::getDouble(L,2);
    o.g1span.y    = Lua::getDouble(L,3);
    o.g2          = Lua::getInt   (L,4);
    o.g2span.x    = Lua::getDouble(L,5);
    o.g2span.y    = Lua::getDouble(L,6);
    o.temperature = Lua::getDouble(L,7);
    //const char * matn = Lua::getString(L,4);
    //auto it = materials.find(matn);
    //if( it != materials.end() ){r.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id   = theSpaceCraft->radiators.size();
    printf( "Radiator %i(%.2f,%.2f) %i(%.2f,%.2f) %f -> %i\n", o.g1, o.g1span.x, o.g1span.y,   o.g2, o.g2span.x, o.g2span.y,  o.temperature, o.id );
    theSpaceCraft->radiators.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Shield (lua_State * L){
    Shield o;
    o.g1          = Lua::getInt   (L,1);
    o.g1span.x    = Lua::getDouble(L,2);
    o.g1span.y    = Lua::getDouble(L,3);
    o.g2          = Lua::getInt   (L,4);
    o.g2span.x    = Lua::getDouble(L,5);
    o.g2span.y    = Lua::getDouble(L,6);
    //const char * matn = Lua::getString(L,4);
    //auto it = materials.find(matn);
    //if( it != materials.end() ){r.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id   = theSpaceCraft->shields.size();
    printf( "Radiator %i(%.2f,%.2f) %i(%.2f,%.2f) %f -> %i\n", o.g1, o.g1span.x, o.g1span.y,   o.g2, o.g2span.x, o.g2span.y, o.id );
    theSpaceCraft->shields.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Tank (lua_State * L){
    Tank o;
    Lua::getVec3(L, 1, o.pose.pos   );
    Lua::getVec3(L, 2, o.pose.rot.c );
    o.pose.rot.c.normalize();
    o.pose.rot.c.getSomeOrtho( o.pose.rot.b, o.pose.rot.a );
    Lua::getVec3(L, 3, o.span );
    const char * matn = Lua::getString(L,4);
    o.id   = theSpaceCraft->tanks.size();
    printf( "Tank pos (%f,%f,%f) dir (%f,%f,%f) l %f r %f -> %i\n",  o.pose.pos.x,o.pose.pos.x,o.pose.pos.x,  o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,   o.span.b, o.span.c, o.id );
    theSpaceCraft->tanks.push_back( o );

    lua_pushnumber(L, o.id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Thruster (lua_State * L){
    Thruster o;
    Lua::getVec3(L, 1, o.pose.pos   );
    Lua::getVec3(L, 2, o.pose.rot.c );
    o.pose.rot.c.normalize();
    o.pose.rot.c.getSomeOrtho( o.pose.rot.b, o.pose.rot.a );
    Lua::getVec3(L, 3, o.span );
    const char * skind = Lua::getString(L,4);
    o.id   = theSpaceCraft->thrusters.size();
    printf( "Thruster pos (%f,%f,%f) dir (%f,%f,%f) l %f r %f -> %i\n",  o.pose.pos.x,o.pose.pos.x,o.pose.pos.x,  o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,   o.span.b, o.span.c, o.id );
    theSpaceCraft->thrusters.push_back( o );

    lua_pushnumber(L, o.id);
    return 1;
};

int l_Gun     (lua_State * L){
    Gun o;
    o.suppId     = Lua::getInt   (L,1);
    o.suppSpan.x  = Lua::getDouble(L,2);
    o.suppSpan.y  = Lua::getDouble(L,3);
    const char * skind = Lua::getString(L,4);
    //auto it = materials.find(matn);
    //if( it != materials.end() ){r.material = it->second;}else{ printf( "Material `%s` not found\n", matn ); }
    o.id   = theSpaceCraft->guns.size();
    printf( "Gun %i(%.2f,%.2f) %s -> %i\n", o.suppId, o.suppSpan.x, o.suppSpan.y,   skind, o.id );
    theSpaceCraft->guns.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};


int l_Rock     (lua_State * L){
    Rock o;
    Vec3d dir,up;
    Lua::getVec3(L, 1, o.pose.pos   );

    Lua::getVec3(L, 2, dir );
    Lua::getVec3(L, 3, up  );
    o.pose.rot.fromDirUp(dir, up);
    Lua::getVec3(L, 4, o.span );

    o.id   = theSpaceCraft->rocks.size();
    //printf( "Gun %i(%.2f,%.2f) %s -> %i\n", o., o.suppSpan.x, o.suppSpan.y,   skind, o.id );
    theSpaceCraft->rocks.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

int l_Balloon     (lua_State * L){
    Balloon o;
    Vec3d dir,up;
    Lua::getVec3(L, 1, o.pose.pos   );

    Lua::getVec3(L, 2, dir );
    Lua::getVec3(L, 3, up  );
    o.pose.rot.fromDirUp(dir, up);
    Lua::getVec3(L, 4, o.span );

    o.id   = theSpaceCraft->balloons.size();
    //printf( "Gun %i(%.2f,%.2f) %s -> %i\n", o., o.suppSpan.x, o.suppSpan.y,   skind, o.id );
    theSpaceCraft->balloons.push_back( o );
    lua_pushnumber(L, o.id);
    return 1;
};

//int l_Truster (lua_State * L){ return 0; };
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
    lua_register(L, "Thruster", l_Thruster );
    lua_register(L, "Tank",     l_Tank     );
    lua_register(L, "Radiator", l_Radiator );
    lua_register(L, "Shield",   l_Shield   );
    lua_register(L, "Balloon",  l_Balloon  );
    lua_register(L, "Rock",     l_Rock     );
}

}; // namespace SpaceCrafting

#endif
