
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
//int fontTex = 0;
lua_State  * theLua=0;

SpaceCraft * theSpaceCraft=0;

SpaceCraftWorkshop workshop;

Radiosity radiositySolver;

static bool bUseOriginalLuaWrappers = false; 


constexpr int NBIT_kind = 8;
constexpr int NBIT_id   = sizeof(int)*8-NBIT_kind;
constexpr int MASK_id   = (1<<NBIT_id)-1;

inline void i2idKind( int i, int& kind, int& id ){ kind=(i>>NBIT_id);  id=(i&&MASK_id); }
inline int  idKind2i( int kind, int id ){ return (kind<<NBIT_id)|(MASK_id&&id); }

int l_Material (lua_State * L){
    Material *o = new Material();
    //Lua::dumpStack(L);   // DEBUG
    strcpy( o->name,  Lua::getStringField(L, "name"    ) );
    o->density      = Lua::getDoubleField(L, "density" );
    o->Spull        = Lua::getDoubleField(L, "Spull"   );
    o->Spush        = Lua::getDoubleField(L, "Spush"   );
    o->Kpull        = Lua::getDoubleField(L, "Kpull"   );
    o->Kpush        = Lua::getDoubleField(L, "Kpush"   );
    o->reflectivity = Lua::getDoubleField(L, "reflectivity" );
    o->Tmelt        = Lua::getDoubleField(L, "Tmelt"   );
    if( workshop.materials.add(o) && (verbosity>0) ) printf( "Material(%s) replaced\n", o->name );
    o->id =  workshop.materials.getId( o->name );
    if(verbosity>1) o->print();
    lua_pushnumber(L, o->id);
    return 1;
};

/*
int l_PanelMaterial (lua_State * L){
    PanelMaterials *mat = new PanelMaterials();
    //Lua::dumpStack(L); DEBUG
    strcpy( mat->name,  Lua::getStringField(L, "name"    ) );
    auto it = materials.find( mat->name );
    if( it == materials.end() ){ materials.insert({mat->name,mat}); }else{ printf( "Material `%s` replaced\n", mat->name ); delete it->second; it->second = mat;  }
    printf( "Material %s dens %g Strength (%g,%g) Stiffness (%g,%g) Refl %g Tmelt %g \n", mat->name, mat->density, mat->Kpull,mat->Kpush,   mat->Spull,mat->Spush,  mat->reflectivity, mat->Tmelt );
    return 0;
};
*/

int l_StickMaterial(lua_State * L){
    //printf( "l_StickMaterial()\n" );
    char mat_name[NAME_LEN];
    StickMaterial *o = new StickMaterial();
    //Lua::dumpStack(L); //DEBUG
    strcpy( o->name,   Lua::getString(L, 1 ) );
    strcpy( mat_name,  Lua::getString(L, 2 ) );
    o->diameter      = Lua::getDouble(L, 3 ) * 1e-3; // [mm]->[m]
    o->wallThickness = Lua::getDouble(L, 4 ) * 1e-3; // [mm]->[m]
    o->preStrain     = Lua::getDouble(L, 5 );
    //o->area        =  Lua::getDouble(L, 4 );
    o->materialId = workshop.materials.getId( mat_name );
    o->update( workshop.materials.vec[o->materialId] );
    if( workshop.stickMaterials.add(o) && (verbosity>0) ) printf( "StickMaterial(%s) replaced\n", o->name );
    o->id = workshop.stickMaterials.getId( o->name );
    if(verbosity>1) o->print();
    lua_pushnumber(L, o->id);
    return 1;
};

// ToDo: we need to be able generate nodes at some point along girder or rope. Slider should be child of Node
int l_Node    (lua_State * L){
    Node* o = new Node();
    //Vec3d pos;
    o->pos = Lua::getVec3(L, 1);
    o->id =  theSpaceCraft->nodes.size();
    theSpaceCraft->nodes.push_back( o  );
    if(verbosity>1) o->print();
    lua_pushnumber(L, o->id);
    return 1;
};

// n1 = Node( g1, 0.5, p0 )
int l_BoundNode(lua_State * L){
    //printf( "!!!! l_BoundNode()\n" );
    Node* o = new Node();
    //Vec3d pos;
    int gid   = Lua::getInt(L, 1);       // index of girder
    int kind  = Lua::getInt(L, 2);       // type of component
    o->calong = Lua::getDouble(L, 3);   // position along girder
    Vec3d p0  = Lua::getVec3(L, 4);
    o->boundTo = theSpaceCraft->getStructuralComponent( gid, kind);
    //o->along.x = o->boundTo->pointAlong( o->calong, -1, &o->pos, p0 );
    o->updateBound( p0 );
    o->id =  theSpaceCraft->nodes.size();
    theSpaceCraft->nodes.push_back( o  );
    //if(verbosity>1) 
    o->print();
    printf(" -- boundTo: "); o->boundTo->print();
    lua_pushnumber(L, o->id);
    return 1;
};

int l_Slider(lua_State * L){
    //Lua::dumpStack(L);
    Slider* o = new Slider();
    int gid   = Lua::getInt   (L,1);
    int kind  = Lua::getInt   (L,2);
    o->calong = Lua::getDouble(L, 3);   // position along girder
    Vec3d p0  = Lua::getVec3(L, 4);
    o->boundTo = theSpaceCraft->getStructuralComponent( gid, kind);
    o->updateBound( p0 );
    //o->updatePath();  // path vertices are not created yet
    o->id   = theSpaceCraft->sliders.size();
    //if(verbosity>1) 
    o->print();
    printf(" - boundTo:"); o->boundTo->print();
    theSpaceCraft->sliders.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

/*
int l_Slider(lua_State * L){
    //Lua::dumpStack(L);
    Slider* o = new Slider();
    int id1  = Lua::getInt   (L,1);
    int id2  = Lua::getInt   (L,2);
    int it1  = Lua::getInt   (L,3);
    int it2  = Lua::getInt   (L,4);
    int is1  = Lua::getInt   (L,5);
    int is2  = Lua::getInt   (L,6);
    double c1 = Lua::getDouble(L,7);
    double c2 = Lua::getDouble(L,8);
    o->comp1 = theSpaceCraft->getStructuralComponent( id1, it1);
    o->comp2 = theSpaceCraft->getStructuralComponent( id2, it2);
    o->along=Vec2d{c1,c2};
    o->sides=Vec2i{is1,is2};
    o->id   = theSpaceCraft->sliders.size();
    //if(verbosity>1) 
    o->print();
    printf(" - compi1:"); o->comp1->print();
    printf(" - compi2:"); o->comp2->print();
    theSpaceCraft->sliders.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};
*/

int l_Rope    (lua_State * L){
    printf("l_Rope(): Calling make_Rope() C++ method.\n");
    int    node1_id  = Lua::getInt(L,1);
    int    node2_id  = Lua::getInt(L,2);
    double thick_mm  = Lua::getDouble(L,3);
    int    nseg      = Lua::getInt(L,4);
    double preStrain = Lua::getDouble(L,5);
    const char *mat_name = Lua::getString(L,6);
    int id = theSpaceCraft->make_Rope(node1_id, node2_id, thick_mm, nseg, preStrain, mat_name);
    if(id<0) return 0;
    lua_pushnumber(L, id);
    return 1;
};

int l_Rope2 (lua_State * L){
    //printf( "####### l_Rope2()\n" );
    //ToDo: Ropes should be pre-strained (pre-tensioned) to avoid slack, we should set pre-strain force for each rope, the leght should be calculated from the rope material properties and pre-strain force
    Rope* o = new Rope();
    int   gs[2];
    float cs[2];
    Lua::getLuaArr( L, 2, gs, 1 );
    Lua::getLuaArr( L, 2, cs, 2 );
    Node* nd[2];
    Vec3d p0 = (theSpaceCraft->girders[gs[0]]->nodes.x->pos + theSpaceCraft->girders[gs[0]]->nodes.y->pos 
             +  theSpaceCraft->girders[gs[1]]->nodes.x->pos + theSpaceCraft->girders[gs[1]]->nodes.y->pos)*0.25;
    //printf( "l_Rope2() p0 (%g,%g,%g)\n", p0.x, p0.y, p0.z );
    for(int i=0; i<2; i++){ 
        if       (cs[i]<0){  // start point
            //printf( "l_Rope2() node[%i] is start \n", i  );
            nd[i] = theSpaceCraft->girders[gs[i]]->nodes.x;
        }else if (cs[i]>1){  // end point
            //printf( "l_Rope2() node[%i] is end \n", i  );
            nd[i] = theSpaceCraft->girders[gs[i]]->nodes.y;
        }else{   // along girder
            //printf( "l_Rope2() node[%i] is bound \n", i  );
            nd[i] = new Node();
            nd[i]->boundTo = theSpaceCraft->getStructuralComponent( gs[i], (int)ComponetKind::Girder );
            nd[i]->calong = cs[i];
            nd[i]->updateBound( p0 );
            nd[i]->id = theSpaceCraft->nodes.size(); 
            theSpaceCraft->nodes.push_back( nd[i] ); 
        }
        ((Node**)&(o->nodes))[i] = nd[i];
    }
    o->thick = Lua::getDouble(L,3) * 1e-3; // [mm]->[m]
    o->nseg  = Lua::getInt(L,4);
    double preStrain = Lua::getDouble(L,5);
    const char * matn = Lua::getString(L,6);
    //o->face_mat = workshop.panelMaterials.getId( matn );
    //o->face_mat = workshop.stickMaterials.getId( matn );
    o->id   = theSpaceCraft->ropes.size();
    { // define new StickMaterial
        StickMaterial *mat = new StickMaterial();
        //char mat_name[NAME_LEN];
        sprintf( mat->name, "rope%i", o->id ); //printf( "l_Rope() mat->name(%s)\n", mat->name );
        mat->diameter      =  o->thick;
        mat->wallThickness =  o->thick*0.5;
        mat->preStrain     =  preStrain;
        //mat->area          =  o->thick*o->thick*0.25*M_PI;
        mat->materialId    = workshop.materials.getId( matn );
        mat->update( workshop.materials.vec[mat->materialId] );
        if( workshop.stickMaterials.add(mat) && (verbosity>0) ) printf( "StickMaterial(%s) replaced\n", mat->name );
        mat->id = workshop.stickMaterials.getId( mat->name );
        if(verbosity>1) mat->print();
        o->face_mat = mat->id;
    }
    //if(verbosity>1) o->print();
    theSpaceCraft->ropes.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

int l_Girder  (lua_State * L){
    printf("l_Girder(): Calling make_Girder() C++ method.\n");
    int    node1_id = Lua::getInt(L,1);
    int    node2_id = Lua::getInt(L,2);
    Vec3d  up       = Lua::getVec3(L,3);
    int    nseg     = Lua::getInt(L,4);
    int    mseg     = Lua::getInt(L,5);
    Vec2d  wh       = Lua::getVec2(L,6);
    const char* mat_name = Lua::getString(L,7);
    Quat4i st       = Lua::getVec4i(L,8);
    int id = theSpaceCraft->make_Girder(node1_id, node2_id, up, nseg, mseg, wh, mat_name, st);
    if(id<0) return 0;
    lua_pushnumber(L, id);
    return 1;
};

// ToDo: Make l_Ring function which takes as input 3-4 girders to which the wheel should be attached, that way we can avoid specification of axis, and others

int l_Ring    (lua_State * L){
    printf("l_Ring(): Calling make_Ring() C++ method.\n");
    Vec3d  pos      = Lua::getVec3(L,1);
    Vec3d  dir      = Lua::getVec3(L,2);
    Vec3d  up       = Lua::getVec3(L,3);
    int    nseg     = Lua::getInt(L,4);
    double R        = Lua::getDouble(L,5);
    Vec2d  wh       = Lua::getVec2(L,6);
    const char* mat_name = Lua::getString(L,7);
    Quat4i st       = Lua::getVec4i(L,8);
    int id = theSpaceCraft->make_Ring(pos, dir, up, R, nseg, wh, mat_name, st);
    if(id<0) return 0;
    lua_pushnumber(L, id);
    return 1;
};

// this ring is attached to girders by 3-4 nodes, 3 points define a circle
// Ring( [g1,g2,g3,g4], [c1,c1,c1,c1], p0, nseg, wh={4.0,4.0}, "Steel", g1_st )
int l_Ring2    (lua_State * L){
    int   gs[4];
    float cs[4] ;
    Lua::getLuaArr( L, 4, gs, 1 );
    Lua::getLuaArr( L, 4, cs, 2 );
    Vec3d p0        = Lua::getVec3(L,3);
    int    nseg     = Lua::getInt(L,4);
    Vec2d  wh       = Lua::getVec2(L,5);
    const char* mat_name = Lua::getString(L,6);
    Quat4i st       = Lua::getVec4i(L,7);
    int icontrol    = Lua::getInt(L,8);
    int id = theSpaceCraft->make_Ring2(gs, cs, p0, nseg, wh, mat_name, st, icontrol);
    if(id<0) return 0;
    lua_pushnumber(L, id);
    return 1;
};

int l_Weld    (lua_State * L){
    //Lua::dumpStack(L);
    Weld* o = new Weld();
    int ic1     = Lua::getInt   (L,1);
    int ic2     = Lua::getInt   (L,2);
    o->Rmax     = Lua::getDouble(L,3);
    o->face_mat = Lua::getInt   (L,4);
    o->comps.x = theSpaceCraft->getStructuralComponent( ic1, (int)ComponetKind::Girder );
    o->comps.y = theSpaceCraft->getStructuralComponent( ic2, (int)ComponetKind::Girder );
    o->id     = theSpaceCraft->welds.size();
    if(verbosity>1) o->print();
    theSpaceCraft->welds.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Radiator (lua_State * L){
    Radiator* o = new Radiator();
    o->g1          = Lua::getInt   (L,1);
    o->g1span.x    = Lua::getDouble(L,2);
    o->g1span.y    = Lua::getDouble(L,3);
    o->g2          = Lua::getInt   (L,4);
    o->g2span.x    = Lua::getDouble(L,5);
    o->g2span.y    = Lua::getDouble(L,6);
    o->temperature = Lua::getDouble(L,7);
    o->id   = theSpaceCraft->radiators.size();
    if(verbosity>1) o->print();
    theSpaceCraft->radiators.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Shield (lua_State * L){
    Shield* o = new Shield();
    o->g1          = Lua::getInt   (L,1);
    o->g1span.x    = Lua::getDouble(L,2);
    o->g1span.y    = Lua::getDouble(L,3);
    o->g2          = Lua::getInt   (L,4);
    o->g2span.x    = Lua::getDouble(L,5);
    o->g2span.y    = Lua::getDouble(L,6);
    o->id   = theSpaceCraft->shields.size();
    if(verbosity>1) o->print();
    theSpaceCraft->shields.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Tank (lua_State * L){
    Tank* o = new Tank();
    o->pose.pos   = Lua::getVec3(L, 1);
    o->pose.rot.c = Lua::getVec3(L, 2);
    o->pose.rot.c.normalize();
    o->pose.rot.c.getSomeOrtho( o->pose.rot.b, o->pose.rot.a );
    o->span = Lua::getVec3(L, 3);
    const char * matn = Lua::getString(L,4);
    o->id   = theSpaceCraft->tanks.size();
    if(verbosity>1) o->print();
    theSpaceCraft->tanks.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

// Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
int l_Thruster (lua_State * L){
    Thruster* o = new Thruster();
    o->pose.pos   = Lua::getVec3(L, 1);
    o->pose.rot.c = Lua::getVec3(L, 2);
    o->pose.rot.c.normalize();
    o->pose.rot.c.getSomeOrtho( o->pose.rot.b, o->pose.rot.a );
    o->span       = Lua::getVec3(L, 3);
    const char * skind = Lua::getString(L,4);
    o->id   = theSpaceCraft->thrusters.size();
    if(verbosity>1) o->print();
    theSpaceCraft->thrusters.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

int l_Gun     (lua_State * L){
    Gun* o = new Gun();
    o->suppId          = Lua::getInt   (L,1);
    o->suppSpan.x      = Lua::getDouble(L,2);
    o->suppSpan.y      = Lua::getDouble(L,3);
    const char * skind = Lua::getString(L,4);
    o->id   = theSpaceCraft->guns.size();
    if(verbosity>1) o->print();
    theSpaceCraft->guns.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

int l_Rock     (lua_State * L){
    Rock* o = new Rock();
    Vec3d dir,up;
    o->pose.pos = Lua::getVec3(L, 1);
    dir     = Lua::getVec3(L, 2);
    up      = Lua::getVec3(L, 3);
    o->pose.rot.fromDirUp(dir, up);
    o->span = Lua::getVec3(L, 4);
    o->id   = theSpaceCraft->rocks.size();
    if(verbosity>1) o->print();
    theSpaceCraft->rocks.push_back( o );
    lua_pushnumber(L, o->id);
    return 1;
};

int l_Balloon     (lua_State * L){
    Balloon* o = new Balloon();
    Vec3d dir,up;
    o->pose.pos = Lua::getVec3(L, 1);
    dir     = Lua::getVec3(L, 2);
    up      = Lua::getVec3(L, 3);
    o->pose.rot.fromDirUp(dir, up);
    o->span = Lua::getVec3(L, 4);
    o->id   = theSpaceCraft->balloons.size();
    if(verbosity>1) o->print();
    theSpaceCraft->balloons.push_back( o );
    lua_pushnumber(L, o->id);
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
    lua_register(L, "StickMaterial", l_StickMaterial );
    theSpaceCraft->workshop = &workshop; // Initialize the workshop pointer in SpaceCraft
    //lua_register(L, "PanelMaterial", l_PanelMaterial );
    lua_register(L, "Node",     l_Node     );
    lua_register(L, "BoundNode",l_BoundNode);
    lua_register(L, "Rope",     l_Rope     );
    lua_register(L, "Rope2",    l_Rope2    );
    lua_register(L, "Girder",   l_Girder   );
    lua_register(L, "Ring",     l_Ring     );
    lua_register(L, "Ring2",    l_Ring2    );
    lua_register(L, "Slider",   l_Slider   );
    lua_register(L, "Weld",     l_Weld     ); 
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
