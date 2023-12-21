
#ifndef  SpaceCraftComponents_h
#define  SpaceCraftComponents_h

#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <unordered_map>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"

namespace SpaceCrafting{

const int NAME_LEN = 16;

class BodyPose{ public:
    Vec3d pos;
    Mat3d rot;
    //Quat4d  rot;
};


/*
//class BBox : public BodyPose{ public:
class BBox : public BodyPose{ public:
    //int type;     //  elipsoide, box, ... ?
    Vec3d pmin;
    Vec3d pmax;
};
*/

// === SpaceShip Component Types

enum class ComponetKind:int{ Node, Rope, Girder, Ring, Thruster, Gun, Radiator, Shield, Tank, Pipe, Balloon, Rock };

// ==== Materials

class CatalogItem{ public:
	int  id;
	int  kind;
	char name[NAME_LEN] = "\n";

    virtual void print(){ printf("CatalogItem(%i|%s) kidn=%i \n", id, name, kind ); }
    virtual ~CatalogItem(){}   // this is necessary to avoid error see. https://stackoverflow.com/questions/12994920/how-to-delete-an-object-of-a-polymorphic-class-type-that-has-no-virtual-destruct
};

class Material : public CatalogItem { public:
	double density;      // [kg/m3]
	double Spull,Spush;  // [Pa] Strenght
	double Kpull,Kpush;  // [Pa] elastic modulus
	double reflectivity; // [1] reflectivity
	double Tmelt;        // [T] temperature of failure

    virtual void print(){ printf( "Material(%i|%-16s) dens %g Strength (%g,%g) Stiffness (%g,%g) Refl %g Tmelt %g \n", id, name, density, Kpull,Kpush, Spull,Spush, reflectivity, Tmelt ); };
 	// What about heat, electricity etc. ?
};

class Commodity : public CatalogItem { public:
	double density;      // [kg/m3]
	// What about heat, electricity etc. ?

    virtual void print(){ printf( "Commodity(%i|%s) %s dens %g \n", id, name, density ); };
};

class FuelType : public Commodity { public:
    double EnergyDesity;   // [J/kg]

    virtual void print(){ printf( "FuelType(%i|%s) %s dens %g E=%g[J/kg] \n", id, name, density, EnergyDesity ); };
};


class PanelLayer{ public:
    double thickness;    // [m]
    int    materiallId;  // index of material in catalog
};

class StickMaterial : public CatalogItem { public:
    int    materiallId;    // index of material in catalog
    double diameter;       // [m]
    double area;           // [m2]
    // -- auxuliary
    double linearDensity;  // [kg/m]
    double reflectivity;   // [1] reflectivity
    double Tmelt;          // [T] temperature of failure  
    double damping;
    double Kpull,Kpush;    // [N/m] elastic modulus
    double Spull,Spush;    // [N]   Strenght

    void update( const Material* mat ){
        //area          = M_PI*diameter*diameter*0.25;
        //const Material* mat = mats[materiallId];
        linearDensity = area*mat->density;
        Kpull         = area*mat->Kpull;
        Kpush         = area*mat->Kpush;
        Spull         = area*mat->Spull;
        Spush         = area*mat->Spush;
        reflectivity  = mat->reflectivity;
        Tmelt         = mat->Tmelt;
    }

    virtual void print(){ printf( "StickMaterial(%i|%-16s) d=%g[m] A=%g[m2] ldens=%g[kg/m] S(%g,%g)[N] K(%g,%g)[N/m] Refl=%g Tmelt=%g[K] \n", id, name, diameter, area, linearDensity, Kpull,Kpush, Spull,Spush, reflectivity, Tmelt ); };
};

class PanelMaterial : public CatalogItem { public:
    std::vector<PanelLayer> layers; 
    double areaDensity;  // [kg/m2]
    int stickMaterialId; // index of material in catalog  ... because each plate can have to be supported by some frame
     //double reflectivity; // [1] reflectivity
 	// ToDo: how to calculate response to impact of projectile ? Or irradiation by laser ?

    void evalAreaDensity(){
        areaDensity = 0.0;
        for( PanelLayer& l : layers ){ areaDensity += l.thickness * l.materiallId; }
    }

    virtual void print(){ printf( "PanelMaterial(%i|%s) dens=%g[kg/m2] nlayer=%i \n", id, name, areaDensity, layers.size() ); };
};

class ThrusterType : public CatalogItem { public:
	double efficiency;    // how much power is used for acceleration of propelant
	double veMin;         // minimal exhaust velocity [m/s]
	double veMax;         // maximal exhaust velocity [m/s]
	bool   exhaustFuel;   // if true the burned fuel is added to propellant mass
	FuelType  * fuel      = NULL;
	Commodity  * Propelant = NULL;

    virtual void print(){ printf( "ThrusterType(%i|%s) eff=%g ve(%g,%g)[m/s] fuel=%s propelant=%s \n", id, name, efficiency, veMin, veMax, fuel->name, Propelant->name ); };
};

class GunType : public CatalogItem { public:
	double recoil;
	// scaling laws - how performace (power, accuracy, penetration, time of flight ...) scales with size ?

    virtual void print(){ printf( "GunType(%i|%s) recoil=%g \n", id, name, recoil ); };
};

// ==== Components

class Node{ public:
    Vec3d pos;
    std::vector<Vec2i> components; // {kind,index}  // TODO: is this still valid ?
    int id;
    //Node(Vec3d pos):pos(pos){};

    virtual void print(){ printf("Node(id=%i) kidn=%i face_mat=%i \n", id, pos.x,pos.y,pos.z ); };
};

class ShipComponent{ public:
    int    id;
    int    kind;
    int    shape;
    int    face_mat=-1;
    //int    edge_mat=-1;
    //int    p0; // anchor node
    // char name[NAME_LEN];
	double mass;           // [kg]
	//RigidBody pose;
    Vec2i poitRange;  // index of start and end in Truss
    Vec2i stickRange; // --,,--

    virtual void print(){ printf("ShipComponent(id=%i) kidn=%i face_mat=%i \n", id, kind, face_mat ); };
};

class Modul: public ShipComponent { public:
    BodyPose pose;
    Box      bbox;
    Vec3d    span;
    double   volume;

    void pick(const Vec3d& ro, const Vec3d& rd){

    }

    virtual void print() override { printf("Modul(id=%i) kind=%i face_mat=%i \n", id, kind, face_mat ); };
};

class Tank : public Modul { public:
    int commodityId;
	//Commodity * typ;
	//double radius;
	//double length;
	        // [m^3]
	double filled;         // [1]

    virtual void print() override { printf("Tank(id=%i) %g [m^3] of Commodity[%i] \n", id, filled, commodityId ); };
};


class Balloon : public Modul { public:
    virtual void print() override { printf("Balloon(id=%i) kind=%i face_mat=%i \n", id, kind, face_mat ); };
};

class Rock : public Modul { public:
    virtual void print() override { printf("Rock(id=%i) kind=%i face_mat=%i \n", id, kind, face_mat ); };

};

//class GirderType : public CatalogItem { public:
//    int mseg;
//    Vec3d wh;
//};


class NodeLinker : public ShipComponent { public:
    int p0,p1;
    double length;

    virtual void print() override { printf("NodeLinker(id=%i) between(%i,%i) L=%g \n", id, p0,p1, length ); };

    void ray( const Vec3d& ro, const Vec3d& rd ){

    }
};

class Girder : public NodeLinker { public:
    //int p1; // anchor node; p0 inherate
    //double length;
    int nseg;
    int mseg;
    Vec2d wh;  // [m] width and height
    Vec3d up;
    Quat4i st;
    //double SPull,SPush;
    //double kPull,kPush;
    //Material * material;
    //Vec2i poitRange;  // index of start and end in Truss
    //Vec2i stickRange; // --,,---
    //GirderType * type = NULL;
    
    virtual void print() override {
        printf( "Girder(%i) ps(%i,%i) up(%g,%g,%g) nm(%i,%i) wh(%g,%g) st(%i,%i,%i,%i)\n", id, p0, p1, up.x, up.y, up.z, nseg, mseg, wh.x, wh.y, st.x,st.y,st.z,st.w );
    }

};

//class Ring : public NodeLinker { public:
class Ring : public ShipComponent { public:
    //int p1; // anchor node; p0 inherate
    //double length;
    BodyPose pose;
    int nseg;
    double R;
    Vec2d wh;
    Quat4i st;
    //Material * material;
    //Vec2i poitRange;  // index of start and end in Truss
    //Vec2i stickRange; // --,,---
    //GirderType * type = NULL;

    virtual void print() override {
        printf( "Ring(%i) n(%i) pos(%g,%g,%g) rot(%g,%g,%g)(%g,%g,%g)(%g,%g,%g) R=%g wh(%g,%g) st(%i,%i,%i,%i)\n", id, nseg,
        pose.pos.x,   pose.pos.y,   pose.pos.z,
        pose.rot.a.x, pose.rot.a.y, pose.rot.a.z,
        pose.rot.b.x, pose.rot.b.y, pose.rot.b.z,
        pose.rot.c.x, pose.rot.c.y, pose.rot.c.z,
        R, wh.x, wh.y, st.x,st.y,st.z,st.w );
    }

};

class Rope : public NodeLinker { public:
    double thick;
    //Material * material;

    virtual void print() override { printf("Rope(id=%i) between(%i,%i) L=%g \n", id, p0,p1, length ); };
};

class Pipe : public ShipComponent { public:
    double maxFlow;   // units depend on commodity
	ShipComponent * a;
	ShipComponent * b;
};

//class Hub : public ShipComponent { public:
//	int   npipes;      // number of pipes going to hub;
//	Pipe * pipes;      //
//}

// ==== Motors

class Plate : public ShipComponent { public:
    double area;
    int g1,g2;    // anchor girders
    Vec2d g1span; // pos along girdes
    Vec2d g2span;
    int plate_mat;
    //Vec3d normal;
	//int ntris;
	//int * tris;  // triangles from points of spaceship

    virtual void print() override { printf("Plate(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f))  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y ); };

};

class Radiator : public Plate{ public:
    double temperature;
    virtual void print() override { printf("Radiator(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f)) T=%g  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y, temperature ); };
};

class Shield : public Plate{ public:
    virtual void print() override { printf("Plate(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f))  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y ); };
};

class Collector : public Plate{ public:
    virtual void print() override { printf("Collector(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f))  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y ); };
};

// ==== Motors

class Thruster : public Modul { public:
	//ThrusterType * typ = NULL;
    int    type;
	double thrust;
	double power;
	double consumption;

    virtual void print() override { printf("Thruster(id=%i) type=%i kidn=%i face_mat=%i P=%g[W] F=%g C=%g \n", id, type, kind, face_mat, power, thrust, consumption ); };
};

class Rotor : public ShipComponent { public:
    // Move ship componenets with respect to each other in inner manuevers
    double Radius;
    double power;    //  [kg.m^2]
    double torque;   //  [kg.m^2]
    double Inertia;  //  [kg.m^2]  moment of inertia

    virtual void print() override { printf("Rotor(id=%i) kidn=%i face_mat=%i I=%g P=%g[W] tq=%g \n", id, kind, face_mat, Inertia, power, torque ); };
};

class Slider : public ShipComponent { public:
    // allow slide a node over a girder
    int  girder;
    double power;    //  [kg.m^2]
    double Force;

    virtual void print() override { printf("Slider(id=%i) kidn=%i face_mat=%i girder=%i P=%g[W] F=%g \n", id, kind, face_mat, girder, power, Force ); };
};

// === Guns


// Also Accelerator?

class Accelerator : public ShipComponent{ public:
    // TODO: can be also attached to Ring ?
    //       This can be perhaps determined from type

    int   suppType; // support type - ring or girder
    int   suppId;   // support Id   - anchor which girder or ring ?
    Vec2d suppSpan; // pos along the support  (girder or ring)

    double lenght;        // [m]
    double PowerPeak;     // [W]
    double PulseEnergy;   // [J]
    double PulseDuration; // [s]
    double PulsePeriod;   // [s]

    virtual void print() override { printf("Accelerator(id=%i) kidn=%i face_mat=%i supp(%i(%4.2f,%4.2f)) L=%g P=%g[W] E=%g[J] ts(%g,%g)[s] \n", id, kind, face_mat, suppType,suppSpan.x,suppSpan.y, lenght, PowerPeak, PulseEnergy, PulseDuration, PulsePeriod); };

    //std::vector<int> anchors; // anchor points
};


class Gun : public Accelerator{ public:
	int  gunId;             // index of gun type in catalog
    double Aperture;        // [m^2]
    double divergence;      // [1] tangens of angle
    // attached to girder?
    // incorporated in girder?

    virtual void print() override { printf("Gun(id=%i) kidn=%i face_mat=%i supp(%i(%4.2f,%4.2f)) A=%g[m^2] D=%g[1] L=%g[m] P=%g[W] E=%g[J] ts(%g,%g)[s] \n", id, kind, face_mat, suppType,suppSpan.x,suppSpan.y, Aperture, divergence, lenght, PowerPeak, PulseEnergy, PulseDuration, PulsePeriod); };

};

template<typename T>
class Dict{ public:
    std::unordered_map<std::string,int> map;
    std::vector<T*>                     vec;

    int getId( const char* name )const{
        auto it = map.find(name);
        if( it != map.end() ){ return it->second; }else{ return -1; }
    }

    T* get( const char* name )const{
        auto it = map.find(name);
        if( it != map.end() ){ return vec[it->second]; }else{ return 0; }
    }

    int add( T* mat, bool bDel=true ){
        auto it = map.find( mat->name );
        if( it == map.end() ){ int i=vec.size(); vec.push_back(mat); map.insert({mat->name,i}); return i; }else{ if(bDel)delete vec[it->second]; vec[it->second] = mat; return -1; }
    }
};


// === SpaceCraftWorkshop

class SpaceCraftWorkshop{ public:
    bool bPrint = true;
    Dict<Material>      materials;
    Dict<Commodity>     commodities;
    Dict<FuelType>      fuels;
    Dict<PanelMaterial> panelMaterials;
    Dict<StickMaterial> stickMaterials;
    Dict<ThrusterType>  thrusterTypes;
    Dict<GunType>       gunTypes;

    int add_Material ( const char* name, double density, double Spull, double Spush, double Kpull, double Kpush, double reflectivity, double Tmelt ){
        Material *mat = new Material();
        //auto it = materials.find( mat->name );
        //if( it == materials.end() ){ materials.insert({mat->name,mat}); }else{ printf( "Material `%s` replaced\n", mat->name ); delete it->second; it->second = mat;  }
        strcpy( mat->name, name );
        mat->density     = density;
        mat->Spull       = Spull;
        mat->Spush       = Spush;
        mat->Kpull       = Kpull;
        mat->Kpush       = Kpush;
        mat->reflectivity= reflectivity;
        mat->Tmelt       = Tmelt;
        materials.add( mat );
        if(bPrint)mat->print();
        return 0;
    };

    /*
    int l_PanelMaterial (lua_State * L){
        PanelMaterials *mat = new PanelMaterials();
        Lua::dumpStack(L);
        strcpy( mat->name,  Lua::getStringField(L, "name"    ) );
        auto it = materials.find( mat->name );
        if( it == materials.end() ){ materials.insert({mat->name,mat}); }else{ printf( "Material `%s` replaced\n", mat->name ); delete it->second; it->second = mat;  }
        printf( "Material %s dens %g Strength (%g,%g) Stiffness (%g,%g) Refl %g Tmelt %g \n", mat->name, mat->density, mat->Kpull,mat->Kpush,   mat->Spull,mat->Spush,  mat->reflectivity, mat->Tmelt );
        return 0;
    };
    */

}; // class SpaceCraftWorkshop

} // namespace SpaceCrafting

#endif
