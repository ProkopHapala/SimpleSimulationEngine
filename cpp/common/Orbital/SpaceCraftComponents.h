
#ifndef  SpaceCraftComponents_h
#define  SpaceCraftComponents_h

#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <unordered_map>

#include "datatypes.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "macroUtils.h"

#include "containers.h"

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

//enum class ComponetKind:int{ Node, Rope, Girder, Ring, Thruster, Gun, Radiator, Shield, Tank, Pipe, Balloon, Rock };
enum class ComponetKind:int{Node=0,ShipComponent=1,StructuralComponent=2,NodeLinker=3,Girder=4,Ring=5,Rope=6,Pipe=7,Plate=8,Radiator=9,Shield=10,Collector=11,Thruster=12,Rotor=13,Slider=14,Accelerator=15,Gun=16,Modul=17,Tank=18,Balloon=19,Rock=20};
// ==== Materials

class CatalogItem{ public:
	int  id;
	int  kind;
	char name[NAME_LEN] = "\n";

    virtual void print() const { printf("CatalogItem(%i|%s) kidn=%i \n", id, name, kind ); }
    virtual ~CatalogItem(){}   // this is necessary to avoid error see. https://stackoverflow.com/questions/12994920/how-to-delete-an-object-of-a-polymorphic-class-type-that-has-no-virtual-destruct
};

class Material : public CatalogItem { public:
	double density;      // [kg/m3]
	double Spull,Spush;  // [Pa] Strenght
	double Kpull,Kpush;  // [Pa] elastic modulus
	double reflectivity; // [1] reflectivity
	double Tmelt;        // [T] temperature of failure

    virtual void print()const override{ printf( "Material(%i|%-16s) dens=%5.2e[kg/m3] S(%5.2e,%5.2e)[Pa] K(%5.2e,%5.2e)[Pa] Refl=%5.2e Tmelt %4.0f[K]\n", id, name, density, Spull,Spush,  Kpull,Kpush, reflectivity, Tmelt ); };
 	// What about heat, electricity etc. ?
};

class Commodity : public CatalogItem { public:
	double density;      // [kg/m3]
	// What about heat, electricity etc. ?

    virtual void print()const override{ printf( "Commodity(%i|%s) %s dens %g \n", id, name, density ); };
};

class FuelType : public Commodity { public:
    double EnergyDesity;   // [J/kg]

    virtual void print()const override{ printf( "FuelType(%i|%s) %s dens %g E=%g[J/kg] \n", id, name, density, EnergyDesity ); };
};


class PanelLayer{ public:
    double thickness;    // [m]
    int    materialId;  // index of material in catalog
};

class StickMaterial : public CatalogItem { public:
    int    materialId;     // index of material in catalog
    double diameter;       // [m]
    double wallThickness;  // [m]
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
        //const Material* mat = mats[materialId];
        double din    = diameter - 2*wallThickness;
        area          = M_PI*(diameter*diameter - din*din)*0.25; 
        linearDensity = area*mat->density;
        Kpull         = area*mat->Kpull;
        Kpush         = area*mat->Kpush;
        Spull         = area*mat->Spull;
        Spush         = area*mat->Spush;
        reflectivity  = mat->reflectivity;
        Tmelt         = mat->Tmelt;
    }

    virtual void print()const override{ printf( "StickMaterial(%3i|%-16s) d=%8.5f[m] w=%8.5f[mm] A=%8.5e[m2] ldens=%5.2e[kg/m] S(%5.2e,%5.2e)[N] K(%5.2e,%5.2e)[N] Refl=%5.2e Tmelt=%4.0f[K] \n", id, name, diameter, wallThickness*1e+3, area, linearDensity, Spull,Spush, Kpull,Kpush, reflectivity, Tmelt ); };
};

class PanelMaterial : public CatalogItem { public:
    std::vector<PanelLayer> layers; 
    double areaDensity;  // [kg/m2]
    int stickMaterialId; // index of material in catalog  ... because each plate can have to be supported by some frame
     //double reflectivity; // [1] reflectivity
 	// ToDo: how to calculate response to impact of projectile ? Or irradiation by laser ?

    void evalAreaDensity(){
        areaDensity = 0.0;
        for( PanelLayer& l : layers ){ areaDensity += l.thickness * l.materialId; }
    }

    virtual void print()const override{ printf( "PanelMaterial(%i|%s) dens=%g[kg/m2] nlayer=%i \n", id, name, areaDensity, layers.size() ); };
};

class ThrusterType : public CatalogItem { public:
	double efficiency;    // how much power is used for acceleration of propelant
	double veMin;         // minimal exhaust velocity [m/s]
	double veMax;         // maximal exhaust velocity [m/s]
	bool   exhaustFuel;   // if true the burned fuel is added to propellant mass
	FuelType  * fuel      = NULL;
	Commodity  * Propelant = NULL;

    virtual void print()const override{ printf( "ThrusterType(%i|%s) eff=%g ve(%g,%g)[m/s] fuel=%s propelant=%s \n", id, name, efficiency, veMin, veMax, fuel->name, Propelant->name ); };
};

class GunType : public CatalogItem { public:
	double recoil;
	// scaling laws - how performace (power, accuracy, penetration, time of flight ...) scales with size ?

    virtual void print()const override{ printf( "GunType(%i|%s) recoil=%g \n", id, name, recoil ); };
};


// ====

class Path{ public:
    int    n;    // number of points in the path
    int*   ps;   // indexes of points in the path
    double cur;  // position along the path
    bool closed=false; // is it line point-to-point or closed loop? 

    inline void   realloc(int n_){ printf("Path::realoc(%i)\n",n_);   n=n_; _realloc(ps,n); }
    inline int    icur   (  ){ return (int)cur;          }
    inline double dcur   (  ){ return cur-(int)cur;      }
};

// ====================
// ==== Components
// ====================

class Object{ public:
    int    id;
    int    kind;

    virtual ~Object(){  };
    virtual void print(bool bShort=false)const{ printf("Object(id=%i,kind=%i)",id,kind); }
    virtual int component_kind(){ return -1; }
    //virtual Object (){ kind = component_kind(); }
};

class ShipComponent : public Object { public:
    //int    id;
    //int    kind;
    //int    compKind = (int)ComponetKind::ShipComponent;
    int    shape;
    int    face_mat=-1;
    //int    edge_mat=-1;
    //int    p0; // anchor node
    // char name[NAME_LEN];
	double mass;           // [kg]
	//RigidBody pose;
    Vec2i pointRange;  // index of start and end in Truss
    Vec2i stickRange; // --,,--

    virtual ~ShipComponent(){};
    virtual void print(bool bShort=false)const{ if(bShort){printf("ShipComponent(id=%i)",id);}else{
        printf("ShipComponent(id=%i) kidn=%i face_mat=%i \n", id, kind, face_mat );} 
    }
    virtual int component_kind(){ return (int)ComponetKind::ShipComponent; }; 
};


//class GirderType : public CatalogItem { public:
//    int mseg;
//    Vec3d wh;
//};


/*
class NodeLinker : public ShipComponent { public:
    int p0,p1;
    double length;
    
    virtual void print() override { printf("NodeLinker(id=%i) between(%i,%i) L=%g \n", id, p0,p1, length ); };
    void ray( const Vec3d& ro, const Vec3d& rd ){}
};
*/

class Node; // forward declaration

class StructuralComponent : public ShipComponent { public:
    //Quat4i nodes;
    vec4<Node*> nodes;
    //double length;
    virtual void print(bool bShort=false)const override;
    void ray( const Vec3d& ro, const Vec3d& rd ){}
    virtual int component_kind(){ return (int)ComponetKind::StructuralComponent; };

    virtual double rotMat( Mat3d& rot)const = 0;
    virtual int nearSide  ( Vec3d p, const Mat3d* rot=0 ) const = 0;
    virtual int pointAlong( double c, int side, Vec3d* pout=0 ) const = 0;
    virtual int sideToPath( int side, int* inds ) const =0;

};

class Node : public Object{ public:
    int ivert=-1; // to which vertex int the mesh/truss it corresponds
    Vec3d pos;
    //std::vector<Vec2i> components; // {kind,index}  // TODO: is this still valid ?
    //int id;
    StructuralComponent* boundTo=0; // node can be bound to a girder, rope or ring. if boundTo==0 then node is free in space
    double calong;           // position along the bound component
    Vec2i  along{-1,-1};     // index of
    double length;           // length of the bound component
    //Node(Vec3d pos):pos(pos){};

    virtual ~Node(){};
    int updateBound(Vec3d p0=Vec3dZero){ if(boundTo){ if(along.y<0)along.y=boundTo->nearSide(p0); along.x = boundTo->pointAlong( calong, along.y, &pos); }else{ along.x=-1; } return along.x; } 
    virtual void print(bool bShort=false)const{  if(bShort){printf("Node(id=%i)",id);}else{ 
        printf("Node(id=%i) iv=%i pos(%g,%g,%g) \n", id, ivert, pos.x,pos.y,pos.z ); if(boundTo){printf(" -- boundTo(along.x=%i calong=%g ", along.x, calong ); boundTo->print(true); printf(")\n");} } 
    }
};

void StructuralComponent::print(bool bShort)const{ 
    if(bShort){ printf("StructuralComponent(id=%i)",id); }else{
        printf("StructuralComponent(id=%i) nodes(%i,%i,%i,%i) \n", id, nodes.x->id,nodes.y->id,nodes.z->id,nodes.w->id );
    } 
}


class NodeLinker : public StructuralComponent { public:
    double length;
    virtual void print(bool bShort=false)const override {  if(bShort){printf("NodeLinker(id=%i)",id);}else{
        printf("NodeLinker(id=%i) nodes(%i,%i) L=%g \n", id, nodes.x,nodes.y, length ); }
    }
    void ray( const Vec3d& ro, const Vec3d& rd ){}
    virtual int component_kind(){ return (int)ComponetKind::NodeLinker; };

    //virtual double rotMat( Mat3d& rot)const override = 0;
    //virtual int nearSide  ( Vec3d p, const Mat3d* rot=0) const override = 0;
    //virtual int pointAlong( double c, int side=-1, Vec3d* pout=0, Vec3d p0=Vec3dZero ) const override = 0;
};


class Girder : public NodeLinker { public:
    //int p1; // anchor node; p0 inherate
    //double length;
    int nseg;
    int mseg;
    //Vec2i nseg;
    Vec2d wh;  // [m] width and height
    Vec3d up;
    Quat4i st;
    //double SPull,SPush;
    //double kPull,kPush;
    //Material * material;
    //Vec2i pointRange;  // index of start and end in Truss
    //Vec2i stickRange; // --,,---
    //GirderType * type = NULL;
    
    virtual void print(bool bShort=false)const override { if(bShort){printf("Girder(id=%i)",id);}else{
        printf( "Girder(%i) nodes(%i,%i) up(%g,%g,%g) nm(%i,%i) wh(%g,%g) st(%i,%i,%i,%i) pointRange(%i,%i)\n", id, nodes.x->id, nodes.y->id, up.x, up.y, up.z, nseg, mseg, wh.x, wh.y, st.x,st.y,st.z,st.w, pointRange.x,pointRange.y ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Girder; };

    virtual double rotMat( Mat3d& rot)const override{
        rot.c = nodes.y->pos - nodes.x->pos;
        double length = rot.c.normalize();
        rot.b = up;
        rot.b.makeOrthoU( rot.c );
        rot.a.set_cross( rot.c, rot.b );
        rot.a.normalize();
        return length;
    }
    virtual int nearSide( Vec3d p, const Mat3d* rot=0)const override{
        Mat3d rot_;
        if(rot==0){ rotMat(rot_); rot=&rot_; }
        Vec3d d = p - nodes.x->pos;
        double ca = rot->a.dot(d);
        double cb = rot->b.dot(d);
        //printf( "Girder::nearSide() ca,cb(%g,%g) d(%g,%g,%g) rot(%g,%g,%g)(%g,%g,%g)(%g,%g,%g)\n", ca, cb,  d.x,d.y,d.z,   rot->a.x,rot->a.y,rot->a.z, rot->b.x,rot->b.y,rot->b.z, rot->c.x,rot->c.y,rot->c.z );
        int side;
        if    ( fabs(ca)>fabs(cb) ){ if( ca<0 ){ side=0; }else{ side=1; } }  // height
        else                       { if( cb<0 ){ side=2; }else{ side=3; } }  // width
        return side;
    };
    virtual int pointAlong( double c, int side, Vec3d* pout=0 )const override{
        int i;
        if(side>1){ i = (int)(c*(nseg+0.5)-0.5); }
        else      { i = (int)(c*(nseg+0.5));     }
        if(i<0){i=0;}else if(i>=nseg){i=nseg-1;}
        int ip=i*4+side;
        if(pout){
            Mat3d rot;
            rotMat(rot);
            Vec3d d = (nodes.y->pos - nodes.x->pos)*(1./(nseg*2+1));
             *pout = nodes.x->pos + d*(i*2.0+1);
            if(side>1){ 
                pout->add( d );
                pout->add_mul( rot.b, (side-2.5)*2*wh.y );
            }else{     
                pout->add_mul( rot.a, (side-0.5)*2*wh.x );
            }            
        }
        //if(side_out){ *side_out=side; }
        //printf( "Girder::pointAlong(c=%g) i=%i side=%i nseg=%i (nv/4)=%i \n", c, i, side, nseg, pointRange.y-pointRange.x );
        return ip; 
    };
    virtual int sideToPath( int side, int* inds ) const override{
        int i0 = pointRange.x;
        int n  = pointRange.y-i0;
        //if((n>1000)||(n<=0)){printf( "Girder::sideToPath() n=%i seems wrong\n", n ); exit(0);   }
        for(int i=0; i<n; i++){ inds[i] = i0+4*i+side; }
        return n;
    }

};


//class Ring : public NodeLinker { public:
class Ring : public StructuralComponent { public:
    //int p1; // anchor node; p0 inherate
    //double length;
    BodyPose pose;
    int nseg;
    double R;
    Vec2d wh;
    Quat4i st;
    //Material * material;
    //Vec2i pointRange;  // index of start and end in Truss
    //Vec2i stickRange; // --,,---
    //GirderType * type = NULL;

    virtual void print(bool bShort=false)const override { if(bShort){printf("Ring(id=%i)",id);}else{
        printf( "Ring(%i) n(%i) pos(%g,%g,%g) rot(%g,%g,%g)(%g,%g,%g)(%g,%g,%g) R=%g wh(%g,%g) st(%i,%i,%i,%i)\n", id, nseg,
        pose.pos.x,   pose.pos.y,   pose.pos.z,
        pose.rot.a.x, pose.rot.a.y, pose.rot.a.z,
        pose.rot.b.x, pose.rot.b.y, pose.rot.b.z,
        pose.rot.c.x, pose.rot.c.y, pose.rot.c.z,
        R, wh.x, wh.y, st.x,st.y,st.z,st.w );}
    }
    virtual int component_kind(){ return (int)ComponetKind::Ring; };

    virtual double rotMat( Mat3d& rot)const override{ rot=Mat3dIdentity; return 0; };
    virtual int nearSide  ( Vec3d p, const Mat3d* rot=0 )const override{ return -1; };
    virtual int pointAlong( double c, int side, Vec3d* pout=0 )const override{ return -1; };

    virtual int sideToPath( int side, int* inds ) const override{
        int i0 = pointRange.x;
        int n  = pointRange.y-i0;
        //if((n>1000)||(n<=0)){printf( "Girder::sideToPath() n=%i seems wrong\n", n ); exit(0);   }
        for(int i=0; i<n; i++){ inds[i] = i0+i; }
        return n;
    }
};

class Rope : public NodeLinker { public:
    double thick;
    int nseg;
    //Material * material;

    virtual void print(bool bShort=false)const override { if(bShort){printf("Rope(id=%i)",id);}else{
        printf("Rope(id=%i) nodes(%i,%i) L=%g \n", id, nodes.x->id,nodes.y->id, length ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Rope; };

    virtual double rotMat( Mat3d& rot)const override{ rot.c = nodes.y->pos - nodes.x->pos; double l=rot.c.normalize(); return l; };
    virtual int nearSide  ( Vec3d p, const Mat3d* rot=0 ) const override{ return 0; };
    virtual int pointAlong( double c, int side, Vec3d* pout=0 )const override{ 
        int i = (int)(c*nseg+0.5);
        if(pout){
            Vec3d d = (nodes.y->pos - nodes.x->pos)*(1./nseg);
            *pout = nodes.x->pos + d*i;
        }
        return i; 
    };
    virtual int sideToPath( int side, int* inds ) const override{
        int i0 = pointRange.x;
        int n  = pointRange.y-i0;
        //if((n>1000)||(n<=0)){printf( "Rope::sideToPath() n=%i seems wrong\n", n ); exit(0);   }
        for(int i=0; i<n; i++){ inds[i] = i0+4*i+side; }
        return n;
    }
};

class Modul: public ShipComponent { public:
    BodyPose pose;
    Box      bbox;
    Vec3d    span;
    double   volume;

    void pick(const Vec3d& ro, const Vec3d& rd){}
    virtual void print(bool bShort=false)const override{ if(bShort){printf("Modul(id=%i)",id);}else{
        printf("Modul(id=%i) kind=%i face_mat=%i \n", id, kind, face_mat ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Modul; }; 
};

class Tank : public Modul { public:
    int commodityId;
	//Commodity * typ;
	//double radius;
	//double length;
	        // [m^3]
	double filled;         // [1]

    virtual void print(bool bShort=false)const override { if(bShort){printf("Tank(id=%i)",id);}else{
        printf("Tank(id=%i) %g [m^3] of Commodity[%i] \n", id, filled, commodityId ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Tank; }; 
};


class Balloon : public Modul { public:
    virtual void print(bool bShort=false)const  override { if(bShort){printf("Balloon(id=%i)",id);}else{
        printf("Balloon(id=%i) kind=%i face_mat=%i \n", id, kind, face_mat ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Balloon; };
};

class Rock : public Modul { public:
    virtual void print(bool bShort=false)const override { if(bShort){printf("Rock(id=%i)",id);}else{
        printf("Rock(id=%i) kind=%i face_mat=%i \n", id, kind, face_mat ); };
    }
    virtual int component_kind(){ return (int)ComponetKind::Rock; };
};

class Pipe : public ShipComponent { public:
    double maxFlow;   // units depend on commodity
    //int  npath; // number of vertex along path
    //int* path;  // indexes of vertexes along the path
	Path path;
    ShipComponent * a;
	ShipComponent * b;
    virtual void print(bool bShort=false)const override{ if(bShort){printf("Pipe(id=%i)",id);}else{
        printf("Pipe(id=%i) Comp_a(kind=%i,id=%i) Comp_b(kind=%i,id=%i) nvert=%i maxFlow=%g \n", id, kind, a->kind,a->id, b->kind,b->id, path.n, maxFlow ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Rope; };
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

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Plate(id=%i)",id);}else{
        printf("Plate(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f))  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Plate; };
};

class Radiator : public Plate{ public:
    double temperature;
    virtual void print(bool bShort=false)const override{ if(bShort){printf("Radiator(id=%i)",id);}else{
        printf("Radiator(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f)) T=%g  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y, temperature ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Radiator; };
};

class Shield : public Plate{ public:
    virtual void print(bool bShort=false)const override{ if(bShort){printf("Shield(id=%i)",id);}else{
        printf("Shield(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f))  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Shield; };
};

class Collector : public Plate{ public:
    virtual void print(bool bShort=false)const override{ if(bShort){printf("Collector(id=%i)",id);}else{
        printf("Collector(id=%i) kidn=%i face_mat=%i Girders(%i(%4.2f,%4.2f),%i(%4.2f,%4.2f))  \n", id, kind, face_mat, g1,g1span.x,g1span.y, g2,g2span.x,g2span.y ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Collector; };
};

// ==== Motors

class Thruster : public Modul { public:
	//ThrusterType * typ = NULL;
    int    type;
	double thrust;
	double power;
	double consumption;

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Thruster(id=%i)",id);}else{
        printf("Thruster(id=%i) type=%i kidn=%i face_mat=%i P=%g[W] F=%g C=%g \n", id, type, kind, face_mat, power, thrust, consumption ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Thruster; };
};

class Rotor : public ShipComponent { public:
    // Move ship componenets with respect to each other in inner manuevers
    double Radius;
    double power;    //  [kg.m^2]
    double torque;   //  [kg.m^2]
    double Inertia;  //  [kg.m^2]  moment of inertia

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Rotor(id=%i)",id);}else{
        printf("Rotor(id=%i) kidn=%i face_mat=%i I=%g P=%g[W] tq=%g \n", id, kind, face_mat, Inertia, power, torque ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Rotor; };
};

/*
class Slider : public ShipComponent { public:
    // allow slide a node over a girder
    int  girder;
    double power;    //  [kg.m^2]
    double Force;

    virtual void print()const override{ printf("Slider(id=%i) kidn=%i face_mat=%i girder=%i P=%g[W] F=%g \n", id, kind, face_mat, girder, power, Force ); };
    virtual int component_kind(){ return (int)ComponetKind::Slider; };
};
*/

//ToDo : We need to move wheel with respect to girder


//class Slider : public Node { public:    // If we make Slider child of Node we can easily generate it somewhere along a girder  
class Slider : public Node { public:
    //int  ifix;    // to which vertex it is anchored
    Path path;
    double range;    // how much deflection of the slider perpedicular to the edge on loop this slider allows? 
    double forceMax; // max force which ca ne exerted by this slider
    double powerMax;

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Slider(id=%i)",id);}else{
        printf("Slider(id=%i) kind=%i girder=%i F=%g[N] P=%g[W] \n", id, kind, forceMax, powerMax ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Slider; };


    void updatePath(){
        printf("Slider::updatePath()\n");
        print();
        printf(" - boundTo:"); boundTo->print();

        /*
        //o->ifix = getVertAlong(o->comp2, o->along.y, o->sides.y );
        //int t1 = o->comp1->component_kind();
        printf( "updateSliderPaths[%i] t=%i | Girder=%i Ring=%i Rope=%i\n", io, t1, (int)ComponetKind::Girder, (int)ComponetKind::Ring, (int)ComponetKind::Rope );
        int i0 = o->comp1->pointRange.x;
        int n  = o->comp1->pointRange.y-i0;
        if((n>1000)||(n<=0)){printf( "updateSliderPaths() n=%i seems wrong\n", n ); exit(0);   }
        printf("SpaceCraft::updateSliderPaths() i0=%i n=%i\n", i0, n );
        if( (t1 == (int)ComponetKind::Girder) || (t1 == (int)ComponetKind::Ring) ){
            if((t1 == (int)ComponetKind::Ring)) o->path.closed=true;
            n/=4;
            o->path.realloc(n);
            for(int i=0; i<n; i++){ o->path.ps[i] = i0+4*i+o->sides.x; }
        }else if (t1 == (int)ComponetKind::Rope){
            o->path.realloc(n);
            for(int i=0; i<n; i++){ o->path.ps[i] = i0+i; }
        }
        */
    }

};


/*
class Slider_old : public Node { public:
    // allow slide a node over a girder
    // this slider moves one vertex (fixed point) which respect to vertex-loop (e.g. girder,wheel,rope). It will interpolate the position along current edge on the vertex loop. It will apply force to the two vertexes of the current edge.
    StructuralComponent* comp1;  // Warrning - this becomes invalid when arrays are re-allocated !!!!
    StructuralComponent* comp2;
    Vec2d along;
    Vec2i sides;
    int  ifix;    // to which vertex it is anchored
    //int  iloop;   // index along vetex loop
    //int  npath;   // number of points int he vetrex loop
    //int* path;   // vertex loop - typically along girder, wheel or rope
    //Vec3d dir;  // orientation of slinding
    Path path;
    double range;    // how much deflection of the slider perpedicular to the edge on loop this slider allows? 
    double forceMax; // max force which ca ne exerted by this slider
    double powerMax;

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Slider(id=%i)",id);}else{
        printf("Slider(id=%i) kidn=%i girder=%i P=%g[W] F=%g \n", id, kind, ifix, forceMax, powerMax ); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Slider; };
};
*/

// === Guns


// Also Accelerator?

class Accelerator : public ShipComponent{ public:
    // TODO: can be also attached to Ring ?
    //       This can be perhaps determined from type

    int   suppType; // support type - ring or girder
    int   suppId;   // support Id   - anchor which girder or ring ?
    Vec2d suppSpan; // pos along the support  (girder or ring)

    Path path; // along which vertexes it goes?

    double lenght;        // [m]
    double PowerPeak;     // [W]
    double PulseEnergy;   // [J]
    double PulseDuration; // [s]
    double PulsePeriod;   // [s]

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Accelerator(id=%i)",id);}else{
        printf("Accelerator(id=%i) kidn=%i face_mat=%i supp(%i(%4.2f,%4.2f)) L=%g P=%g[W] E=%g[J] ts(%g,%g)[s] \n", id, kind, face_mat, suppType,suppSpan.x,suppSpan.y, lenght, PowerPeak, PulseEnergy, PulseDuration, PulsePeriod); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Accelerator; };
    //std::vector<int> anchors; // anchor points
};


class Gun : public Accelerator{ public:
	int  gunId;             // index of gun type in catalog
    double Aperture;        // [m^2]
    double divergence;      // [1] tangens of angle
    // attached to girder?
    // incorporated in girder?

    virtual void print(bool bShort=false)const override{ if(bShort){printf("Gun(id=%i)",id);}else{
        printf("Gun(id=%i) kidn=%i face_mat=%i supp(%i(%4.2f,%4.2f)) A=%g[m^2] D=%g[1] L=%g[m] P=%g[W] E=%g[J] ts(%g,%g)[s] \n", id, kind, face_mat, suppType,suppSpan.x,suppSpan.y, Aperture, divergence, lenght, PowerPeak, PulseEnergy, PulseDuration, PulsePeriod); }
    }
    virtual int component_kind(){ return (int)ComponetKind::Gun; };
};


// ============= SpaceCraftWorkshop

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
        //Lua::dumpStack(L);
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
