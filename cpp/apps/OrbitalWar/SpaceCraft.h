
#ifndef  SpaceCraft_h
#define  SpaceCraft_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"

#include "Truss.h"

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

// ==== Materials

class CatalogItem{ public:
	int  id;
	int  kind;
	char name[NAME_LEN] = "\n";
};

class Material : public CatalogItem { public:
	double density;      // [kg/m3]
	double Spull,Spush;  // [Pa] Strenght
	double Kpull,Kpush;  // [Pa] elastic modulus
	double reflectivity; // [1] reflectivity
	double Tmelt;        // [T] temperature of failure
 	// What about heat, electricity etc. ?
};

class Commodity : public CatalogItem { public:
	double density;      // [kg/m3]
	// What about heat, electricity etc. ?
};

class FuelType : public Commodity { public:
    double EnergyDesity;   // [J/kg]
};

// ==== Components

class Node{ public:
    Vec3d pos;
    std::vector<Vec2i> components; // {kind,index}
    Node(Vec3d pos):pos(pos){};
};

class ShipComponent{ public:
    int    id;
    int    kind;
    int    shape;
    //int    p0; // anchor node
    // char name[NAME_LEN];
	double mass;           // [kg]
	//RigidBody pose;
};

class Modul: public ShipComponent { public:
    BodyPose pose;
    Box      bbox;
    Vec3d    span;
    double   volume;

    void pick(const Vec3d& ro, const Vec3d& rd){

    }

};

class Tank : public Modul { public:
	Commodity * typ;
	//double radius;
	//double length;
	        // [m^3]
	double filled;         // [1]
};


class Balloon : public Modul { public:

};

class Rock : public Modul { public:

};



//class GirderType : public CatalogItem { public:
//    int mseg;
//    Vec3d wh;
//};


class NodeLinker : public ShipComponent { public:
    int p0,p1;
    double length;
};

class Girder : public NodeLinker { public:
    //int p1; // anchor node; p0 inherate
    //double length;
    int nseg;
    int mseg;
    Vec2d wh;  // [m] width and height
    Vec3d up;
    //double SPull,SPush;
    //double kPull,kPush;
    Material * material;
    Vec2i poitRange;  // index of start and end in Truss
    Vec2i stickRange; // --,,---
    //GirderType * type = NULL;
};

//class Ring : public NodeLinker { public:
class Ring : public ShipComponent { public:
    //int p1; // anchor node; p0 inherate
    //double length;
    BodyPose pose;
    int nseg;
    double R;
    Vec2d wh;
    Material * material;
    Vec2i poitRange;  // index of start and end in Truss
    Vec2i stickRange; // --,,---
    //GirderType * type = NULL;
};

class Rope : public NodeLinker { public:
    double thick;
    Material * material;
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
    Vec2f g1span; // pos along girdes
    Vec2f g2span;
    //Vec3d normal;
	//int ntris;
	//int * tris;  // triangles from points of spaceship
};

class Radiator : public Plate{ public:
    double temperature;
};

class Shield : public Plate{ public:
};

class Collector : public Plate{ public:
};

// ==== Motors

class ThrusterType : public CatalogItem { public:
	double efficiency;    // how much power is used for acceleration of propelant
	double veMin;         // minimal exhaust velocity [m/s]
	double veMax;         // maximal exhaust velocity [m/s]
	bool   exhaustFuel;   // if true the burned fuel is added to propellant mass
	FuelType  * fuel      = NULL;
	Commodity  * Propelant = NULL;
};

class Thruster : public Modul { public:
	ThrusterType * typ = NULL;
	double thrust;
	double power;
	double consumption;
};

class Rotor : public ShipComponent { public:
    // Move ship componenets with respect to each other in inner manuevers
    double Radius;
    double power;    //  [kg.m^2]
    double torque;   //  [kg.m^2]
    double Inertia;  //  [kg.m^2]  moment of inertia
};

class Slider : public ShipComponent { public:
    // allow slide a note over a girder
    int  girder;
    double power;    //  [kg.m^2]
    double Force;
};

// === Guns

class GunType : public CatalogItem { public:
	double recoil;
	// scaling laws - how performace (power, accuracy, penetration, time of flight ...) scales with size ?
};

// Also Accelerator?

class Accelerator : public ShipComponent{ public:
    // TODO: can be also attached to Ring ?
    //       This can be perhaps determined from type

    int   suppType; // ring or girder
    int   suppId;   // anchor girders
    Vec2f suppSpan; // pos along girdes

    double lenght;        // [m]
    double PowerPeak;     // [W]
    double PulseEnergy;   // [J]
    double PulseDuration; // [s]
    double PulsePerios;   // [s]

    //std::vector<int> anchors; // anchor points
};


class Gun : public Accelerator{ public:
	GunType * typ = NULL;

    double Aperture;      // [m^2]
    double divergence;    // [1] tangens of angle
    // attached to girder?
    // incorporated in girder?
};

// === SpaceShip

class SpaceCraft : public CatalogItem { public:

    std::vector<int>       LODs;
    Truss truss;
    std::vector<Node>      nodes;
    std::vector<Rope>      ropes;
    std::vector<Girder>    girders;
    std::vector<Ring>      rings;
	std::vector<Thruster>  thrusters;
	std::vector<Gun>       guns;
	std::vector<Radiator>  radiators;
	std::vector<Shield>    shields;
	std::vector<Tank>      tanks;
	std::vector<Pipe>      pipes;

	std::vector<Balloon>  balloons;
	std::vector<Rock>     rocks;
	// Truss * coarse = NULL;
	// Truss * fine   = NULL;

	void clear(){
        nodes.clear(); ropes.clear(); girders.clear(); rings.clear(); thrusters.clear(); guns.clear(); radiators.clear(); shields.clear(); tanks.clear(); pipes.clear();
        rocks.clear(); balloons.clear();
        truss.clear();
	};

	void pick(Vec3d ro, Vec3d rd){

	}

	void toTruss(Truss& truss){
        int i=0;
        truss.newBlock();
        int ip0 = truss.points.size();
        for(Node o: nodes){
            truss.points.push_back( o.pos );
        }
        for(Rope o: ropes){
            truss.edges.push_back( (TrussEdge){o.p0,o.p1,0} );
        }
        for(Girder o: girders){
            //printf("DEBUG toTruss : girder #%i \n", i);
            truss.girder1( nodes[o.p0].pos, nodes[o.p1].pos, o.up, o.nseg, o.wh.a );
            truss.girder1_caps( o.p0, o.p1, 1 );

            Vec2i& bak = truss.blocks.back();
            o.poitRange  = {bak.x,truss.points.size()};
            o.stickRange = {bak.y,truss.edges .size()};
            i++;
        }
        i=0;
        for(Ring o: rings){
            printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o.nseg, o.wh.a );
            truss.wheel( o.pose.pos, o.pose.pos+o.pose.rot.b*o.R, o.pose.rot.c, o.nseg, o.wh.a );
            Vec2i& bak = truss.blocks.back();
            o.poitRange  = {bak.x,truss.points.size()};
            o.stickRange = {bak.y,truss.edges .size()};
            i++;
        }
        printf( "npoint %i nstick %i nblocks %i \n", truss.points.size(), truss.edges.size(), truss.blocks.size()  );
	}
	void toTruss(){ truss.clear(); toTruss(truss); };

};

} // namespace SpaceCrafting

#endif
