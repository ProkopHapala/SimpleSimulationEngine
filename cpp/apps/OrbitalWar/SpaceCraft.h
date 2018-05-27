
#ifndef  SpaceCraft_h
#define  SpaceCraft_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

namespace SpaceCrafting{

const int NAME_LEN = 16;

class BodyPose{ public:
    Vec3d pos;
    //Quat4d  rot;
    Mat3d rot;
};

class BBox : public BodyPose{ public:
    //int type;     //  elipsoide, box, ... ?
    Vec3d pmin;
    Vec3d pmax;
};

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
    double volume;
};

class Tank : public Modul { public:
	Commodity * typ;
	double radius;
	double length;
	        // [m^3]
	double filled;         // [1]
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

class TrusterType : public CatalogItem { public:
	double efficiency;    // how much power is used for acceleration of propelant
	double veMin;         // minimal exhaust velocity [m/s]
	double veMax;         // maximal exhaust velocity [m/s]
	bool   exhaustFuel;   // if true the burned fuel is added to propellant mass
	FuelType  * fuel      = NULL;
	Commodity  * Propelant = NULL;
};

class Truster : public ShipComponent { public:
	TrusterType * typ = NULL;
	double thrust;
	double power;
	double consumption;
	double mass;
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
};

class Gun : public ShipComponent{ public:
	GunType * typ = NULL;
    double Aperture;      // [m^2]
    double lenght;        // [m]
    double PowerPeak;     // [W]
    double PulseEnergy;   // [J]
    double PulseDuration; // [s]

    std::vector<int> anchors; // anchor points
};

// === SpaceShip

class SpaceCraft : public CatalogItem { public:
    std::vector<Node>      nodes;
    std::vector<Rope>      ropes;
    std::vector<Girder>    girders;
	std::vector<Truster>   thrustes;
	std::vector<Gun>       guns;
	std::vector<Radiator>  radiators;
	std::vector<Shield>    shields;
	std::vector<Tank>      tanks;
	std::vector<Pipe>      pipes;
	// Truss * coarse = NULL;
	// Truss * fine   = NULL;

	void clear(){
        nodes.clear(); ropes.clear(); girders.clear(); thrustes.clear(); guns.clear(); radiators.clear(); shields.clear(); tanks.clear(); pipes.clear();
	};
};

} // namespace SpaceCrafting

#endif
