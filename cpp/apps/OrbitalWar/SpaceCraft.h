
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
    Vec2d g1span; // pos along girdes
    Vec2d g2span;
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


enum class ComponetKind:int{ Node, Rope, Girder, Ring, Thruster, Gun, Radiator, Shield, Tank, Pipe, Balloon, Rock };

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

	// aux
	int pickedTyp = -1;




	void clear(){
        nodes.clear(); ropes.clear(); girders.clear(); rings.clear(); thrusters.clear(); guns.clear(); radiators.clear(); shields.clear(); tanks.clear(); pipes.clear();
        rocks.clear(); balloons.clear();
        truss.clear();
	};

	inline void linker2line( const NodeLinker& o, Line3d& l ){ l.a=nodes[o.p0].pos; l.b=nodes[o.p1].pos; }
	inline void plate2quad ( const Plate& o, Quad3d& qd ){
        //qd.l1.a=nodes[ girders[o.g1].p0 ].pos;
        //qd.l1.b=nodes[ girders[o.g1].p1 ].pos;
        //qd.l2.a=nodes[ girders[o.g2].p0 ].pos;
        //qd.l2.b=nodes[ girders[o.g2].p1 ].pos;
        //Vec3d d;
        //d = p01-p00; qd.p01=qd.p00+d*o.g1span.x;  qd.p00.add_mul( d,o.g1span.y);
        //d = p11-p10; qd.p11=qd.p10+d*o.g2span.x;  qd.p10.add_mul( d,o.g2span.y);
        qd.l1.fromSubLine( nodes[ girders[o.g1].p0 ].pos, nodes[ girders[o.g1].p1 ].pos, o.g1span.x, o.g1span.y );
        qd.l2.fromSubLine( nodes[ girders[o.g2].p0 ].pos, nodes[ girders[o.g2].p1 ].pos, o.g2span.x, o.g2span.y );
	}

    inline double rayPlate( const Plate& plate, const Vec3d& ro, const Vec3d& rd, Vec3d& normal, const Vec3d& hX, const Vec3d& hY ){
        Quad3d qd;
        //getQuad(qd);
        plate2quad(plate, qd);
        return qd.ray(ro, rd, normal, hX, hY );
	}

    inline double rayLinkLine( const NodeLinker& o, const Vec3d& ro, const Vec3d& rd, double rmax ){
        Vec3d lp0=nodes[o.p0].pos;
        Vec3d lpd=nodes[o.p1].pos-lp0;
        double lmax = lpd.normalize();
        double t,l;
        double r = rayLine(ro, rd, lp0, lpd, t, l );
        //printf( "rayLinkLine t,l,r %g %g %g \n", t,l,r );
        if( (r>rmax) || (l>lmax) || (l<0) ) return t_inf;
        return t;
	}

	int pick(Vec3d ro, Vec3d rd, double rmax = 1.0 ){
        Quad3d qd;
        Vec3d  nr;
        Vec3d hX,hY;
        rd.getSomeOrtho(hX,hY);
        int imin=-1;
        pickedTyp=-1;
        double tmin = t_inf;
        for(int i=0; i<radiators.size(); i++){
            //plate2quad( plate, qd );
            double t = rayPlate( radiators[i], ro, rd, nr, hX, hY );
            //printf( "pick.radiators[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Radiator; }
        };
        for(int i=0; i<shields.size(); i++){
            double t = rayPlate( shields[i], ro, rd, nr, hX, hY );
            //printf( "pick.shields[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Shield; }
        }
        //double rmax = 5.0;
        for(int i=0; i<girders.size(); i++){
            double t = rayLinkLine( girders[i], ro, rd, rmax );
            //printf( "pick.girders[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Girder; }
        }
        for(int i=0; i<ropes.size(); i++){
            double t = rayLinkLine( ropes[i], ro, rd, rmax );
            //printf( "pick.girders[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Rope; }
        }
        //return tmin<0.9*t_inf;
        return imin;
	}

	ShipComponent* getPicked( int picked ){
	    switch( (ComponetKind)pickedTyp){
            case ComponetKind::Radiator : return &radiators[picked]; break;
            case ComponetKind::Shield   : return &shields  [picked]; break;
            case ComponetKind::Girder   : return &girders  [picked]; break;
            case ComponetKind::Rope     : return &ropes    [picked]; break;
        }
        return 0;
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
            //printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o.nseg, o.wh.a );
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
