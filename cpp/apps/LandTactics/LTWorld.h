

#ifndef LTWorld_h
#define LTWorld_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
//#include "HashMap2D.h"
#include "Body2D.h"

#include "SimplexRuler.h"
#include "Ruler2DFast.h"
#include "TerrainHydraulics.h"


//#include "TileBuffer2D.h"
//#include "Ruler2DFast.h"
//#include "TerrainCubic.h"

#include "LTUnit.h"
#include "LTShelter.h"

class LTMapSquare{ public:
    std::vector<LTStaticObject*> objects;
    //vector<LTUnit>  objects;
};

class LTWorld{
	public:
    //std::vector<BattleLine*> battleLines;
    std::vector<LTUnit*>        units;
    std::vector<LTFaction*>     factions;
    std::vector<LTObjectType*> objectTypes;
    //std::vector<LTShelter>     shelters;

    SimplexRuler       ruler;
    Ruler2DFast        square_ruler;
    TerrainHydraulics  hydraulics;
    double * ground       = NULL;

    //LTStaticObject * objects = NULL;
    std::vector<LTStaticObject>  objects;
    LTMapSquare    * squares    = NULL;

    Vec2d map_center;
    double maxHeight = 500.0;

    static const int   ntg = 4;
    static const int nAngles = 64;
    double  tgs    [ntg] = {-0.2,-0.1,-0.0,0.1};
    double  Ttgs   [ntg];
    Vec2d   hitcontours[nAngles*ntg];

    static constexpr int nTypesMax = 16;
    LTUnitType      unitTypes[nTypesMax];
    int             nTypes = 0;

    double RmaxInteract = 1.5;

    //Ruler2DFast colruler;
    //TerrainCubic   terrain;

    int nSoldierInteractions = 0;
    int nSoldiers = 0;

	//HashMap2D<Particle2D> map;

    double dt_frame  = 0.1;
    int    per_frame = 1;
    double damping   = 0.2;
    double damp;
    double dt;

    void init();
    void initStaticObject();
    void update();
    void simulationStep( double dt );
    int  getUnitAt( const Vec2d& p, LTFaction * faction );

    inline void setSimParams( double dt_frame_, double per_frame_, double damping_ ){
        dt_frame  = dt_frame_;
        per_frame = per_frame_;
        damping   = damping_;
        evalAuxSimParams();
        //printf( " dt_frame, per_frame,  dt, damp " );
    }

    inline void evalAuxSimParams(){
        const double dampMin = 0.5;
        dt   = dt_frame / per_frame;
        damp = ( 1 - damping * dt );
        if( damp < dampMin ){ damp = dampMin; }
    };

};

#endif
